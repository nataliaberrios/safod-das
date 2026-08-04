#!/usr/bin/env python3
"""Aggregate completed seasonal F-K chunks without pooling raw days first."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "ambient_transfer"
MANIFEST = OUT / "seasonal_day_selection.json"
RNG = np.random.default_rng(20260803)

def velocity_scores(top, lags, dist):
    vv = np.linspace(1200, 6000, 193)
    scores = np.array([np.nanmedian([r[np.argmin(abs(lags-d/v))]
                                     for r, d in zip(top, dist)]) for v in vv])
    return vv, scores

def season(date):
    m = int(date[5:7])
    return "DJF" if m in (12,1,2) else "MAM" if m in (3,4,5) else "JJA" if m in (6,7,8) else "SON"

def combine_day(date, n):
    sums = {m: None for m in ("negative", "positive", "both")}
    total = 0; lags = dist = vv = None
    for start in range(0, n, 10):
        count = min(10, n-start)
        js = OUT / f"fk_controls_{date}_start{start}_n{count}.json"
        ns = OUT / f"fk_controls_{date}_start{start}_n{count}.npz"
        if not js.exists() or not ns.exists():
            continue
        meta = json.loads(js.read_text())
        used = int(meta[next(iter(meta))]["used"])
        z = np.load(ns)
        if lags is None: lags, dist, vv = z["lags"], z["dist"], z["vv"]
        for m in sums:
            a = z[f"{m}_top"]
            sums[m] = a * used if sums[m] is None else sums[m] + a * used
        total += used
    if total == 0: return None
    results = {"date": date, "season": season(date), "used_files": total,
               "norm_seconds": 5.0, "passband_hz": [5.0, 20.0], "modes": {}}
    for m, s in sums.items():
        top = s / total; vv, scores = velocity_scores(top, lags, dist)
        null = np.array([np.nanmax(velocity_scores(top, lags, RNG.permutation(dist))[1]) for _ in range(2000)])
        ip = int(np.nanargmax(scores)); i32 = int(np.argmin(abs(vv-3200)))
        results["modes"][m] = {"peak_v_mps": float(vv[ip]), "peak_score": float(scores[ip]),
            "score_3200": float(scores[i32]), "null95": float(np.quantile(null,.95)),
            "p_peak": float((1 + np.sum(null >= scores[ip])) / (len(null)+1))}
        np.savez(OUT / f"fk_seasonal_{date}_{m}.npz", top=top, lags=lags, dist=dist, vv=vv, scores=scores)
    results["files_complete"] = total >= 1400
    (OUT / f"fk_seasonal_{date}.json").write_text(json.dumps(results, indent=2))
    return results, {m: s/total for m,s in sums.items()}, lags, dist, vv

def main():
    days = json.loads(MANIFEST.read_text())["days"]
    reports=[]; tops=[]; weights=[]; lags=dist=vv=None
    for item in days:
        out = combine_day(item["date"], int(item["nfiles"]))
        if out is None: continue
        r, t, lags, dist, vv = out; reports.append(r); tops.append(t); weights.append(r["used_files"])
    (OUT/"seasonal_fk_day_reports.json").write_text(json.dumps(reports, indent=2))
    if not tops: raise SystemExit("No completed seasonal days")
    # Across-day aggregate is computed only after day-level products exist.
    aggregate={m:sum(t[m]*w for t,w in zip(tops,weights))/sum(weights) for m in tops[0]}
    agg_report={"days": [r["date"] for r in reports], "n_days": len(reports),
                "weighted_files": int(sum(weights)), "norm_seconds": 5.0,
                "passband_hz": [5.0,20.0], "modes": {}}
    for m, top in aggregate.items():
        vv, scores=velocity_scores(top,lags,dist); null=np.array([np.nanmax(velocity_scores(top,lags,RNG.permutation(dist))[1]) for _ in range(5000)])
        ip=int(np.nanargmax(scores)); i32=int(np.argmin(abs(vv-3200)))
        agg_report["modes"][m]={"peak_v_mps":float(vv[ip]),"peak_score":float(scores[ip]),"score_3200":float(scores[i32]),"null95":float(np.quantile(null,.95)),"p_peak":float((1+np.sum(null>=scores[ip]))/(len(null)+1))}
        np.savez(OUT/f"fk_seasonal_aggregate_{m}.npz",top=top,lags=lags,dist=dist,vv=vv,scores=scores)
    (OUT/"fk_seasonal_aggregate.json").write_text(json.dumps(agg_report, indent=2))
    seasons=[r["season"] for r in reports]; order=["DJF","MAM","JJA","SON"]
    fig,axs=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    x=np.arange(len(reports)); colors=[{"DJF":"#2c7bb6","MAM":"#1a9850","JJA":"#fdae61","SON":"#d7191c"}[s] for s in seasons]
    axs[0].scatter(x,[r["modes"]["negative"]["peak_v_mps"]/1000 for r in reports],c=colors,s=65)
    axs[0].axhline(3.2,color="k",ls="--",lw=1,label="3.2 km/s reference"); axs[0].set_ylabel("Peak apparent velocity (km/s)"); axs[0].set_xticks(x, [r["date"] for r in reports], rotation=55, ha="right"); axs[0].legend(frameon=False)
    axs[1].scatter(x,[r["modes"]["negative"]["score_3200"] for r in reports],c=colors,s=65,label="negative")
    axs[1].scatter(x,[r["modes"]["positive"]["score_3200"] for r in reports],facecolors="none",edgecolors=colors,s=65,label="positive")
    axs[1].axhline(0,color="0.5",lw=.8); axs[1].set_ylabel("Median normalized correlation at 3.2 km/s"); axs[1].set_xticks(x, [r["date"] for r in reports], rotation=55, ha="right"); axs[1].legend(frameon=False)
    fig.savefig(OUT/"fk_seasonal_day_comparison.png",dpi=400); fig.savefig(OUT/"fk_seasonal_day_comparison.pdf")
    print(json.dumps(agg_report, indent=2))

if __name__ == "__main__": main()
