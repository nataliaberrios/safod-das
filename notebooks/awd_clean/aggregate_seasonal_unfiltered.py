#!/usr/bin/env python3
"""Aggregate exact-resolution unfiltered seasonal transfer chunks."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent; OUT=ROOT/"ambient_transfer"; MANIFEST=OUT/"seasonal_day_selection.json"
RNG=np.random.default_rng(20260803)

def score(top,lags,dist):
    vv=np.linspace(1200,6000,193)
    s=np.array([np.nanmedian([r[np.argmin(abs(lags-d/v))] for r,d in zip(top,dist)]) for v in vv])
    return vv,s

def main():
    items=json.loads(MANIFEST.read_text())["days"]; reports=[]; tops=[]; weights=[]; lags=dist=None
    for item in items:
        date=item["date"]; n=int(item["nfiles"]); total=0; sm=None
        for start in range(0,n,10):
            count=min(10,n-start); js=OUT/f"transfer_{date}_start{start}_n{count}.json"; ns=OUT/f"transfer_{date}_start{start}_n{count}.npz"
            if not (js.exists() and ns.exists()): continue
            meta=json.loads(js.read_text()); used=int(meta["used_files"]); z=np.load(ns)
            if lags is None: lags,dist=z["lags"],z["distances"]
            sm=z["top_stack"]*used if sm is None else sm+z["top_stack"]*used; total+=used
        if total==0: continue
        top=sm/total; vv,s=score(top,lags,dist); null=np.array([np.nanmax(score(top,lags,RNG.permutation(dist))[1]) for _ in range(2000)])
        ip=int(np.nanargmax(s)); i32=int(np.argmin(abs(vv-3200)))
        r={"date":date,"used_files":total,"files_complete":total>=1400,"norm_seconds":5.0,"passband_hz":[5.0,20.0],"peak_v_mps":float(vv[ip]),"peak_score":float(s[ip]),"score_3200":float(s[i32]),"null95":float(np.quantile(null,.95)),"p_peak":float((1+np.sum(null>=s[ip]))/(len(null)+1))}
        reports.append(r); tops.append(top); weights.append(total)
        np.savez(OUT/f"transfer_seasonal_{date}.npz",top_stack=top,lags=lags,distances=dist,vv=vv,scores=s)
    (OUT/"seasonal_unfiltered_day_reports.json").write_text(json.dumps(reports,indent=2))
    if not reports: raise SystemExit("No complete unfiltered day products")
    top=sum(t*w for t,w in zip(tops,weights))/sum(weights); vv,s=score(top,lags,dist); null=np.array([np.nanmax(score(top,lags,RNG.permutation(dist))[1]) for _ in range(5000)])
    ip=int(np.nanargmax(s)); i32=int(np.argmin(abs(vv-3200)))
    agg={"days":[r["date"] for r in reports],"n_days":len(reports),"weighted_files":int(sum(weights)),"norm_seconds":5.0,"passband_hz":[5.0,20.0],"peak_v_mps":float(vv[ip]),"peak_score":float(s[ip]),"score_3200":float(s[i32]),"null95":float(np.quantile(null,.95)),"p_peak":float((1+np.sum(null>=s[ip]))/(len(null)+1))}
    (OUT/"seasonal_unfiltered_aggregate.json").write_text(json.dumps(agg,indent=2)); np.savez(OUT/"seasonal_unfiltered_aggregate.npz",top_stack=top,lags=lags,distances=dist,vv=vv,scores=s)
    print(json.dumps(agg,indent=2))

if __name__=="__main__": main()
