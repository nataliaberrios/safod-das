"""Conditional test of the Lellouch 500--520 m anomaly hypothesis.

Working assumption: current channel 0 and 1.020952-m spacing inherit the
Lellouch same-fiber position/depth convention.  This script uses available
ten-minute negative-branch F-K control products, measures the row-maximum
moveout, compares early and late apparent velocities, and scans for a
descriptive two-segment breakpoint.  It is not a lithology inversion.
"""
from pathlib import Path
import json
import re

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import theilslopes


HERE = Path(__file__).resolve().parent
PRODUCTS = HERE / "ambient_transfer"
OUT_FIG = HERE / "fk_geology_anomaly_test.png"
OUT_JSON = HERE / "fk_geology_anomaly_test.json"
EARLY = (50.0, 350.0)
LATE = (450.0, 650.0)
BREAKPOINTS = np.arange(250.0, 551.0, 25.0)


def fit_velocity(distance, lag, lo, hi):
    q = (distance >= lo) & (distance <= hi) & np.isfinite(lag)
    if q.sum() < 4:
        return {"velocity_m_s": None, "intercept_ms": None, "n": int(q.sum())}
    fit = theilslopes(lag[q], distance[q])
    return {
        "velocity_m_s": float(1.0 / fit.slope),
        "intercept_ms": float(1000.0 * fit.intercept),
        "n": int(q.sum()),
    }


def main():
    records = []
    for path in sorted(PRODUCTS.glob("fk_controls_*.npz")):
        match = re.search(r"fk_controls_(\d{4}-\d{2}-\d{2})_start(\d+)_n10", path.name)
        if match is None:
            continue
        z = np.load(path)
        distance = np.asarray(z["dist"], float)
        lag = np.asarray(z["lags"], float)[np.argmax(z["negative_top"], axis=1)]
        q = (distance <= 650.0) & np.isfinite(lag)
        if q.sum() < 10:
            continue
        d, t = distance[q], lag[q]
        early = fit_velocity(d, t, *EARLY)
        late = fit_velocity(d, t, *LATE)

        # Descriptive breakpoint scan, using least-squares residual reduction.
        one = np.polyfit(d, t, 1)
        one_sse = float(np.sum((t - np.polyval(one, d)) ** 2))
        improvements = []
        for breakpoint in BREAKPOINTS:
            left, right = d <= breakpoint, d > breakpoint
            if left.sum() < 4 or right.sum() < 4:
                improvements.append(np.nan)
                continue
            p_left, p_right = np.polyfit(d[left], t[left], 1), np.polyfit(d[right], t[right], 1)
            two_sse = np.sum((t[left] - np.polyval(p_left, d[left])) ** 2)
            two_sse += np.sum((t[right] - np.polyval(p_right, d[right])) ** 2)
            improvements.append(float(1.0 - two_sse / one_sse))
        best = int(np.nanargmax(improvements))
        records.append({
            "date": match.group(1),
            "start_minute": int(match.group(2)),
            "early": early,
            "late": late,
            "best_breakpoint_m": float(BREAKPOINTS[best]),
            "two_segment_sse_reduction": float(improvements[best]),
            "ridge_lag_ms": (1000.0 * t).tolist(),
            "distance_m": d.tolist(),
        })

    dates = sorted({r["date"] for r in records})
    by_date = {}
    for date in dates:
        subset = [r for r in records if r["date"] == date]
        early = np.array([r["early"]["velocity_m_s"] for r in subset], float)
        late = np.array([r["late"]["velocity_m_s"] for r in subset], float)
        bp = np.array([r["best_breakpoint_m"] for r in subset], float)
        by_date[date] = {
            "n_blocks": len(subset),
            "early_median_m_s": float(np.nanmedian(early)),
            "late_median_m_s": float(np.nanmedian(late)),
            "late_minus_early_median_m_s": float(np.nanmedian(late - early)),
            "breakpoint_median_m": float(np.nanmedian(bp)),
            "breakpoint_q16_q84_m": [float(np.nanquantile(bp, 0.16)), float(np.nanquantile(bp, 0.84))],
        }

    summary = {
        "assumption": "Current channel 0 and 1.020952-m spacing inherit the Lellouch same-cemented-fiber position/depth convention.",
        "hypothesis_interval_m": [500.0, 520.0],
        "products": str(PRODUCTS),
        "n_products": len(records),
        "dates": dates,
        "by_date": by_date,
        "interpretation": "Repeatable early/late curvature supports a common transition-like feature, but breakpoint estimates are broader and shallower than 500--520 m; this does not identify the Lellouch anomaly uniquely.",
        "limitations": [
            "Negative F-K branch and row maxima are conditional on the processing selection.",
            "Ten-minute products are not independent samples; date-level summaries are descriptive.",
            "The breakpoint scan is not a lithology or Vp inversion.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), constrained_layout=True)
    colors = plt.get_cmap("tab10").colors
    for i, date in enumerate(dates):
        subset = [r for r in records if r["date"] == date]
        lag_stack = np.array([r["ridge_lag_ms"] for r in subset])
        d = np.array(subset[0]["distance_m"])
        axes[0].plot(d, np.nanmedian(lag_stack, axis=0), color=colors[i % 10], lw=2, label=f"{date} (n={len(subset)})")
        early = [r["early"]["velocity_m_s"] / 1000 for r in subset]
        late = [r["late"]["velocity_m_s"] / 1000 for r in subset]
        axes[1].scatter(np.full(len(early), i - 0.13), early, color=colors[i % 10], s=12, alpha=0.35)
        axes[1].scatter(np.full(len(late), i + 0.13), late, color=colors[i % 10], s=12, alpha=0.35, marker="s")
        axes[1].plot([i - 0.13, i + 0.13], [np.median(early), np.median(late)], color=colors[i % 10], lw=2)
        bp = [r["best_breakpoint_m"] for r in subset]
        axes[2].scatter(np.full(len(bp), i), bp, color=colors[i % 10], s=12, alpha=0.35)
        axes[2].plot(i, np.median(bp), "o", color=colors[i % 10], ms=8)

    axes[0].axvspan(500, 520, color="tab:red", alpha=0.18, label="Lellouch anomaly hypothesis")
    axes[0].set(xlabel="Inherited measured-depth coordinate (m)", ylabel="Median ridge lag (ms)", title="Repeated ridge curvature")
    axes[0].legend(fontsize=7)
    axes[1].axhline(3.2, color="k", ls="--", lw=1, label="3.2 km/s reference")
    axes[1].set_xticks(range(len(dates)), [d.replace("-", "\n") for d in dates], fontsize=8)
    axes[1].set(xlabel="Date", ylabel="Apparent velocity (km/s)", title="Early vs late moveout")
    axes[1].legend(fontsize=8)
    axes[2].axhspan(500, 520, color="tab:red", alpha=0.18)
    axes[2].set_xticks(range(len(dates)), [d.replace("-", "\n") for d in dates], fontsize=8)
    axes[2].set(xlabel="Date", ylabel="Best two-segment breakpoint (m)", title="Breakpoint scan")
    fig.suptitle("Conditional test of the SAFOD 500–520 m anomaly hypothesis")
    fig.savefig(OUT_FIG, dpi=250)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
