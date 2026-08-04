"""Fit apparent moveout to the five-hour negative-branch F-K interim product.

This is a diagnostic for the Lellouch-style wiggle plot.  It reports the
row-maximum ridge, an intercept-aware constant-velocity fit, and a simple
two-segment fit.  The result is conditional on the F-K wedge and is not an
absolute Vp estimate.
"""
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import theilslopes


HERE = Path(__file__).resolve().parent
INPUT = HERE / "ambient_transfer" / "interim_2024-12-20_n300_fk_negative.npz"
OUT_FIG = HERE / "fk_velocity_moveout_diagnostic.png"
OUT_JSON = HERE / "fk_velocity_moveout_diagnostic.json"


def main():
    z = np.load(INPUT)
    distance_m = np.asarray(z["dist"], float)
    lag_s = np.asarray(z["lags"], float)
    panel = np.asarray(z["top"], float)
    ridge_idx = np.argmax(panel, axis=1)
    ridge_lag_s = lag_s[ridge_idx]

    # The 700-m row has an abrupt lag jump and is excluded from the primary
    # fit; it is retained in the diagnostic plot for transparency.
    usable = distance_m <= 650.0
    d = distance_m[usable]
    t = ridge_lag_s[usable]
    robust = theilslopes(t, d)
    v_robust = 1.0 / robust.slope
    # Least-squares intercept-aware fit, useful for comparison with candidate
    # reference velocities.
    ols_slope, ols_intercept = np.polyfit(d, t, 1)
    v_ols = 1.0 / ols_slope

    # A breakpoint near the visible bend is descriptive only.
    breakpoint_m = 350.0
    left = d <= breakpoint_m
    right = d > breakpoint_m
    left_fit = np.polyfit(d[left], t[left], 1)
    right_fit = np.polyfit(d[right], t[right], 1)
    v_left = 1.0 / left_fit[0]
    v_right = 1.0 / right_fit[0]

    candidate = np.arange(2400.0, 3301.0, 25.0)
    rms_ms = []
    intercept_ms = []
    for v in candidate:
        t0 = np.mean(t - d / v)
        intercept_ms.append(1000.0 * t0)
        rms_ms.append(1000.0 * np.sqrt(np.mean((t - (t0 + d / v)) ** 2)))
    best_i = int(np.argmin(rms_ms))

    summary = {
        "input": str(INPUT),
        "selection": "distance <= 650 m; 700 m row retained in plot but excluded from fit",
        "ridge_distance_m": distance_m.tolist(),
        "ridge_lag_ms": (1000.0 * ridge_lag_s).tolist(),
        "theil_sen_velocity_m_s": float(v_robust),
        "theil_sen_intercept_ms": float(1000.0 * robust.intercept),
        "ols_velocity_m_s": float(v_ols),
        "ols_intercept_ms": float(1000.0 * ols_intercept),
        "descriptive_breakpoint_m": breakpoint_m,
        "left_velocity_m_s": float(v_left),
        "right_velocity_m_s": float(v_right),
        "best_constant_velocity_m_s": float(candidate[best_i]),
        "best_constant_intercept_ms": float(intercept_ms[best_i]),
        "best_constant_rms_ms": float(rms_ms[best_i]),
        "caveat": "F-K-selected five-hour interim diagnostic; not an absolute formation Vp estimate.",
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    norm = np.max(np.abs(panel), axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    im = ax.imshow(
        panel / norm,
        extent=[1000 * lag_s[0], 1000 * lag_s[-1], distance_m[-1], distance_m[0]],
        aspect="auto",
        cmap="gray",
        vmin=-1,
        vmax=1,
    )
    ax.scatter(1000 * ridge_lag_s, distance_m, s=20, c="tab:red", label="row maximum")
    dline = np.linspace(0, 720, 300)
    ax.plot(1000 * dline / 3200, dline, "--", lw=1.5, label="3.2 km/s; zero intercept")
    ax.plot(1000 * (robust.intercept + dline / v_robust), dline, "-", lw=1.8,
            label=f"robust fit: {v_robust/1000:.2f} km/s")
    ax.set(xlabel="Lag (ms)", ylabel="Separation (m)", title="Negative F-K branch")
    ax.legend(fontsize=8, loc="lower right")

    bx.plot(candidate / 1000, rms_ms, "k-")
    bx.axvline(3.2, color="tab:blue", ls="--", label="Lellouch reference")
    bx.axvline(candidate[best_i] / 1000, color="tab:red", ls="-",
                label=f"best constant: {candidate[best_i]/1000:.2f} km/s")
    bx.set(xlabel="Candidate apparent velocity (km/s)", ylabel="RMS moveout residual (ms)",
           title="Intercept-aware velocity scan")
    bx.legend(fontsize=8)
    fig.suptitle("Preliminary velocity diagnostic: 2024-12-20, 300 one-minute files (5 h)")
    fig.savefig(OUT_FIG, dpi=220)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
