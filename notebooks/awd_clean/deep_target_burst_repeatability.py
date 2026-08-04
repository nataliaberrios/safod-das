"""Burst-level repeatability of the targeted Deep strand-depth observable.

Odd non-empty epochs select the fixed trajectory.  Every non-empty burst stack
is then evaluated independently at the selected slowness/time, separately for
the outbound and reversed-return legs.  This tests whether the target-depth
result is present burst-to-burst rather than only after epoch stacking.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binomtest

from deep_target_scan import (
    BANDS,
    CDZ_M,
    DX_M,
    PRE_S,
    SDZ_M,
    TIME_DECIMATION,
    TIME_WINDOW,
    TURNAROUND_CH,
    fixed_score,
    leg_section,
    prepare,
    scan,
    target_to_leg_coordinate,
    window_for_target,
)
from fk_dispersion import weighted_stack


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
TARGETS = (SDZ_M, CDZ_M)
CONTROL_DEPTHS = (2800.0, 3000.0)
N_BOOT = 5000
SEED = 20260802

OUT_CSV = HERE / "deep_target_burst_repeatability.csv"
OUT_NPZ = HERE / "deep_target_burst_repeatability.npz"
OUT_PNG = HERE / "deep_target_burst_repeatability.png"
OUT_PDF = HERE / "deep_target_burst_repeatability.pdf"
OUT_TXT = HERE / "deep_target_burst_repeatability.txt"
OUT_JSON = HERE / "deep_target_burst_repeatability.json"


def trajectory_rms(section, coordinate, fs, p, target_time):
    from deep_target_scan import trajectory_rms as _rms
    return _rms(section, coordinate, fs, p, target_time)


def select_fixed(stacks, counts, fs0, i0, i1, leg, depth, band):
    epochs = np.flatnonzero(counts > 0)
    odd = epochs[epochs % 2 == 0]
    odd_counts = counts.copy()
    odd_counts[epochs[epochs % 2 == 1]] = 0
    discovery = weighted_stack(stacks, odd_counts)
    section, _ = leg_section(discovery, leg, i0, i1)
    target_s = target_to_leg_coordinate(depth, leg)
    start, stop, width = window_for_target(section.shape[0], target_s, DX_M)
    coordinate = np.arange(start, stop, 2) * DX_M
    dsec = prepare(section[start:stop], fs0, band)[::2]
    result = scan(dsec, coordinate, fs0 / TIME_DECIMATION,
                  np.arange(0, i1 - i0, TIME_DECIMATION) / fs0 + TIME_WINDOW[0])
    return {
        "start": start, "stop": stop, "width": width,
        "window_start_m": start * DX_M, "window_stop_m": (stop - 1) * DX_M,
        "positive_p": result["positive"][1], "positive_t": result["positive"][2],
        "negative_p": result["negative"][1], "negative_t": result["negative"][2],
    }


def bootstrap_difference(values, rng):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, np.nan
    samples = rng.choice(values, size=(N_BOOT, len(values)), replace=True).mean(axis=1)
    return float(np.median(values)), float(np.percentile(samples, 2.5)), float(np.percentile(samples, 97.5))


def main():
    rng = np.random.default_rng(SEED)
    with np.load(STACKS) as data:
        stacks = data["deep_stacks"]
        counts = data["n_common"]
        fs0 = float(data["fs"])
        epochs = np.flatnonzero(counts > 0)
    i0 = int(round((PRE_S + TIME_WINDOW[0]) * fs0))
    i1 = int(round((PRE_S + TIME_WINDOW[1]) * fs0))
    time = np.arange(0, i1 - i0, TIME_DECIMATION) / fs0 + TIME_WINDOW[0]

    rows = []
    summaries = []
    fixed_by_key = {}
    depths = TARGETS + CONTROL_DEPTHS
    for depth in depths:
        for leg in ("outbound", "return"):
            for band in BANDS:
                fixed = select_fixed(stacks, counts, fs0, i0, i1, leg, depth, band)
                fixed_by_key[(depth, leg, band[0])] = fixed
                for epoch in epochs:
                    section, _ = leg_section(stacks[epoch], leg, i0, i1)
                    start, stop = fixed["start"], fixed["stop"]
                    coordinate = np.arange(start, stop, 2) * DX_M
                    filtered = prepare(section[start:stop], fs0, band)[::2]
                    pp = fixed_score(filtered, coordinate, fs0 / TIME_DECIMATION,
                                     time, fixed["positive_p"], fixed["positive_t"])
                    pn = fixed_score(filtered, coordinate, fs0 / TIME_DECIMATION,
                                     time, fixed["negative_p"], fixed["negative_t"])
                    rows.append({
                        "epoch": int(epoch), "n_common_drops": int(counts[epoch]),
                        "split": "discovery" if epoch % 2 == 0 else "validation",
                        "assumed_depth_m": float(depth), "leg": leg,
                        "band_low_hz": band[0], "band_high_hz": band[1],
                        "positive_semblance": float(pp), "negative_semblance": float(pn),
                        "positive_minus_negative": float(pp - pn),
                        "positive_slowness_ms_per_m": float(fixed["positive_p"] * 1e3),
                        "positive_velocity_mps": float(1.0 / fixed["positive_p"]),
                        "positive_time_s": float(fixed["positive_t"]),
                        "n_common_drops": int(counts[epoch]),
                    })

                q = [r for r in rows if r["assumed_depth_m"] == depth and r["leg"] == leg and r["band_low_hz"] == band[0]]
                for split in ("discovery", "validation", "all"):
                    z = [r for r in q if split == "all" or r["split"] == split]
                    diffs = np.asarray([r["positive_minus_negative"] for r in z])
                    positive_count = int(np.sum(diffs > 0))
                    p_sign = float(binomtest(positive_count, len(diffs), 0.5, alternative="greater").pvalue) if len(diffs) else np.nan
                    median, lo, hi = bootstrap_difference(diffs, rng)
                    summaries.append({
                        "assumed_depth_m": float(depth), "leg": leg,
                        "band_low_hz": band[0], "band_high_hz": band[1], "split": split,
                        "n_bursts": len(diffs), "n_positive_minus_negative": positive_count,
                        "positive_fraction": positive_count / len(diffs) if len(diffs) else np.nan,
                        "sign_test_p": p_sign, "median_difference": median,
                        "bootstrap_ci_low": lo, "bootstrap_ci_high": hi,
                    })

    with OUT_CSV.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    summary_csv = HERE / "deep_target_burst_repeatability_summary.csv"
    with summary_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0])); writer.writeheader(); writer.writerows(summaries)

    # Figure: each panel is one target depth and band; lines compare legs.
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True, sharey=True, constrained_layout=True)
    panel_specs = ((SDZ_M, (15.0, 30.0)), (CDZ_M, (15.0, 30.0)),
                   (SDZ_M, (3.0, 15.0)), (CDZ_M, (3.0, 15.0)))
    for ax, (depth, band) in zip(axes.flat, panel_specs):
        for leg, color, ls in (("outbound", "#1565c0", "-"), ("return", "#d95f02", "--")):
            q = [r for r in rows if r["assumed_depth_m"] == depth and r["leg"] == leg and r["band_low_hz"] == band[0]]
            q.sort(key=lambda r: r["epoch"])
            x = np.arange(len(q))
            y = np.asarray([r["positive_minus_negative"] for r in q])
            ax.plot(x, y, marker="o", ms=3.5, lw=1.2, ls=ls, color=color, label=leg)
            ax.scatter(x[np.array([r["split"] == "validation" for r in q])], y[np.array([r["split"] == "validation" for r in q])],
                       s=28, facecolors="white", edgecolors=color, zorder=4)
        ax.axhline(0, color="0.25", lw=0.8)
        ax.set_title(f"Assumed depth {depth:.0f} m; {band[0]:.0f}–{band[1]:.0f} Hz")
        ax.set_ylabel("positive − negative semblance")
        ax.grid(alpha=0.2)
        ax.legend(frameon=False, fontsize=8)
    axes[1, 0].set_xlabel("non-empty burst epoch")
    axes[1, 1].set_xlabel("non-empty burst epoch")
    fig.suptitle("Deep target-window burst repeatability under provisional depth mapping", fontsize=15, fontweight="bold")
    fig.savefig(OUT_PNG, dpi=300); fig.savefig(OUT_PDF); plt.close(fig)

    report = [
        "SAFOD AWD Deep target-window burst repeatability",
        "Fixed trajectories selected from even-numbered discovery epochs; all 46 non-empty burst stacks evaluated independently.",
        "Positive-minus-negative semblance is a directional repeatability statistic, not a calibrated amplitude or creep measurement.",
        "Primary band: 15–30 Hz; 3–15 Hz is the control. Assumed depths remain conditional on the provisional fiber mapping.",
        "",
        "Summary (validation split):",
    ]
    for s in summaries:
        if s["split"] == "validation":
            report.append(f"{s['assumed_depth_m']:.0f} m {s['leg']} {s['band_low_hz']:.0f}–{s['band_high_hz']:.0f} Hz: n={s['n_bursts']}, positive fraction={s['positive_fraction']:.3f}, sign p={s['sign_test_p']:.4g}, median difference={s['median_difference']:.4f} [{s['bootstrap_ci_low']:.4f},{s['bootstrap_ci_high']:.4f}]")
    report.extend(["", "Interpretation: burst-level support strengthens the active-source baseline; absence of burst-level support would limit the creep-monitoring framing. Neither outcome measures a creep rate from this 24-hour experiment."])
    OUT_TXT.write_text("\n".join(report) + "\n")
    OUT_JSON.write_text(json.dumps({"status": "burst_repeatability_conditional_depth", "n_epochs": int(len(epochs)), "target_depths_m": list(TARGETS), "control_depths_m": list(CONTROL_DEPTHS), "bands_hz": [list(b) for b in BANDS], "summary": summaries}, indent=2) + "\n")
    np.savez(OUT_NPZ, epoch=np.asarray([r["epoch"] for r in rows]), assumed_depth_m=np.asarray([r["assumed_depth_m"] for r in rows]), positive_minus_negative=np.asarray([r["positive_minus_negative"] for r in rows]), positive_semblance=np.asarray([r["positive_semblance"] for r in rows]), negative_semblance=np.asarray([r["negative_semblance"] for r in rows]))
    print(OUT_TXT.read_text())


if __name__ == "__main__":
    main()
