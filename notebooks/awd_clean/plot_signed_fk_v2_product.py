#!/usr/bin/env python3
"""Plot one corrected signed-lag chunk with branch-specific null tests."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from aggregate_seasonal_signed_fk_v2 import BRANCH_SIGN, VELOCITIES, metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stem", type=Path, help="Product stem without .npz/.json")
    args = parser.parse_args()
    product = np.load(args.stem.with_suffix(".npz"))
    metadata = json.loads(args.stem.with_suffix(".json").read_text())
    lags = product["lags"]
    distance = product["distance"]
    tops = {mode: product[f"{mode}_top"] for mode in BRANCH_SIGN}
    results = {
        mode: metrics(tops[mode], lags, distance, sign, permutations=2000)
        for mode, sign in BRANCH_SIGN.items()
    }
    common = np.concatenate([top.ravel() for top in tops.values()])
    limit = float(np.nanpercentile(np.abs(common), 98.5))

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), constrained_layout=True)
    for ax, mode in zip(axes[0], ("negative", "positive")):
        image = ax.imshow(
            tops[mode],
            extent=[lags[0], lags[-1], distance[-1], distance[0]],
            aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit,
            interpolation="nearest",
        )
        reference = np.r_[0.0, distance]
        ax.plot(BRANCH_SIGN[mode] * reference / 3200.0, reference, "k--", lw=1.25)
        ax.scatter(0.0, 0.0, marker="*", s=110, color="#ffd84d", edgecolor="black",
                   linewidth=0.8, zorder=5, clip_on=False)
        result = results[mode]
        ax.set_title(
            f"{mode.capitalize()} branch: {result['peak_velocity_m_s']/1000:.3f} km/s, "
            f"p={result['p_peak']:.4g}"
        )
        ax.set_xlim(-0.25, 0.25)
        ax.set_ylim(700, -15)
        ax.set_xlabel("Correlation lag (s)")
        ax.set_ylabel("Receiver offset from channel 0 (m)")
        ax.grid(alpha=0.15, lw=0.4)
    colorbar = fig.colorbar(image, ax=axes[0].tolist(), shrink=0.9, pad=0.02)
    colorbar.set_label("Normalized correlation (common display scale)")

    ax = axes[1, 0]
    for mode, color in (("negative", "tab:blue"), ("positive", "tab:orange")):
        result = results[mode]
        ax.plot(VELOCITIES / 1000.0, np.abs(result["physical_scores"]), color=color,
                lw=2, label=f"{mode}: physical lag")
        ax.axhline(result["null95"], color=color, ls="--", lw=1,
                   label=f"{mode}: null 95%")
    ax.axvline(3.2, color="0.25", ls=":", lw=1.2)
    ax.set_xlabel("Trial apparent velocity (km/s)")
    ax.set_ylabel("Absolute median correlation")
    ax.set_title("Signed-lag velocity scores")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    labels = ["Negative\nphysical", "Negative\nleakage", "Positive\nphysical", "Positive\nleakage"]
    values = [
        results["negative"]["absolute_score_3200"],
        results["negative"]["absolute_opposite_lag_leakage_3200"],
        results["positive"]["absolute_score_3200"],
        results["positive"]["absolute_opposite_lag_leakage_3200"],
    ]
    ax.bar(np.arange(4), values, color=["tab:blue", "#9ecae1", "tab:orange", "#fdd0a2"])
    ax.axhline(results["negative"]["null95"], color="tab:blue", ls="--", lw=1)
    ax.axhline(results["positive"]["null95"], color="tab:orange", ls="--", lw=1)
    ax.set_xticks(np.arange(4), labels)
    ax.set_ylabel("Absolute median correlation at 3.2 km/s")
    ax.set_title("Physical branches versus opposite-lag leakage")
    ax.grid(axis="y", alpha=0.2)

    fig.suptitle(
        f"Corrected signed-lag pilot: {metadata['date']}, {metadata['used_files']} files",
        fontsize=14,
    )
    output = args.stem.with_name(args.stem.name + "_diagnostic.png")
    fig.savefig(output, dpi=350)
    summary = {
        mode: {
            key: value for key, value in result.items()
            if key not in {"physical_scores", "leakage_scores", "null"}
        }
        for mode, result in results.items()
    }
    args.stem.with_name(args.stem.name + "_diagnostic.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
