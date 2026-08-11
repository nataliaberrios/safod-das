#!/usr/bin/env python3
"""Aggregate corrected seasonal signed-lag chunks and make comparison figures."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
LEGACY = ROOT / "ambient_transfer"
OUT = LEGACY / "signed_lag_v2"
MANIFEST = LEGACY / "seasonal_day_selection.json"
VELOCITIES = np.linspace(1200.0, 6000.0, 193)
BRANCH_SIGN = {"negative": 1.0, "positive": -1.0}
RNG = np.random.default_rng(20260805)


def velocity_scores(top, lags, distance, lag_sign):
    return np.array([
        np.nanmedian([
            row[np.argmin(np.abs(lags - lag_sign * offset / velocity))]
            for row, offset in zip(top, distance)
        ])
        for velocity in VELOCITIES
    ])


def metrics(top, lags, distance, lag_sign, permutations=5000):
    physical = velocity_scores(top, lags, distance, lag_sign)
    leakage = velocity_scores(top, lags, distance, -lag_sign)
    null = np.array([
        np.nanmax(np.abs(velocity_scores(top, lags, RNG.permutation(distance), lag_sign)))
        for _ in range(permutations)
    ])
    peak = int(np.nanargmax(np.abs(physical)))
    index_3200 = int(np.argmin(abs(VELOCITIES - 3200.0)))
    return {
        "physical_lag_sign": int(lag_sign),
        "peak_velocity_m_s": float(VELOCITIES[peak]),
        "peak_signed_score": float(physical[peak]),
        "peak_absolute_score": float(abs(physical[peak])),
        "score_3200": float(physical[index_3200]),
        "absolute_score_3200": float(abs(physical[index_3200])),
        "opposite_lag_leakage_3200": float(leakage[index_3200]),
        "absolute_opposite_lag_leakage_3200": float(abs(leakage[index_3200])),
        "null95": float(np.quantile(null, 0.95)),
        "p_peak": float((1 + np.sum(null >= abs(physical[peak]))) / (len(null) + 1)),
        "physical_scores": physical,
        "leakage_scores": leakage,
        "null": null,
    }


def combine_day(date, nfiles):
    sums = {mode: None for mode in BRANCH_SIGN}
    total = 0
    lags = distance = None
    for start in range(0, nfiles, 10):
        count = min(10, nfiles - start)
        stem = OUT / f"signed_fk_v2_{date}_start{start}_n{count}"
        if not stem.with_suffix(".json").exists() or not stem.with_suffix(".npz").exists():
            continue
        metadata = json.loads(stem.with_suffix(".json").read_text())
        used = int(metadata["used_files"])
        product = np.load(stem.with_suffix(".npz"))
        if lags is None:
            lags = product["lags"]
            distance = product["distance"]
        for mode in sums:
            contribution = product[f"{mode}_top"] * used
            sums[mode] = contribution if sums[mode] is None else sums[mode] + contribution
        total += used
    if total == 0:
        return None
    tops = {mode: value / total for mode, value in sums.items()}
    report = {"date": date, "used_files": total, "branches": {}}
    arrays = {"lags": lags, "distance": distance, "velocities_m_s": VELOCITIES}
    for mode, lag_sign in BRANCH_SIGN.items():
        result = metrics(tops[mode], lags, distance, lag_sign, permutations=1000)
        report["branches"][mode] = {
            key: value for key, value in result.items()
            if key not in {"physical_scores", "leakage_scores", "null"}
        }
        arrays[f"{mode}_top"] = tops[mode]
        arrays[f"{mode}_physical_scores"] = result["physical_scores"]
        arrays[f"{mode}_leakage_scores"] = result["leakage_scores"]
    np.savez_compressed(OUT / f"seasonal_signed_fk_v2_{date}.npz", **arrays)
    (OUT / f"seasonal_signed_fk_v2_{date}.json").write_text(json.dumps(report, indent=2))
    return report, tops, lags, distance


def plot_aggregate(tops, lags, distance, results, weighted_files):
    common = np.concatenate([tops["negative"].ravel(), tops["positive"].ravel()])
    limit = float(np.nanpercentile(np.abs(common), 98.5))
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), constrained_layout=True)
    for ax, mode in zip(axes[0], ("negative", "positive")):
        sign = BRANCH_SIGN[mode]
        image = ax.imshow(
            tops[mode],
            extent=[lags[0], lags[-1], distance[-1], distance[0]],
            aspect="auto",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            interpolation="nearest",
        )
        reference_distance = np.r_[0.0, distance]
        peak_velocity = results[mode]["peak_velocity_m_s"]
        ax.plot(
            sign * reference_distance / peak_velocity,
            reference_distance,
            "k--",
            lw=1.4,
            label=f"score peak: {peak_velocity / 1000:.3f} km/s",
        )
        ax.plot(
            sign * reference_distance / 3200.0,
            reference_distance,
            color="0.4",
            ls=":",
            lw=1.0,
            label="fixed reference: 3.2 km/s",
        )
        ax.scatter(0.0, 0.0, marker="*", s=110, color="#ffd84d", edgecolor="black",
                   linewidth=0.8, zorder=5, clip_on=False)
        ax.annotate(
            "virtual source (channel 0)",
            (0.0, 0.0),
            xytext=(8, -14),
            textcoords="offset points",
            fontsize=8,
            va="top",
        )
        physical_lag = "positive" if sign > 0 else "negative"
        branch_symbol = "<" if mode == "negative" else ">"
        ax.set_title(f"F K {branch_symbol} 0 evaluated at {physical_lag} lag")
        legend_location = "lower left" if sign > 0 else "lower right"
        ax.legend(loc=legend_location, frameon=True, framealpha=0.88, fontsize=7)
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
                lw=2, label=f"{mode}: physical signed lag")
        ax.axhline(result["null95"], color=color, ls="--", lw=1,
                   label=f"{mode}: null 95%")
    ax.axvline(3.2, color="0.25", ls=":", lw=1.2)
    ax.set_xlabel("Trial apparent velocity (km/s)")
    ax.set_ylabel("Absolute median correlation")
    ax.set_title("Branch-specific signed-lag velocity scores")
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
    colors = ["tab:blue", "#9ecae1", "tab:orange", "#fdd0a2"]
    ax.bar(np.arange(4), values, color=colors)
    ax.axhline(results["negative"]["null95"], color="tab:blue", ls="--", lw=1)
    ax.axhline(results["positive"]["null95"], color="tab:orange", ls="--", lw=1)
    ax.set_xticks(np.arange(4), labels)
    ax.set_ylabel("Absolute median correlation at 3.2 km/s")
    ax.set_title("Physical-lag signals and opposite-lag leakage")
    ax.grid(axis="y", alpha=0.2)

    fig.suptitle(
        f"Corrected seasonal signed-lag ambient F-K comparison ({weighted_files:,} files)",
        fontsize=14,
        y=1.025,
    )
    fig.savefig(OUT / "seasonal_signed_fk_v2_aggregate.png", dpi=350,
                bbox_inches="tight")
    fig.savefig(OUT / "seasonal_signed_fk_v2_aggregate.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_days(day_products):
    all_values = np.concatenate([
        top.ravel() for _, tops, _, _ in day_products for top in tops.values()
    ])
    limit = float(np.nanpercentile(np.abs(all_values), 98.5))
    fig, axes = plt.subplots(8, 2, figsize=(11, 22), constrained_layout=True,
                             sharex=True, sharey=True)
    for row, (report, tops, lags, distance) in enumerate(day_products):
        for column, mode in enumerate(("negative", "positive")):
            ax = axes[row, column]
            ax.imshow(
                tops[mode],
                extent=[lags[0], lags[-1], distance[-1], distance[0]],
                aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit,
                interpolation="nearest",
            )
            ax.plot(BRANCH_SIGN[mode] * distance / 3200.0, distance, "k--", lw=0.9)
            ax.set_title(f"{report['date']} — {mode}", fontsize=9)
            ax.set_xlim(-0.25, 0.25)
            ax.set_ylim(700, 0)
            if column == 0:
                ax.set_ylabel("Offset (m)")
            if row == 7:
                ax.set_xlabel("Lag (s)")
    fig.suptitle("Corrected day-level causal and anti-causal signed F-K sections", fontsize=14)
    fig.savefig(OUT / "seasonal_signed_fk_v2_days.png", dpi=300)
    fig.savefig(OUT / "seasonal_signed_fk_v2_days.pdf")


def main():
    days = json.loads(MANIFEST.read_text())["days"]
    day_products = []
    for item in days:
        product = combine_day(item["date"], int(item["nfiles"]))
        if product is not None:
            day_products.append(product)
    if not day_products:
        raise SystemExit("No corrected seasonal chunks available")

    weighted_files = sum(report["used_files"] for report, _, _, _ in day_products)
    aggregate = {
        mode: sum(tops[mode] * report["used_files"] for report, tops, _, _ in day_products)
        / weighted_files
        for mode in BRANCH_SIGN
    }
    lags = day_products[0][2]
    distance = day_products[0][3]
    full_results = {
        mode: metrics(aggregate[mode], lags, distance, sign, permutations=5000)
        for mode, sign in BRANCH_SIGN.items()
    }
    report = {
        "workflow_version": "signed_lag_v2",
        "days": [entry[0]["date"] for entry in day_products],
        "weighted_files": weighted_files,
        "branch_lag_sign": BRANCH_SIGN,
        "branches": {
            mode: {
                key: value for key, value in result.items()
                if key not in {"physical_scores", "leakage_scores", "null"}
            }
            for mode, result in full_results.items()
        },
    }
    (OUT / "seasonal_signed_fk_v2_aggregate.json").write_text(json.dumps(report, indent=2))
    np.savez_compressed(
        OUT / "seasonal_signed_fk_v2_aggregate.npz",
        lags=lags,
        distance=distance,
        velocities_m_s=VELOCITIES,
        negative_top=aggregate["negative"],
        positive_top=aggregate["positive"],
        negative_physical_scores=full_results["negative"]["physical_scores"],
        positive_physical_scores=full_results["positive"]["physical_scores"],
        negative_leakage_scores=full_results["negative"]["leakage_scores"],
        positive_leakage_scores=full_results["positive"]["leakage_scores"],
    )
    plot_aggregate(aggregate, lags, distance, full_results, weighted_files)
    plot_days(day_products)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
