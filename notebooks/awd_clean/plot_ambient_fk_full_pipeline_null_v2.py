#!/usr/bin/env python3
"""Plot observed statistics against pre-F-K full-pipeline null surrogates.

The figure is intentionally diagnostic rather than declarative.  Each panel
shows the empirical cumulative distribution of a coherence-destroying null,
the observed statistic, and the one-sided empirical exceedance probability.
All panels share an x scale so branch and null-method differences remain
visually comparable.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


METHOD_LABELS = {
    "channel_permutation": "Channel permutation",
    "circular_time_shift": "Independent circular time shift",
}
STATISTICS = (
    ("negative", "Negative signed branch"),
    ("positive", "Positive signed branch"),
    ("familywise", "Maximum across both branches"),
)
METHOD_COLORS = {
    "channel_permutation": "#2C6E9B",
    "circular_time_shift": "#C96B32",
}


def empirical_p(observed: float, surrogate: np.ndarray) -> float:
    """Return the plus-one one-sided exceedance probability."""
    return float((1 + np.count_nonzero(surrogate >= observed)) / (surrogate.size + 1))


def p_text(value: float, n: int) -> str:
    """Format p without implying more precision than the null ensemble allows."""
    minimum = 1.0 / (n + 1)
    if np.isclose(value, minimum):
        return f"p = {value:.3f} (minimum)"
    return f"p = {value:.3f}"


def load_inputs(json_path: Path, npz_path: Path) -> tuple[dict, dict[str, np.ndarray]]:
    """Load and validate a paired aggregate report and numerical product."""
    report = json.loads(json_path.read_text())
    if report.get("workflow_version") != "ambient_fk_full_pipeline_null_v2_aggregate":
        raise ValueError(
            "expected workflow_version ambient_fk_full_pipeline_null_v2_aggregate"
        )

    with np.load(npz_path) as product:
        arrays = {name: product[name] for name in product.files}

    n_report = int(report["null_realizations"])
    ids = np.asarray(arrays["null_realization_ids"])
    if ids.size != n_report or np.unique(ids).size != n_report:
        raise ValueError("null realization IDs are missing, duplicated, or inconsistent")

    for method in report["null_methods"]:
        for statistic, _ in STATISTICS:
            key = (
                f"null_{method}_familywise_statistics"
                if statistic == "familywise"
                else f"null_{method}_{statistic}_selection_statistics"
            )
            if key not in arrays:
                raise KeyError(f"missing aggregate array: {key}")
            if np.asarray(arrays[key]).shape != (n_report,):
                raise ValueError(f"unexpected shape for {key}: {arrays[key].shape}")
    return report, arrays


def statistic_values(
    report: dict,
    arrays: dict[str, np.ndarray],
    method: str,
    statistic: str,
) -> tuple[float, np.ndarray, float]:
    """Return observed statistic, null ensemble, and recorded empirical p."""
    if statistic == "familywise":
        observed = float(report["observed_familywise_score"])
        null = np.asarray(arrays[f"null_{method}_familywise_statistics"], dtype=float)
        recorded_p = float(report["null_results"][method]["familywise"]["p"])
    else:
        branch = report["null_results"][method]["branches"][statistic]
        observed = float(branch["observed_peak_absolute_score_in_wedge"])
        null = np.asarray(
            arrays[f"null_{method}_{statistic}_selection_statistics"], dtype=float
        )
        recorded_p = float(branch["p_peak_in_wedge"])

    calculated_p = empirical_p(observed, null)
    if not np.isclose(calculated_p, recorded_p):
        raise ValueError(
            f"empirical p mismatch for {method}/{statistic}: "
            f"JSON={recorded_p}, recomputed={calculated_p}"
        )
    return observed, null, calculated_p


def plot_panel(
    ax: plt.Axes,
    observed: float,
    null: np.ndarray,
    p_value: float,
    color: str,
) -> None:
    """Draw an ECDF with individual surrogate values and observed statistic."""
    ordered = np.sort(null)
    y = np.arange(1, ordered.size + 1, dtype=float) / ordered.size
    q95 = float(np.quantile(ordered, 0.95))

    ax.step(
        np.r_[ordered[0], ordered],
        np.r_[0.0, y],
        where="post",
        color=color,
        linewidth=1.6,
        zorder=2,
    )
    ax.scatter(
        ordered,
        y,
        s=20,
        facecolor="white",
        edgecolor=color,
        linewidth=0.9,
        zorder=3,
    )
    ax.axvline(q95, color="#777777", linestyle=(0, (2, 2)), linewidth=1.1, zorder=1)
    ax.axvline(observed, color="#111111", linewidth=1.8, zorder=4)
    ax.scatter(
        [observed],
        [1.035],
        marker="v",
        s=42,
        color="#111111",
        clip_on=False,
        zorder=5,
    )
    ax.text(
        0.04,
        0.94,
        f"{p_text(p_value, null.size)}\nn = {null.size}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        linespacing=1.25,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.5},
        zorder=6,
    )
    ax.set_ylim(0.0, 1.08)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", length=3)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot an aggregate pre-F-K full-pipeline null test."
    )
    parser.add_argument("--json", type=Path, required=True, help="aggregate JSON report")
    parser.add_argument(
        "--npz",
        type=Path,
        help="paired aggregate NPZ (default: JSON path with .npz suffix)",
    )
    parser.add_argument(
        "--output-stem",
        type=Path,
        required=True,
        help="output path without extension",
    )
    args = parser.parse_args()

    npz_path = args.npz if args.npz is not None else args.json.with_suffix(".npz")
    report, arrays = load_inputs(args.json, npz_path)
    methods = list(report["null_methods"])
    unknown = [method for method in methods if method not in METHOD_LABELS]
    if unknown:
        raise ValueError(f"unsupported null methods: {unknown}")

    values: dict[tuple[str, str], tuple[float, np.ndarray, float]] = {}
    scale_values: list[np.ndarray] = []
    for method in methods:
        for statistic, _ in STATISTICS:
            item = statistic_values(report, arrays, method, statistic)
            values[(method, statistic)] = item
            scale_values.extend((item[1], np.asarray([item[0]])))

    all_values = np.concatenate(scale_values)
    span = float(np.ptp(all_values))
    pad = 0.06 * span if span > 0 else 0.05
    x_limits = (max(0.0, float(np.min(all_values)) - pad), float(np.max(all_values)) + pad)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(
        len(methods),
        len(STATISTICS),
        figsize=(11.2, 6.5),
        sharex=True,
        sharey=True,
        constrained_layout=False,
    )
    axes = np.atleast_2d(axes)

    for row, method in enumerate(methods):
        color = METHOD_COLORS[method]
        for col, (statistic, title) in enumerate(STATISTICS):
            ax = axes[row, col]
            observed, null, p_value = values[(method, statistic)]
            plot_panel(ax, observed, null, p_value, color)
            ax.set_xlim(*x_limits)
            if row == 0:
                ax.set_title(title, pad=8, fontweight="semibold")
            if col == 0:
                ax.set_ylabel(
                    f"{METHOD_LABELS[method]}\nEmpirical cumulative probability",
                    labelpad=8,
                )
            if row == len(methods) - 1:
                ax.set_xlabel("Peak absolute median-correlation statistic")

    date = report.get("date", "unknown date")
    used_files = report.get("used_files", "?")
    n_null = int(report["null_realizations"])
    fig.suptitle(
        "Pre–F–K full-pipeline null test",
        x=0.075,
        y=0.985,
        ha="left",
        fontsize=15,
        fontweight="semibold",
    )
    fig.text(
        0.075,
        0.946,
        (
            "Frozen 2.5–4.5 km s⁻¹ velocity wedge; all surrogates pass through "
            "the same F–K filtering, correlation, stacking, and velocity scan"
        ),
        ha="left",
        va="top",
        fontsize=9.5,
        color="#333333",
    )
    fig.text(
        0.925,
        0.982,
        f"{date}  |  {used_files} input files  |  {n_null} surrogates per method",
        ha="right",
        va="top",
        fontsize=8.5,
        color="#555555",
    )

    legend_items = [
        Line2D([0], [0], color="#555555", marker="o", markerfacecolor="white", lw=1.5,
               markersize=5, label="Surrogate ECDF and realizations"),
        Line2D([0], [0], color="#111111", lw=1.8, label="Observed statistic"),
        Line2D([0], [0], color="#777777", lw=1.1, linestyle=(0, (2, 2)),
               label="Surrogate 95th percentile"),
    ]
    fig.legend(
        handles=legend_items,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.053),
        ncol=3,
        frameon=False,
        fontsize=8.5,
        handlelength=2.5,
    )
    fig.text(
        0.075,
        0.012,
        (
            "One-sided empirical p counts surrogate statistics at least as large as observed, "
            "with a +1 correction. Higher statistics indicate stronger selected coherence. "
            f"With n = {n_null}, p-value resolution is 1/(n+1) = {1/(n_null+1):.3f}."
        ),
        ha="left",
        va="bottom",
        fontsize=8.1,
        color="#4A4A4A",
    )
    fig.subplots_adjust(left=0.105, right=0.975, top=0.875, bottom=0.17, wspace=0.16, hspace=0.22)

    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = args.output_stem.with_suffix(".png")
    pdf_path = args.output_stem.with_suffix(".pdf")
    fig.savefig(png_path, dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(png_path)
    print(pdf_path)


if __name__ == "__main__":
    main()
