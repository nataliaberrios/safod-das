#!/usr/bin/env python3
"""Aggregate exact-stack ordered and matched-null products.

The decision statistic is the maximum fixed-top simple-stack envelope score
over the declared 1.5--5.0 km/s scan.  Taking the maximum in every null
realization accounts for selecting velocity after inspecting Figure 7c.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
FULL = ROOT / "ambient_transfer" / "lellouch2019_exact_stack"
NULLS = ROOT / "ambient_transfer" / "lellouch2019_exact_stack_nulls"
OUT = ROOT / "ambient_transfer" / "lellouch2019_exact_stack_summary"
MODES = ("unfiltered", "downgoing_2p5_4p5")
NULL_METHODS = ("white_noise", "channel_permutation", "circular_time_shift")


def records(directory: Path) -> list[dict]:
    output = []
    for path in sorted(directory.glob("lellouch_exact_*.json")):
        report = json.loads(path.read_text())
        report["_json"] = path
        report["_npz"] = path.with_suffix(".npz")
        if report["_npz"].is_file():
            output.append(report)
    return output


def empirical_p(observed: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= observed)) / (null.size + 1))


def main() -> None:
    available = records(FULL) + records(NULLS)
    five_hour = [item for item in available if item["used_files"] == 300]
    summary = {
        "workflow_version": "summarize_lellouch2019_exact_stack_v1",
        "decision_statistic": (
            "maximum fixed-top simple-stack envelope score over 1.5-5.0 km/s; "
            "velocity selection repeated independently for every null"
        ),
        "modes": {},
    }
    loaded = {}
    for mode in MODES:
        ordered_records = [
            item for item in five_hour
            if item["fk_mode"] == mode
            and item["null_method"] == "ordered"
            and item["spectral_mode"] == "cross_correlation"
            and not item["common_mode_removal"]
        ]
        if len(ordered_records) != 1:
            raise RuntimeError(
                f"expected one five-hour ordered {mode} product, found {len(ordered_records)}"
            )
        ordered_report = ordered_records[0]
        ordered = np.load(ordered_report["_npz"])
        observed_stat = float(np.max(ordered["fixed_simple_velocity_scores"]))
        observed_i = int(np.argmax(ordered["fixed_simple_velocity_scores"]))
        mode_report = {
            "ordered_product": str(ordered_report["_npz"]),
            "ordered_used_windows": ordered_report["used_windows"],
            "ordered_selected_velocity_m_s": float(ordered["velocity_grid_m_s"][observed_i]),
            "ordered_maximum_score": observed_stat,
            "nulls": {},
        }
        loaded[(mode, "ordered")] = ordered
        for method in NULL_METHODS:
            candidates = [
                item for item in five_hour
                if item["fk_mode"] == mode
                and item["null_method"] == method
                and item["spectral_mode"] == "cross_correlation"
                and not item["common_mode_removal"]
            ]
            candidates.sort(key=lambda item: item["null_realization"])
            if len(candidates) != 20:
                raise RuntimeError(
                    f"expected 20 {mode}/{method} products, found {len(candidates)}"
                )
            statistics = []
            selected_velocities = []
            curves = []
            for item in candidates:
                if item["used_windows"] != ordered_report["used_windows"]:
                    raise RuntimeError(f"window mismatch in {item['_json']}")
                product = np.load(item["_npz"])
                curve = product["fixed_simple_velocity_scores"]
                index = int(np.argmax(curve))
                statistics.append(float(curve[index]))
                selected_velocities.append(float(product["velocity_grid_m_s"][index]))
                curves.append(curve)
            statistics = np.asarray(statistics)
            mode_report["nulls"][method] = {
                "realizations": len(candidates),
                "null95_maximum_score": float(np.quantile(statistics, 0.95)),
                "empirical_p": empirical_p(observed_stat, statistics),
                "selected_velocity_m_s_median": float(np.median(selected_velocities)),
            }
            loaded[(mode, method)] = {
                "statistics": statistics,
                "curves": np.asarray(curves),
            }
        summary["modes"][mode] = mode_report

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "lellouch2019_exact_stack_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    fig, axes = plt.subplots(3, 2, figsize=(12.5, 12), constrained_layout=True)
    colors = {
        "white_noise": "tab:green",
        "channel_permutation": "tab:blue",
        "circular_time_shift": "tab:orange",
    }
    for column, mode in enumerate(MODES):
        ordered = loaded[(mode, "ordered")]
        lags = ordered["lags_s"]
        distance = ordered["fixed_coordinates_m"]
        section = ordered["fixed_simple"]
        limit = float(np.percentile(np.abs(section), 99.0))
        axes[0, column].imshow(
            section,
            extent=[lags[0], lags[-1], distance[-1], distance[0]],
            aspect="auto",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            interpolation="nearest",
        )
        axes[0, column].plot(distance / 3200.0, distance, "k--", lw=1.2)
        axes[0, column].set(
            xlim=(0.0, 0.35),
            xlabel="Correlation lag (s)",
            ylabel="Receiver coordinate (m)",
            title=f"{mode}: ordered five-hour Figure 7c stack",
        )
        velocity = ordered["velocity_grid_m_s"] / 1000.0
        axes[1, column].plot(
            velocity,
            ordered["fixed_simple_velocity_scores"],
            color="black",
            lw=2,
            label="ordered",
        )
        for method in NULL_METHODS:
            curves = loaded[(mode, method)]["curves"]
            lo, hi = np.quantile(curves, [0.05, 0.95], axis=0)
            axes[1, column].fill_between(
                velocity, lo, hi, color=colors[method], alpha=0.16, label=method
            )
        axes[1, column].axvline(3.2, color="0.35", ls=":")
        axes[1, column].set(
            xlabel="Trial apparent velocity (km/s)",
            ylabel="Envelope moveout score",
            title="Ordered curve and matched-null 5--95% ranges",
        )
        axes[1, column].legend(frameon=False, fontsize=8)
        observed = summary["modes"][mode]["ordered_maximum_score"]
        bins = np.linspace(
            0.0,
            max(
                observed,
                max(
                    np.max(loaded[(mode, method)]["statistics"])
                    for method in NULL_METHODS
                ),
            )
            * 1.08,
            18,
        )
        for method in NULL_METHODS:
            axes[2, column].hist(
                loaded[(mode, method)]["statistics"],
                bins=bins,
                histtype="step",
                lw=1.7,
                color=colors[method],
                label=method,
            )
        axes[2, column].axvline(observed, color="black", lw=2, label="ordered")
        axes[2, column].set(
            xlabel="Maximum score selected over velocity",
            ylabel="Null count",
            title="Familywise matched-null decision",
        )
        axes[2, column].legend(frameon=False, fontsize=8)
    fig.savefig(OUT / "lellouch2019_exact_stack_summary.png", dpi=350)
    fig.savefig(OUT / "lellouch2019_exact_stack_summary.pdf")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
