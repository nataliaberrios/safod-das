#!/usr/bin/env python3
"""Exploratory filter-kernel residual diagnostic for the ambient F-K result.

The earlier full-pipeline surrogates inherit a velocity-corridor-compatible
correlation kernel from the selected F-K operator. This script estimates that
kernel as the mean surrogate correlation section, subtracts it from the
observed section, and uses leave-one-out subtraction for each surrogate.

This is a secondary sensitivity diagnostic, not an independent confirmation:
the mean artificial-surrogate kernel need not equal the observation's filter
kernel. Channel-permutation and independent-time-shift results are therefore
reported separately, boundary peaks are flagged, and disagreement prevents
promotion as physical-wave evidence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ambient_signed_fk_v2 import BRANCH_LAG_SIGN, VELOCITIES, velocity_scores


METHODS = ("channel_permutation", "circular_time_shift")
BRANCHES = ("negative", "positive")
WEDGE_M_S = (2500.0, 4500.0)


def empirical_upper_p(observed: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= observed)) / (null.size + 1))


def load_products(input_dir: Path) -> tuple[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]:
    records = []
    for path in sorted(input_dir.glob("fk_full_pipeline_null_v2_*_r*-*.npz")):
        z = np.load(path, allow_pickle=False)
        rid = int(z["null_realization_ids"][0])
        records.append((rid, {key: np.asarray(z[key]) for key in z.files}))
        z.close()
    records.sort(key=lambda item: item[0])
    ids = [item[0] for item in records]
    if len(records) < 3 or ids != list(range(len(records))):
        raise ValueError(f"expected contiguous null IDs from zero; found {ids}")

    first = records[0][1]
    common = {
        "lags_s": first["lags_s"],
        "receiver_offsets_m": first["receiver_offsets_m"],
        "velocities_m_s": first["velocities_m_s"],
    }
    if not np.allclose(common["velocities_m_s"], VELOCITIES):
        raise ValueError("saved velocities differ from ambient_signed_fk_v2.VELOCITIES")
    observed_samples = {
        branch: np.stack(
            [product[f"observed_{branch}_top"] for _, product in records]
        )
        for branch in BRANCHES
    }
    observed = {
        branch: np.mean(observed_samples[branch], axis=0) for branch in BRANCHES
    }
    observed_max_difference = max(
        float(np.max(np.ptp(observed_samples[branch], axis=0)))
        for branch in BRANCHES
    )
    if observed_max_difference > 1e-4:
        raise ValueError(
            "distributed observed sections differ by more than 1e-4: "
            f"{observed_max_difference}"
        )
    nulls = {
        method: {
            branch: np.stack(
                [
                    product[f"null_{method}_{branch}_top_first"]
                    for _, product in records
                ]
            )
            for branch in BRANCHES
        }
        for method in METHODS
    }
    for _, product in records[1:]:
        for key, value in common.items():
            if not np.allclose(product[key], value):
                raise ValueError(f"distributed products disagree for {key}")
    common["null_realization_ids"] = np.asarray(ids, dtype=int)
    common["maximum_distributed_observed_top_difference"] = np.asarray(
        observed_max_difference
    )
    common.update({f"observed_{key}_top": value for key, value in observed.items()})
    return common, nulls


def branch_scores(top: np.ndarray, lags: np.ndarray, offsets: np.ndarray, branch: str) -> np.ndarray:
    return velocity_scores(top, lags, offsets, BRANCH_LAG_SIGN[branch])


def positive_peak(scores: np.ndarray) -> tuple[float, float, bool]:
    inside = (VELOCITIES >= WEDGE_M_S[0]) & (VELOCITIES <= WEDGE_M_S[1])
    indices = np.flatnonzero(inside)
    peak_index = int(indices[np.nanargmax(scores[inside])])
    near_boundary = peak_index <= indices[2] or peak_index >= indices[-3]
    return float(scores[peak_index]), float(VELOCITIES[peak_index]), bool(near_boundary)


def run(input_dir: Path, output_stem: Path) -> dict:
    common, null_tops = load_products(input_dir)
    lags = common["lags_s"]
    offsets = common["receiver_offsets_m"]
    index_3200 = int(np.argmin(abs(VELOCITIES - 3200.0)))
    arrays = dict(common)
    report = {
        "workflow_version": "ambient_fk_filter_kernel_residual_v1",
        "input_dir": str(input_dir),
        "null_realizations": int(common["null_realization_ids"].size),
        "maximum_distributed_observed_top_difference": float(
            common["maximum_distributed_observed_top_difference"]
        ),
        "residual_definition": (
            "observed minus mean surrogate top; each surrogate minus the "
            "leave-one-out mean of the remaining surrogate tops"
        ),
        "selection_statistic": (
            "maximum signed-positive median correlation inside 2.5-4.5 km/s; "
            "familywise maximum across both signed branches"
        ),
        "methods": {},
        "decision": (
            "Exploratory only. Null-family disagreement or a near-boundary peak "
            "prevents promotion as physical-wave evidence."
        ),
    }

    for method in METHODS:
        method_report = {"branches": {}}
        observed_family = -np.inf
        null_family = np.full(common["null_realization_ids"].size, -np.inf)
        total = {
            branch: np.sum(null_tops[method][branch], axis=0)
            for branch in BRANCHES
        }
        for branch in BRANCHES:
            null_stack = null_tops[method][branch]
            kernel = total[branch] / null_stack.shape[0]
            observed_residual = common[f"observed_{branch}_top"] - kernel
            observed_scores = branch_scores(observed_residual, lags, offsets, branch)
            observed_peak, observed_velocity, boundary = positive_peak(observed_scores)
            observed_family = max(observed_family, observed_peak)

            null_scores = []
            null_peaks = []
            for i in range(null_stack.shape[0]):
                leave_one_out = (total[branch] - null_stack[i]) / (null_stack.shape[0] - 1)
                residual = null_stack[i] - leave_one_out
                scores = branch_scores(residual, lags, offsets, branch)
                null_scores.append(scores)
                null_peaks.append(positive_peak(scores)[0])
            null_scores = np.asarray(null_scores)
            null_peaks = np.asarray(null_peaks)
            null_family = np.maximum(null_family, null_peaks)

            arrays[f"{method}_{branch}_kernel_top"] = kernel
            arrays[f"{method}_{branch}_observed_residual_top"] = observed_residual
            arrays[f"{method}_{branch}_observed_scores"] = observed_scores
            arrays[f"{method}_{branch}_null_scores"] = null_scores
            arrays[f"{method}_{branch}_null_peaks"] = null_peaks
            method_report["branches"][branch] = {
                "observed_positive_peak": observed_peak,
                "observed_peak_velocity_m_s": observed_velocity,
                "peak_is_near_search_boundary": boundary,
                "observed_score_3200": float(observed_scores[index_3200]),
                "null95_positive_peak": float(np.quantile(null_peaks, 0.95)),
                "p_upper": empirical_upper_p(observed_peak, null_peaks),
            }

        arrays[f"{method}_null_familywise"] = null_family
        method_report["familywise"] = {
            "observed_maximum_across_branches": float(observed_family),
            "null95": float(np.quantile(null_family, 0.95)),
            "p_upper": empirical_upper_p(observed_family, null_family),
        }
        report["methods"][method] = method_report

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    output_stem.with_suffix(".json").write_text(json.dumps(report, indent=2))
    np.savez_compressed(output_stem.with_suffix(".npz"), **arrays)

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.5), sharex=True, constrained_layout=True)
    for row, method in enumerate(METHODS):
        for col, branch in enumerate(BRANCHES):
            ax = axes[row, col]
            observed = arrays[f"{method}_{branch}_observed_scores"]
            null_scores = arrays[f"{method}_{branch}_null_scores"]
            ax.fill_between(
                VELOCITIES / 1000.0,
                np.quantile(null_scores, 0.05, axis=0),
                np.quantile(null_scores, 0.95, axis=0),
                color="0.85",
                label="leave-one-out null 5-95%",
            )
            ax.plot(VELOCITIES / 1000.0, observed, color="tab:blue", lw=2.0, label="observed residual")
            ax.axvspan(2.5, 4.5, color="tab:orange", alpha=0.08)
            ax.axvline(3.2, color="0.2", ls=":", lw=1.4)
            values = report["methods"][method]["branches"][branch]
            ax.set_title(
                f"{method.replace('_', ' ')}, {branch}: "
                f"peak {values['observed_peak_velocity_m_s']/1000:.2f} km/s, "
                f"p={values['p_upper']:.3f}"
            )
            ax.set_ylabel("Residual median correlation")
            ax.grid(alpha=0.2)
    for ax in axes[-1]:
        ax.set_xlabel("Trial apparent velocity (km s$^{-1}$)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Exploratory filter-kernel residual diagnostic")
    fig.savefig(output_stem.with_name(output_stem.name + "_figure").with_suffix(".png"), dpi=320)
    fig.savefig(output_stem.with_name(output_stem.name + "_figure").with_suffix(".pdf"))
    plt.close(fig)
    print(json.dumps(report, indent=2))
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path(__file__).resolve().parent
        / "ambient_transfer"
        / "fk_full_pipeline_null_v2_n300_r20",
    )
    parser.add_argument(
        "--output-stem",
        type=Path,
        default=Path(__file__).resolve().parent
        / "ambient_transfer"
        / "fk_filter_kernel_residual_v1"
        / "fk_filter_kernel_residual_v1",
    )
    args = parser.parse_args()
    run(args.input_dir, args.output_stem)


if __name__ == "__main__":
    main()
