#!/usr/bin/env python3
"""Independently audit the real-noise ambient F-K injection recovery.

The audit reloads every per-file scenario array, reconstructs pre-filter and
post-filter statistics, repeats paired bootstraps with an independent random
seed, and checks the aggregate classifications.  It deliberately does not
trust the aggregate JSON as the source of the measured values.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from aggregate_ambient_fk_injection_recovery_v1 import (
    DIRECTIONS,
    VELOCITIES,
    VELOCITY_TOLERANCE_M_S,
    amplitude_tag,
    null_thresholds,
    paired_bootstrap_mean,
    paired_bootstrap_stack_score,
    physical_trajectory_values,
    scenario_paths,
)
from ambient_signed_fk_v2 import BRANCH_LAG_SIGN, velocity_scores


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "ambient_transfer" / "fk_injection_recovery_v1_n300"
DEFAULT_NULL = (
    ROOT
    / "ambient_transfer"
    / "fk_prefilter_energy_v1_n300_r20"
    / "fk_prefilter_energy_v1_aggregate.json"
)
EXPECTED_AMPLITUDES = [0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]


def close(left: float, right: float, atol: float = 1e-7) -> bool:
    return bool(np.isclose(left, right, rtol=1e-6, atol=atol))


def run(args: argparse.Namespace) -> Path:
    aggregate_path = args.input_dir / "ambient_fk_injection_recovery_v1_aggregate.json"
    if not aggregate_path.is_file():
        raise FileNotFoundError(aggregate_path)
    aggregate = json.loads(aggregate_path.read_text())
    aggregate_by_key = {
        (float(item["velocity_m_s"]), int(item["direction"])): item
        for item in aggregate["scenarios"]
    }
    expected_keys = {
        (velocity, direction) for velocity in VELOCITIES for direction in DIRECTIONS
    }
    thresholds = null_thresholds(args.null_json)
    checks = {
        "scenario_grid_complete": set(aggregate_by_key) == expected_keys,
        "all_file_lists_identical": True,
        "all_files_unique": True,
        "all_arrays_finite": True,
        "zero_injection_identical": True,
        "zero_injection_not_recovered": True,
        "saved_prefilter_metrics_reproduced": True,
        "saved_postfilter_metrics_reproduced": True,
        "saved_peak_velocities_reproduced": True,
        "independent_bootstrap_classifications_stable": True,
        "reported_minima_reproduced": True,
    }
    reference_files = None
    reference_zero = None
    audited_scenarios = []
    for velocity, direction, json_path, npz_path in scenario_paths(
        args.input_dir, args.date, args.start, args.nfiles
    ):
        if not json_path.is_file() or not npz_path.is_file():
            raise FileNotFoundError(f"missing {json_path} or {npz_path}")
        source_report = json.loads(json_path.read_text())
        product = np.load(npz_path)
        files = product["used_files"].astype(str)
        checks["all_files_unique"] &= len(set(files)) == args.expected_files
        if reference_files is None:
            reference_files = files
        else:
            checks["all_file_lists_identical"] &= np.array_equal(reference_files, files)
        if files.size != args.expected_files:
            checks["all_file_lists_identical"] = False
        amplitudes = product["amplitude_ratios"].astype(float).tolist()
        if amplitudes != EXPECTED_AMPLITUDES:
            raise ValueError(f"amplitude grid mismatch for {(velocity, direction)}")
        branch = "negative" if direction == 1 else "positive"
        zero_tag = amplitude_tag(0.0)
        zero_vectors = np.concatenate([
            product[f"{zero_tag}_{name}_prefilter_per_file"]
            for name in ("negative", "positive")
        ])
        if reference_zero is None:
            reference_zero = zero_vectors
        else:
            checks["zero_injection_identical"] &= np.allclose(
                reference_zero, zero_vectors, rtol=5e-5, atol=1e-8
            )
        baseline_pre = product[f"{zero_tag}_{branch}_prefilter_per_file"]
        baseline_corr = product[f"{zero_tag}_{branch}_correlation_per_file"]
        baseline_trajectory = physical_trajectory_values(
            baseline_corr,
            product["lags_s"],
            product["distance_m"],
            velocity,
            BRANCH_LAG_SIGN[branch],
        )
        saved = aggregate_by_key[(velocity, direction)]
        audited_points = []
        for point_index, amplitude in enumerate(amplitudes):
            tag = amplitude_tag(amplitude)
            pre = product[f"{tag}_{branch}_prefilter_per_file"]
            corr = product[f"{tag}_{branch}_correlation_per_file"]
            for array in (pre, corr):
                checks["all_arrays_finite"] &= bool(np.all(np.isfinite(array)))
            trajectory = physical_trajectory_values(
                corr,
                product["lags_s"],
                product["distance_m"],
                velocity,
                BRANCH_LAG_SIGN[branch],
            )
            rng_pre = np.random.default_rng(
                args.seed + int(velocity) * 10 + (0 if direction == 1 else 1) + point_index
            )
            rng_post = np.random.default_rng(
                args.seed + int(velocity) * 20 + (0 if direction == 1 else 1) + point_index
            )
            pre_uplift, pre_low, pre_high = paired_bootstrap_mean(
                pre, baseline_pre, rng_pre, args.bootstraps
            )
            post_uplift, post_low, post_high = paired_bootstrap_stack_score(
                trajectory, baseline_trajectory, rng_post, args.bootstraps
            )
            top = np.mean(corr, axis=0)
            physical = velocity_scores(
                top,
                product["lags_s"],
                product["distance_m"],
                BRANCH_LAG_SIGN[branch],
            )
            peak_velocity = float(product["velocities_m_s"][np.nanargmax(np.abs(physical))])
            pre_mean = float(np.mean(pre))
            exact_score = float(np.median(np.mean(trajectory, axis=0)))
            independent = pre_mean > thresholds[branch] and pre_low > 0.0
            postfilter = (
                post_low > 0.0
                and 2500.0 <= velocity <= 4500.0
                and abs(peak_velocity - velocity) <= VELOCITY_TOLERANCE_M_S
            )
            saved_point = saved["amplitudes"][point_index]
            source_point = source_report["amplitudes"][str(float(amplitude))]
            checks["saved_prefilter_metrics_reproduced"] &= (
                close(pre_mean, saved_point["prefilter_mean_log_enrichment"])
                and close(
                    pre_mean,
                    source_point["prefilter"][branch]["mean_log_enrichment"],
                )
            )
            checks["saved_postfilter_metrics_reproduced"] &= (
                close(exact_score, saved_point["postfilter_score_at_injected_velocity"])
                and close(post_uplift, saved_point["postfilter_score_uplift"])
            )
            checks["saved_peak_velocities_reproduced"] &= (
                close(peak_velocity, saved_point["postfilter_peak_velocity_m_s"])
                and close(
                    peak_velocity,
                    source_point["postfilter"][branch]["peak_velocity_m_s"],
                )
            )
            checks["independent_bootstrap_classifications_stable"] &= (
                independent == bool(saved_point["independent_recovery"])
                and postfilter == bool(saved_point["postfilter_recovery"])
            )
            if amplitude == 0.0:
                checks["zero_injection_not_recovered"] &= not independent and not postfilter
            audited_points.append({
                "amplitude_ratio": amplitude,
                "prefilter_mean_log_enrichment": pre_mean,
                "prefilter_uplift": pre_uplift,
                "prefilter_uplift_ci95_independent_seed": [pre_low, pre_high],
                "postfilter_score_at_injected_velocity": exact_score,
                "postfilter_score_uplift": post_uplift,
                "postfilter_score_uplift_ci95_independent_seed": [post_low, post_high],
                "postfilter_peak_velocity_m_s": peak_velocity,
                "independent_recovery": independent,
                "postfilter_recovery": postfilter,
            })
        positive = [point for point in audited_points if point["amplitude_ratio"] > 0]
        independent_levels = [point["amplitude_ratio"] for point in positive if point["independent_recovery"]]
        post_levels = [point["amplitude_ratio"] for point in positive if point["postfilter_recovery"]]
        independent_minimum = min(independent_levels) if independent_levels else None
        post_minimum = min(post_levels) if post_levels else None
        checks["reported_minima_reproduced"] &= (
            independent_minimum == saved["minimum_independent_recovery_ratio"]
            and post_minimum == saved["minimum_postfilter_recovery_ratio"]
        )
        audited_scenarios.append({
            "velocity_m_s": velocity,
            "direction": direction,
            "injected_branch": branch,
            "minimum_independent_recovery_ratio": independent_minimum,
            "minimum_postfilter_recovery_ratio": post_minimum,
            "amplitudes": audited_points,
        })

    report = {
        "workflow_version": "ambient_fk_injection_recovery_v1_completion_audit",
        "source_workflow_version": aggregate.get("workflow_version"),
        "expected_files": args.expected_files,
        "fixed_prefilter_null95": thresholds,
        "independent_bootstrap_realizations": args.bootstraps,
        "independent_bootstrap_seed": args.seed,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "scenarios": audited_scenarios,
        "interpretive_boundary": (
            "The audit validates saved computations, direction-specific trajectories, "
            "and threshold stability. It does not establish Green's-function convergence "
            "or the physical identity of the uninjected filtered ridge."
        ),
    }
    output = args.input_dir / "ambient_fk_injection_recovery_v1_completion_audit.json"
    output.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if not report["all_checks_pass"]:
        raise RuntimeError(f"completion audit failed: {checks}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--null-json", type=Path, default=DEFAULT_NULL)
    parser.add_argument("--date", default="2024-12-20")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=300)
    parser.add_argument("--expected-files", type=int, default=300)
    parser.add_argument("--bootstraps", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
