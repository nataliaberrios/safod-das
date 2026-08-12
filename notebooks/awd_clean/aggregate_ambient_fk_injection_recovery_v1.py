#!/usr/bin/env python3
"""Aggregate real-noise ambient F-K injection-recovery scenarios."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ambient_fk_injection_recovery_v1 import exact_velocity_score
from ambient_signed_fk_v2 import BRANCH_LAG_SIGN


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "ambient_transfer" / "fk_injection_recovery_v1_n300"
DEFAULT_NULL = (
    ROOT
    / "ambient_transfer"
    / "fk_prefilter_energy_v1_n300_r20"
    / "fk_prefilter_energy_v1_aggregate.json"
)
VELOCITIES = (1800.0, 2750.0, 3200.0, 4000.0)
DIRECTIONS = (1, -1)
VELOCITY_TOLERANCE_M_S = 250.0


def amplitude_tag(amplitude: float) -> str:
    return f"a{amplitude:.6g}".replace(".", "p")


def physical_trajectory_values(
    correlations: np.ndarray,
    lags: np.ndarray,
    distance: np.ndarray,
    velocity_m_s: float,
    lag_sign: float,
) -> np.ndarray:
    """Sample per-file, per-receiver correlations on one fixed trajectory."""
    indices = np.asarray([
        np.argmin(np.abs(lags - lag_sign * offset / velocity_m_s))
        for offset in distance
    ])
    return np.stack([
        correlations[:, receiver, index]
        for receiver, index in enumerate(indices)
    ], axis=1)


def paired_bootstrap_mean(
    injected: np.ndarray,
    baseline: np.ndarray,
    rng: np.random.Generator,
    nboot: int,
) -> tuple[float, float, float]:
    """Return estimate and 95% paired bootstrap interval for mean uplift."""
    difference = np.asarray(injected) - np.asarray(baseline)
    estimate = float(np.mean(difference))
    indices = rng.integers(0, difference.size, size=(nboot, difference.size))
    draws = np.mean(difference[indices], axis=1)
    low, high = np.quantile(draws, (0.025, 0.975))
    return estimate, float(low), float(high)


def paired_bootstrap_stack_score(
    injected: np.ndarray,
    baseline: np.ndarray,
    rng: np.random.Generator,
    nboot: int,
) -> tuple[float, float, float]:
    """Bootstrap the exact difference of median-after-stacking scores.

    The production score first averages files at every receiver and then takes
    the median across receivers. Subtraction must occur after both scores are
    formed; a median of receiver-wise differences is not generally equivalent.
    Paired resamples use the same file indices for injected and baseline data.
    """
    injected = np.asarray(injected)
    baseline = np.asarray(baseline)
    if injected.shape != baseline.shape or injected.ndim != 2:
        raise ValueError("injected and baseline must be equal file-by-receiver arrays")
    estimate = float(
        np.median(np.mean(injected, axis=0))
        - np.median(np.mean(baseline, axis=0))
    )
    indices = rng.integers(0, injected.shape[0], size=(nboot, injected.shape[0]))
    injected_draws = np.median(np.mean(injected[indices], axis=1), axis=1)
    baseline_draws = np.median(np.mean(baseline[indices], axis=1), axis=1)
    draws = injected_draws - baseline_draws
    low, high = np.quantile(draws, (0.025, 0.975))
    return estimate, float(low), float(high)


def scenario_paths(input_dir: Path, date: str, start: int, nfiles: int):
    for velocity in VELOCITIES:
        for direction in DIRECTIONS:
            direction_tag = "inc" if direction == 1 else "dec"
            stem = (
                f"fk_injection_v1_{date}_start{start}_n{nfiles}_"
                f"v{int(velocity)}_{direction_tag}"
            )
            yield velocity, direction, input_dir / f"{stem}.json", input_dir / f"{stem}.npz"


def null_thresholds(path: Path) -> dict[str, float]:
    """Return the stricter branch threshold across the two independent nulls."""
    report = json.loads(path.read_text())
    thresholds = {}
    for branch in ("negative", "positive"):
        thresholds[branch] = max(
            float(report["null_results"][method]["branches"][branch]["null95"])
            for method in ("channel_permutation", "circular_time_shift")
        )
    return thresholds


def plot_summary(summary: dict, output: Path) -> None:
    """Make the publication-style sensitivity overview."""
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.2), constrained_layout=True)
    colors = {1800.0: "0.45", 2750.0: "#0072B2", 3200.0: "#D55E00", 4000.0: "#009E73"}
    styles = {1: "-", -1: "--"}
    direction_labels = {1: "increasing coordinate", -1: "decreasing coordinate"}
    for scenario in summary["scenarios"]:
        velocity = float(scenario["velocity_m_s"])
        direction = int(scenario["direction"])
        points = scenario["amplitudes"]
        amplitude = np.asarray([point["amplitude_ratio"] for point in points])
        label = f"{velocity / 1000:g} km/s, {direction_labels[direction]}"
        color, style = colors[velocity], styles[direction]

        pre = np.asarray([point["prefilter_mean_log_enrichment"] for point in points])
        axes[0, 0].plot(amplitude, pre, style, color=color, marker="o", ms=3, label=label)

        pre_delta = np.asarray([point["prefilter_uplift"] for point in points])
        pre_low = np.asarray([point["prefilter_uplift_ci95"][0] for point in points])
        pre_high = np.asarray([point["prefilter_uplift_ci95"][1] for point in points])
        axes[0, 1].plot(amplitude, pre_delta, style, color=color, marker="o", ms=3)
        axes[0, 1].fill_between(amplitude, pre_low, pre_high, color=color, alpha=0.08)

        post_delta = np.asarray([point["postfilter_score_uplift"] for point in points])
        post_low = np.asarray([point["postfilter_score_uplift_ci95"][0] for point in points])
        post_high = np.asarray([point["postfilter_score_uplift_ci95"][1] for point in points])
        axes[1, 0].plot(amplitude, post_delta, style, color=color, marker="o", ms=3)
        axes[1, 0].fill_between(amplitude, post_low, post_high, color=color, alpha=0.08)

        recovered = np.asarray([point["postfilter_peak_velocity_m_s"] / 1000 for point in points])
        axes[1, 1].plot(amplitude, recovered, style, color=color, marker="o", ms=3)
        axes[1, 1].axhline(velocity / 1000, color=color, lw=0.5, alpha=0.2)

    positive_x = [
        point["amplitude_ratio"]
        for scenario in summary["scenarios"]
        for point in scenario["amplitudes"]
        if point["amplitude_ratio"] > 0
    ]
    xmin = min(positive_x) / 1.4
    for ax in axes.flat:
        ax.set_xscale("symlog", linthresh=xmin, linscale=0.8)
        ax.set_xlabel("Injected RMS / median real 5–20 Hz RMS")
        ax.grid(alpha=0.2)
    axes[0, 0].set_ylabel("Mean pre-filter log target/reference power")
    axes[0, 0].set_title("(a) Independent pre-filter response")
    axes[0, 0].legend(frameon=False, fontsize=7, ncol=2)
    axes[0, 1].axhline(0.0, color="0.2", lw=0.8)
    axes[0, 1].set_ylabel("Paired pre-filter uplift")
    axes[0, 1].set_title("(b) Change relative to identical uninjected files")
    axes[1, 0].axhline(0.0, color="0.2", lw=0.8)
    axes[1, 0].set_ylabel("Paired signed-correlation score uplift")
    axes[1, 0].set_title("(c) Production F–K plus correlation response")
    axes[1, 1].set_ylabel("Recovered conditional peak velocity (km/s)")
    axes[1, 1].set_title("(d) Conditional peak after production wedge")
    fig.suptitle(
        f"Ambient F–K injection–recovery in real noise ({summary['used_files']} files)",
        fontsize=14,
    )
    fig.savefig(output.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace) -> tuple[Path, Path]:
    expected = list(scenario_paths(args.input_dir, args.date, args.start, args.nfiles))
    missing = [str(path) for _, _, json_path, npz_path in expected for path in (json_path, npz_path) if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing scenario products:\n" + "\n".join(missing))
    thresholds = null_thresholds(args.null_json)
    rng = np.random.default_rng(args.seed)
    summary: dict[str, object] = {
        "workflow_version": "ambient_fk_injection_recovery_v1_aggregate",
        "date": args.date,
        "start": args.start,
        "requested_files": args.nfiles,
        "used_files": None,
        "amplitude_definition": "injected broadband RMS / median real-channel 5-20 Hz RMS",
        "prefilter_fixed_null95": thresholds,
        "bootstrap_realizations": args.bootstraps,
        "velocity_tolerance_m_s": VELOCITY_TOLERANCE_M_S,
        "decision_rule": (
            "Independent recovery requires the correct-branch mean pre-filter statistic "
            "to exceed the stricter original channel-permutation/time-shift 95% threshold "
            "and the paired-bootstrap pre-filter uplift lower bound to exceed zero. "
            "Post-filter recovery additionally requires positive paired score uplift and a "
            "conditional peak within 250 m/s of the injected in-wedge velocity."
        ),
        "scenarios": [],
    }
    reference_files = None
    reference_zero = None
    for velocity, direction, json_path, npz_path in expected:
        report = json.loads(json_path.read_text())
        product = np.load(npz_path)
        files = product["used_files"].astype(str)
        if reference_files is None:
            reference_files = files
            summary["used_files"] = int(files.size)
        elif not np.array_equal(reference_files, files):
            raise RuntimeError("scenario file lists differ")
        injected_branch = "negative" if direction == 1 else "positive"
        amplitudes = product["amplitude_ratios"].astype(float)
        baseline_tag = amplitude_tag(0.0)
        baseline_pre = product[f"{baseline_tag}_{injected_branch}_prefilter_per_file"]
        baseline_corr = product[f"{baseline_tag}_{injected_branch}_correlation_per_file"]
        baseline_trajectory = physical_trajectory_values(
            baseline_corr,
            product["lags_s"],
            product["distance_m"],
            velocity,
            BRANCH_LAG_SIGN[injected_branch],
        )
        current_zero = np.concatenate([
            product[f"{baseline_tag}_{branch}_prefilter_per_file"]
            for branch in ("negative", "positive")
        ])
        if reference_zero is None:
            reference_zero = current_zero
        elif not np.allclose(reference_zero, current_zero, rtol=5e-5, atol=1e-8):
            raise RuntimeError("zero-injection baselines differ across scenarios")
        scenario = {
            "velocity_m_s": velocity,
            "direction": direction,
            "injected_branch": injected_branch,
            "inside_production_wedge": 2500.0 <= velocity <= 4500.0,
            "amplitudes": [],
        }
        for amplitude in amplitudes:
            tag = amplitude_tag(amplitude)
            pre = product[f"{tag}_{injected_branch}_prefilter_per_file"]
            pre_uplift, pre_low, pre_high = paired_bootstrap_mean(
                pre, baseline_pre, rng, args.bootstraps
            )
            corr = product[f"{tag}_{injected_branch}_correlation_per_file"]
            trajectory = physical_trajectory_values(
                corr,
                product["lags_s"],
                product["distance_m"],
                velocity,
                BRANCH_LAG_SIGN[injected_branch],
            )
            post_uplift, post_low, post_high = paired_bootstrap_stack_score(
                trajectory, baseline_trajectory, rng, args.bootstraps
            )
            entry = report["amplitudes"][str(float(amplitude))]
            peak_velocity = float(entry["postfilter"][injected_branch]["peak_velocity_m_s"])
            independent = (
                float(np.mean(pre)) > thresholds[injected_branch]
                and pre_low > 0.0
            )
            postfilter = (
                post_low > 0.0
                and abs(peak_velocity - velocity) <= VELOCITY_TOLERANCE_M_S
                and scenario["inside_production_wedge"]
            )
            scenario["amplitudes"].append({
                "amplitude_ratio": float(amplitude),
                "prefilter_mean_log_enrichment": float(np.mean(pre)),
                "prefilter_uplift": pre_uplift,
                "prefilter_uplift_ci95": [pre_low, pre_high],
                "prefilter_passes_fixed_null95": bool(float(np.mean(pre)) > thresholds[injected_branch]),
                "independent_recovery": bool(independent),
                "postfilter_score_at_injected_velocity": exact_velocity_score(
                    np.mean(corr, axis=0),
                    product["lags_s"],
                    product["distance_m"],
                    velocity,
                    BRANCH_LAG_SIGN[injected_branch],
                ),
                "postfilter_score_uplift": post_uplift,
                "postfilter_score_uplift_ci95": [post_low, post_high],
                "postfilter_peak_velocity_m_s": peak_velocity,
                "postfilter_recovery": bool(postfilter),
            })
        positive = [item for item in scenario["amplitudes"] if item["amplitude_ratio"] > 0]
        independent_levels = [item["amplitude_ratio"] for item in positive if item["independent_recovery"]]
        post_levels = [item["amplitude_ratio"] for item in positive if item["postfilter_recovery"]]
        scenario["minimum_independent_recovery_ratio"] = min(independent_levels) if independent_levels else None
        scenario["minimum_postfilter_recovery_ratio"] = min(post_levels) if post_levels else None
        summary["scenarios"].append(scenario)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.output_dir / "ambient_fk_injection_recovery_v1_aggregate"
    json_path = stem.with_suffix(".json")
    json_path.write_text(json.dumps(summary, indent=2))
    plot_summary(summary, stem)
    print(json.dumps(summary, indent=2))
    return json_path, stem.with_suffix(".png")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--null-json", type=Path, default=DEFAULT_NULL)
    parser.add_argument("--date", default="2024-12-20")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=300)
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260811)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
