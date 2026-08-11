#!/usr/bin/env python3
"""Pre-F-K surrogate nulls for the corrected ambient signed-lag workflow.

This script asks a narrower and more rigorous question than the existing
receiver-distance permutation: does the observed velocity-score statistic
require a physically ordered, phase-coherent wavefield *before* the selected
F-K wedge is applied?

For every input segment, the observed and surrogate wavefields receive the
same 5--20 Hz preprocessing, signed F-K filter, channel-0 virtual-source
correlation, file stack, and apparent-velocity scan.  The surrogate transform
is applied after channel-wise preprocessing and before ``fk_filter``.  Channel
permutation is equivalent to permuting the raw channels before preprocessing
because the production preprocessing operates independently on each channel.

The two supported surrogate families are deliberately reported separately:

``channel_permutation``
    Keep the channel-0 trace at coordinate zero, but randomly reassign every
    other trace to a fiber coordinate.  This preserves individual trace
    spectra and amplitudes while destroying ordered moveout along the fiber.

``circular_time_shift``
    Keep channel 0 fixed and independently circularly shift every other trace.
    This preserves each trace's samples and spectrum while destroying its
    phase relationship to channel 0 and neighboring channels.

The familywise statistic is the largest absolute median correlation over both
signed branches and all trial velocities inside the predeclared 2.5--4.5 km/s
F-K wedge.  It therefore accounts for selecting the larger branch and velocity
within that frozen search region.  Results do not establish Green's-function
recovery, wave type, or formation Vp, and are conditional on the wedge having
been specified independently of the evaluated data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ambient_fk_transfer_test import fk_filter
from ambient_signed_fk_v2 import BRANCH_LAG_SIGN, VELOCITIES, velocity_scores
from ambient_transfer_test import (
    CSV,
    corrected_path,
    load_segment,
    normalized_corr_pairs,
    preprocess,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "ambient_transfer" / "fk_full_pipeline_null_v2"
SUPPORTED_NULLS = ("channel_permutation", "circular_time_shift")
VELOCITY_WEDGE_M_S = (2500.0, 4500.0)
PASSBAND_HZ = (5.0, 20.0)
MAX_LAG_S = 0.35


def stable_rng(seed: int, method: str, realization: int, file_key: str) -> np.random.Generator:
    """Return a reproducible generator independent of Python's hash seed."""
    token = f"{seed}|{method}|{realization}|{file_key}".encode("utf-8")
    digest = hashlib.sha256(token).digest()
    child_seed = int.from_bytes(digest[:8], byteorder="little", signed=False)
    return np.random.default_rng(child_seed)


def channel_permutation(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Destroy spatial ordering while keeping the virtual-source trace fixed."""
    if x.shape[0] < 3:
        raise ValueError("channel permutation requires at least three channels")
    order = np.arange(x.shape[0])
    order[1:] = rng.permutation(order[1:])
    return x[order]


def circular_time_shift(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Destroy interchannel phase coherence while keeping channel 0 fixed."""
    if x.shape[1] < 2:
        raise ValueError("circular time shifts require at least two samples")
    out = np.empty_like(x)
    out[0] = x[0]
    shifts = rng.integers(0, x.shape[1], size=x.shape[0] - 1)
    for channel, shift in enumerate(shifts, start=1):
        out[channel] = np.roll(x[channel], int(shift))
    return out


def make_surrogate(
    x: np.ndarray,
    method: str,
    rng: np.random.Generator,
) -> np.ndarray:
    if method == "channel_permutation":
        return channel_permutation(x, rng)
    if method == "circular_time_shift":
        return circular_time_shift(x, rng)
    raise ValueError(f"unsupported null method: {method}")


def branch_correlations(
    x: np.ndarray,
    fs: float,
    dx: float,
    mode: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply the production F-K filter and channel-0 correlation geometry."""
    filtered, filtered_fs, filtered_dx = fk_filter(x, fs, dx, mode)
    targets = [
        int(round(50.0 * index / filtered_dx))
        for index in range(1, 15)
    ]
    targets = [target for target in targets if target < filtered.shape[0]]
    if not targets:
        raise RuntimeError("no receiver targets fit inside the filtered array")
    lags, correlation = normalized_corr_pairs(
        filtered,
        [(0, target) for target in targets],
        filtered_fs,
        max_lag=MAX_LAG_S,
    )
    distance = np.asarray(targets, dtype=float) * filtered_dx
    return lags, distance, correlation


def score_stack(top: np.ndarray, lags: np.ndarray, distance: np.ndarray, mode: str) -> np.ndarray:
    return velocity_scores(top, lags, distance, BRANCH_LAG_SIGN[mode])


def selection_statistic(scores: np.ndarray) -> tuple[float, float]:
    """Maximum absolute score and its velocity inside the frozen F-K wedge."""
    inside = (
        (VELOCITIES >= VELOCITY_WEDGE_M_S[0])
        & (VELOCITIES <= VELOCITY_WEDGE_M_S[1])
    )
    local = np.flatnonzero(inside)
    peak = int(local[np.nanargmax(np.abs(scores[inside]))])
    return float(abs(scores[peak])), float(VELOCITIES[peak])


def empirical_p(observed: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= observed)) / (null.size + 1))


def load_rows(date: str, start: int, nfiles: int) -> pd.DataFrame:
    database = pd.read_csv(CSV, sep=r"\s+")
    database = database[database.nSamples > 0].copy()
    database["time"] = pd.to_datetime(database.startTime, utc=True, errors="coerce")
    rows = database[database.time.dt.strftime("%Y-%m-%d") == date].sort_values("time")
    return rows.iloc[start:start + nfiles]


def run(args: argparse.Namespace) -> tuple[Path, Path]:
    methods = tuple(item.strip() for item in args.null_methods.split(",") if item.strip())
    invalid = sorted(set(methods) - set(SUPPORTED_NULLS))
    if invalid:
        raise ValueError(f"unsupported null methods: {invalid}; choose from {SUPPORTED_NULLS}")
    if not methods:
        raise ValueError("at least one null method is required")
    if args.nulls < 1:
        raise ValueError("--nulls must be at least one")

    rows = load_rows(args.date, args.start, args.nfiles)
    observed_sums = {mode: None for mode in BRANCH_LAG_SIGN}
    null_sums = {
        method: {mode: None for mode in BRANCH_LAG_SIGN}
        for method in methods
    }
    used_files: list[str] = []
    lags = distance = None

    for file_index, row in enumerate(rows.itertuples(index=False)):
        path = corrected_path(row.file)
        if not path.exists():
            continue
        raw, fs, dx = load_segment(path)
        processed = preprocess(raw, fs, norm_seconds=5.0)

        for mode in BRANCH_LAG_SIGN:
            branch_lags, branch_distance, correlation = branch_correlations(
                processed, fs, dx, mode
            )
            if observed_sums[mode] is None:
                observed_sums[mode] = np.zeros_like(correlation)
            observed_sums[mode] += correlation
            if lags is None:
                lags = branch_lags
                distance = branch_distance

        file_key = str(path)
        for method in methods:
            for realization in range(args.null_start, args.null_start + args.nulls):
                rng = stable_rng(args.seed, method, realization, file_key)
                surrogate = make_surrogate(processed, method, rng)
                null_index = realization - args.null_start
                for mode in BRANCH_LAG_SIGN:
                    _, surrogate_distance, correlation = branch_correlations(
                        surrogate, fs, dx, mode
                    )
                    if not np.allclose(surrogate_distance, distance):
                        raise RuntimeError("receiver geometry changed inside null pipeline")
                    if null_sums[method][mode] is None:
                        null_sums[method][mode] = np.zeros(
                            (args.nulls,) + correlation.shape,
                            dtype=float,
                        )
                    null_sums[method][mode][null_index] += correlation

        used_files.append(file_key)
        print(
            f"processed {len(used_files)}/{len(rows)} "
            f"(nulls={args.nulls}, methods={len(methods)})",
            flush=True,
        )

    if not used_files:
        raise RuntimeError(
            f"no usable files for {args.date}, start={args.start}, nfiles={args.nfiles}"
        )

    count = len(used_files)
    observed_top = {mode: value / count for mode, value in observed_sums.items()}
    observed_scores = {
        mode: score_stack(observed_top[mode], lags, distance, mode)
        for mode in BRANCH_LAG_SIGN
    }
    observed_stats = {
        mode: selection_statistic(observed_scores[mode])
        for mode in BRANCH_LAG_SIGN
    }
    observed_familywise = max(value[0] for value in observed_stats.values())
    velocity_3200_index = int(np.argmin(abs(VELOCITIES - 3200.0)))

    report: dict[str, object] = {
        "workflow_version": "ambient_fk_full_pipeline_null_v2",
        "scientific_question": (
            "Does the signed-lag velocity statistic require spatially ordered, "
            "phase-coherent data before the frozen F-K projection?"
        ),
        "date": args.date,
        "start": args.start,
        "requested_files": args.nfiles,
        "used_files": count,
        "null_start": args.null_start,
        "null_realizations": args.nulls,
        "random_seed": args.seed,
        "null_methods": list(methods),
        "coordinate_assumption": (
            "channel 0 is the top virtual source; increasing channel coordinate "
            "follows the downhole fiber direction"
        ),
        "receiver_offsets_m": distance.tolist(),
        "preprocessing": "5-20 Hz bandpass followed by 5-s running-absolute-mean normalization",
        "surrogate_stage": "after channel-wise preprocessing and before F-K filtering",
        "fk_velocity_wedge_m_s": list(VELOCITY_WEDGE_M_S),
        "fk_signed_branches": {
            "negative": "F*K < 0 evaluated at positive lag",
            "positive": "F*K > 0 evaluated at negative lag",
        },
        "selection_statistic": (
            "maximum absolute median normalized correlation over 2.5-4.5 km/s; "
            "familywise maximum also spans both signed branches"
        ),
        "observed": {
            mode: {
                "peak_absolute_score_in_wedge": observed_stats[mode][0],
                "peak_velocity_m_s_in_wedge": observed_stats[mode][1],
                "signed_score_3200": float(observed_scores[mode][velocity_3200_index]),
                "absolute_score_3200": float(abs(observed_scores[mode][velocity_3200_index])),
            }
            for mode in BRANCH_LAG_SIGN
        },
        "observed_familywise_score": observed_familywise,
        "null_results": {},
        "supported_claim": (
            "A small empirical p value supports that the filtered/correlated ridge "
            "depends on pre-F-K spatial/phase coherence and is not produced by the "
            "fixed wedge acting on these coherence-destroyed surrogates."
        ),
        "limitations": [
            "The inference is conditional on the 2.5-4.5 km/s wedge being frozen before evaluation.",
            "These artificial surrogates do not reproduce every form of spatially correlated ambient noise.",
            "The test does not establish Green's-function convergence, wave type, causality, or formation Vp.",
            "Circular shifts introduce one wrap boundary per shifted trace.",
            "The empirical p-value resolution is limited to 1/(number of nulls + 1).",
        ],
    }

    arrays: dict[str, np.ndarray] = {
        "null_realization_ids": np.arange(
            args.null_start, args.null_start + args.nulls, dtype=int
        ),
        "lags_s": lags,
        "receiver_offsets_m": distance,
        "velocities_m_s": VELOCITIES,
        "used_files": np.asarray(used_files),
    }
    for mode in BRANCH_LAG_SIGN:
        arrays[f"observed_{mode}_top"] = observed_top[mode]
        arrays[f"observed_{mode}_scores"] = observed_scores[mode]

    for method in methods:
        method_report: dict[str, object] = {"branches": {}}
        familywise = np.zeros(args.nulls, dtype=float)
        for mode in BRANCH_LAG_SIGN:
            tops = null_sums[method][mode] / count
            score_curves = np.asarray([
                score_stack(top, lags, distance, mode)
                for top in tops
            ])
            statistics = np.asarray([
                selection_statistic(curve)[0]
                for curve in score_curves
            ])
            familywise = np.maximum(familywise, statistics)
            fixed_3200 = np.abs(score_curves[:, velocity_3200_index])
            observed_value, observed_velocity = observed_stats[mode]
            method_report["branches"][mode] = {
                "observed_peak_absolute_score_in_wedge": observed_value,
                "observed_peak_velocity_m_s_in_wedge": observed_velocity,
                "null95_peak_score_in_wedge": float(np.quantile(statistics, 0.95)),
                "p_peak_in_wedge": empirical_p(observed_value, statistics),
                "null95_absolute_score_3200": float(np.quantile(fixed_3200, 0.95)),
                "p_absolute_score_3200": empirical_p(
                    float(abs(observed_scores[mode][velocity_3200_index])),
                    fixed_3200,
                ),
            }
            arrays[f"null_{method}_{mode}_scores"] = score_curves
            arrays[f"null_{method}_{mode}_selection_statistics"] = statistics
            arrays[f"null_{method}_{mode}_top_first"] = tops[0]
        method_report["familywise"] = {
            "observed_maximum_across_branches_and_wedge": observed_familywise,
            "null95": float(np.quantile(familywise, 0.95)),
            "p": empirical_p(observed_familywise, familywise),
        }
        report["null_results"][method] = method_report
        arrays[f"null_{method}_familywise_statistics"] = familywise

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.stem or (
        f"fk_full_pipeline_null_v2_{args.date}_start{args.start}_"
        f"n{args.nfiles}_r{args.null_start}-{args.null_start + args.nulls - 1}"
    )
    json_path = args.output_dir / f"{stem}.json"
    npz_path = args.output_dir / f"{stem}.npz"
    if not args.overwrite and (json_path.exists() or npz_path.exists()):
        raise FileExistsError(
            f"output exists for {stem}; choose a new --stem or pass --overwrite"
        )
    np.savez_compressed(npz_path, **arrays)
    json_path.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return json_path, npz_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run coherence-destroying surrogates before the signed F-K pipeline."
    )
    parser.add_argument("--date", required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=1)
    parser.add_argument("--null-start", type=int, default=0)
    parser.add_argument("--nulls", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260805)
    parser.add_argument(
        "--null-methods",
        default=",".join(SUPPORTED_NULLS),
        help="comma-separated subset of channel_permutation,circular_time_shift",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--stem")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
