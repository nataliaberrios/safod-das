#!/usr/bin/env python3
"""Run the reported Lellouch et al. (2019) ambient stack as one auditable operator.

The baseline contains only operations reported in section 4.1 of the paper:
time differentiation, running-absolute-mean normalization, 30 s correlation
windows with 15 s overlap, the fixed-top and 50 m constant-offset geometries,
R-10:R+10 correlation sums, a first unshifted stack, travel-time-shifted
restacking, a final 5--20 Hz bandpass, and three-sample quadratic peak picks.

F-K filtering, common-mode removal, and stabilized source-power division are
available only as explicitly labelled sensitivity branches.  The paper gives
the Green-function proportionality containing |S(omega)|^2 but does not state
how that spectrum was estimated or stabilized, so no spectral-deconvolution
branch is called an exact reproduction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import hilbert

from ambient_lellouch2019_reproduction_v1 import (
    MAX_LAG_SECONDS,
    OUTPUT_BAND_HZ,
    RAM_SECONDS,
    STEP_SECONDS,
    WINDOW_SECONDS,
    acquisition_metadata,
    day_rows,
    final_bandpass,
    release_allocator_memory,
    strain_rate_and_ram,
)
from ambient_lellouch_fk_matched_v2 import MODE_SPECS, fk_filter
from ambient_transfer_test import corrected_path


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "ambient_transfer" / "lellouch2019_exact_stack"
TARGET_SPACING_M = 50.0
LOCAL_HALF_WIDTH_CHANNELS = 10
VELOCITY_GRID_M_S = np.arange(1500.0, 5000.1, 25.0)
NULL_METHODS = ("ordered", "white_noise", "channel_permutation", "circular_time_shift")


def stable_rng(seed: int, method: str, realization: int, key: str) -> np.random.Generator:
    token = f"{seed}|{method}|{realization}|{key}".encode()
    child = int.from_bytes(hashlib.sha256(token).digest()[:8], "little")
    return np.random.default_rng(child)


def build_geometries(dx: float, n_channels: int) -> dict[str, np.ndarray]:
    """Return pair arrays for the two geometries shown in paper Figure 7a-b."""
    step = int(round(TARGET_SPACING_M / dx))
    fixed_targets = np.arange(step, n_channels, step, dtype=int)
    fixed_targets = fixed_targets[
        fixed_targets + LOCAL_HALF_WIDTH_CHANNELS < n_channels
    ]
    fixed_targets = fixed_targets[:16]
    fixed_sources = np.zeros((fixed_targets.size, 21), dtype=int)
    fixed_receivers = np.stack([
        np.arange(target - 10, target + 11, dtype=int)
        for target in fixed_targets
    ])

    constant_sources_center = np.arange(0, n_channels, step, dtype=int)
    constant_targets = constant_sources_center + step
    usable = constant_targets + LOCAL_HALF_WIDTH_CHANNELS < n_channels
    constant_sources_center = constant_sources_center[usable][:15]
    constant_targets = constant_targets[usable][:15]
    constant_sources = np.repeat(constant_sources_center[:, None], 21, axis=1)
    constant_receivers = np.stack([
        np.arange(target - 10, target + 11, dtype=int)
        for target in constant_targets
    ])

    maximum = int(max(fixed_receivers.max(), constant_receivers.max()))
    return {
        "fixed_sources": fixed_sources,
        "fixed_receivers": fixed_receivers,
        "fixed_targets": fixed_targets,
        "fixed_coordinates_m": fixed_targets.astype(float) * dx,
        "constant_sources": constant_sources,
        "constant_receivers": constant_receivers,
        "constant_source_coordinates_m": constant_sources_center.astype(float) * dx,
        "constant_target_coordinates_m": constant_targets.astype(float) * dx,
        "constant_midpoint_coordinates_m": (
            (constant_sources_center + constant_targets).astype(float) * dx / 2.0
        ),
        "maximum_channel": np.asarray(maximum),
    }


def load_contiguous(path: Path, maximum_channel: int) -> tuple[np.ndarray, float]:
    with h5py.File(path, "r") as handle:
        raw_group = handle["Acquisition/Raw[0]"]
        data = raw_group["RawData"]
        fs = float(raw_group.attrs.get("OutputDataRate", data.attrs.get("OutputDataRate", 500.0)))
        return np.asarray(data[:, : maximum_channel + 1], dtype=np.float32).T, fs


def transform_input(
    normalized: np.ndarray,
    method: str,
    rng: np.random.Generator,
) -> np.ndarray:
    """Apply a null before any optional F-K operation."""
    if method == "ordered":
        return normalized
    if method == "white_noise":
        return normalized
    if method == "channel_permutation":
        return normalized[rng.permutation(normalized.shape[0])]
    if method == "circular_time_shift":
        out = np.empty_like(normalized)
        shifts = rng.integers(0, normalized.shape[1], size=normalized.shape[0])
        for channel, shift in enumerate(shifts):
            out[channel] = np.roll(normalized[channel], int(shift))
        return out
    raise ValueError(method)


def remove_common_mode(data: np.ndarray) -> np.ndarray:
    """Optional, non-published sensitivity operation."""
    return data - np.median(data, axis=0, keepdims=True)


def pair_correlations(
    window: np.ndarray,
    sources: np.ndarray,
    receivers: np.ndarray,
    fs: float,
    spectral_mode: str,
    waterlevel: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Correlate all flattened pairs and restore geometry x neighbor axes."""
    n_time = window.shape[1]
    n_fft = 1 << int(np.ceil(np.log2(2 * n_time - 1)))
    maximum_lag = int(round(MAX_LAG_SECONDS * fs))
    spectra = np.fft.rfft(window, n=n_fft, axis=1)
    source_spectra = spectra[sources.ravel()]
    receiver_spectra = spectra[receivers.ravel()]
    cross = np.conj(source_spectra) * receiver_spectra
    if spectral_mode == "source_power_stabilized":
        power = np.abs(source_spectra) ** 2
        floor = waterlevel * np.max(power, axis=1, keepdims=True)
        cross /= np.maximum(power, np.maximum(floor, np.finfo(float).eps))
    elif spectral_mode != "cross_correlation":
        raise ValueError(spectral_mode)
    correlation = np.fft.irfft(cross, n=n_fft, axis=1)
    correlation = np.concatenate(
        (correlation[:, -maximum_lag:], correlation[:, : maximum_lag + 1]),
        axis=1,
    )
    correlation /= float(n_time)
    correlation = correlation.reshape(sources.shape + (correlation.shape[-1],))
    lags = np.arange(-maximum_lag, maximum_lag + 1, dtype=float) / fs
    return lags, correlation


def add_windows(
    data: np.ndarray,
    fs: float,
    geometries: dict[str, np.ndarray],
    accumulators: dict[str, np.ndarray | None],
    spectral_mode: str,
    waterlevel: float,
) -> tuple[dict[str, np.ndarray], np.ndarray, int]:
    window_samples = int(round(WINDOW_SECONDS * fs))
    step_samples = int(round(STEP_SECONDS * fs))
    used = 0
    lags = None
    for start in range(0, data.shape[1] - window_samples + 1, step_samples):
        window = data[:, start : start + window_samples]
        for name in ("fixed", "constant"):
            lags, current = pair_correlations(
                window,
                geometries[f"{name}_sources"],
                geometries[f"{name}_receivers"],
                fs,
                spectral_mode,
                waterlevel,
            )
            if accumulators[name] is None:
                accumulators[name] = np.zeros_like(current)
            accumulators[name] += current
        used += 1
    return accumulators, lags, used


def sum_neighbors(correlations: np.ndarray) -> np.ndarray:
    """Implement C_S,R = sum from Z=R-10 through R+10."""
    return np.sum(correlations, axis=1)


def shifted_neighbor_sum(
    correlations: np.ndarray,
    lags: np.ndarray,
    receivers: np.ndarray,
    targets: np.ndarray,
    dx: float,
    velocity_m_s: float,
) -> np.ndarray:
    output = np.zeros((correlations.shape[0], correlations.shape[-1]), dtype=float)
    for row in range(correlations.shape[0]):
        for neighbor, receiver in enumerate(receivers[row]):
            delta = (float(receiver) - float(targets[row])) * dx / velocity_m_s
            output[row] += np.interp(
                lags + delta,
                lags,
                correlations[row, neighbor],
                left=0.0,
                right=0.0,
            )
    return output


def fixed_top_velocity_scan(
    section: np.ndarray,
    lags: np.ndarray,
    distances: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Descriptive envelope scan used to expose, not hide, first-pass ambiguity."""
    envelope = np.abs(hilbert(section, axis=1))
    scale = np.max(envelope, axis=1, keepdims=True)
    envelope = envelope / np.maximum(scale, np.finfo(float).eps)
    scores = []
    for velocity in VELOCITY_GRID_M_S:
        values = [
            trace[int(np.argmin(np.abs(lags - distance / velocity)))]
            for trace, distance in zip(envelope, distances)
        ]
        scores.append(float(np.median(values)))
    scores = np.asarray(scores)
    best = float(VELOCITY_GRID_M_S[int(np.argmax(scores))])
    return VELOCITY_GRID_M_S.copy(), scores, best


def quadratic_peak(
    trace: np.ndarray,
    lags: np.ndarray,
    lag_min: float = 0.005,
    lag_max: float = 0.060,
) -> float:
    """Three-sample quadratic interpolation around the largest causal sample."""
    eligible = np.flatnonzero((lags >= lag_min) & (lags <= lag_max))
    index = int(eligible[np.argmax(trace[eligible])])
    if index == 0 or index == trace.size - 1:
        return float(lags[index])
    left, center, right = trace[index - 1 : index + 2]
    denominator = left - 2.0 * center + right
    fraction = 0.0 if abs(denominator) < np.finfo(float).eps else 0.5 * (left - right) / denominator
    fraction = float(np.clip(fraction, -1.0, 1.0))
    return float(lags[index] + fraction * (lags[1] - lags[0]))


def plot_outputs(
    output_path: Path,
    lags: np.ndarray,
    fixed_simple: np.ndarray,
    fixed_distances: np.ndarray,
    constant_aligned: np.ndarray,
    midpoints: np.ndarray,
    picks: np.ndarray,
    alignment_velocity: float,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 6.5), constrained_layout=True)
    for trace, distance in zip(fixed_simple, fixed_distances):
        scale = max(np.max(np.abs(trace)), np.finfo(float).eps)
        axes[0].plot(lags, distance + 15.0 * trace / scale, color="black", lw=0.55)
    axes[0].plot(fixed_distances / alignment_velocity, fixed_distances, "r--", lw=1.4)
    axes[0].set(
        xlim=(0.0, 0.35),
        ylim=(fixed_distances[-1] + 30.0, fixed_distances[0] - 30.0),
        xlabel="Correlation lag (s)",
        ylabel="Receiver coordinate relative to channel 0 (m)",
        title="a  Fixed-top simple R±10 sum",
    )
    for trace, midpoint, pick in zip(constant_aligned, midpoints, picks):
        scale = max(np.max(np.abs(trace)), np.finfo(float).eps)
        axes[1].plot(lags, midpoint + 15.0 * trace / scale, color="black", lw=0.55)
        axes[1].plot(pick, midpoint, "o", ms=3.5, color="tab:red")
    axes[1].set(
        xlim=(0.0, 0.08),
        ylim=(midpoints[-1] + 30.0, midpoints[0] - 30.0),
        xlabel="Correlation lag (s)",
        ylabel="Pair midpoint coordinate (m)",
        title="b  Constant 50 m offset, aligned and picked",
    )
    fig.savefig(output_path.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    if args.null_method not in NULL_METHODS:
        raise ValueError(f"--null-method must be one of {NULL_METHODS}")
    rows = day_rows(args.date, args.start, args.nfiles)
    first = next(
        (corrected_path(row.file) for row in rows.itertuples(index=False)
         if corrected_path(row.file).is_file()),
        None,
    )
    if first is None:
        raise RuntimeError("no mounted input files")
    fs, dx, n_channels, duration = acquisition_metadata(first)
    if not np.isclose(duration, 60.0, atol=0.1):
        raise ValueError(f"expected 60 s records, found {duration}")
    geometries = build_geometries(dx, min(n_channels, 850))
    maximum_channel = int(geometries["maximum_channel"])
    accumulators: dict[str, np.ndarray | None] = {"fixed": None, "constant": None}
    used_files: list[str] = []
    used_windows = 0
    lags = None
    previous_tail = None
    previous_end = None
    step_samples = int(round(STEP_SECONDS * fs))

    for number, row in enumerate(rows.itertuples(index=False), 1):
        path = corrected_path(row.file)
        if not path.is_file():
            continue
        raw, current_fs = load_contiguous(path, maximum_channel)
        if not np.isclose(current_fs, fs):
            raise ValueError(f"sample-rate change at {path}")
        rng = stable_rng(args.seed, args.null_method, args.realization, str(path))
        if args.null_method == "white_noise":
            raw_input = rng.standard_normal(raw.shape).astype(np.float32)
            normalized = strain_rate_and_ram(raw_input, fs)
            del raw_input
        else:
            normalized = strain_rate_and_ram(raw, fs)
        wavefield = transform_input(normalized, args.null_method, rng)
        if args.common_mode:
            wavefield = remove_common_mode(wavefield)
        if args.fk_mode != "unfiltered":
            wavefield = fk_filter(wavefield, fs, dx, MODE_SPECS[args.fk_mode])
        row_time = pd.to_datetime(row.startTime, utc=True)
        contiguous = (
            previous_tail is not None
            and previous_end is not None
            and abs((row_time - previous_end).total_seconds()) < 0.02
        )
        if contiguous:
            bridge = np.concatenate(
                (previous_tail, wavefield[:, :step_samples]), axis=1
            )
            accumulators, lags, count = add_windows(
                bridge,
                fs,
                geometries,
                accumulators,
                args.spectral_mode,
                args.waterlevel,
            )
            used_windows += count
        accumulators, lags, count = add_windows(
            wavefield,
            fs,
            geometries,
            accumulators,
            args.spectral_mode,
            args.waterlevel,
        )
        used_windows += count
        previous_tail = wavefield[:, -step_samples:].copy()
        previous_end = row_time + pd.to_timedelta(
            wavefield.shape[1] / fs, unit="s"
        )
        used_files.append(str(path))
        del raw, normalized, wavefield
        release_allocator_memory()
        if number % 10 == 0:
            print(f"processed {number}/{len(rows)}; windows={used_windows}", flush=True)

    if not used_files or lags is None:
        raise RuntimeError("no usable records")
    averaged = {key: value / used_windows for key, value in accumulators.items()}
    fixed_simple = final_bandpass(sum_neighbors(averaged["fixed"]), fs)
    velocities, velocity_scores, data_velocity = fixed_top_velocity_scan(
        fixed_simple,
        lags,
        geometries["fixed_coordinates_m"],
    )
    alignment_velocity = (
        data_velocity if args.alignment_velocity == "data"
        else float(args.alignment_velocity)
    )
    fixed_aligned = final_bandpass(
        shifted_neighbor_sum(
            averaged["fixed"], lags, geometries["fixed_receivers"],
            geometries["fixed_targets"], dx, alignment_velocity
        ),
        fs,
    )
    constant_simple = final_bandpass(sum_neighbors(averaged["constant"]), fs)
    constant_targets = geometries["constant_receivers"][:, 10]
    constant_aligned = final_bandpass(
        shifted_neighbor_sum(
            averaged["constant"], lags, geometries["constant_receivers"],
            constant_targets, dx, alignment_velocity
        ),
        fs,
    )
    picks = np.asarray([quadratic_peak(trace, lags) for trace in constant_aligned])
    separation = (
        geometries["constant_target_coordinates_m"]
        - geometries["constant_source_coordinates_m"]
    )
    local_velocity = separation / picks

    args.output_dir.mkdir(parents=True, exist_ok=True)
    label = (
        f"lellouch_exact_{args.date}_start{args.start}_n{len(used_files)}_"
        f"{args.null_method}_r{args.realization}_{args.fk_mode}_"
        f"{args.spectral_mode}_align{alignment_velocity:.0f}"
    )
    stem = args.output_dir / label
    np.savez_compressed(
        stem.with_suffix(".npz"),
        lags_s=lags,
        fixed_coordinates_m=geometries["fixed_coordinates_m"],
        fixed_simple=fixed_simple,
        fixed_aligned=fixed_aligned,
        constant_midpoint_coordinates_m=geometries["constant_midpoint_coordinates_m"],
        constant_simple=constant_simple,
        constant_aligned=constant_aligned,
        constant_picks_s=picks,
        local_velocity_m_s=local_velocity,
        velocity_grid_m_s=velocities,
        fixed_simple_velocity_scores=velocity_scores,
        used_files=np.asarray(used_files),
    )
    report = {
        "workflow_version": "ambient_lellouch2019_exact_stack_v1",
        "citation": "Lellouch et al. (2019), doi:10.1029/2019JB017533",
        "date": args.date,
        "used_files": len(used_files),
        "used_windows": used_windows,
        "fs_hz": fs,
        "dx_m": dx,
        "null_method": args.null_method,
        "null_realization": args.realization,
        "fk_mode": args.fk_mode,
        "common_mode_removal": args.common_mode,
        "spectral_mode": args.spectral_mode,
        "source_power_waterlevel": args.waterlevel,
        "first_pass_data_velocity_m_s": data_velocity,
        "alignment_velocity_m_s": alignment_velocity,
        "alignment_velocity_argument": args.alignment_velocity,
        "mean_local_velocity_m_s": float(np.nanmean(local_velocity)),
        "reported_steps": {
            "time_derivative": True,
            "running_absolute_mean": True,
            "ram_seconds_project_choice": RAM_SECONDS,
            "window_seconds": WINDOW_SECONDS,
            "overlap_seconds": WINDOW_SECONDS - STEP_SECONDS,
            "fixed_top_geometry": True,
            "constant_offset_m": TARGET_SPACING_M,
            "neighbor_sum": "R-10:R+10, 21 correlations",
            "two_pass": "simple fixed-top velocity scan, then shifted neighbor restack",
            "output_band_hz": list(OUTPUT_BAND_HZ),
            "quadratic_pick_samples": 3,
        },
        "departures": [
            "The archive is 2024-2025, not the 2017 recording.",
            "The paper does not report the RAM duration; 5 s is a project choice.",
            "The paper does not specify the algorithm used to read average velocity from Figure 7c; the data option uses a declared envelope scan.",
            "Common-mode removal is not reported in section 4.1 and is off by default.",
            "F-K filtering is not reported for ambient interferometry and is off by default.",
            "Stabilized source-power division is a sensitivity branch because the estimator and stabilization of |S(omega)|^2 are not reported.",
        ],
    }
    stem.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    plot_outputs(
        stem,
        lags,
        fixed_simple,
        geometries["fixed_coordinates_m"],
        constant_aligned,
        geometries["constant_midpoint_coordinates_m"],
        picks,
        alignment_velocity,
    )
    print(json.dumps(report, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", default="2024-12-20")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=10)
    parser.add_argument("--null-method", default="ordered", choices=NULL_METHODS)
    parser.add_argument("--realization", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--fk-mode", default="unfiltered", choices=tuple(MODE_SPECS))
    parser.add_argument(
        "--spectral-mode",
        default="cross_correlation",
        choices=("cross_correlation", "source_power_stabilized"),
    )
    parser.add_argument("--waterlevel", type=float, default=1e-3)
    parser.add_argument("--common-mode", action="store_true")
    parser.add_argument("--alignment-velocity", default="data")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
