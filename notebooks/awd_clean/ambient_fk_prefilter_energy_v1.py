#!/usr/bin/env python3
"""Test for ordered target-slowness energy before ambient F-K filtering.

The production ambient workflow uses a signed 2.5--4.5 km/s F-K wedge to
isolate weak directional energy.  Scoring a correlation after applying that
same wedge is useful for describing the selected output, but is not an
independent test that the input wavefield contained target-slowness energy.

This module therefore measures spectral concentration directly in the
pre-filter F-K power spectrum.  No F-K mask is applied before the statistic is
computed.  Channel permutations and independent circular time shifts are
performed before the spatial FFT and compared with the identical statistic.

The test establishes only whether the normalized input wavefield contains
more spatially ordered energy in the predeclared target slowness interval than
the stated coherence-destroying surrogates.  It does not establish ambient
Green's-function convergence, wave type, causality, or formation velocity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal.windows import tukey

from ambient_transfer_test import CSV, corrected_path, load_segment, preprocess


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "ambient_transfer" / "fk_prefilter_energy_v1"
PASSBAND_HZ = (5.0, 20.0)
TARGET_VELOCITY_M_S = (2500.0, 4500.0)
REFERENCE_VELOCITY_M_S = (1200.0, 6500.0)
BRANCHES = ("negative", "positive")
NULL_METHODS = ("channel_permutation", "circular_time_shift")
SPACE_DECIMATION = 2
TIME_DECIMATION = 2


def stable_rng(seed: int, method: str, realization: int, file_key: str) -> np.random.Generator:
    """Return a deterministic generator independent of Python hash state."""
    token = f"{seed}|{method}|{realization}|{file_key}".encode("utf-8")
    digest = hashlib.sha256(token).digest()
    return np.random.default_rng(
        int.from_bytes(digest[:8], byteorder="little", signed=False)
    )


def time_spectrum(x: np.ndarray, fs: float) -> tuple[np.ndarray, np.ndarray]:
    """Return tapered channel spectra in the passband before any spatial FFT."""
    if x.ndim != 2:
        raise ValueError("expected channel-by-time input")
    time_window = tukey(x.shape[1], alpha=0.1).astype(x.dtype, copy=False)
    frequency = np.fft.rfftfreq(x.shape[1], 1.0 / fs)
    selected = (
        (frequency >= PASSBAND_HZ[0])
        & (frequency <= PASSBAND_HZ[1])
    )
    return frequency[selected], np.fft.rfft(x * time_window, axis=1)[:, selected]


def spatial_power(
    channel_spectrum: np.ndarray,
    dx: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return wavenumber coordinates and tapered pre-filter F-K power."""
    space_window = tukey(channel_spectrum.shape[0], alpha=0.1)
    spectrum = np.fft.fft(
        channel_spectrum * space_window[:, None],
        axis=0,
    )
    wavenumber = np.fft.fftfreq(channel_spectrum.shape[0], dx)
    return wavenumber, np.abs(spectrum) ** 2


def supports(
    frequency: np.ndarray,
    wavenumber: np.ndarray,
    branch: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return target and broad-reference supports without filtering the data."""
    k_grid, f_grid = np.meshgrid(wavenumber, frequency, indexing="ij")
    slowness = np.abs(k_grid) / np.maximum(f_grid, np.finfo(float).tiny)
    target = (
        (slowness >= 1.0 / TARGET_VELOCITY_M_S[1])
        & (slowness <= 1.0 / TARGET_VELOCITY_M_S[0])
    )
    reference = (
        (slowness >= 1.0 / REFERENCE_VELOCITY_M_S[1])
        & (slowness <= 1.0 / REFERENCE_VELOCITY_M_S[0])
    )
    if branch == "negative":
        direction = k_grid < 0.0
    elif branch == "positive":
        direction = k_grid > 0.0
    else:
        raise ValueError(f"unknown branch: {branch}")
    return target & direction, reference & direction


def enrichment(power: np.ndarray, target: np.ndarray, reference: np.ndarray) -> float:
    """Return log target/reference mean-power concentration per F-K cell."""
    if not np.any(target) or not np.any(reference):
        raise ValueError("empty F-K support")
    target_mean = float(np.mean(power[target]))
    reference_mean = float(np.mean(power[reference]))
    tiny = np.finfo(float).tiny
    return float(np.log((target_mean + tiny) / (reference_mean + tiny)))


def channel_permutation(
    channel_spectrum: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Destroy coordinate ordering while retaining channel zero and trace spectra."""
    order = np.arange(channel_spectrum.shape[0])
    order[1:] = rng.permutation(order[1:])
    return channel_spectrum[order]


def circular_time_shift_spectrum(
    channel_spectrum: np.ndarray,
    frequency: np.ndarray,
    fs: float,
    nsamples: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Apply independent time shifts as exact phase ramps in the passband."""
    shifts = np.zeros(channel_spectrum.shape[0], dtype=int)
    shifts[1:] = rng.integers(0, nsamples, size=channel_spectrum.shape[0] - 1)
    phase = np.exp(
        -2j * np.pi * shifts[:, None] * frequency[None, :] / fs
    )
    return channel_spectrum * phase


def surrogate_spectrum(
    channel_spectrum: np.ndarray,
    frequency: np.ndarray,
    fs: float,
    nsamples: int,
    method: str,
    rng: np.random.Generator,
) -> np.ndarray:
    """Construct one pre-spatial-FFT surrogate."""
    if method == "channel_permutation":
        return channel_permutation(channel_spectrum, rng)
    if method == "circular_time_shift":
        return circular_time_shift_spectrum(
            channel_spectrum, frequency, fs, nsamples, rng
        )
    raise ValueError(f"unknown null method: {method}")


def empirical_upper_p(observed: float, null: np.ndarray) -> float:
    """Return a plus-one upper-tail empirical probability."""
    return float((1 + np.sum(null >= observed)) / (null.size + 1))


def prepare_segment(raw: np.ndarray, fs: float) -> tuple[np.ndarray, float]:
    """Match production temporal preprocessing before safe stride decimation."""
    processed = preprocess(raw, fs, norm_seconds=5.0)
    decimated = processed[::SPACE_DECIMATION, ::TIME_DECIMATION]
    return decimated, fs / TIME_DECIMATION


def segment_statistics(
    x: np.ndarray,
    fs: float,
    dx: float,
    file_key: str,
    seed: int,
    null_start: int,
    nulls: int,
) -> tuple[dict[str, float], dict[str, dict[str, np.ndarray]], dict[str, np.ndarray]]:
    """Compute observed and surrogate pre-filter enrichment for one segment."""
    frequency, channel_spectrum = time_spectrum(x, fs)
    wavenumber, observed_power = spatial_power(channel_spectrum, dx)
    masks = {
        branch: supports(frequency, wavenumber, branch)
        for branch in BRANCHES
    }
    observed = {
        branch: enrichment(observed_power, *masks[branch])
        for branch in BRANCHES
    }
    null_values = {
        method: {branch: np.empty(nulls, dtype=float) for branch in BRANCHES}
        for method in NULL_METHODS
    }
    for method in NULL_METHODS:
        for local_index, realization in enumerate(
            range(null_start, null_start + nulls)
        ):
            rng = stable_rng(seed, method, realization, file_key)
            surrogate = surrogate_spectrum(
                channel_spectrum,
                frequency,
                fs,
                x.shape[1],
                method,
                rng,
            )
            _, power = spatial_power(surrogate, dx)
            for branch in BRANCHES:
                null_values[method][branch][local_index] = enrichment(
                    power, *masks[branch]
                )
    coordinates = {
        "frequency_hz": frequency,
        "wavenumber_cycles_m": wavenumber,
        "observed_power": observed_power,
        "negative_target_support": masks["negative"][0],
        "positive_target_support": masks["positive"][0],
        "negative_reference_support": masks["negative"][1],
        "positive_reference_support": masks["positive"][1],
    }
    return observed, null_values, coordinates


def load_rows(date: str, start: int, nfiles: int) -> pd.DataFrame:
    """Load an ordered slice of the SAFOD file catalog."""
    database = pd.read_csv(CSV, sep=r"\s+")
    database = database[database.nSamples > 0].copy()
    database["time"] = pd.to_datetime(
        database.startTime, utc=True, errors="coerce"
    )
    rows = database[
        database.time.dt.strftime("%Y-%m-%d") == date
    ].sort_values("time")
    return rows.iloc[start:start + nfiles]


def run(args: argparse.Namespace) -> tuple[Path, Path]:
    """Run and save one resumable pre-filter-energy batch."""
    if args.nulls < 1:
        raise ValueError("--nulls must be positive")
    rows = load_rows(args.date, args.start, args.nfiles)
    observed_values = {branch: [] for branch in BRANCHES}
    null_values = {
        method: {branch: [] for branch in BRANCHES}
        for method in NULL_METHODS
    }
    used_files: list[str] = []
    coordinates = None

    for row in rows.itertuples(index=False):
        path = corrected_path(row.file)
        if not path.exists():
            continue
        raw, fs, dx = load_segment(path)
        x, filtered_fs = prepare_segment(raw, fs)
        filtered_dx = dx * SPACE_DECIMATION
        observed, null, current_coordinates = segment_statistics(
            x,
            filtered_fs,
            filtered_dx,
            str(path),
            args.seed,
            args.null_start,
            args.nulls,
        )
        for branch in BRANCHES:
            observed_values[branch].append(observed[branch])
            for method in NULL_METHODS:
                null_values[method][branch].append(null[method][branch])
        used_files.append(str(path))
        if coordinates is None:
            coordinates = current_coordinates
        else:
            if not (
                np.array_equal(
                    coordinates["frequency_hz"],
                    current_coordinates["frequency_hz"],
                )
                and np.array_equal(
                    coordinates["wavenumber_cycles_m"],
                    current_coordinates["wavenumber_cycles_m"],
                )
            ):
                raise RuntimeError("F-K coordinates changed between input files")
            coordinates["observed_power"] += current_coordinates["observed_power"]
        print(f"processed {len(used_files)}/{len(rows)}", flush=True)

    if not used_files or coordinates is None:
        raise RuntimeError("no usable files")

    coordinates["observed_power"] /= len(used_files)

    observed_means = {
        branch: float(np.mean(observed_values[branch]))
        for branch in BRANCHES
    }
    arrays: dict[str, np.ndarray] = {
        "used_files": np.asarray(used_files),
        "null_realization_ids": np.arange(
            args.null_start, args.null_start + args.nulls, dtype=int
        ),
        **coordinates,
    }
    for branch in BRANCHES:
        arrays[f"observed_{branch}_per_file"] = np.asarray(
            observed_values[branch]
        )

    report: dict[str, object] = {
        "workflow_version": "ambient_fk_prefilter_energy_v1",
        "date": args.date,
        "start": args.start,
        "requested_files": args.nfiles,
        "used_files": len(used_files),
        "null_start": args.null_start,
        "null_realizations": args.nulls,
        "random_seed": args.seed,
        "passband_hz": list(PASSBAND_HZ),
        "target_velocity_m_s": list(TARGET_VELOCITY_M_S),
        "reference_velocity_m_s": list(REFERENCE_VELOCITY_M_S),
        "statistic": (
            "mean across files of log(mean target F-K power per cell / "
            "mean broad signed-reference F-K power per cell), measured before filtering"
        ),
        "saved_power_map": "mean pre-filter F-K power across all used files",
        "preprocessing": (
            "linear detrend, 5-20 Hz Butterworth bandpass, 5-s running-absolute-mean "
            "normalization, 2x temporal and spatial decimation, 10% Tukey tapers"
        ),
        "observed": {},
        "null_results": {},
        "supported_claim": (
            "A small upper-tail probability supports excess ordered input-wavefield "
            "energy in the predeclared target slowness interval before F-K filtering."
        ),
        "limitations": [
            "The statistic tests raw directional spectral concentration, not Green's-function convergence.",
            "The velocity interval is a prior from the Lellouch comparison and the existing production mask.",
            "The broad-reference normalization does not remove all colored spatial-noise structure.",
            "Channel permutation and independent time shifts are artificial nulls and are reported separately.",
        ],
    }

    for branch in BRANCHES:
        report["observed"][branch] = {
            "mean_log_enrichment": observed_means[branch],
            "geometric_enrichment_ratio": float(np.exp(observed_means[branch])),
        }
    familywise_observed = max(observed_means.values())

    for method in NULL_METHODS:
        method_report: dict[str, object] = {"branches": {}}
        familywise = np.full(args.nulls, -np.inf, dtype=float)
        for branch in BRANCHES:
            values = np.asarray(null_values[method][branch])
            mean_by_realization = np.mean(values, axis=0)
            arrays[f"null_{method}_{branch}_per_file"] = values
            arrays[f"null_{method}_{branch}_mean"] = mean_by_realization
            familywise = np.maximum(familywise, mean_by_realization)
            method_report["branches"][branch] = {
                "observed_mean_log_enrichment": observed_means[branch],
                "null95": float(np.quantile(mean_by_realization, 0.95)),
                "p_upper": empirical_upper_p(
                    observed_means[branch], mean_by_realization
                ),
            }
        arrays[f"null_{method}_familywise"] = familywise
        method_report["familywise"] = {
            "observed_maximum_across_branches": familywise_observed,
            "null95": float(np.quantile(familywise, 0.95)),
            "p_upper": empirical_upper_p(familywise_observed, familywise),
        }
        report["null_results"][method] = method_report

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.stem or (
        f"fk_prefilter_energy_v1_{args.date}_start{args.start}_n{args.nfiles}_"
        f"r{args.null_start}-{args.null_start + args.nulls - 1}"
    )
    json_path = args.output_dir / f"{stem}.json"
    npz_path = args.output_dir / f"{stem}.npz"
    if not args.overwrite and (json_path.exists() or npz_path.exists()):
        raise FileExistsError(f"output exists for {stem}")
    np.savez_compressed(npz_path, **arrays)
    json_path.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return json_path, npz_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Measure target-slowness energy before ambient F-K filtering."
    )
    parser.add_argument("--date", default="2024-12-20")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=1)
    parser.add_argument("--null-start", type=int, default=0)
    parser.add_argument("--nulls", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stem")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
