#!/usr/bin/env python3
"""Inject broadband plane waves into real ambient data before preprocessing.

This sensitivity experiment separates two questions that are otherwise easy to
confound in the SAFOD ambient-noise workflow:

1. At what input amplitude does the independent, pre-filter F-K energy
   statistic recover an ordered plane wave?
2. At what input amplitude does the production signed F-K plus correlation
   workflow recover the same wave?

Synthetic waves are injected through a precision-safe linear decomposition:
real data and synthetic data are detrended and bandpassed separately, then
added in float64 before running-absolute-mean temporal normalization,
decimation, F-K selection, or correlation. This is mathematically equivalent
to adding before those linear operators without losing sub-count signals
against the large raw offset. Amplitude is expressed as injected broadband RMS
divided by the median 5--20 Hz RMS of the real channels in each file. Outputs
are scenario-specific and resumable; aggregation is performed separately.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, detrend, sosfiltfilt
from scipy.signal.windows import tukey

from ambient_fk_prefilter_energy_v1 import (
    BRANCHES,
    enrichment,
    load_rows,
    spatial_power,
    supports,
    time_spectrum,
)
from ambient_fk_transfer_test import fk_filter
from ambient_signed_fk_v2 import BRANCH_LAG_SIGN, VELOCITIES, velocity_scores
from ambient_transfer_test import (
    corrected_path,
    load_segment,
    normalized_corr_pairs,
    preprocess,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "ambient_transfer" / "fk_injection_recovery_v2"
PASSBAND_HZ = (5.0, 20.0)
SPACE_DECIMATION = 2
TIME_DECIMATION = 2
DEFAULT_AMPLITUDES = (0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0)


def stable_seed(seed: int, *tokens: object) -> int:
    """Return a deterministic 64-bit seed independent of Python hash state."""
    payload = "|".join([str(seed), *(str(token) for token in tokens)]).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def broadband_plane_wave(
    nchannels: int,
    nsamples: int,
    fs: float,
    dx: float,
    velocity_m_s: float,
    direction: int,
    seed: int,
) -> np.ndarray:
    """Return a unit-RMS, 5--20 Hz plane wave in the project convention.

    ``direction=+1`` propagates toward increasing fiber coordinate and is
    retained by F*K<0. ``direction=-1`` propagates toward decreasing fiber
    coordinate and is retained by F*K>0. The same random broadband source
    spectrum is phase shifted across channels; a Tukey envelope suppresses
    circular wraparound at the one-minute record edges.
    """
    if direction not in (-1, 1):
        raise ValueError("direction must be -1 or +1")
    if velocity_m_s <= 0:
        raise ValueError("velocity_m_s must be positive")
    rng = np.random.default_rng(seed)
    frequency = np.fft.rfftfreq(nsamples, 1.0 / fs)
    source = np.fft.rfft(rng.standard_normal(nsamples))
    source[(frequency < PASSBAND_HZ[0]) | (frequency > PASSBAND_HZ[1])] = 0.0
    coordinate = np.arange(nchannels, dtype=float)[:, None] * dx
    phase = np.exp(
        -1j * 2.0 * np.pi * frequency[None, :] * direction
        * coordinate / velocity_m_s
    )
    wave = np.fft.irfft(source[None, :] * phase, n=nsamples, axis=1)
    wave *= tukey(nsamples, alpha=0.1)[None, :]
    rms = np.sqrt(np.mean(wave ** 2, axis=1, keepdims=True))
    return (wave / np.maximum(rms, np.finfo(float).tiny)).astype(np.float32)


def real_band_rms(raw: np.ndarray, fs: float) -> float:
    """Return the median channel RMS in the injection passband."""
    x = detrend(raw, axis=1, type="linear")
    sos = butter(4, PASSBAND_HZ, btype="bandpass", fs=fs, output="sos")
    x = sosfiltfilt(sos, x, axis=1)
    channel_rms = np.sqrt(np.mean(x ** 2, axis=1))
    value = float(np.median(channel_rms[np.isfinite(channel_rms)]))
    if not np.isfinite(value) or value <= 0:
        raise ValueError("real-data band RMS is not positive and finite")
    return value


def preprocess_injected(
    raw: np.ndarray,
    wave: np.ndarray,
    amplitude_ratio: float,
    real_band_scale: float,
    fs: float,
    norm_seconds: float = 5.0,
) -> np.ndarray:
    """Apply a pre-detrend injection without float32 raw-offset rounding.

    Detrending and bandpass filtering are linear. We therefore apply both
    operations separately to the production float32 real data and unit-RMS
    synthetic wave, then add the two float64 filtered arrays before the
    nonlinear running-absolute-mean normalization. This retains the exact
    zero-injection baseline of the production preprocessing while preventing
    sub-count injections from being rounded away against the large raw integer
    offset.
    """
    if raw.shape != wave.shape:
        raise ValueError("raw and wave shapes differ")
    real_detrended = detrend(raw, axis=1, type="linear")
    wave_detrended = detrend(wave, axis=1, type="linear")
    sos = butter(4, PASSBAND_HZ, btype="bandpass", fs=fs, output="sos")
    real_band = sosfiltfilt(sos, real_detrended, axis=1)
    wave_band = sosfiltfilt(sos, wave_detrended, axis=1)
    x = real_band + float(amplitude_ratio * real_band_scale) * wave_band
    nwin = max(3, int(norm_seconds * fs))
    running = uniform_filter1d(np.abs(x), size=nwin, axis=1, mode="nearest")
    floor = np.percentile(running, 5, axis=1, keepdims=True) * 0.1 + 1e-12
    return x / np.maximum(running, floor)


def prefilter_enrichment(processed: np.ndarray, fs: float, dx: float) -> dict[str, float]:
    """Evaluate the independent statistic before any F-K selection."""
    x = processed[::SPACE_DECIMATION, ::TIME_DECIMATION]
    filtered_fs = fs / TIME_DECIMATION
    filtered_dx = dx * SPACE_DECIMATION
    frequency, channel_spectrum = time_spectrum(x, filtered_fs)
    wavenumber, power = spatial_power(channel_spectrum, filtered_dx)
    return {
        branch: enrichment(power, *supports(frequency, wavenumber, branch))
        for branch in BRANCHES
    }


def branch_correlations(
    processed: np.ndarray,
    fs: float,
    dx: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Run both production signed wedges and top-channel correlations."""
    products: dict[str, np.ndarray] = {}
    lags = distance = None
    for branch in BRANCHES:
        filtered, filtered_fs, filtered_dx = fk_filter(processed, fs, dx, branch)
        targets = [
            int(round(50.0 * index / filtered_dx))
            for index in range(1, 15)
        ]
        targets = [target for target in targets if target < filtered.shape[0]]
        current_lags, correlation = normalized_corr_pairs(
            filtered,
            [(0, target) for target in targets],
            filtered_fs,
            max_lag=0.35,
        )
        current_distance = np.asarray(targets, dtype=float) * filtered_dx
        if lags is None:
            lags, distance = current_lags, current_distance
        elif not (
            np.array_equal(lags, current_lags)
            and np.array_equal(distance, current_distance)
        ):
            raise RuntimeError("signed-branch correlation coordinates changed")
        products[branch] = correlation.astype(np.float32)
    assert lags is not None and distance is not None
    return lags, distance, products


def exact_velocity_score(
    top: np.ndarray,
    lags: np.ndarray,
    distance: np.ndarray,
    velocity_m_s: float,
    lag_sign: float,
) -> float:
    """Return the median signed correlation along one fixed trajectory."""
    values = [
        row[np.argmin(np.abs(lags - lag_sign * offset / velocity_m_s))]
        for row, offset in zip(top, distance)
    ]
    return float(np.median(values))


def parse_amplitudes(value: str) -> tuple[float, ...]:
    amplitudes = tuple(float(item) for item in value.split(","))
    if not amplitudes or any(item < 0 for item in amplitudes):
        raise ValueError("amplitudes must be a non-empty nonnegative list")
    if 0.0 not in amplitudes:
        raise ValueError("amplitudes must include zero for the paired baseline")
    if len(set(amplitudes)) != len(amplitudes):
        raise ValueError("amplitudes must be unique")
    return tuple(sorted(amplitudes))


def run(args: argparse.Namespace) -> tuple[Path, Path]:
    """Run one velocity-direction scenario over an ordered set of real files."""
    amplitudes = parse_amplitudes(args.amplitudes)
    if args.direction not in (-1, 1):
        raise ValueError("--direction must be -1 or +1")
    rows = load_rows(args.date, args.start, args.nfiles)
    injected_branch = "negative" if args.direction == 1 else "positive"
    opposite_branch = "positive" if injected_branch == "negative" else "negative"
    per_file_prefilter = {
        amplitude: {branch: [] for branch in BRANCHES}
        for amplitude in amplitudes
    }
    per_file_correlations = {
        amplitude: {branch: [] for branch in BRANCHES}
        for amplitude in amplitudes
    }
    per_file_band_rms: list[float] = []
    used_files: list[str] = []
    lags = distance = None

    for row in rows.itertuples(index=False):
        path = corrected_path(row.file)
        if not path.exists():
            continue
        raw, fs, dx = load_segment(path)
        scale = real_band_rms(raw, fs)
        wave = broadband_plane_wave(
            raw.shape[0],
            raw.shape[1],
            fs,
            dx,
            args.velocity,
            args.direction,
            stable_seed(args.seed, path, args.velocity, args.direction),
        )
        for amplitude in amplitudes:
            processed = preprocess_injected(raw, wave, amplitude, scale, fs)
            prefilter = prefilter_enrichment(processed, fs, dx)
            current_lags, current_distance, correlations = branch_correlations(
                processed, fs, dx
            )
            if lags is None:
                lags, distance = current_lags, current_distance
            elif not (
                np.array_equal(lags, current_lags)
                and np.array_equal(distance, current_distance)
            ):
                raise RuntimeError("coordinates changed between files")
            for branch in BRANCHES:
                per_file_prefilter[amplitude][branch].append(prefilter[branch])
                per_file_correlations[amplitude][branch].append(correlations[branch])
        per_file_band_rms.append(scale)
        used_files.append(str(path))
        print(f"processed {len(used_files)}/{len(rows)}", flush=True)

    if not used_files or lags is None or distance is None:
        raise RuntimeError("no usable real-data files")

    arrays: dict[str, np.ndarray] = {
        "amplitude_ratios": np.asarray(amplitudes),
        "lags_s": lags,
        "distance_m": distance,
        "velocities_m_s": VELOCITIES,
        "used_files": np.asarray(used_files),
        "real_band_rms_per_file": np.asarray(per_file_band_rms),
    }
    report: dict[str, object] = {
        "workflow_version": "ambient_fk_injection_recovery_v2",
        "date": args.date,
        "start": args.start,
        "requested_files": args.nfiles,
        "used_files": len(used_files),
        "velocity_m_s": args.velocity,
        "direction": args.direction,
        "injected_branch": injected_branch,
        "opposite_branch": opposite_branch,
        "amplitude_definition": (
            "injected broadband plane-wave RMS / median real-channel 5-20 Hz RMS, "
            "computed separately for each one-minute file"
        ),
        "injection_point": (
            "mathematically before the linear detrend and 5-20 Hz bandpass; real and "
            "synthetic arrays are detrended and bandpassed separately, then added in "
            "float64 before 5-s normalization to avoid rounding against the large raw "
            "integer offset, decimation, F-K filtering, and correlation"
        ),
        "numerical_precision": {
            "raw_load_dtype": "float32 (unchanged production loader)",
            "linear_operator_output_dtype": "float64",
            "injection_addition_dtype": "float64",
            "zero_baseline": (
                "identical production detrend, 5-20 Hz bandpass, and 5-s "
                "running-absolute-mean normalization at amplitude zero"
            ),
        },
        "signal": "deterministic random broadband 5-20 Hz plane wave with 10% Tukey time taper",
        "amplitudes": {},
        "interpretive_boundary": (
            "This is a processing sensitivity calibration in real noise. It does not "
            "prove that the uninjected ambient ridge is a physical Green's function."
        ),
    }

    for amplitude in amplitudes:
        tag = f"a{amplitude:.6g}".replace(".", "p")
        report["amplitudes"][str(amplitude)] = {"prefilter": {}, "postfilter": {}}
        for branch in BRANCHES:
            pre_values = np.asarray(per_file_prefilter[amplitude][branch])
            correlations = np.asarray(per_file_correlations[amplitude][branch])
            top = np.mean(correlations, axis=0)
            physical = velocity_scores(
                top, lags, distance, BRANCH_LAG_SIGN[branch]
            )
            leakage = velocity_scores(
                top, lags, distance, -BRANCH_LAG_SIGN[branch]
            )
            peak = int(np.nanargmax(np.abs(physical)))
            exact = exact_velocity_score(
                top,
                lags,
                distance,
                args.velocity,
                BRANCH_LAG_SIGN[branch],
            )
            arrays[f"{tag}_{branch}_prefilter_per_file"] = pre_values
            arrays[f"{tag}_{branch}_correlation_per_file"] = correlations
            arrays[f"{tag}_{branch}_top"] = top
            arrays[f"{tag}_{branch}_physical_scores"] = physical
            arrays[f"{tag}_{branch}_leakage_scores"] = leakage
            report["amplitudes"][str(amplitude)]["prefilter"][branch] = {
                "mean_log_enrichment": float(np.mean(pre_values)),
                "geometric_enrichment_ratio": float(np.exp(np.mean(pre_values))),
            }
            report["amplitudes"][str(amplitude)]["postfilter"][branch] = {
                "peak_velocity_m_s": float(VELOCITIES[peak]),
                "peak_absolute_score": float(abs(physical[peak])),
                "score_at_injected_velocity": exact,
            }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    direction_tag = "inc" if args.direction == 1 else "dec"
    stem = args.stem or (
        f"fk_injection_v2_{args.date}_start{args.start}_n{args.nfiles}_"
        f"v{int(round(args.velocity))}_{direction_tag}"
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
        description="Inject plane waves into real noise before ambient F-K processing."
    )
    parser.add_argument("--date", default="2024-12-20")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=10)
    parser.add_argument("--velocity", type=float, required=True)
    parser.add_argument("--direction", type=int, choices=(-1, 1), required=True)
    parser.add_argument(
        "--amplitudes",
        default=",".join(str(value) for value in DEFAULT_AMPLITUDES),
    )
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--stem")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
