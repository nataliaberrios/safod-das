#!/usr/bin/env python3
"""Reproduce the published Lellouch et al. (2019) SAFOD ambient workflow.

This script is deliberately separate from the F-K-assisted production branch.
It implements the published top-virtual-source geometry, 30 s windows with
15 s overlap, running-absolute-mean temporal normalization, correlations to
receivers every 50 m, stacking over receiver R +/- 10 channels, a 3.2 km/s
alignment control, and a final 5--20 Hz correlation filter.  The paper does
not report the running-mean duration; the established project value of 5 s is
therefore retained and declared rather than treated as a published constant.

Raw 2024--2025 interrogator phase is proportional to strain.  Following the
paper, a time derivative is applied before temporal normalization, producing a
strain-rate proxy.  Only the channels required by the published geometry are
read, so a complete day can be processed without loading the full archive.
"""
from __future__ import annotations

import argparse
import ctypes
import gc
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, detrend, sosfiltfilt

from ambient_transfer_test import CSV, corrected_path


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "ambient_transfer" / "lellouch2019_reproduction_v1"
REFERENCE_VELOCITY_M_S = 3200.0
TARGET_SPACING_M = 50.0
LOCAL_HALF_WIDTH_CHANNELS = 10
WINDOW_SECONDS = 30.0
STEP_SECONDS = 15.0
RAM_SECONDS = 5.0
MAX_LAG_SECONDS = 0.35
OUTPUT_BAND_HZ = (5.0, 20.0)

try:
    LIBC = ctypes.CDLL("libc.so.6")
except OSError:  # pragma: no cover - non-Linux fallback
    LIBC = None


def release_allocator_memory() -> None:
    """Release per-record FFT/HDF5 temporaries during long batch runs."""
    gc.collect()
    if LIBC is not None:
        LIBC.malloc_trim(0)


def day_rows(date: str, start: int, nfiles: int | None) -> pd.DataFrame:
    """Return chronologically ordered, usable manifest rows for one UTC day."""
    database = pd.read_csv(CSV, sep=r"\s+")
    database = database[database.nSamples > 0].copy()
    database["time"] = pd.to_datetime(
        database.startTime, utc=True, errors="coerce"
    )
    selected = database[
        database.time.dt.strftime("%Y-%m-%d") == date
    ].sort_values("time")
    if nfiles is None:
        return selected.iloc[start:].reset_index(drop=True)
    return selected.iloc[start : start + nfiles].reset_index(drop=True)


def acquisition_metadata(path: Path) -> tuple[float, float, int, float]:
    """Return sample rate, channel spacing, channel count, and duration."""
    with h5py.File(path, "r") as handle:
        raw_group = handle["Acquisition/Raw[0]"]
        data = raw_group["RawData"]
        acquisition_attributes = handle["Acquisition"].attrs
        fs = float(raw_group.attrs.get("OutputDataRate", data.attrs.get("OutputDataRate", 500.0)))
        dx = float(acquisition_attributes.get("SpatialSamplingInterval", 1.0))
        duration = data.shape[0] / fs
        return fs, dx, int(data.shape[1]), float(duration)


def geometry(dx: float, n_channels: int) -> tuple[np.ndarray, np.ndarray]:
    """Return 50 m target channels and all required receiver channels."""
    targets = np.asarray(
        [
            int(round(TARGET_SPACING_M * index / dx))
            for index in range(1, 17)
        ],
        dtype=int,
    )
    targets = targets[targets + LOCAL_HALF_WIDTH_CHANNELS < n_channels]
    neighborhoods = np.asarray(
        [
            np.arange(
                target - LOCAL_HALF_WIDTH_CHANNELS,
                target + LOCAL_HALF_WIDTH_CHANNELS + 1,
                dtype=int,
            )
            for target in targets
        ]
    )
    required = np.unique(np.concatenate(([0], neighborhoods.ravel())))
    return targets, required


def load_required_channels(path: Path, channels: np.ndarray) -> tuple[np.ndarray, float]:
    """Read selected channels as channel x time float32 phase records."""
    with h5py.File(path, "r") as handle:
        raw_group = handle["Acquisition/Raw[0]"]
        data = raw_group["RawData"]
        fs = float(raw_group.attrs.get("OutputDataRate", data.attrs.get("OutputDataRate", 500.0)))
        selected = np.asarray(data[:, channels], dtype=np.float32).T
    return selected, fs


def strain_rate_and_ram(data: np.ndarray, fs: float) -> np.ndarray:
    """Convert phase/strain to a temporally normalized strain-rate proxy."""
    data = detrend(data, axis=1, type="linear").astype(np.float32, copy=False)
    rate = np.empty_like(data)
    rate[:, 0] = 0.0
    np.subtract(data[:, 1:], data[:, :-1], out=rate[:, 1:])
    rate *= np.float32(fs)
    window = max(3, int(round(RAM_SECONDS * fs)))
    amplitude = uniform_filter1d(
        np.abs(rate), size=window, axis=1, mode="nearest"
    )
    floor = np.percentile(amplitude, 5.0, axis=1, keepdims=True) * 0.1
    floor = np.maximum(floor, np.finfo(np.float32).eps)
    rate /= np.maximum(amplitude, floor)
    return rate


def correlation_window(
    window: np.ndarray,
    source_index: int,
    fs: float,
    max_lag_seconds: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Correlate one virtual source with every selected receiver."""
    n_time = window.shape[1]
    n_fft = 1 << int(np.ceil(np.log2(2 * n_time - 1)))
    maximum_lag = int(round(max_lag_seconds * fs))
    spectra = np.fft.rfft(window, n=n_fft, axis=1)
    cross_spectra = np.conj(spectra[source_index]) * spectra
    correlation = np.fft.irfft(cross_spectra, n=n_fft, axis=1)
    correlation = np.concatenate(
        (correlation[:, -maximum_lag:], correlation[:, : maximum_lag + 1]),
        axis=1,
    )
    correlation /= float(n_time)
    lags = np.arange(-maximum_lag, maximum_lag + 1, dtype=float) / fs
    return lags, correlation


def final_bandpass(correlation: np.ndarray, fs: float) -> np.ndarray:
    """Apply the paper's reported 5--20 Hz band to stacked correlations."""
    sos = butter(4, OUTPUT_BAND_HZ, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, correlation, axis=1)


def align_and_stack(
    correlations: np.ndarray,
    lags: np.ndarray,
    required_channels: np.ndarray,
    targets: np.ndarray,
    dx: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build isolated, local simple, and 3.2 km/s aligned receiver stacks."""
    lookup = {int(channel): index for index, channel in enumerate(required_channels)}
    isolated = []
    local_simple = []
    local_aligned = []
    for target in targets:
        nearby = np.arange(
            target - LOCAL_HALF_WIDTH_CHANNELS,
            target + LOCAL_HALF_WIDTH_CHANNELS + 1,
            dtype=int,
        )
        traces = np.asarray([correlations[lookup[int(ch)]] for ch in nearby])
        isolated.append(correlations[lookup[int(target)]])
        local_simple.append(np.mean(traces, axis=0))
        shifted = []
        for trace, channel in zip(traces, nearby):
            delta = (float(channel) - float(target)) * dx / REFERENCE_VELOCITY_M_S
            shifted.append(
                np.interp(lags + delta, lags, trace, left=0.0, right=0.0)
            )
        local_aligned.append(np.mean(shifted, axis=0))
    return np.asarray(isolated), np.asarray(local_simple), np.asarray(local_aligned)


def velocity_curve(
    section: np.ndarray,
    lags: np.ndarray,
    distances: np.ndarray,
    lag_sign: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a descriptive median moveout score over trial velocities."""
    velocities = np.arange(2000.0, 5000.1, 25.0)
    scale = np.max(np.abs(section), axis=1, keepdims=True)
    normalized = section / np.maximum(scale, np.finfo(float).eps)
    scores = []
    for velocity in velocities:
        samples = [
            trace[
                int(
                    np.argmin(
                        np.abs(lags - lag_sign * distance / velocity)
                    )
                )
            ]
            for trace, distance in zip(normalized, distances)
        ]
        scores.append(float(np.median(samples)))
    return velocities, np.asarray(scores)


def ridge_metrics(
    section: np.ndarray, lags: np.ndarray, distances: np.ndarray
) -> dict[str, float]:
    """Summarize causal and anti-causal fixed-velocity coherence."""
    velocities, causal = velocity_curve(section, lags, distances, 1.0)
    _, anti_causal = velocity_curve(section, lags, distances, -1.0)
    peak = int(np.argmax(np.abs(causal)))
    reference = int(np.argmin(np.abs(velocities - REFERENCE_VELOCITY_M_S)))
    return {
        "descriptive_causal_peak_velocity_m_s": float(velocities[peak]),
        "descriptive_causal_peak_signed_score": float(causal[peak]),
        "causal_score_at_3200_m_s": float(causal[reference]),
        "anti_causal_score_at_3200_m_s": float(anti_causal[reference]),
        "abs_causal_to_anti_causal_ratio_at_3200": float(
            abs(causal[reference])
            / max(abs(anti_causal[reference]), np.finfo(float).eps)
        ),
    }


def wiggle_panel(
    axis: plt.Axes,
    section: np.ndarray,
    lags: np.ndarray,
    distances: np.ndarray,
    title: str,
) -> None:
    """Draw a Lellouch-Figure-7c-style black wiggle section."""
    causal = lags >= 0.0
    time = lags[causal]
    traces = section[:, causal]
    trace_scale = np.percentile(np.max(np.abs(traces), axis=1), 80.0)
    trace_scale = max(float(trace_scale), np.finfo(float).eps)
    spacing = TARGET_SPACING_M * 0.35 / trace_scale
    for trace, distance in zip(traces, distances):
        plotted = distance + spacing * trace
        axis.plot(time, plotted, color="black", lw=0.65)
        axis.fill_between(
            time,
            distance,
            plotted,
            where=plotted >= distance,
            color="black",
            linewidth=0.0,
        )
    axis.plot(
        distances / REFERENCE_VELOCITY_M_S,
        distances,
        color="tab:red",
        lw=1.6,
        label="3.2 km s$^{-1}$ reference",
    )
    axis.set_xlim(0.0, MAX_LAG_SECONDS)
    axis.set_ylim(distances[-1] + 35.0, distances[0] - 35.0)
    axis.set_xlabel("Positive correlation lag (s)")
    axis.set_ylabel("Receiver position relative to channel 0 (m)")
    axis.set_title(title)
    axis.legend(frameon=False, fontsize=8, loc="lower right")


def run(args: argparse.Namespace) -> None:
    rows = day_rows(args.date, args.start, args.nfiles)
    if rows.empty:
        raise RuntimeError(f"No manifest rows for {args.date}")

    first_path = corrected_path(rows.iloc[0].file)
    fs, dx, n_channels, duration = acquisition_metadata(first_path)
    if not np.isclose(duration, 60.0, atol=0.1):
        raise ValueError(f"Expected one-minute records, found {duration:.3f} s")
    targets, required = geometry(dx, n_channels)
    source_index = int(np.flatnonzero(required == 0)[0])

    correlation_sum = None
    lags = None
    used_files: list[str] = []
    n_windows = 0
    previous_tail = None
    previous_end = None
    window_samples = int(round(WINDOW_SECONDS * fs))
    step_samples = int(round(STEP_SECONDS * fs))

    for row_number, row in enumerate(rows.itertuples(index=False), start=1):
        path = corrected_path(row.file)
        if not path.is_file():
            continue
        data, record_fs = load_required_channels(path, required)
        if not np.isclose(record_fs, fs):
            raise ValueError(f"Mixed sample rates: {fs} and {record_fs}")
        processed = strain_rate_and_ram(data, fs)
        row_time = pd.to_datetime(row.startTime, utc=True)
        contiguous = (
            previous_tail is not None
            and previous_end is not None
            and abs((row_time - previous_end).total_seconds()) < 0.02
        )
        windows = []
        if contiguous:
            windows.append(np.concatenate((previous_tail, processed[:, :step_samples]), axis=1))
        for start_sample in range(0, processed.shape[1] - window_samples + 1, step_samples):
            windows.append(processed[:, start_sample : start_sample + window_samples])

        for window in windows:
            window_lags, correlations = correlation_window(
                window, source_index, fs, MAX_LAG_SECONDS
            )
            if correlation_sum is None:
                correlation_sum = correlations
                lags = window_lags
            else:
                correlation_sum += correlations
            n_windows += 1
        previous_tail = processed[:, -step_samples:].copy()
        previous_end = row_time + pd.to_timedelta(processed.shape[1] / fs, unit="s")
        used_files.append(str(path))
        del data, processed, windows, window, correlations
        if row_number % 10 == 0:
            release_allocator_memory()
        if row_number % 25 == 0 or row_number == len(rows):
            print(
                f"processed {row_number}/{len(rows)} files; "
                f"usable={len(used_files)} windows={n_windows}",
                flush=True,
            )

    if correlation_sum is None or not used_files:
        raise RuntimeError("No usable records")
    average_correlation = correlation_sum / float(n_windows)
    filtered_correlation = final_bandpass(average_correlation, fs)
    isolated, simple, aligned = align_and_stack(
        filtered_correlation, lags, required, targets, dx
    )
    distances = targets.astype(float) * dx
    products = {
        "isolated_receiver": isolated,
        "local_21_channel_simple": simple,
        "local_21_channel_aligned_3200": aligned,
    }
    metrics = {
        name: ridge_metrics(section, lags, distances)
        for name, section in products.items()
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    requested_label = "all" if args.nfiles is None else str(args.nfiles)
    stem = (
        f"lellouch2019_{args.date}_start{args.start}_"
        f"requested{requested_label}_used{len(used_files)}"
    )
    np.savez_compressed(
        args.output_dir / f"{stem}.npz",
        lags_s=lags,
        distances_m=distances,
        target_channels=targets,
        required_channels=required,
        isolated_receiver=isolated,
        local_21_channel_simple=simple,
        local_21_channel_aligned_3200=aligned,
        used_files=np.asarray(used_files),
        fs_hz=fs,
        dx_m=dx,
        n_windows=n_windows,
    )
    report = {
        "workflow_version": "ambient_lellouch2019_reproduction_v1",
        "citation": "Lellouch et al. (2019), doi:10.1029/2019JB017533",
        "date": args.date,
        "start_file_index": args.start,
        "requested_files": None if args.nfiles is None else args.nfiles,
        "used_files": len(used_files),
        "used_30_s_windows": n_windows,
        "sample_rate_hz": fs,
        "channel_spacing_m": dx,
        "virtual_source_channel": 0,
        "virtual_source_coordinate_m": 0.0,
        "target_spacing_m": TARGET_SPACING_M,
        "target_channels": targets.tolist(),
        "target_coordinates_m": distances.tolist(),
        "published_steps_reproduced": {
            "time_derivative_to_strain_rate_proxy": True,
            "running_absolute_mean_normalization": True,
            "window_seconds": WINDOW_SECONDS,
            "window_overlap_seconds": WINDOW_SECONDS - STEP_SECONDS,
            "receiver_stack": "R-10 through R+10, same channel-0 source",
            "travel_time_alignment_velocity_m_s": REFERENCE_VELOCITY_M_S,
            "correlation_output_band_hz": list(OUTPUT_BAND_HZ),
        },
        "declared_project_choice_not_reported_by_paper": {
            "running_absolute_mean_window_seconds": RAM_SECONDS,
        },
        "important_departures_or_boundaries": [
            "The input is the 2024--2025 acquisition, not the June--July 2017 records used by Lellouch et al.",
            "The 2024--2025 channel-0 origin and measured spacing are assumed to inherit the same cemented-fiber coordinate convention.",
            "The 3.2 km/s aligned stack is a signal-enhancement operator fixed from the published average velocity, not an independent velocity estimate.",
            "No F-K filter is applied in this reproduction because Lellouch et al. did not report F-K filtering in their ambient-interferometry subsection.",
        ],
        "metrics": metrics,
    }
    (args.output_dir / f"{stem}.json").write_text(json.dumps(report, indent=2) + "\n")

    figure, axes = plt.subplots(1, 3, figsize=(16.5, 6.8), constrained_layout=True)
    wiggle_panel(axes[0], isolated, lags, distances, "a  Isolated 50 m receivers")
    wiggle_panel(axes[1], simple, lags, distances, "b  Published local stack: R $\\pm$ 10")
    wiggle_panel(
        axes[2],
        aligned,
        lags,
        distances,
        "c  Local stack aligned at 3.2 km s$^{-1}$",
    )
    figure.suptitle(
        f"Lellouch et al. (2019) ambient-workflow reproduction: {args.date}, "
        f"{len(used_files)} one-minute files ({n_windows} windows)",
        fontsize=13,
    )
    figure.savefig(args.output_dir / f"{stem}.png", dpi=350, bbox_inches="tight")
    figure.savefig(args.output_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(report, indent=2), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", default="2024-12-20")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument(
        "--nfiles",
        type=int,
        default=None,
        help="Number of one-minute files; omit for the remaining complete day.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
