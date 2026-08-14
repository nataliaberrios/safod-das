#!/usr/bin/env python3
"""Paper-faithful Lellouch et al. (2019) Figure 7c reproduction.

The baseline implements only operations reported in section 4.1:

1. time derivative from recorded strain/phase to a strain-rate proxy;
2. running-absolute-mean (RAM) temporal normalization;
3. 30 s correlation windows with 15 s overlap over a continuous day;
4. fixed top virtual source and receiver centers every 50 m;
5. C_S,R = sum_{Z=R-10}^{R+10} C_S,Z, with the same source channel;
6. simple (unshifted) stacking for Figure 7c; and
7. a 5--20 Hz bandpass applied to the resulting stacked correlations.

The paper does not report a RAM duration.  The default 0.1 s follows Bensen et
al. (2007), who recommend approximately half the maximum period of the analysis
band (0.5 * 1/5 Hz = 0.1 s).  RAM duration is therefore a declared sensitivity,
not a claimed Lellouch constant.

Common-mode subtraction and stabilized source-power division are available only
as labelled sensitivity branches.  F-K filtering is intentionally absent: it is
not reported in the paper's ambient-interferometry section.  Hourly chunk files
store summed cross spectra, allowing exact full-day aggregation and convergence
analysis without rebuilding raw-data correlations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, hilbert, sosfiltfilt


HERE = Path(__file__).resolve().parent
CSV = Path(
    "/oak/stanford/groups/ettore88/data/SAFOD/"
    "SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv"
)
DEFAULT_OUT = HERE / "ambient_transfer" / "lellouch2019_exact_stack"
CITATION = "Lellouch et al. (2019), doi:10.1029/2019JB017533"
BENSEN_CITATION = "Bensen et al. (2007), doi:10.1111/j.1365-246X.2007.03374.x"
WINDOW_SECONDS = 30.0
STEP_SECONDS = 15.0
MAX_LAG_SECONDS = 0.35
OUTPUT_BAND_HZ = (5.0, 20.0)
TARGET_OFFSETS_M = np.arange(50.0, 700.1, 50.0)
NEIGHBOR_HALF_WIDTH = 10
REFERENCE_VELOCITY_M_S = 3200.0
VELOCITY_GRID_M_S = np.arange(1500.0, 6000.1, 25.0)
VALID_NULLS = ("ordered", "white_noise", "channel_permutation")
VALID_SPECTRAL_MODES = ("cross_correlation", "source_power_stabilized")
WORKFLOW_VERSION = "lellouch2019_exact_stack_v2"


def corrected_path(value: str | Path) -> Path:
    """Translate the original transfer-manifest prefix to the mounted archive."""
    return Path(
        str(value).replace(
            "/data/SAFODAS1-harddrive-transfer",
            "/data/SAFOD/SAFODAS1-harddrive-transfer",
        )
    )


def day_rows(date: str) -> pd.DataFrame:
    """Return every usable manifest row for a UTC date in chronological order."""
    database = pd.read_csv(CSV, sep=r"\s+")
    database = database[database.nSamples > 0].copy()
    database["time"] = pd.to_datetime(
        database.startTime, utc=True, errors="coerce", format="mixed"
    )
    rows = database[
        database.time.dt.strftime("%Y-%m-%d") == date
    ].sort_values("time").reset_index(drop=True)
    if rows.empty:
        raise RuntimeError(f"No manifest rows for {date}")
    return rows


def acquisition_metadata(path: Path) -> dict[str, float | int | str]:
    """Read the acquisition quantities that define sampling and data provenance."""
    with h5py.File(path, "r") as handle:
        raw_group = handle["Acquisition/Raw[0]"]
        data = raw_group["RawData"]
        acquisition = handle["Acquisition"].attrs
        fs = float(
            raw_group.attrs.get(
                "OutputDataRate", data.attrs.get("OutputDataRate", 500.0)
            )
        )
        return {
            "sample_rate_hz": fs,
            "channel_spacing_m": float(acquisition["SpatialSamplingInterval"]),
            "gauge_length_m": float(acquisition["GaugeLength"]),
            "start_locus_index": int(acquisition["StartLocusIndex"]),
            "number_of_channels": int(data.shape[1]),
            "record_samples": int(data.shape[0]),
            "record_seconds": float(data.shape[0] / fs),
            "raw_data_unit": raw_group.attrs["RawDataUnit"].decode(),
            "raw_description": raw_group.attrs["RawDescription"].decode(),
        }


def validate_continuity(rows: pd.DataFrame, record_seconds: float) -> None:
    """Require the selected day to be continuous at the manifest time precision."""
    if len(rows) < 2:
        return
    differences = np.diff(rows.time.astype("int64").to_numpy()) / 1.0e9
    bad = np.flatnonzero(np.abs(differences - record_seconds) > 0.02)
    if bad.size:
        examples = [
            (int(index), float(differences[index])) for index in bad[:8]
        ]
        raise RuntimeError(f"Manifest is not continuous: {examples}")


def geometry(
    source_channel: int,
    dx: float,
    number_of_channels: int,
) -> dict[str, np.ndarray]:
    """Build the fixed-top Figure 7c receiver centers and R±10 neighborhoods."""
    centers = source_channel + np.rint(TARGET_OFFSETS_M / dx).astype(int)
    keep = centers + NEIGHBOR_HALF_WIDTH < number_of_channels
    centers = centers[keep]
    offsets = (centers - source_channel).astype(float) * dx
    neighborhoods = np.stack(
        [
            np.arange(
                center - NEIGHBOR_HALF_WIDTH,
                center + NEIGHBOR_HALF_WIDTH + 1,
                dtype=int,
            )
            for center in centers
        ]
    )
    required_channels = np.unique(
        np.concatenate(([source_channel], centers, neighborhoods.ravel()))
    )
    if required_channels.min() < 0 or required_channels.max() >= number_of_channels:
        raise ValueError("Figure 7c geometry extends outside the recorded array")
    lookup = {int(channel): index for index, channel in enumerate(required_channels)}
    return {
        "source_channel": np.asarray(source_channel),
        "centers": centers,
        "offsets_m": offsets,
        "neighborhoods": neighborhoods,
        "required_channels": required_channels,
        "source_local": np.asarray(lookup[source_channel]),
        "center_local": np.asarray([lookup[int(value)] for value in centers]),
        "neighbor_local": np.asarray(
            [[lookup[int(value)] for value in row] for row in neighborhoods]
        ),
    }


def stable_seed(seed: int, *parts: object) -> int:
    token = "|".join([str(seed), *map(str, parts)]).encode()
    return int.from_bytes(hashlib.sha256(token).digest()[:8], "little")


def receiver_mapping(
    required_channels: np.ndarray,
    source_channel: int,
    null_method: str,
    seed: int,
    realization: int,
) -> np.ndarray:
    """Map analysis coordinates to raw channels for a stable receiver scramble."""
    mapped = required_channels.copy()
    if null_method != "channel_permutation":
        return mapped
    receiver_mask = mapped != source_channel
    receiver_values = mapped[receiver_mask].copy()
    rng = np.random.default_rng(
        stable_seed(seed, "receiver-permutation", source_channel, realization)
    )
    mapped[receiver_mask] = rng.permutation(receiver_values)
    return mapped


def read_file_channels(path: Path, mapped_channels: np.ndarray) -> np.ndarray:
    """Read channels in analysis order, even if a null permutes physical channels."""
    unique_channels = np.unique(mapped_channels)
    with h5py.File(path, "r") as handle:
        data = handle["Acquisition/Raw[0]/RawData"]
        selected = np.asarray(data[:, unique_channels], dtype=np.float32).T
    lookup = {int(channel): index for index, channel in enumerate(unique_channels)}
    order = np.asarray([lookup[int(channel)] for channel in mapped_channels])
    return selected[order]


def white_noise_file(
    shape: tuple[int, int],
    seed: int,
    date: str,
    file_index: int,
    realization: int,
) -> np.ndarray:
    """Generate reproducible iid broadband raw-input noise file by file."""
    rng = np.random.default_rng(
        stable_seed(seed, "white-noise", date, file_index, realization)
    )
    return rng.standard_normal(shape, dtype=np.float32)


def odd_ram_samples(ram_seconds: float, fs: float) -> int:
    samples = max(1, int(round(ram_seconds * fs)))
    return samples if samples % 2 else samples + 1


def strain_rate_and_ram_continuous(
    raw: np.ndarray,
    fs: float,
    ram_seconds: float,
    common_mode: bool,
) -> tuple[np.ndarray, float]:
    """Differentiate and apply centered Bensen running-absolute-mean weights.

    No detrend, bandpass, taper, percentile floor, or per-window normalization is
    applied.  A machine-scale denominator floor prevents division by exact zero
    and is reported as a fraction of weights.
    """
    rate = np.empty_like(raw, dtype=np.float32)
    rate[:, 0] = 0.0
    np.subtract(raw[:, 1:], raw[:, :-1], out=rate[:, 1:])
    rate *= np.float32(fs)
    del raw
    if common_mode:
        rate -= np.median(rate, axis=0, keepdims=True).astype(np.float32)
    absolute = np.abs(rate)
    weights = uniform_filter1d(
        absolute,
        size=odd_ram_samples(ram_seconds, fs),
        axis=1,
        mode="nearest",
    )
    del absolute
    scale = float(np.nanmedian(weights))
    floor = max(np.finfo(np.float32).tiny, np.finfo(np.float32).eps * scale)
    floored_fraction = float(np.mean(weights <= floor))
    rate /= np.maximum(weights, np.float32(floor))
    return rate, floored_fraction


def cross_spectra_for_starts(
    normalized: np.ndarray,
    local_starts: np.ndarray,
    fs: float,
    geom: dict[str, np.ndarray],
    spectral_mode: str,
    waterlevel: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Accumulate central-receiver and published R±10 cross spectra."""
    window_samples = int(round(WINDOW_SECONDS * fs))
    n_fft = 1 << int(math.ceil(math.log2(2 * window_samples - 1)))
    n_frequency = n_fft // 2 + 1
    n_targets = len(geom["centers"])
    central_sum = np.zeros((n_targets, n_frequency), dtype=np.complex128)
    neighbor_sum = np.zeros_like(central_sum)
    source_local = int(geom["source_local"])
    center_local = geom["center_local"]
    neighbor_local = geom["neighbor_local"]
    for start in local_starts:
        stop = int(start) + window_samples
        source = normalized[source_local, int(start):stop]
        central = normalized[center_local, int(start):stop]
        neighbors = normalized[neighbor_local, int(start):stop].sum(axis=1)
        source_spectrum = np.fft.rfft(source, n=n_fft)
        central_spectrum = np.fft.rfft(central, n=n_fft, axis=1)
        neighbor_spectrum = np.fft.rfft(neighbors, n=n_fft, axis=1)
        cross_central = np.conj(source_spectrum)[None, :] * central_spectrum
        cross_neighbors = np.conj(source_spectrum)[None, :] * neighbor_spectrum
        if spectral_mode == "source_power_stabilized":
            power = np.abs(source_spectrum) ** 2
            denominator = np.maximum(
                power, waterlevel * max(float(np.max(power)), np.finfo(float).eps)
            )
            cross_central /= denominator[None, :]
            cross_neighbors /= denominator[None, :]
        central_sum += cross_central
        neighbor_sum += cross_neighbors
    return central_sum, neighbor_sum, n_fft


def chunk_tag(args: argparse.Namespace) -> str:
    ram = ("%.6g" % args.ram_seconds).replace(".", "p")
    return (
        f"{args.date}_src{args.source_channel}_ram{ram}_"
        f"{args.spectral_mode}_{args.null_method}_r{args.realization}"
        f"{'_cm' if args.common_mode else ''}"
    )


def chunk_path(args: argparse.Namespace, start: int, nfiles: int) -> Path:
    return args.output_dir / (
        f"chunk_{chunk_tag(args)}_start{start:04d}_n{nfiles:04d}.npz"
    )


def run_chunk(args: argparse.Namespace) -> None:
    """Process a core file range with context halos and save summed spectra."""
    rows = day_rows(args.date)
    first_path = corrected_path(rows.iloc[0].file)
    metadata = acquisition_metadata(first_path)
    fs = float(metadata["sample_rate_hz"])
    dx = float(metadata["channel_spacing_m"])
    record_samples = int(metadata["record_samples"])
    record_seconds = float(metadata["record_seconds"])
    validate_continuity(rows, record_seconds)
    if args.ram_seconds > record_seconds:
        raise ValueError("RAM duration exceeds one-file context halo")
    core_first = int(args.start)
    core_stop = min(core_first + int(args.nfiles), len(rows))
    if core_first < 0 or core_first >= core_stop:
        raise ValueError("Empty core range")
    geom = geometry(args.source_channel, dx, int(metadata["number_of_channels"]))
    mapped_channels = receiver_mapping(
        geom["required_channels"],
        args.source_channel,
        args.null_method,
        args.seed,
        args.realization,
    )
    central_total = None
    neighbor_total = None
    total_windows = 0
    total_floored = 0.0
    blocks_used = 0
    n_fft = None
    total_day_samples = len(rows) * record_samples
    window_samples = int(round(WINDOW_SECONDS * fs))
    step_samples = int(round(STEP_SECONDS * fs))

    for block_first in range(core_first, core_stop, args.block_files):
        block_stop = min(block_first + args.block_files, core_stop)
        context_first = max(0, block_first - 1)
        context_stop = min(len(rows), block_stop + 1)
        pieces = []
        used_paths = []
        for file_index in range(context_first, context_stop):
            path = corrected_path(rows.iloc[file_index].file)
            if not path.is_file():
                raise FileNotFoundError(path)
            if args.null_method == "white_noise":
                piece = white_noise_file(
                    (len(mapped_channels), record_samples),
                    args.seed,
                    args.date,
                    file_index,
                    args.realization,
                )
            else:
                piece = read_file_channels(path, mapped_channels)
            if piece.shape[1] != record_samples:
                raise ValueError(f"Unexpected record shape for {path}: {piece.shape}")
            pieces.append(piece)
            used_paths.append(str(path))
        raw = np.concatenate(pieces, axis=1)
        del pieces
        normalized, floored_fraction = strain_rate_and_ram_continuous(
            raw, fs, args.ram_seconds, args.common_mode
        )
        core_sample_first = block_first * record_samples
        core_sample_stop = block_stop * record_samples
        global_starts = np.arange(
            core_sample_first, core_sample_stop, step_samples, dtype=np.int64
        )
        global_starts = global_starts[
            global_starts + window_samples <= total_day_samples
        ]
        local_starts = global_starts - context_first * record_samples
        if np.any(local_starts < 0) or np.any(
            local_starts + window_samples > normalized.shape[1]
        ):
            raise RuntimeError("Context halo is insufficient for assigned windows")
        central, neighbors, this_n_fft = cross_spectra_for_starts(
            normalized,
            local_starts,
            fs,
            geom,
            args.spectral_mode,
            args.waterlevel,
        )
        if central_total is None:
            central_total = np.zeros_like(central)
            neighbor_total = np.zeros_like(neighbors)
            n_fft = this_n_fft
        if this_n_fft != n_fft:
            raise RuntimeError("FFT length changed within chunk")
        central_total += central
        neighbor_total += neighbors
        total_windows += len(local_starts)
        total_floored += floored_fraction
        blocks_used += 1
        print(
            f"block {block_first}:{block_stop}; windows={len(local_starts)}; "
            f"cumulative={total_windows}",
            flush=True,
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = chunk_path(args, core_first, core_stop - core_first)
    np.savez_compressed(
        output,
        central_cross_spectrum_sum=central_total,
        neighbor_cross_spectrum_sum=neighbor_total,
        n_windows=np.asarray(total_windows),
        n_fft=np.asarray(n_fft),
        fs_hz=np.asarray(fs),
        dx_m=np.asarray(dx),
        source_channel=np.asarray(args.source_channel),
        center_channels=geom["centers"],
        offsets_m=geom["offsets_m"],
        required_channels=geom["required_channels"],
        mapped_channels=mapped_channels,
        start_file=np.asarray(core_first),
        used_core_files=np.asarray(core_stop - core_first),
        ram_seconds=np.asarray(args.ram_seconds),
        ram_samples=np.asarray(odd_ram_samples(args.ram_seconds, fs)),
        denominator_floored_fraction=np.asarray(total_floored / blocks_used),
        common_mode=np.asarray(args.common_mode),
        spectral_mode=np.asarray(args.spectral_mode),
        null_method=np.asarray(args.null_method),
        realization=np.asarray(args.realization),
    )
    print(f"wrote {output}")
    print(f"core files={core_stop-core_first}; windows={total_windows}")


def correlation_from_spectrum(
    cross_spectrum: np.ndarray,
    n_windows: int,
    n_fft: int,
    fs: float,
) -> tuple[np.ndarray, np.ndarray]:
    average = cross_spectrum / float(n_windows)
    circular = np.fft.irfft(average, n=n_fft, axis=1)
    maximum_lag = int(round(MAX_LAG_SECONDS * fs))
    correlation = np.concatenate(
        (circular[:, -maximum_lag:], circular[:, : maximum_lag + 1]), axis=1
    )
    lags = np.arange(-maximum_lag, maximum_lag + 1, dtype=float) / fs
    return lags, correlation


def final_bandpass(correlation: np.ndarray, fs: float) -> np.ndarray:
    sos = butter(4, OUTPUT_BAND_HZ, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, correlation, axis=1)


def normalized_envelope(section: np.ndarray) -> np.ndarray:
    envelope = np.abs(hilbert(section, axis=1))
    scale = np.median(envelope, axis=1, keepdims=True)
    return envelope / np.maximum(scale, np.finfo(float).eps)


def moveout_scores(
    section: np.ndarray,
    lags: np.ndarray,
    offsets: np.ndarray,
    velocities: np.ndarray = VELOCITY_GRID_M_S,
    sign: float = 1.0,
    half_window_s: float = 0.012,
) -> np.ndarray:
    """Median normalized envelope in a fixed gate along each trial moveout."""
    envelope = normalized_envelope(section)
    scores = []
    for velocity in velocities:
        values = []
        for trace, offset in zip(envelope, offsets):
            center = sign * offset / velocity
            mask = np.abs(lags - center) <= half_window_s
            values.append(float(np.mean(trace[mask])) if mask.any() else np.nan)
        scores.append(float(np.nanmedian(values)))
    return np.asarray(scores)


def receiver_order_null(
    section: np.ndarray,
    lags: np.ndarray,
    offsets: np.ndarray,
    observed_peak: float,
    seed: int,
    count: int,
) -> tuple[np.ndarray, float]:
    """Repeat the full velocity scan after permuting receiver-coordinate groups."""
    rng = np.random.default_rng(stable_seed(seed, "aggregate-row-null"))
    maxima = np.empty(count, dtype=float)
    for index in range(count):
        permuted = section[rng.permutation(section.shape[0])]
        maxima[index] = np.max(moveout_scores(permuted, lags, offsets))
    p_value = (1.0 + np.sum(maxima >= observed_peak)) / (count + 1.0)
    return maxima, float(p_value)


def trace_normalize_for_display(section: np.ndarray) -> np.ndarray:
    scale = np.percentile(np.abs(section), 99.0, axis=1, keepdims=True)
    return section / np.maximum(scale, np.finfo(float).eps)


def plot_wiggles(
    axis: plt.Axes,
    section: np.ndarray,
    lags: np.ndarray,
    offsets: np.ndarray,
    title: str,
    reference_velocity: float = REFERENCE_VELOCITY_M_S,
) -> None:
    shown = trace_normalize_for_display(section)
    amplitude = 14.0
    for trace, offset in zip(shown, offsets):
        curve = offset + amplitude * trace
        axis.plot(lags, curve, color="black", linewidth=0.55)
    axis.plot(
        offsets / reference_velocity,
        offsets,
        color="tab:red",
        linestyle="--",
        linewidth=1.5,
        label=f"{reference_velocity/1000:.1f} km s$^{{-1}}$",
    )
    axis.set(
        xlim=(0.0, MAX_LAG_SECONDS),
        ylim=(offsets[-1] + 30.0, offsets[0] - 30.0),
        xlabel="Positive correlation lag (s)",
        ylabel="Receiver offset from virtual source (m)",
        title=title,
    )
    axis.legend(frameon=False, fontsize=8, loc="lower right")


def aggregate(args: argparse.Namespace) -> None:
    rows = day_rows(args.date)
    expected_total = min(args.total_files, len(rows))
    starts = list(range(0, expected_total, args.nfiles))
    files = [
        chunk_path(args, start, min(args.nfiles, expected_total - start))
        for start in starts
    ]
    missing = [str(path) for path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing exact-stack chunks:\n" + "\n".join(missing[:20])
        )
    central_sum = None
    neighbor_sum = None
    n_windows = 0
    core_files = 0
    metadata_reference = None
    for path in files:
        with np.load(path, allow_pickle=False) as data:
            metadata = {
                "fs_hz": float(data["fs_hz"]),
                "dx_m": float(data["dx_m"]),
                "n_fft": int(data["n_fft"]),
                "source_channel": int(data["source_channel"]),
                "center_channels": data["center_channels"].tolist(),
                "offsets_m": data["offsets_m"].tolist(),
                "ram_seconds": float(data["ram_seconds"]),
                "ram_samples": int(data["ram_samples"]),
                "common_mode": bool(data["common_mode"]),
                "spectral_mode": str(data["spectral_mode"]),
                "null_method": str(data["null_method"]),
                "realization": int(data["realization"]),
            }
            if metadata_reference is None:
                metadata_reference = metadata
                central_sum = np.zeros_like(data["central_cross_spectrum_sum"])
                neighbor_sum = np.zeros_like(data["neighbor_cross_spectrum_sum"])
            elif metadata != metadata_reference:
                raise ValueError(f"Chunk metadata mismatch at {path}")
            central_sum += data["central_cross_spectrum_sum"]
            neighbor_sum += data["neighbor_cross_spectrum_sum"]
            n_windows += int(data["n_windows"])
            core_files += int(data["used_core_files"])
    fs = metadata_reference["fs_hz"]
    n_fft = metadata_reference["n_fft"]
    offsets = np.asarray(metadata_reference["offsets_m"], dtype=float)
    lags, central = correlation_from_spectrum(central_sum, n_windows, n_fft, fs)
    _, neighbors = correlation_from_spectrum(neighbor_sum, n_windows, n_fft, fs)
    central_filtered = final_bandpass(central, fs)
    neighbors_filtered = final_bandpass(neighbors, fs)
    causal_scores = moveout_scores(neighbors_filtered, lags, offsets, sign=1.0)
    acausal_scores = moveout_scores(neighbors_filtered, lags, offsets, sign=-1.0)
    peak_index = int(np.argmax(causal_scores))
    observed_peak = float(causal_scores[peak_index])
    null_maxima, p_value = receiver_order_null(
        neighbors_filtered,
        lags,
        offsets,
        observed_peak,
        args.seed,
        args.null_count,
    )
    reference_index = int(
        np.argmin(np.abs(VELOCITY_GRID_M_S - REFERENCE_VELOCITY_M_S))
    )
    record_seconds = float(
        acquisition_metadata(corrected_path(rows.iloc[0].file))["record_seconds"]
    )
    if core_files == len(rows):
        expected_windows = int(
            math.floor((core_files * record_seconds - WINDOW_SECONDS) / STEP_SECONDS) + 1
        )
    else:
        # Each partial aggregate owns all window starts inside its core duration;
        # its final start may extend into the following context file.
        expected_windows = int(round(core_files * record_seconds / STEP_SECONDS))
    if core_files == len(rows) and n_windows != expected_windows:
        raise RuntimeError(
            f"Full-day window count mismatch: {n_windows} != {expected_windows}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.output_dir / f"aggregate_{chunk_tag(args)}"
    np.savez_compressed(
        stem.with_suffix(".npz"),
        lags_s=lags,
        offsets_m=offsets,
        center_channels=np.asarray(metadata_reference["center_channels"]),
        central_receiver_correlation=central_filtered,
        r_plus_minus_10_correlation=neighbors_filtered,
        velocity_grid_m_s=VELOCITY_GRID_M_S,
        causal_moveout_scores=causal_scores,
        acausal_moveout_scores=acausal_scores,
        receiver_order_null_maxima=null_maxima,
        n_windows=np.asarray(n_windows),
        n_files=np.asarray(core_files),
    )
    report = {
        "workflow_version": WORKFLOW_VERSION,
        "citation": CITATION,
        "bensen_citation": BENSEN_CITATION,
        "date": args.date,
        "raw_manifest": str(CSV),
        "source_channel": args.source_channel,
        "source_channel_interpretation": (
            "provisional wellhead from G0 noise transition"
            if args.source_channel == 23
            else "channel-origin sensitivity; not independently the wellhead"
        ),
        "files": core_files,
        "windows_30_s_15_s_step": n_windows,
        "expected_windows_if_contiguous": expected_windows,
        "ram_seconds": args.ram_seconds,
        "ram_duration_status": (
            "project choice following Bensen half-maximum-period guidance; "
            "not specified by Lellouch et al."
        ),
        "ram_samples": metadata_reference["ram_samples"],
        "common_mode_removal": args.common_mode,
        "common_mode_status": "not reported; sensitivity only",
        "spectral_mode": args.spectral_mode,
        "spectral_mode_status": (
            "reported cross-correlation baseline"
            if args.spectral_mode == "cross_correlation"
            else "stabilized Equation-6 sensitivity; estimator not specified in paper"
        ),
        "fk_filter": False,
        "linear_detrend": False,
        "input_bandpass_before_correlation": False,
        "per_window_correlation_normalization": False,
        "published_neighbor_stack": "sum of receiver correlations R-10 through R+10",
        "final_output_band_hz": list(OUTPUT_BAND_HZ),
        "best_causal_velocity_m_s": float(VELOCITY_GRID_M_S[peak_index]),
        "best_causal_score": observed_peak,
        "causal_score_at_3200": float(causal_scores[reference_index]),
        "acausal_score_at_3200": float(acausal_scores[reference_index]),
        "receiver_order_scan_max_null_p": p_value,
        "receiver_order_null_count": args.null_count,
        "opinion_audit": {
            "common_mode_removal": (
                "not a published ambient step; must not be in the exact baseline"
            ),
            "equation_6_power_division": (
                "theory is published, but estimator/stabilization/application are not; "
                "tested separately, not called exact"
            ),
            "r_plus_minus_10_stack": (
                "explicitly published and described as required"
            ),
            "twenty_one_x_snr": (
                "incorrect as a general SNR statement; ideal independent-noise "
                "amplitude SNR gain is sqrt(21), and correlated noise gives less"
            ),
        },
    }
    stem.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")

    figure, axes = plt.subplots(1, 3, figsize=(15.5, 6.3), constrained_layout=True)
    plot_wiggles(
        axes[0], central_filtered, lags, offsets,
        "a  Single center receivers\n(no R±10 enhancement)"
    )
    plot_wiggles(
        axes[1], neighbors_filtered, lags, offsets,
        "b  Published Figure 7c operator\n(simple R±10 sum)"
    )
    axes[2].plot(VELOCITY_GRID_M_S / 1000.0, causal_scores, color="tab:blue", label="causal")
    axes[2].plot(VELOCITY_GRID_M_S / 1000.0, acausal_scores, color="tab:orange", label="acausal")
    axes[2].axvline(3.2, color="black", linestyle="--", linewidth=1.2, label="3.2 km s$^{-1}$")
    axes[2].axhline(np.percentile(null_maxima, 95.0), color="0.45", linestyle=":", label="row-null 95% max")
    axes[2].set(
        xlabel="Trial apparent velocity (km s$^{-1}$)",
        ylabel="Median normalized-envelope gate score",
        title="c  Predeclared moveout scan\n(selection repeated in null)",
    )
    axes[2].grid(alpha=0.25)
    axes[2].legend(frameon=False, fontsize=8)
    figure.suptitle(
        f"Lellouch et al. (2019) Figure 7c reproduction — {args.date}; "
        f"source ch {args.source_channel}; {core_files} files; {n_windows} windows",
        fontsize=12,
    )
    figure.savefig(stem.with_suffix(".png"), dpi=350, bbox_inches="tight")
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(report, indent=2))


def synthetic_validation(args: argparse.Namespace) -> None:
    """Validate lag sign, R±10 linearity, window ownership, and white noise."""
    fs = 100.0
    dx = 1.0
    source_channel = 0
    n_channels = 180
    geom = geometry(source_channel, dx, n_channels)
    duration = 120.0
    time = np.arange(int(duration * fs)) / fs
    velocity = 500.0
    frequency = 8.0
    raw = np.empty((len(geom["required_channels"]), time.size), dtype=np.float32)
    for local, channel in enumerate(geom["required_channels"]):
        delay = (channel - source_channel) * dx / velocity
        raw[local] = np.sin(2.0 * np.pi * frequency * (time - delay))
    normalized, _ = strain_rate_and_ram_continuous(raw, fs, 0.1, False)
    starts = np.arange(
        0,
        normalized.shape[1] - int(WINDOW_SECONDS * fs) + 1,
        int(STEP_SECONDS * fs),
    )
    central, neighbors, n_fft = cross_spectra_for_starts(
        normalized, starts, fs, geom, "cross_correlation", 1.0e-3
    )
    lags, section = correlation_from_spectrum(neighbors, len(starts), n_fft, fs)
    expected = geom["offsets_m"] / velocity
    picked = []
    for trace, target in zip(section, expected):
        gate = np.abs(lags - target) <= 0.03
        picked.append(float(lags[gate][np.argmax(trace[gate])]))
    lag_error = float(np.median(np.abs(np.asarray(picked) - expected)))

    # Correlation is linear: correlating S with the sum of 21 receivers must
    # equal the sum of the 21 individual S-Z correlations used in the paper.
    window_samples = int(WINDOW_SECONDS * fs)
    one_source = normalized[int(geom["source_local"]), :window_samples]
    one_neighbors = normalized[geom["neighbor_local"][0], :window_samples]
    one_source_spectrum = np.fft.rfft(one_source, n=n_fft)
    individual_spectra = np.fft.rfft(one_neighbors, n=n_fft, axis=1)
    direct_sum = np.sum(
        np.conj(one_source_spectrum)[None, :] * individual_spectra, axis=0
    )
    linear_sum = np.conj(one_source_spectrum) * np.fft.rfft(
        np.sum(one_neighbors, axis=0), n=n_fft
    )
    linearity_error = float(
        np.linalg.norm(direct_sum - linear_sum)
        / max(np.linalg.norm(direct_sum), np.finfo(float).eps)
    )

    expected_day_windows = int((86400.0 - WINDOW_SECONDS) / STEP_SECONDS) + 1
    assigned = []
    record_samples = int(60.0 * fs)
    total_samples = int(86400.0 * fs)
    step_samples = int(STEP_SECONDS * fs)
    win_samples = int(WINDOW_SECONDS * fs)
    for core_first in range(0, 1440, 60):
        core_stop = min(core_first + 60, 1440)
        starts_global = np.arange(
            core_first * record_samples, core_stop * record_samples, step_samples
        )
        assigned.extend(starts_global[starts_global + win_samples <= total_samples])
    assigned = np.asarray(assigned)
    if len(assigned) != expected_day_windows or len(np.unique(assigned)) != len(assigned):
        raise AssertionError("Window ownership does not produce 5,759 unique windows")

    rng = np.random.default_rng(12345)
    white = rng.standard_normal((2, 12000), dtype=np.float32)
    frequencies = np.fft.rfftfreq(white.shape[1], 1.0 / fs)
    power = np.abs(np.fft.rfft(white, axis=1)) ** 2
    mean_power = np.mean(power, axis=0)
    band = (frequencies >= 1.0) & (frequencies <= 0.9 * fs / 2.0)
    log_power_slope = float(
        np.polyfit(frequencies[band], np.log(mean_power[band]), 1)[0]
    )
    white_corr = float(np.corrcoef(white)[0, 1])
    white_n_fft = 1 << int(math.ceil(math.log2(2 * white.shape[1] - 1)))
    white_cross = np.fft.irfft(
        np.conj(np.fft.rfft(white[0], n=white_n_fft))
        * np.fft.rfft(white[1], n=white_n_fft),
        n=white_n_fft,
    )
    white_max_lag = int(fs)
    white_cross = np.concatenate(
        (white_cross[-white_max_lag:], white_cross[: white_max_lag + 1])
    )
    white_cross /= math.sqrt(float(np.sum(white[0] ** 2) * np.sum(white[1] ** 2)))
    white_lags = np.arange(-white_max_lag, white_max_lag + 1) / fs

    report = {
        "synthetic_velocity_m_s": velocity,
        "median_absolute_lag_error_s": lag_error,
        "lag_sign_pass": bool(lag_error <= 1.5 / fs),
        "r_plus_minus_10_linearity_relative_error": linearity_error,
        "r_plus_minus_10_linearity_pass": bool(linearity_error < 1.0e-6),
        "expected_full_day_windows": expected_day_windows,
        "assigned_unique_windows": int(len(assigned)),
        "window_ownership_pass": bool(len(assigned) == 5759),
        "white_noise_raw_cross_correlation": white_corr,
        "white_noise_independence_pass": bool(abs(white_corr) < 0.03),
        "white_noise_log_psd_slope_per_hz": log_power_slope,
        "white_noise_broadband_pass": bool(abs(log_power_slope) < 0.02),
    }
    required_passes = [
        "lag_sign_pass",
        "r_plus_minus_10_linearity_pass",
        "window_ownership_pass",
        "white_noise_independence_pass",
        "white_noise_broadband_pass",
    ]
    if not all(report[name] for name in required_passes):
        raise AssertionError(report)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "synthetic_validation.json"
    output.write_text(json.dumps(report, indent=2) + "\n")

    figure, axes = plt.subplots(1, 3, figsize=(14.2, 4.5), constrained_layout=True)
    axes[0].plot(frequencies[1:], mean_power[1:] / np.median(mean_power[band]), color="black", linewidth=0.55)
    axes[0].axvspan(1.0, 0.9 * fs / 2.0, color="tab:blue", alpha=0.08)
    axes[0].set(
        xlabel="Frequency (Hz)",
        ylabel="Power / in-band median",
        title="a  Raw iid control is broadband",
        yscale="log",
        xlim=(0.0, fs / 2.0),
    )
    axes[1].plot(white_lags, white_cross, color="black", linewidth=0.8)
    axes[1].axvline(0.0, color="0.6", linestyle=":")
    axes[1].set(
        xlabel="Lag (s)",
        ylabel="Normalized cross-correlation",
        title=f"b  Independent channels: r = {white_corr:.4f}",
    )
    plot_wiggles(
        axes[2], section, lags, geom["offsets_m"],
        "c  Known synthetic moveout recovered", reference_velocity=velocity
    )
    figure.savefig(args.output_dir / "synthetic_validation.png", dpi=350, bbox_inches="tight")
    figure.savefig(args.output_dir / "synthetic_validation.pdf", bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(report, indent=2))


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--action", choices=("chunk", "aggregate", "synthetic"), default="chunk")
    result.add_argument("--date", default="2024-12-20")
    result.add_argument("--start", type=int, default=0)
    result.add_argument("--nfiles", type=int, default=60)
    result.add_argument("--total-files", type=int, default=1440)
    result.add_argument("--block-files", type=int, default=10)
    result.add_argument("--source-channel", type=int, default=23)
    result.add_argument("--ram-seconds", type=float, default=0.1)
    result.add_argument("--spectral-mode", choices=VALID_SPECTRAL_MODES, default="cross_correlation")
    result.add_argument("--waterlevel", type=float, default=1.0e-3)
    result.add_argument("--common-mode", action="store_true")
    result.add_argument("--null-method", choices=VALID_NULLS, default="ordered")
    result.add_argument("--realization", type=int, default=0)
    result.add_argument("--seed", type=int, default=20260814)
    result.add_argument("--null-count", type=int, default=5000)
    result.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    return result


def main() -> None:
    args = parser().parse_args()
    if args.action == "chunk":
        run_chunk(args)
    elif args.action == "aggregate":
        aggregate(args)
    else:
        synthetic_validation(args)


if __name__ == "__main__":
    main()
