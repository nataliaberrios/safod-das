"""Publication-quality multiscale AWD record sections for Nano and Deep DAS.

The figure compares one deterministically selected same-drop observation, its
paired-drop burst stack, and the drop-count-weighted full stack.  The selected
drop is the temporally central paired drop in the burst containing the largest
number of common Nano--Deep drops; it is not selected by visual signal quality.

All time axes are stated literally relative to the ``UTC_Date`` value in
``p26.cc9.txt``.  That value is used only as an alignment timestamp here; this
script does not call it a source time, impact time, catalog time, or phase pick.
"""

from __future__ import annotations

import csv
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, detrend, sosfiltfilt


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
MANIFEST = HERE / "awd_manifest.csv"
NANO_DIR = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/Nano")
DEEP_DIR = Path(
    "/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/"
    "01_--_recording_2026-06-15T230629Z_--_active_source"
)

UTILITY_DIRS = (
    Path("/home/groups/edunham/nberrios/safod_das/DAS-utilities/python"),
    Path("/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python"),
)
for directory in UTILITY_DIRS:
    if directory.exists() and str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from DASutils import readFile_HDF, readFile_protobuf  # noqa: E402


PRE_S = 0.5
POST_S = 3.0
FS_EXPECTED = 1000.0
NANO_BAND = (30.0, 60.0)
DEEP_BAND = (3.0, 15.0)
NANO_APERTURE_M = (80.0, 500.0)
NANO_TIME_S = (-0.08, 0.30)
DEEP_TIME_S = (-0.08, 2.85)
DEEP_TURNAROUND_CH = 1702
NOISE_WINDOW_S = (-0.45, -0.15)
DISPLAY_LIMIT_NOISE_RMS = 12.0
NANO_DISPLAY_STRIDE = 2
DEEP_DISPLAY_STRIDE = 4


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    paired = [row for row in rows if int(row["paired_available"]) == 1]
    if not paired:
        raise RuntimeError(f"No paired drops found in {MANIFEST}")
    return paired


def select_representative_drop(
    rows: list[dict[str, str]], counts: np.ndarray
) -> tuple[int, dict[str, str], list[dict[str, str]]]:
    """Select without looking at DAS amplitude or coherence."""
    valid = np.flatnonzero(counts > 0)
    if valid.size == 0:
        raise RuntimeError("Canonical product contains no paired bursts")
    # Nominal manifest pairing can exceed complete canonical extraction near
    # raw-file boundaries. Select the largest burst whose nominal paired rows
    # exactly match the canonical complete-window count, without viewing data.
    for candidate in valid[np.argsort(counts[valid])[::-1]]:
        burst_id = int(candidate)
        burst_rows = sorted(
            (row for row in rows if int(row["burst_id"]) == burst_id),
            key=lambda row: int(row["drop_id"]),
        )
        if len(burst_rows) == int(counts[burst_id]):
            selected = burst_rows[len(burst_rows) // 2]
            return burst_id, selected, burst_rows
    raise RuntimeError("No burst has matching manifest and canonical complete-window counts")


def parse_nano_start(name: str) -> datetime:
    fields = os.path.basename(name).split("_")
    value = f"{fields[1]}T{fields[2].replace('.', ':')}Z"
    return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)


def parse_deep_start(name: str) -> datetime:
    match = re.search(r"_(\d{4}-\d{2}-\d{2}T\d{6}Z)\.h5$", os.path.basename(name))
    if not match:
        raise ValueError(f"Cannot parse Deep start time from {name}")
    return datetime.strptime(match.group(1), "%Y-%m-%dT%H%M%SZ").replace(
        tzinfo=timezone.utc
    )


def extract_window(data: np.ndarray, fs: float, offset_s: float) -> np.ndarray:
    i0 = int(round((offset_s - PRE_S) * fs))
    i1 = int(round((offset_s + POST_S) * fs))
    if i0 < 0 or i1 > data.shape[1]:
        raise RuntimeError(
            f"Requested [{i0}:{i1}] outside raw sample interval [0:{data.shape[1]}]"
        )
    return np.asarray(data[:, i0:i1], dtype=np.float32)


def load_individual(row: dict[str, str]) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Read the same paired drop using the canonical reader settings."""
    nano_path = NANO_DIR / row["nano_file"]
    deep_path = DEEP_DIR / row["deep_file"]
    if not nano_path.exists() or not deep_path.exists():
        raise FileNotFoundError(f"Missing raw file: {nano_path if not nano_path.exists() else deep_path}")

    nano, nano_info = readFile_protobuf(
        [str(nano_path)], fmin=1.0, fmax=100.0, desampling=False
    )
    deep, deep_info = readFile_HDF(
        [str(deep_path)], fmin=1.0, fmax=100.0,
        desampling=False, verbose=False,
    )
    fs_nano = float(nano_info["fs"])
    fs_deep = float(deep_info["fs"])
    if not np.isclose(fs_nano, fs_deep) or not np.isclose(fs_nano, FS_EXPECTED):
        raise RuntimeError(f"Unexpected sampling rates: Nano={fs_nano}, Deep={fs_deep}")

    utc = datetime.fromisoformat(row["utc_time"])
    # Recompute offsets from literal UTC values as an independent manifest check.
    nano_offset = (utc - parse_nano_start(row["nano_file"])).total_seconds()
    deep_offset = (utc - parse_deep_start(row["deep_file"])).total_seconds()
    if not np.isclose(nano_offset, float(row["nano_offset_s"]), atol=5e-4):
        raise RuntimeError("Nano offset does not reproduce the manifest")
    if not np.isclose(deep_offset, float(row["deep_offset_s"]), atol=5e-4):
        raise RuntimeError("Deep offset does not reproduce the manifest")

    return (
        extract_window(nano, fs_nano, nano_offset),
        extract_window(deep, fs_deep, deep_offset),
        fs_nano,
        float(nano_info["dx"]),
    )


def finite_weighted_stack(stacks: np.ndarray, counts: np.ndarray) -> np.ndarray:
    good = counts > 0
    weights = counts[good].astype(np.float64)[:, None, None]
    values = np.asarray(stacks[good], dtype=np.float64)
    finite = np.isfinite(values)
    numerator = np.sum(np.where(finite, values, 0.0) * weights, axis=0)
    denominator = np.sum(finite * weights, axis=0)
    return np.divide(
        numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0
    ).astype(np.float32)


def bandpass(section: np.ndarray, fs: float, band: tuple[float, float]) -> np.ndarray:
    data = np.nan_to_num(np.asarray(section, dtype=np.float64))
    data = detrend(data, axis=-1, type="linear")
    sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, data, axis=-1)


def noise_standardize(section: np.ndarray, time: np.ndarray) -> np.ndarray:
    """Express traces in pre-reference noise-RMS units, preserving relative signal."""
    noise = (time >= NOISE_WINDOW_S[0]) & (time <= NOISE_WINDOW_S[1])
    if noise.sum() < 20:
        raise RuntimeError("Noise window is outside the extracted record")
    rms = np.sqrt(np.mean(section[:, noise] ** 2, axis=1))
    positive = rms[np.isfinite(rms) & (rms > 0)]
    if positive.size == 0:
        raise RuntimeError("No finite pre-reference noise RMS values")
    floor = np.nanmedian(positive) * 0.05
    scale = np.where(np.isfinite(rms) & (rms > floor), rms, np.nan)
    return np.divide(
        section, scale[:, None], out=np.zeros_like(section),
        where=np.isfinite(scale[:, None]),
    )


def block_mean(section: np.ndarray, coordinate: np.ndarray, stride: int):
    n = (section.shape[0] // stride) * stride
    if n < stride:
        return section, coordinate
    reduced = np.nanmean(section[:n].reshape(-1, stride, section.shape[1]), axis=1)
    coord = np.mean(coordinate[:n].reshape(-1, stride), axis=1)
    return reduced, coord


def prepared_panel(
    section: np.ndarray,
    coordinate: np.ndarray,
    time: np.ndarray,
    fs: float,
    band: tuple[float, float],
    aperture: tuple[float, float],
    time_window: tuple[float, float],
    stride: int,
):
    filtered = noise_standardize(bandpass(section, fs, band), time)
    channels = (coordinate >= aperture[0]) & (coordinate <= aperture[1])
    samples = (time >= time_window[0]) & (time <= time_window[1])
    display, distance = block_mean(filtered[channels][:, samples], coordinate[channels], stride)
    return display.astype(np.float32), distance, time[samples]


def draw_row(
    axes, panels, distances, times, titles, band, installation,
    ridge: tuple[float, float, float] | None = None,
):
    image = None
    for column, (axis, panel, distance, time, title) in enumerate(
        zip(axes, panels, distances, times, titles)
    ):
        image = axis.imshow(
            panel,
            aspect="auto",
            origin="upper",
            cmap="RdBu_r",
            interpolation="nearest",
            extent=[time[0], time[-1], distance[-1], distance[0]],
            vmin=-DISPLAY_LIMIT_NOISE_RMS,
            vmax=DISPLAY_LIMIT_NOISE_RMS,
            rasterized=True,
        )
        if ridge is not None:
            x0, intercept, velocity = ridge
            predicted = intercept + (distance - x0) / velocity
            axis.plot(predicted, distance, color="k", ls="--", lw=1.1)
        axis.axvline(0, color="#E69F00", lw=0.9)
        axis.set_title(f"{chr(65 + column)}  {title}" if installation == "Nano" else
                       f"{chr(68 + column)}  {title}", loc="left", fontsize=10)
        axis.set_xlabel("Time relative to p26.cc9.txt UTC_Date (s)")
        if column == 0:
            label = "Distance along Nano fiber (m)" if installation == "Nano" else (
                "Distance along Deep outbound leg from channel 0 (m)"
            )
            axis.set_ylabel(label)
        axis.text(
            0.99, 0.02, f"{installation}, {band[0]:g}–{band[1]:g} Hz",
            transform=axis.transAxes, ha="right", va="bottom", fontsize=8,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.76, pad=1.8),
        )
    return image


def main() -> None:
    rows = read_manifest()
    with np.load(STACKS) as product:
        nano_stacks = product["nano_stacks"]
        deep_stacks = product["deep_stacks"]
        counts = product["n_common"].astype(int)
        fs = float(product["fs"])
        dx_nano = float(product["dx_nano"])
        dx_deep = float(product["dx_deep"])
        burst_id, selected, burst_rows = select_representative_drop(rows, counts)
        nano_burst = np.array(nano_stacks[burst_id], copy=True)
        deep_burst = np.array(deep_stacks[burst_id], copy=True)
        nano_full = finite_weighted_stack(nano_stacks, counts)
        deep_full = finite_weighted_stack(deep_stacks, counts)

    nano_drop, deep_drop, raw_fs, raw_dx_nano = load_individual(selected)
    if not np.isclose(raw_fs, fs) or not np.isclose(raw_dx_nano, dx_nano):
        raise RuntimeError("Raw and canonical Nano sampling metadata disagree")
    expected_samples = int(round((PRE_S + POST_S) * fs))
    for name, section in (
        ("nano_drop", nano_drop), ("deep_drop", deep_drop),
        ("nano_burst", nano_burst), ("deep_burst", deep_burst),
        ("nano_full", nano_full), ("deep_full", deep_full),
    ):
        if section.shape[1] != expected_samples:
            raise RuntimeError(f"{name} has {section.shape[1]} rather than {expected_samples} samples")

    time = np.arange(expected_samples) / fs - PRE_S
    nano_coordinate = np.arange(nano_drop.shape[0]) * dx_nano
    deep_coordinate = np.arange(deep_drop.shape[0]) * dx_deep
    # The outbound leg is plotted because it contains the strongest validated
    # positive-slowness slow mode.  No conversion to physical borehole depth is made.
    deep_stop = min(DEEP_TURNAROUND_CH, deep_drop.shape[0], deep_burst.shape[0], deep_full.shape[0])
    deep_coordinate = deep_coordinate[:deep_stop]

    nano_inputs = (nano_drop, nano_burst, nano_full)
    deep_inputs = (deep_drop[:deep_stop], deep_burst[:deep_stop], deep_full[:deep_stop])
    nano_panels, nano_distances, nano_times = [], [], []
    deep_panels, deep_distances, deep_times = [], [], []
    for section in nano_inputs:
        panel, distance, panel_time = prepared_panel(
            section, nano_coordinate, time, fs, NANO_BAND,
            NANO_APERTURE_M, NANO_TIME_S, NANO_DISPLAY_STRIDE,
        )
        nano_panels.append(panel); nano_distances.append(distance); nano_times.append(panel_time)
    deep_aperture = (0.0, deep_coordinate[-1])
    for section in deep_inputs:
        panel, distance, panel_time = prepared_panel(
            section, deep_coordinate, time, fs, DEEP_BAND,
            deep_aperture, DEEP_TIME_S, DEEP_DISPLAY_STRIDE,
        )
        deep_panels.append(panel); deep_distances.append(distance); deep_times.append(panel_time)

    n_burst = int(counts[burst_id])
    n_full = int(counts.sum())
    utc_text = selected["utc_time"].replace("+00:00", "Z")
    titles = (
        f"Single paired drop\nUTC_Date={utc_text}",
        f"Burst {burst_id} mean ({n_burst} paired drops)",
        f"Full weighted mean ({n_full} paired drops; {np.sum(counts > 0)} bursts)",
    )
    fig, axes = plt.subplots(2, 3, figsize=(16.4, 10.2), constrained_layout=True)
    nano_image = draw_row(
        axes[0], nano_panels, nano_distances, nano_times, titles,
        NANO_BAND, "Nano", ridge=(NANO_APERTURE_M[0], -0.022, 2975.0),
    )
    deep_image = draw_row(
        axes[1], deep_panels, deep_distances, deep_times, titles,
        DEEP_BAND, "Deep", ridge=None,
    )
    nano_bar = fig.colorbar(nano_image, ax=axes[0], shrink=0.88, pad=0.012)
    nano_bar.set_label("Band-passed strain-rate amplitude / pre-reference noise RMS")
    deep_bar = fig.colorbar(deep_image, ax=axes[1], shrink=0.88, pad=0.012)
    deep_bar.set_label("Band-passed strain-rate amplitude / pre-reference noise RMS")
    fig.suptitle(
        "SAFOD AWD multiscale record sections: single drop, burst stack, and full stack",
        fontsize=14,
    )
    fig.savefig(HERE / "awd_multiscale_record_sections.png", dpi=240)
    plt.close(fig)

    np.savez_compressed(
        HERE / "awd_multiscale_record_sections.npz",
        nano_drop=nano_panels[0], nano_burst=nano_panels[1], nano_full=nano_panels[2],
        nano_distance_m=nano_distances[0], nano_time_s=nano_times[0],
        deep_drop=deep_panels[0], deep_burst=deep_panels[1], deep_full=deep_panels[2],
        deep_outbound_distance_m=deep_distances[0], deep_time_s=deep_times[0],
        selected_burst_id=burst_id, selected_drop_id=int(selected["drop_id"]),
        selected_utc_date=selected["utc_time"], selected_nano_file=selected["nano_file"],
        selected_deep_file=selected["deep_file"], n_burst_paired=n_burst,
        n_full_paired=n_full, n_full_bursts=int(np.sum(counts > 0)),
        nano_band_hz=NANO_BAND, deep_band_hz=DEEP_BAND,
        noise_window_s=NOISE_WINDOW_S, display_limit_noise_rms=DISPLAY_LIMIT_NOISE_RMS,
    )

    report = (
        "SAFOD AWD multiscale record sections\n"
        "====================================\n"
        f"Selection rule: temporal center of burst with most paired drops\n"
        f"Selected burst/drop: {burst_id}/{selected['drop_id']}\n"
        f"Selected p26.cc9.txt UTC_Date: {selected['utc_time']}\n"
        f"Selected Nano file: {selected['nano_file']}\n"
        f"Selected Deep file: {selected['deep_file']}\n"
        f"Paired drops in selected burst: {n_burst}\n"
        f"Paired drops in full weighted stack: {n_full}\n"
        f"Nonempty paired bursts: {np.sum(counts > 0)}\n"
        f"Nano display: {NANO_BAND[0]:g}-{NANO_BAND[1]:g} Hz, "
        f"{NANO_APERTURE_M[0]:g}-{NANO_APERTURE_M[1]:g} m, "
        f"{NANO_TIME_S[0]:g}-{NANO_TIME_S[1]:g} s\n"
        f"Deep display: outbound leg only, {DEEP_BAND[0]:g}-{DEEP_BAND[1]:g} Hz, "
        f"0-{deep_coordinate[-1]:.1f} m fiber coordinate, "
        f"{DEEP_TIME_S[0]:g}-{DEEP_TIME_S[1]:g} s\n"
        f"Amplitude scaling: each trace divided by RMS in "
        f"{NOISE_WINDOW_S[0]:g} to {NOISE_WINDOW_S[1]:g} s; fixed display limits "
        f"+/-{DISPLAY_LIMIT_NOISE_RMS:g} noise RMS for every panel\n"
        "Timing definition: all times are relative to the UTC_Date value in "
        "p26.cc9.txt; its physical relation to the AWD impact is not assumed.\n"
        "Deep coordinate definition: distance along the outbound fiber leg from "
        "Deep channel 0, not independently registered physical borehole depth.\n"
    )
    (HERE / "awd_multiscale_record_sections.txt").write_text(report)
    print(report, end="")
    print("Saved awd_multiscale_record_sections.png/.npz/.txt")


if __name__ == "__main__":
    main()
