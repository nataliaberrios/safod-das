#!/usr/bin/env python3
"""Matched Lellouch-style ambient comparison with and without F-K filtering.

Every path uses the same contiguous channels, strain-to-strain-rate conversion,
5-s running-absolute-mean normalization, 30-s windows with 15-s overlap,
channel-0 virtual source, 50-m receivers, R +/- 10 receiver stack, and final
5--20 Hz correlation filter.  The only path-specific operation is the signed
F-K mask applied to the normalized strain-rate wavefield before correlation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from ambient_lellouch2019_reproduction_v1 import (
    MAX_LAG_SECONDS,
    STEP_SECONDS,
    WINDOW_SECONDS,
    acquisition_metadata,
    align_and_stack,
    correlation_window,
    day_rows,
    final_bandpass,
    geometry,
    release_allocator_memory,
    ridge_metrics,
    strain_rate_and_ram,
)
from ambient_transfer_test import corrected_path


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "ambient_transfer" / "lellouch_fk_matched_v2"
MODE_SPECS = {
    "unfiltered": {"branch": None, "velocity_m_s": None},
    "downgoing_2p5_4p5": {"branch": "negative", "velocity_m_s": (2500.0, 4500.0)},
    "upgoing_2p5_4p5": {"branch": "positive", "velocity_m_s": (2500.0, 4500.0)},
    "downgoing_direction_only": {"branch": "negative", "velocity_m_s": None},
    "upgoing_direction_only": {"branch": "positive", "velocity_m_s": None},
}


def load_contiguous(path: Path, stop_channel: int) -> tuple[np.ndarray, float]:
    with h5py.File(path, "r") as handle:
        dataset = handle["Acquisition/Raw[0]/RawData"]
        fs = float(dataset.attrs.get("OutputDataRate", 500.0))
        return np.asarray(dataset[:, :stop_channel], dtype=np.float32).T, fs


def fk_filter(data: np.ndarray, fs: float, dx: float, spec: dict) -> np.ndarray:
    frequency = np.fft.fftfreq(data.shape[1], 1.0 / fs)
    wavenumber = np.fft.fftfreq(data.shape[0], dx)
    k_grid, f_grid = np.meshgrid(wavenumber, frequency, indexing="ij")
    selected = (
        (np.abs(f_grid) >= 5.0)
        & (np.abs(f_grid) <= 20.0)
        & (np.abs(k_grid) > 0.0)
    )
    if spec["velocity_m_s"] is not None:
        apparent_velocity = np.abs(f_grid) / np.maximum(np.abs(k_grid), 1e-12)
        vmin, vmax = spec["velocity_m_s"]
        selected &= (apparent_velocity >= vmin) & (apparent_velocity <= vmax)
    if spec["branch"] == "negative":
        selected &= f_grid * k_grid < 0.0
    elif spec["branch"] == "positive":
        selected &= f_grid * k_grid > 0.0
    else:
        raise ValueError(spec)
    return np.fft.ifft2(np.fft.fft2(data) * selected).real.astype(np.float32)


def add_windows(data, required, source_index, fs, accumulator, lags):
    window_samples = int(round(WINDOW_SECONDS * fs))
    step_samples = int(round(STEP_SECONDS * fs))
    used = 0
    for start in range(0, data.shape[1] - window_samples + 1, step_samples):
        window = data[required, start:start + window_samples]
        current_lags, current = correlation_window(
            window, source_index, fs, MAX_LAG_SECONDS
        )
        accumulator = current if accumulator is None else accumulator + current
        lags = current_lags if lags is None else lags
        used += 1
    return accumulator, lags, used


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", default="2024-12-20")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=300)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    rows = day_rows(args.date, args.start, args.nfiles)
    first = next((corrected_path(row.file) for row in rows.itertuples(index=False)
                  if corrected_path(row.file).exists()), None)
    if first is None:
        raise RuntimeError("no mounted input files")
    fs, dx, n_channels, _ = acquisition_metadata(first)
    targets, required = geometry(dx, n_channels)
    stop_channel = int(required.max()) + 1
    source_index = int(np.flatnonzero(required == 0)[0])
    sums = {mode: None for mode in MODE_SPECS}
    windows = {mode: 0 for mode in MODE_SPECS}
    lags = None
    used_files = []

    for number, row in enumerate(rows.itertuples(index=False), 1):
        path = corrected_path(row.file)
        if not path.exists():
            continue
        raw, current_fs = load_contiguous(path, stop_channel)
        if current_fs != fs:
            raise ValueError(f"sample-rate change at {path}: {current_fs} != {fs}")
        normalized = strain_rate_and_ram(raw, fs)
        for mode, spec in MODE_SPECS.items():
            wavefield = normalized if mode == "unfiltered" else fk_filter(
                normalized, fs, dx, spec
            )
            sums[mode], lags, count = add_windows(
                wavefield, required, source_index, fs, sums[mode], lags
            )
            windows[mode] += count
        used_files.append(str(path))
        del raw, normalized
        release_allocator_memory()
        if number % 10 == 0:
            print(f"processed {number}/{len(rows)}; usable={len(used_files)}", flush=True)

    if not used_files:
        raise RuntimeError("no usable records")
    distances = targets.astype(float) * dx
    arrays = {
        "lags_s": lags,
        "distances_m": distances,
        "target_channels": targets,
        "required_channels": required,
        "used_files": np.asarray(used_files),
    }
    report = {
        "workflow_version": "ambient_lellouch_fk_matched_v2",
        "date": args.date,
        "start": args.start,
        "requested_files": args.nfiles,
        "used_files": len(used_files),
        "fs_hz": fs,
        "dx_m": dx,
        "modes": {},
        "fixed_processing": (
            "strain-rate proxy; 5-s RAM normalization; 30-s windows; 15-s "
            "overlap; channel-0 source; 50-m receivers; R+/-10 stack; final "
            "5-20 Hz correlation filter"
        ),
        "direction_mapping": (
            "channel coordinate increases downhole; synthetic and AWD checks "
            "map F*K<0 at positive lag to downgoing/increasing-coordinate and "
            "F*K>0 at negative lag to upgoing/decreasing-coordinate"
        ),
    }
    for mode in MODE_SPECS:
        stacked = final_bandpass(sums[mode] / windows[mode], fs)
        isolated, simple, aligned = align_and_stack(
            stacked, lags, required, targets, dx
        )
        arrays[f"{mode}__isolated"] = isolated
        arrays[f"{mode}__simple"] = simple
        arrays[f"{mode}__aligned_3200"] = aligned
        report["modes"][mode] = {
            "windows": windows[mode],
            "simple_metrics": ridge_metrics(simple, lags, distances),
            "aligned_metrics": ridge_metrics(aligned, lags, distances),
        }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.output_dir / (
        f"lellouch_fk_matched_{args.date}_start{args.start}_n{len(used_files)}"
    )
    np.savez_compressed(stem.with_suffix(".npz"), **arrays)
    stem.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")

    modes = list(MODE_SPECS)
    fig, axes = plt.subplots(1, len(modes), figsize=(18, 4.6),
                             constrained_layout=True, sharex=True, sharey=True)
    common = np.concatenate([arrays[f"{mode}__simple"].ravel() for mode in modes])
    limit = float(np.nanpercentile(np.abs(common), 98.5))
    for ax, mode in zip(axes, modes):
        section = arrays[f"{mode}__simple"]
        ax.imshow(
            section, extent=[lags[0], lags[-1], distances[-1], distances[0]],
            aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit,
            interpolation="nearest",
        )
        sign = -1.0 if mode.startswith("upgoing") else 1.0
        ax.plot(sign * distances / 3200.0, distances, "k--", lw=1.0)
        ax.set(title=mode.replace("_", " "), xlabel="Correlation lag (s)")
    axes[0].set_ylabel("Receiver offset (m)")
    fig.suptitle(
        f"Matched Lellouch comparison: {len(used_files)} one-minute files",
        fontsize=13,
    )
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    fig.savefig(stem.with_suffix(".pdf"))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
