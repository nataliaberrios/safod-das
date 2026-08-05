#!/usr/bin/env python3
"""Compare direct, anti-aliased, and full-resolution ambient F-K processing.

This is a preprocessing-sensitivity test, not a new signal search.  Every
path uses the same Lellouch-style geometry, temporal normalization, frequency
band, signed velocity wedge, and receiver-permutation null.  The only change
is whether the 2x time/space reduction is direct slicing, anti-aliased
polyphase resampling, or omitted.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import resample_poly

from ambient_transfer_test import CSV, OUT, corrected_path, load_segment, preprocess, normalized_corr_pairs


def fk_mask(nz: int, nt: int, fs: float, dx: float, mode: str):
    f = np.fft.fftfreq(nt, 1.0 / fs)
    k = np.fft.fftfreq(nz, dx)
    K, F = np.meshgrid(k, f, indexing="ij")
    af, ak = np.abs(F), np.abs(K)
    v = af / np.maximum(ak, 1e-12)
    mask = (af >= 5.0) & (af <= 20.0) & (v >= 2500.0) & (v <= 4500.0) & (ak > 0)
    if mode == "negative":
        mask &= F * K < 0
    elif mode == "positive":
        mask &= F * K > 0
    else:
        raise ValueError(mode)
    return mask


def prepare(x: np.ndarray, fs: float, dx: float, path: str):
    """Return the array and sampling for one preprocessing path."""
    if path == "direct":
        return x[::2, ::2], fs / 2.0, dx * 2.0
    if path == "antialiased":
        # Polyphase resampling applies a low-pass filter before each factor-2
        # reduction, including along the spatial/channel axis.
        y = resample_poly(x, up=1, down=2, axis=0, window=("kaiser", 8.6))
        y = resample_poly(y, up=1, down=2, axis=1, window=("kaiser", 8.6))
        return y.astype(np.float32, copy=False), fs / 2.0, dx * 2.0
    if path == "full_resolution":
        return x, fs, dx
    raise ValueError(path)


def score_stack(top: np.ndarray, lags: np.ndarray, distances: np.ndarray, seed: int = 20260804):
    velocities = np.linspace(1200.0, 6000.0, 193)
    scores = np.array([
        np.nanmedian([row[np.argmin(np.abs(lags - distance / velocity))]
                      for row, distance in zip(top, distances)])
        for velocity in velocities
    ])
    rng = np.random.default_rng(seed)
    null = np.array([
        np.nanmax([
            np.nanmedian([row[np.argmin(np.abs(lags - distance / velocity))]
                          for row, distance in zip(top, rng.permutation(distances))])
            for velocity in velocities
        ])
        for _ in range(500)
    ])
    peak = int(np.nanargmax(scores))
    return {
        "peak_velocity_m_s": float(velocities[peak]),
        "peak_score": float(scores[peak]),
        "score_3200": float(scores[np.argmin(np.abs(velocities - 3200.0))]),
        "null95": float(np.quantile(null, 0.95)),
        "p_peak": float((1.0 + np.sum(null >= scores[peak])) / (len(null) + 1.0)),
        "velocities_m_s": velocities,
        "scores": scores,
        "null": null,
    }


def run(args):
    db = pd.read_csv(CSV, sep=r"\s+")
    db = db[db.nSamples > 0].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    day = db[db.t.dt.strftime("%Y-%m-%d") == args.date].sort_values("t")
    day = day.iloc[args.start:args.start + args.nfiles]
    paths = ["direct", "antialiased", "full_resolution"]
    stacks = {path: None for path in paths}
    lags_by_path = {}
    distances_by_path = {}
    used = []

    for row in day.itertuples(index=False):
        source = corrected_path(row.file)
        if not source.exists():
            continue
        x, fs, dx = load_segment(source)
        x = preprocess(x, fs, norm_seconds=5.0)
        for path in paths:
            y, fs_y, dx_y = prepare(x, fs, dx, path)
            targets = [int(round(50.0 * j / dx_y)) for j in range(1, 15)]
            targets = [channel for channel in targets if channel < y.shape[0]]
            lags, corr = normalized_corr_pairs(y, [(0, channel) for channel in targets], fs_y)
            distances = np.asarray(targets, dtype=float) * dx_y
            stacks[path] = corr if stacks[path] is None else stacks[path] + corr
            lags_by_path[path] = lags
            distances_by_path[path] = distances
        used.append(str(source))
        print(f"processed {len(used)}/{len(day)}", flush=True)

    if not used:
        raise RuntimeError("No usable files")

    results = {}
    arrays = {"used_files": np.asarray(used)}
    for path in paths:
        top = stacks[path] / len(used)
        result = score_stack(top, lags_by_path[path], distances_by_path[path])
        results[path] = {key: value for key, value in result.items()
                         if key not in {"velocities_m_s", "scores", "null"}}
        arrays[f"{path}_top"] = top
        arrays[f"{path}_lags"] = lags_by_path[path]
        arrays[f"{path}_distances"] = distances_by_path[path]
        arrays[f"{path}_velocities_m_s"] = result["velocities_m_s"]
        arrays[f"{path}_scores"] = result["scores"]
        arrays[f"{path}_null"] = result["null"]

    OUT.mkdir(exist_ok=True)
    stem = f"alias_sensitivity_{args.date}_start{args.start}_n{len(used)}"
    np.savez_compressed(OUT / f"{stem}.npz", **arrays)
    report = {
        "date": args.date,
        "start": args.start,
        "requested_files": args.nfiles,
        "used_files": len(used),
        "paths": paths,
        "temporal_band_hz": [5.0, 20.0],
        "velocity_wedge_m_s": [2500.0, 4500.0],
        "resampling": {
            "direct": "x[::2, ::2]",
            "antialiased": "scipy.signal.resample_poly down=2 along space and time",
            "full_resolution": "no decimation",
        },
        "results": results,
    }
    (OUT / f"{stem}.json").write_text(json.dumps(report, indent=2))

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for col, path in enumerate(paths):
        top = arrays[f"{path}_top"]
        lags = arrays[f"{path}_lags"]
        distances = arrays[f"{path}_distances"]
        scale = np.percentile(np.abs(top), 98)
        axes[0, col].imshow(top, extent=[lags[0], lags[-1], distances[-1], distances[0]],
                            aspect="auto", cmap="RdBu_r", vmin=-scale, vmax=scale)
        axes[0, col].set_title(f"{path}: p={results[path]['p_peak']:.3f}")
        axes[0, col].set_xlabel("Lag (s)")
        axes[0, col].set_ylabel("Receiver position (m)")
        axes[1, col].plot(arrays[f"{path}_velocities_m_s"] / 1000.0, arrays[f"{path}_scores"], label="data")
        axes[1, col].axhline(results[path]["null95"], color="k", ls="--", lw=0.8, label="null 95%")
        axes[1, col].axvline(3.2, color="tab:red", ls=":", lw=0.9)
        axes[1, col].set_xlabel("Trial velocity (km/s)")
        axes[1, col].set_ylabel("Median normalized correlation")
        axes[1, col].legend(frameon=False, fontsize=8)
    fig.suptitle("F-K preprocessing sensitivity: direct vs anti-aliased vs full resolution")
    fig.savefig(OUT / f"{stem}.png", dpi=250)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", default="2024-12-20")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=30)
    run(parser.parse_args())
