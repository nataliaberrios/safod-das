#!/usr/bin/env python3
"""Versioned signed-lag ambient F-K chunk processor.

This workflow preserves the legacy products and writes corrected correlations
to ``ambient_transfer/signed_lag_v2``.  Under the verified project convention,
F*K < 0 is evaluated at positive correlation lag and F*K > 0 at negative lag.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ambient_fk_transfer_test import fk_filter
from ambient_transfer_test import (
    CSV,
    corrected_path,
    load_segment,
    normalized_corr_pairs,
    preprocess,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "ambient_transfer" / "signed_lag_v2"
BRANCH_LAG_SIGN = {"negative": 1.0, "positive": -1.0}
VELOCITIES = np.linspace(1200.0, 6000.0, 193)


def velocity_scores(top, lags, distance, lag_sign):
    scores = np.array([
        np.nanmedian([
            row[np.argmin(np.abs(lags - lag_sign * offset / velocity))]
            for row, offset in zip(top, distance)
        ])
        for velocity in VELOCITIES
    ])
    return scores


def process_chunk(date, start, nfiles, output_dir):
    database = pd.read_csv(CSV, sep=r"\s+")
    database = database[database.nSamples > 0].copy()
    database["time"] = pd.to_datetime(database.startTime, utc=True, errors="coerce")
    rows = database[database.time.dt.strftime("%Y-%m-%d") == date].sort_values("time")
    rows = rows.iloc[start:start + nfiles]

    stacks = {mode: None for mode in BRANCH_LAG_SIGN}
    lags = distance = None
    used_files = []
    for row in rows.itertuples(index=False):
        path = corrected_path(row.file)
        if not path.exists():
            continue
        raw, fs, dx = load_segment(path)
        raw = preprocess(raw, fs, norm_seconds=5.0)
        for mode in BRANCH_LAG_SIGN:
            filtered, filtered_fs, filtered_dx = fk_filter(raw, fs, dx, mode)
            targets = [
                int(round(50.0 * index / filtered_dx))
                for index in range(1, 15)
            ]
            targets = [target for target in targets if target < filtered.shape[0]]
            branch_lags, correlation = normalized_corr_pairs(
                filtered,
                [(0, target) for target in targets],
                filtered_fs,
                max_lag=0.35,
            )
            stacks[mode] = correlation if stacks[mode] is None else stacks[mode] + correlation
            lags = branch_lags
            distance = np.asarray(targets, dtype=float) * filtered_dx
        used_files.append(str(path))
        print(f"processed {len(used_files)}/{len(rows)}", flush=True)

    if not used_files:
        raise RuntimeError(f"No usable files for {date} start={start}")

    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"signed_fk_v2_{date}_start{start}_n{nfiles}"
    arrays = {
        "lags": lags,
        "distance": distance,
        "used_files": np.asarray(used_files),
    }
    report = {
        "workflow_version": "signed_lag_v2",
        "date": date,
        "start": start,
        "requested_files": nfiles,
        "used_files": len(used_files),
        "virtual_source_channel": 0,
        "receiver_offsets_m": distance.tolist(),
        "passband_hz": [5.0, 20.0],
        "velocity_wedge_m_s": [2500.0, 4500.0],
        "normalization": "5-s running absolute mean",
        "correlation": "conj(source_fft)*receiver_fft; negative lags extracted from padded FFT tail",
        "branches": {},
    }
    for mode, lag_sign in BRANCH_LAG_SIGN.items():
        top = stacks[mode] / len(used_files)
        physical = velocity_scores(top, lags, distance, lag_sign)
        leakage = velocity_scores(top, lags, distance, -lag_sign)
        arrays[f"{mode}_top"] = top
        arrays[f"{mode}_physical_scores"] = physical
        arrays[f"{mode}_leakage_scores"] = leakage
        peak = int(np.nanargmax(np.abs(physical)))
        report["branches"][mode] = {
            "physical_lag_sign": int(lag_sign),
            "physical_peak_velocity_m_s": float(VELOCITIES[peak]),
            "physical_peak_signed_score": float(physical[peak]),
            "physical_peak_absolute_score": float(abs(physical[peak])),
            "physical_score_3200": float(physical[np.argmin(abs(VELOCITIES - 3200.0))]),
            "opposite_lag_leakage_3200": float(leakage[np.argmin(abs(VELOCITIES - 3200.0))]),
        }
    arrays["velocities_m_s"] = VELOCITIES
    np.savez_compressed(output_dir / f"{stem}.npz", **arrays)
    (output_dir / f"{stem}.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=10)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    process_chunk(args.date, args.start, args.nfiles, args.output_dir)


if __name__ == "__main__":
    main()
