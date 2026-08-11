#!/usr/bin/env python3
"""Resumable chunk processor for the frozen ambient F-K mask audit.

This audit changes only the F-K support applied after the established
Lellouch-style 5--20 Hz preprocessing.  It preserves the signed-lag-v2
coordinate convention and deliberately treats the post-filter velocity scan
as a *conditional ridge diagnostic*, not as an unbiased velocity estimate.

The production 2.5--4.5 km/s product is reused when an exact corrected-v2
chunk exists.  Alternative masks require a raw rerun because a correlation
stack cannot be re-filtered into a different F-K support after the fact.
"""
from __future__ import annotations

import argparse
import json
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

from ambient_transfer_test import (
    CSV,
    corrected_path,
    load_segment,
    normalized_corr_pairs,
    preprocess,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "ambient_transfer" / "fk_mask_sensitivity_v2"
PRODUCTION_OUT = ROOT / "ambient_transfer" / "signed_lag_v2"
BRANCH_LAG_SIGN = {"negative": 1, "positive": -1}
PASSBAND_HZ = (5.0, 20.0)
DECIMATION = (2, 2)  # space, time; matches signed-lag-v2 production exactly

# Frozen before looking at held-out results.  These masks test whether the
# result survives narrower, broader, and unrestricted directional support.
MASK_SPECS = OrderedDict(
    [
        (
            "production_2p5_4p5",
            {
                "velocity_m_s": [2500.0, 4500.0],
                "role": "predeclared production wedge",
            },
        ),
        (
            "narrow_2p8_3p8",
            {
                "velocity_m_s": [2800.0, 3800.0],
                "role": "narrow sensitivity wedge around the prior target",
            },
        ),
        (
            "broad_2p0_5p5",
            {
                "velocity_m_s": [2000.0, 5500.0],
                "role": "broad sensitivity wedge",
            },
        ),
        (
            "direction_only",
            {
                "velocity_m_s": None,
                "role": "5--20 Hz signed-direction control with no velocity restriction",
            },
        ),
    ]
)


def fk_mask(nchannels, nsamples, fs, dx, branch, spec):
    """Return conjugate-symmetric signed F-K support for one frozen mask."""
    frequency = np.fft.fftfreq(nsamples, 1.0 / fs)
    wavenumber = np.fft.fftfreq(nchannels, dx)
    k_grid, f_grid = np.meshgrid(wavenumber, frequency, indexing="ij")
    abs_f = np.abs(f_grid)
    abs_k = np.abs(k_grid)
    selected = (
        (abs_f >= PASSBAND_HZ[0])
        & (abs_f <= PASSBAND_HZ[1])
        & (abs_k > 0.0)
    )
    velocity_support = spec["velocity_m_s"]
    if velocity_support is not None:
        apparent_velocity = abs_f / np.maximum(abs_k, 1e-12)
        selected &= (
            (apparent_velocity >= velocity_support[0])
            & (apparent_velocity <= velocity_support[1])
        )
    if branch == "negative":
        selected &= f_grid * k_grid < 0.0
    elif branch == "positive":
        selected &= f_grid * k_grid > 0.0
    else:
        raise ValueError(f"Unknown branch {branch!r}")
    return selected


def production_stem(date, start, nfiles):
    return PRODUCTION_OUT / f"signed_fk_v2_{date}_start{start}_n{nfiles}"


def load_reusable_production(date, start, nfiles):
    """Load an exact corrected production chunk, or return ``None``."""
    stem = production_stem(date, start, nfiles)
    json_path = stem.with_suffix(".json")
    npz_path = stem.with_suffix(".npz")
    if not json_path.exists() or not npz_path.exists():
        return None
    metadata = json.loads(json_path.read_text())
    if metadata.get("workflow_version") != "signed_lag_v2":
        return None
    product = np.load(npz_path, allow_pickle=False)
    required = {"lags", "distance", "negative_top", "positive_top", "used_files"}
    if not required.issubset(product.files):
        return None
    return {
        "metadata": metadata,
        "lags": np.asarray(product["lags"]),
        "distance": np.asarray(product["distance"]),
        "used_files": np.asarray(product["used_files"]).astype(str),
        "negative_top": np.asarray(product["negative_top"]),
        "positive_top": np.asarray(product["positive_top"]),
        "source_npz": str(npz_path),
    }


def process_chunk(date, start, nfiles, output_dir, force_raw_production=False):
    database = pd.read_csv(CSV, sep=r"\s+")
    database = database[database.nSamples > 0].copy()
    database["time"] = pd.to_datetime(
        database.startTime, utc=True, errors="coerce"
    )
    rows = database[
        database.time.dt.strftime("%Y-%m-%d") == date
    ].sort_values("time").iloc[start : start + nfiles]

    reusable = None if force_raw_production else load_reusable_production(
        date, start, nfiles
    )
    raw_masks = [
        name for name in MASK_SPECS
        if name != "production_2p5_4p5" or reusable is None
    ]
    sums = {
        (mask_name, branch): None
        for mask_name in raw_masks
        for branch in BRANCH_LAG_SIGN
    }
    lags = distance = None
    raw_used_files = []

    for row_number, row in enumerate(rows.itertuples(index=False), 1):
        path = corrected_path(row.file)
        if not path.exists():
            continue
        raw, fs, dx = load_segment(path)
        normalized = preprocess(raw, fs, norm_seconds=5.0)
        decimated = normalized[:: DECIMATION[0], :: DECIMATION[1]]
        filtered_fs = fs / DECIMATION[1]
        filtered_dx = dx * DECIMATION[0]
        targets = [
            int(round(50.0 * index / filtered_dx))
            for index in range(1, 15)
        ]
        targets = [target for target in targets if target < decimated.shape[0]]
        pairs = [(0, target) for target in targets]
        spectrum = np.fft.fft2(decimated)

        for mask_name in raw_masks:
            spec = MASK_SPECS[mask_name]
            for branch in BRANCH_LAG_SIGN:
                selected = fk_mask(
                    decimated.shape[0],
                    decimated.shape[1],
                    filtered_fs,
                    filtered_dx,
                    branch,
                    spec,
                )
                filtered = np.fft.ifft2(spectrum * selected).real
                branch_lags, correlation = normalized_corr_pairs(
                    filtered,
                    pairs,
                    filtered_fs,
                    max_lag=0.35,
                )
                key = (mask_name, branch)
                sums[key] = (
                    correlation if sums[key] is None else sums[key] + correlation
                )
                lags = branch_lags
                distance = np.asarray(targets, dtype=float) * filtered_dx
        raw_used_files.append(str(path))
        if row_number % 5 == 0 or row_number == len(rows):
            print(
                f"processed {row_number}/{len(rows)}; usable={len(raw_used_files)}",
                flush=True,
            )

    if raw_masks and not raw_used_files:
        raise RuntimeError(f"No usable raw files for {date} start={start}")

    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / f"fk_mask_v2_{date}_start{start}_n{nfiles}"
    arrays = {}
    report = {
        "workflow_version": "ambient_fk_mask_sensitivity_v2",
        "date": date,
        "start": start,
        "requested_files": nfiles,
        "development_date": "2024-12-20",
        "split_role": "development" if date == "2024-12-20" else "held_out",
        "passband_hz": list(PASSBAND_HZ),
        "preprocessing": "linear detrend, 5--20 Hz Butterworth bandpass, 5-s running-absolute-mean normalization",
        "decimation": {
            "space_factor": DECIMATION[0],
            "time_factor": DECIMATION[1],
            "method": "direct stride after temporal bandpass; matches production signed-lag-v2",
            "limitation": "spatial stride decimation is retained for exact production comparability; interpret alongside the separate anti-alias audit",
        },
        "correlation": "conj(source_fft)*receiver_fft with corrected negative-lag tail extraction",
        "branch_lag_sign": BRANCH_LAG_SIGN,
        "virtual_source_channel": 0,
        "masks": MASK_SPECS,
        "scientific_guardrail": "Any velocity scan is conditional on the selected F-K support and is not an unbiased formation-velocity estimate.",
        "raw_rerun_masks": raw_masks,
        "mask_products": {},
    }

    if reusable is not None:
        if lags is None:
            lags = reusable["lags"]
            distance = reusable["distance"]
        elif not (
            np.array_equal(lags, reusable["lags"])
            and np.allclose(distance, reusable["distance"])
        ):
            raise RuntimeError("Reusable production axes differ from raw-rerun axes")
        for branch in BRANCH_LAG_SIGN:
            arrays[f"production_2p5_4p5__{branch}_top"] = reusable[
                f"{branch}_top"
            ]
        report["mask_products"]["production_2p5_4p5"] = {
            "method": "reused exact corrected signed-lag-v2 chunk",
            "used_files": int(reusable["metadata"]["used_files"]),
            "source_npz": reusable["source_npz"],
        }
        arrays["production_2p5_4p5__used_files"] = reusable["used_files"]

    for mask_name in raw_masks:
        for branch in BRANCH_LAG_SIGN:
            arrays[f"{mask_name}__{branch}_top"] = (
                sums[(mask_name, branch)] / len(raw_used_files)
            )
        arrays[f"{mask_name}__used_files"] = np.asarray(raw_used_files)
        report["mask_products"][mask_name] = {
            "method": "raw rerun with frozen mask",
            "used_files": len(raw_used_files),
        }

    arrays["lags"] = lags
    arrays["distance"] = distance
    np.savez_compressed(stem.with_suffix(".npz"), **arrays)
    stem.with_suffix(".json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=10)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--force-raw-production",
        action="store_true",
        help="Recompute the production mask instead of reusing an exact v2 chunk.",
    )
    args = parser.parse_args()
    process_chunk(
        args.date,
        args.start,
        args.nfiles,
        args.output_dir,
        force_raw_production=args.force_raw_production,
    )


if __name__ == "__main__":
    main()
