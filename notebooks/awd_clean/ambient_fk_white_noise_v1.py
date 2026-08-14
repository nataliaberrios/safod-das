#!/usr/bin/env python3
"""Pass deterministic Gaussian white noise through the ambient F-K pipeline.

The synthetic arrays use the production temporal/spatial sampling, channel-0
correlation geometry, preprocessing order, decimation, signed F-K masks, and
velocity scoring.  Multiple independently seeded ensembles provide a direct
null distribution for every frozen mask and both coordinate directions.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ambient_fk_mask_sensitivity_v2 import MASK_SPECS, fk_mask
from ambient_signed_fk_v2 import BRANCH_LAG_SIGN, VELOCITIES, velocity_scores
from ambient_transfer_test import normalized_corr_pairs, preprocess


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "ambient_transfer" / "fk_qc_notebook_v2"
REAL_NPZ = (
    ROOT / "ambient_transfer" / "fk_mask_sensitivity_v2" /
    "ambient_fk_mask_sensitivity_v2.npz"
)
FS_HZ = 500.0
DX_M = 1.0209523439407349
N_CHANNELS = 800
SECONDS_PER_RECORD = 12.0
ENSEMBLES = 10
RECORDS_PER_ENSEMBLE = 3
SEED = 20260813


def correlation_section(x, fs, dx, mask_name, branch):
    decimated = x[::2, ::2]
    filtered_fs = fs / 2.0
    filtered_dx = dx * 2.0
    support = fk_mask(
        decimated.shape[0], decimated.shape[1], filtered_fs, filtered_dx,
        branch, MASK_SPECS[mask_name],
    )
    filtered = np.fft.ifft2(np.fft.fft2(decimated) * support).real
    targets = [
        int(round(50.0 * index / filtered_dx))
        for index in range(1, 15)
    ]
    targets = [target for target in targets if target < filtered.shape[0]]
    lags, top = normalized_corr_pairs(
        filtered, [(0, target) for target in targets], filtered_fs,
        max_lag=0.35,
    )
    distance = np.asarray(targets) * filtered_dx
    return lags, distance, top


def main():
    if not REAL_NPZ.exists():
        raise FileNotFoundError(REAL_NPZ)
    real = np.load(REAL_NPZ, allow_pickle=False)
    rng = np.random.default_rng(SEED)
    branches = tuple(BRANCH_LAG_SIGN)
    masks = tuple(MASK_SPECS)
    ensemble_scores = {
        (mask_name, branch): []
        for mask_name in masks for branch in branches
    }
    final_sections = {}
    lags = distance = None

    for ensemble in range(ENSEMBLES):
        sums = {(mask_name, branch): None for mask_name in masks for branch in branches}
        for _ in range(RECORDS_PER_ENSEMBLE):
            raw = rng.standard_normal(
                (N_CHANNELS, int(FS_HZ * SECONDS_PER_RECORD)),
                dtype=np.float32,
            )
            processed = preprocess(raw, FS_HZ, norm_seconds=5.0)
            for mask_name in masks:
                for branch in branches:
                    lags, distance, top = correlation_section(
                        processed, FS_HZ, DX_M, mask_name, branch
                    )
                    key = (mask_name, branch)
                    sums[key] = top if sums[key] is None else sums[key] + top
        for mask_name in masks:
            for branch in branches:
                key = (mask_name, branch)
                section = sums[key] / RECORDS_PER_ENSEMBLE
                scores = velocity_scores(
                    section, lags, distance, BRANCH_LAG_SIGN[branch]
                )
                ensemble_scores[key].append(scores)
                if ensemble == 0:
                    final_sections[key] = section

    arrays = {
        "lags": lags,
        "distance": distance,
        "velocities_m_s": VELOCITIES,
    }
    report = {
        "workflow_version": "ambient_fk_white_noise_v1",
        "seed": SEED,
        "ensembles": ENSEMBLES,
        "records_per_ensemble": RECORDS_PER_ENSEMBLE,
        "seconds_per_record": SECONDS_PER_RECORD,
        "sampling": {"fs_hz": FS_HZ, "dx_m": DX_M, "channels": N_CHANNELS},
        "preprocessing": (
            "Gaussian samples -> detrend -> 5-20 Hz bandpass -> 5-s "
            "running-absolute-mean normalization -> 2x time/space decimation "
            "-> frozen signed F-K mask -> channel-0 correlations"
        ),
        "masks": {},
    }
    index_3200 = int(np.argmin(abs(VELOCITIES - 3200.0)))
    for mask_name in masks:
        report["masks"][mask_name] = {}
        for branch in branches:
            key = (mask_name, branch)
            null_scores = np.asarray(ensemble_scores[key])
            real_scores = real[f"held_out__{mask_name}__{branch}__scores"]
            observed = float(real_scores[index_3200])
            white_3200 = np.abs(null_scores[:, index_3200])
            white_peak = np.nanmax(np.abs(null_scores), axis=1)
            report["masks"][mask_name][branch] = {
                "real_signed_score_3200": observed,
                "white_abs_score_3200_median": float(np.median(white_3200)),
                "white_abs_score_3200_95": float(np.quantile(white_3200, 0.95)),
                "white_peak_abs_score_median": float(np.median(white_peak)),
                "white_peak_abs_score_95": float(np.quantile(white_peak, 0.95)),
                "real_exceeds_white_95_at_3200": bool(
                    abs(observed) > np.quantile(white_3200, 0.95)
                ),
            }
            arrays[f"{mask_name}__{branch}__ensemble_scores"] = null_scores
            arrays[f"{mask_name}__{branch}__example_top"] = final_sections[key]

    report["decision"] = (
        "A mask passes this synthetic gate only when the held-out real-data "
        "absolute 3.2 km/s score exceeds the 95th percentile of identically "
        "processed white-noise ensembles. This gate is necessary but not "
        "sufficient; the raw-channel permutation gate is evaluated separately."
    )
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "ambient_fk_white_noise_v1.npz", **arrays)
    (OUT / "ambient_fk_white_noise_v1.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )

    fig, axes = plt.subplots(4, 2, figsize=(12, 12.5), constrained_layout=True,
                             sharex=True, sharey=True)
    for row, mask_name in enumerate(masks):
        for column, branch in enumerate(branches):
            ax = axes[row, column]
            top = final_sections[(mask_name, branch)]
            limit = float(np.nanpercentile(np.abs(top), 98.5))
            ax.imshow(
                top,
                extent=[lags[0], lags[-1], distance[-1], distance[0]],
                aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit,
                interpolation="nearest",
            )
            sign = BRANCH_LAG_SIGN[branch]
            ax.plot(sign * distance / 3200.0, distance, "k--", lw=1.0)
            gate = report["masks"][mask_name][branch][
                "real_exceeds_white_95_at_3200"
            ]
            ax.set_title(f"{mask_name}; {branch}; real > white95: {gate}")
            ax.set(xlabel="Correlation lag (s)", ylabel="Receiver offset (m)")
    fig.suptitle(
        "White-noise input after the complete signed F-K/correlation pipeline",
        fontsize=13,
    )
    fig.savefig(OUT / "ambient_fk_white_noise_v1.png", dpi=300)
    fig.savefig(OUT / "ambient_fk_white_noise_v1.pdf")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
