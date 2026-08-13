#!/usr/bin/env python3
"""Validate the ambient F-K fan with a known-direction surface AWD source.

The surface AWD provides an empirical direction label independent of the F-K
filter.  Nano coordinate increases from the top of the cemented main-hole
fiber into the well, so the direct source arrival must move toward increasing
coordinate.  This script applies the two complementary production fans to the
same time-gated AWD stack at full temporal and spatial sampling, reconstructs
both wavefields, and measures leakage into the observed downgoing moveout.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import detrend, windows

from fk_dispersion import weighted_stack


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "canonical_epoch_stacks_paired_deep_all.npz"
OUT_DIR = ROOT / "ambient_transfer" / "awd_production_operator_validation_v1"
OUT = OUT_DIR / "awd_fk_production_operator_validation_v1"

PRE_SOURCE_SECONDS = 0.5
TIME_WINDOW_SECONDS = (-0.08, 0.45)
APERTURE_M = (80.0, 440.0)
FREQUENCY_BANDS_HZ = ((5.0, 20.0), (25.0, 60.0))
VELOCITY_FAN_M_S = (2500.0, 4500.0)
OBSERVED_VELOCITY_M_S = 2975.0
TUBE_HALF_WIDTH_SECONDS = 0.012


def signed_fk_reconstruct(
    section: np.ndarray,
    fs: float,
    dx: float,
    frequency_band: tuple[float, float],
    fk_sign: int,
) -> tuple[np.ndarray, float]:
    """Apply one hard, Hermitian-complete production fan at full resolution."""
    prepared = detrend(section.astype(float), axis=1, type="linear")
    prepared *= windows.tukey(prepared.shape[0], 0.15)[:, None]
    prepared *= windows.tukey(prepared.shape[1], 0.10)[None, :]
    spectrum = np.fft.fft2(prepared)
    frequency = np.fft.fftfreq(prepared.shape[1], 1.0 / fs)
    wavenumber = np.fft.fftfreq(prepared.shape[0], dx)
    k_grid, f_grid = np.meshgrid(wavenumber, frequency, indexing="ij")
    speed = np.abs(f_grid) / np.maximum(np.abs(k_grid), 1e-30)
    mask = (
        (np.abs(f_grid) >= frequency_band[0])
        & (np.abs(f_grid) <= frequency_band[1])
        & (speed >= VELOCITY_FAN_M_S[0])
        & (speed <= VELOCITY_FAN_M_S[1])
        & (np.abs(k_grid) > 0.0)
        & (np.sign(f_grid * k_grid) == fk_sign)
    )
    retained_power = float(np.sum(np.abs(spectrum[mask]) ** 2))
    reconstructed = np.fft.ifft2(spectrum * mask).real
    return reconstructed, retained_power


def fixed_moveout_intercept(
    section: np.ndarray,
    time: np.ndarray,
    coordinate: np.ndarray,
    velocity: float,
) -> float:
    """Estimate one intercept from the unfiltered active gather."""
    trace_scale = np.sqrt(np.mean(section ** 2, axis=1, keepdims=True))
    normalized = section / np.maximum(trace_scale, np.finfo(float).eps)
    candidates = np.arange(-0.08, 0.0801, 0.001)
    scores = []
    relative_coordinate = coordinate - coordinate[0]
    for intercept in candidates:
        samples = []
        for trace, offset in zip(normalized, relative_coordinate):
            predicted = intercept + offset / velocity
            samples.append(trace[int(np.argmin(np.abs(time - predicted)))])
        scores.append(abs(float(np.median(samples))))
    return float(candidates[int(np.argmax(scores))])


def tube_energy(
    section: np.ndarray,
    time: np.ndarray,
    coordinate: np.ndarray,
    intercept: float,
    velocity: float,
) -> float:
    """Return mean squared amplitude in a fixed moveout tube."""
    values = []
    for trace, offset in zip(section, coordinate - coordinate[0]):
        predicted = intercept + offset / velocity
        selected = np.abs(time - predicted) <= TUBE_HALF_WIDTH_SECONDS
        values.append(np.mean(trace[selected] ** 2))
    return float(np.mean(values))


def main() -> None:
    if not INPUT.is_file():
        raise FileNotFoundError(INPUT)
    with np.load(INPUT, allow_pickle=False) as product:
        fs = float(product["fs"])
        dx = float(product["dx_nano"])
        nano = weighted_stack(product["nano_stacks"], product["n_common"])

    channel_start = int(round(APERTURE_M[0] / dx))
    channel_stop = int(round(APERTURE_M[1] / dx)) + 1
    sample_start = int(round((PRE_SOURCE_SECONDS + TIME_WINDOW_SECONDS[0]) * fs))
    sample_stop = int(round((PRE_SOURCE_SECONDS + TIME_WINDOW_SECONDS[1]) * fs))
    section = np.nan_to_num(
        nano[channel_start:channel_stop, sample_start:sample_stop],
        copy=True,
    )
    coordinate = np.arange(channel_start, channel_stop, dtype=float) * dx
    time = np.arange(sample_start, sample_stop, dtype=float) / fs - PRE_SOURCE_SECONDS
    intercept = fixed_moveout_intercept(
        section, time, coordinate, OBSERVED_VELOCITY_M_S
    )

    reconstructed: dict[str, dict[str, np.ndarray]] = {}
    metrics: dict[str, dict[str, float]] = {}
    for fmin, fmax in FREQUENCY_BANDS_HZ:
        band_name = f"{fmin:g}_{fmax:g}_hz"
        expected, expected_power = signed_fk_reconstruct(
            section, fs, dx, (fmin, fmax), fk_sign=-1
        )
        opposite, opposite_power = signed_fk_reconstruct(
            section, fs, dx, (fmin, fmax), fk_sign=1
        )
        expected_tube = tube_energy(
            expected, time, coordinate, intercept, OBSERVED_VELOCITY_M_S
        )
        opposite_tube = tube_energy(
            opposite, time, coordinate, intercept, OBSERVED_VELOCITY_M_S
        )
        reconstructed[band_name] = {
            "expected_fk_negative": expected,
            "opposite_fk_positive": opposite,
        }
        metrics[band_name] = {
            "expected_to_opposite_retained_fourier_power_ratio": float(
                expected_power / opposite_power
            ),
            "expected_to_opposite_energy_in_fixed_downgoing_tube_ratio": float(
                expected_tube / opposite_tube
            ),
            "expected_fraction_of_two_branch_tube_energy": float(
                expected_tube / (expected_tube + opposite_tube)
            ),
        }

    report = {
        "workflow_version": "awd_fk_production_operator_validation_v1",
        "question": (
            "Does the exact signed 2.5--4.5 km/s fan retain the known "
            "increasing-coordinate AWD arrival and reject its complementary branch?"
        ),
        "input": str(INPUT.relative_to(ROOT)),
        "ground_truth": {
            "source": "surface accelerated weight drop",
            "coordinate": "Nano coordinate increases from the top into the well",
            "expected_motion": "toward increasing coordinate (downgoing)",
            "expected_fft_sign": "F*K < 0",
        },
        "operator": {
            "velocity_fan_m_s": list(VELOCITY_FAN_M_S),
            "sampling": "full temporal and spatial sampling; no decimation",
            "mask": "hard complementary Hermitian-complete signed fans",
        },
        "active_gate": {
            "aperture_m": list(APERTURE_M),
            "time_relative_to_catalog_reference_s": list(TIME_WINDOW_SECONDS),
            "independent_observed_velocity_m_s": OBSERVED_VELOCITY_M_S,
            "fixed_moveout_intercept_s": intercept,
            "tube_half_width_s": TUBE_HALF_WIDTH_SECONDS,
        },
        "metrics": metrics,
        "decision": {
            "signed_fan_selects_known_downgoing_direction_at_5_20_hz": bool(
                metrics["5_20_hz"][
                    "expected_to_opposite_energy_in_fixed_downgoing_tube_ratio"
                ]
                > 2.0
            ),
            "signed_fan_selects_known_downgoing_direction_at_25_60_hz": bool(
                metrics["25_60_hz"][
                    "expected_to_opposite_energy_in_fixed_downgoing_tube_ratio"
                ]
                > 10.0
            ),
            "opposite_ambient_branch_can_be_assigned_to_sign_bug": False,
            "interpretation": (
                "The known downgoing AWD moveout strongly calibrates the "
                "physical sign at 25--60 Hz but is nearly balanced at 5--20 Hz. "
                "The low-frequency active gather is therefore not a decisive "
                "empirical direction standard for the ambient band. Exact "
                "5--20 Hz one-way synthetics and real-background injections "
                "remain the appropriate implementation tests; physical ambient "
                "labels must remain coordinate directions unless independently "
                "resolved."
            ),
        },
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    arrays = {
        "time_s": time,
        "coordinate_m": coordinate,
        "unfiltered": section,
    }
    for band_name, branches in reconstructed.items():
        for branch_name, branch in branches.items():
            arrays[f"{band_name}__{branch_name}"] = branch
    np.savez_compressed(OUT.with_suffix(".npz"), **arrays)
    OUT.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")

    figure, axes = plt.subplots(2, 2, figsize=(13.0, 9.0), constrained_layout=True)
    extent = [time[0], time[-1], coordinate[-1], coordinate[0]]
    for row, (fmin, fmax) in enumerate(FREQUENCY_BANDS_HZ):
        band_name = f"{fmin:g}_{fmax:g}_hz"
        panels = (
            ("expected_fk_negative", "expected: $F K<0$"),
            ("opposite_fk_positive", "complement: $F K>0$"),
        )
        common = np.concatenate(
            [reconstructed[band_name][name].ravel() for name, _ in panels]
        )
        limit = float(np.percentile(np.abs(common), 99.0))
        for column, (name, label) in enumerate(panels):
            axis = axes[row, column]
            axis.imshow(
                reconstructed[band_name][name],
                extent=extent,
                aspect="auto",
                cmap="RdBu_r",
                vmin=-limit,
                vmax=limit,
                interpolation="nearest",
            )
            prediction = intercept + (
                coordinate - coordinate[0]
            ) / OBSERVED_VELOCITY_M_S
            axis.plot(prediction, coordinate, "k--", lw=1.2)
            ratio = metrics[band_name][
                "expected_to_opposite_energy_in_fixed_downgoing_tube_ratio"
            ]
            axis.set_title(
                f"{chr(97 + row * 2 + column)}  {fmin:g}--{fmax:g} Hz, {label}"
                + (f"\nexpected/opposite tube energy = {ratio:.2f}" if column == 0 else "")
            )
            axis.set_xlabel("Time relative to catalog reference (s)")
            axis.set_ylabel("Nano coordinate (m)")
    figure.suptitle(
        "Production signed F-K fan tested against the known downgoing AWD wavefield",
        fontsize=14,
    )
    figure.savefig(OUT.with_suffix(".png"), dpi=350, bbox_inches="tight")
    figure.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
