"""Injection--recovery sensitivity test for the Deep guided-mode observable.

Question
--------
What is the smallest fractional change in the apparent speed of the repeatable
Deep guided mode that can be recovered reliably from the 2026 SAFOD AWD survey,
and is it smaller than the corresponding Nano limit?

The recovered quantity is a fractional change in the *apparent along-fiber speed
of the selected Deep guided mode*.  It is not a formation Vp or Vs change, not a
stress or pore-pressure change, and not a depth profile.  Guided-mode physics and
forward modelling would be required to convert it into any of those.

Why the Deep mode might win
---------------------------
The estimator regresses aperture delay against reference travel time,

    dt(s) = a - eps * T0(s),

so the precision of ``eps`` scales inversely with the lever arm in ``T0``.  The
Nano observable spans 80--440 m at 2975 m/s, a lever arm of 0.121 s.  The Deep
outbound leg spans 200--3000 m at ~1.58 km/s, a lever arm of ~1.8 s: roughly 17x
longer.  Whether that converts into 17x better sensitivity depends on how well
the Deep beams correlate burst to burst, which is exactly what this test
measures rather than assumes.

Design decisions, frozen before any recovery was run
----------------------------------------------------
* Held-out half.  ``deep_tube_validation.py`` selects its trajectory on epochs
  with ``epoch % 2 == 1`` and validates on ``epoch % 2 == 0``.  This script
  inherits that convention, so the 23 primary bursts here are strictly held out
  from the trajectory that defines the observable.  Note that
  ``deep_target_burst_repeatability.py`` uses the opposite parity while calling
  it by the same name; that file is not involved here.
* Primary result uses the 23 held-out bursts.  A secondary all-46-burst run is
  also reported and is explicitly labelled trajectory-contaminated; it exists
  only to show that the primary threshold is not an artefact of small N.
* One frozen trajectory per leg and band, obtained by a semblance grid search
  over (intercept, slowness) on the discovery half alone, with slowness confined
  to the already-validated 1300--1800 m/s window.
* Injection grid matches ``nano_dvv_injection_recovery.py`` exactly, so the two
  sensitivity numbers are directly comparable.
* The delay estimator and the robust regression are *imported* from the Nano
  script rather than reimplemented, so no part of the Nano/Deep difference can
  come from a difference in estimator code.

Stages, so recovery never reads the injected values
---------------------------------------------------
1. ``freeze``     grid-searches and writes the trajectory + frozen parameters.
2. ``inject``     writes randomised perturbed gathers and a sealed truth CSV.
3. ``recover``    reads only the blinded gathers and estimates the delay gradient.
4. ``summarize``  joins to truth and measures bias, scatter, null and thresholds.

``python deep_dvv_injection_recovery.py --stage all`` runs all four.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt

from fk_dispersion import weighted_stack

# Identical estimator machinery to the Nano analysis.  Importing rather than
# copying is deliberate: it makes "same estimator, different observable" a fact
# about the code and not a claim in the text.
from nano_dvv_injection_recovery import (
    INJECTED_DVV,
    _normalized_correlation_lag,
    _robust_line,
)


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"

FROZEN = HERE / "deep_dvv_frozen_trajectory.json"
TRUTH = HERE / "deep_dvv_blind_truth.csv"
RECOVERY = HERE / "deep_dvv_recovery.csv"
SUMMARY = HERE / "deep_dvv_summary.csv"
CONTROLS = HERE / "deep_dvv_controls.csv"
COMPARISON = HERE / "deep_dvv_nano_comparison.csv"
RESULTS = HERE / "deep_dvv_injection_recovery.npz"
FIGURE = HERE / "deep_dvv_injection_recovery.png"
REPORT = HERE / "deep_dvv_injection_recovery.txt"

# ---------------------------------------------------------------------------
# Frozen acquisition and geometry constants, inherited from deep_tube_validation
# ---------------------------------------------------------------------------
PRE_S = 0.5
TURNAROUND_CH = 1702
# Analysis range common to both legs, set by the span of the already-validated
# candidate apertures in deep_tube_candidates.csv (outbound 200-3000 m, return
# 0-3000 m).  Frozen to 200-3000 m so the two legs share one coordinate range.
COORD_RANGE_M = (200.0, 3000.0)
# Validated slow-mode window from deep_tube_validation.TUBE_POSITIVE.
SLOWNESS_RANGE = (1.0 / 1800.0, 1.0 / 1300.0)
CHANNEL_STRIDE = 6

APERTURE_M = 400.0
STEP_M = 200.0
APERTURE_SENSITIVITY_M = (200.0, 400.0, 800.0)

# Grid search resolution for the frozen trajectory.
N_SLOWNESS = 101
INTERCEPT_GRID_S = (-0.10, 0.40)
INTERCEPT_STEP_S = 0.002

# Quality control, frozen a priori and applied identically to every pass.
MIN_APERTURE_CORRELATION = 0.30
MIN_APERTURES = 6
BOUNDARY_FRACTION = 0.90

# Control magnitudes.
CONSTANT_SHIFT_S = 5.0e-3
# Source common-mode timing scatter measured for this survey (docs/paper1).
JITTER_SIGMA_S = 0.30e-3
DISPLACED_INTERCEPT_S = 0.12
SLOWNESS_OFFSET_S_PER_M = 5.0e-5

SEED = 20260806

LEGS = ("outbound", "return")

# Per-band configuration.  The 60-120 Hz entry is a wrong-observable control:
# no coherent mode was validated there, so recovery is expected to fail.
BAND_CONFIG = {
    "15_30": {
        "band": (15.0, 30.0),
        "extract_q_s": (-0.28, 0.40),
        "recovery_q_s": (-0.08, 0.16),
        "max_lag_s": 0.030,
        "semblance_half_s": 0.035,
        "role": "primary",
    },
    "3_15": {
        "band": (3.0, 15.0),
        "extract_q_s": (-0.45, 0.65),
        "recovery_q_s": (-0.15, 0.30),
        "max_lag_s": 0.060,
        "semblance_half_s": 0.085,
        "role": "secondary",
    },
    "60_120": {
        "band": (60.0, 120.0),
        "extract_q_s": (-0.12, 0.18),
        "recovery_q_s": (-0.04, 0.08),
        "max_lag_s": 0.012,
        "semblance_half_s": 0.012,
        "role": "control_band",
    },
}
PRIMARY_BAND = "15_30"


# Blinded per-channel gathers are bulky (~4.5 GB total) and are pure
# intermediates, so they go to $SCRATCH per cluster policy.  Every product that
# is actually read afterwards -- JSON, CSV, PNG, TXT -- stays beside the script
# with the rest of the awd_clean outputs.
BLIND_DIR = Path(os.environ.get("SCRATCH", str(HERE))) / "deep_dvv_blind"


def _blind_path(leg: str, band_tag: str, population: str) -> Path:
    return BLIND_DIR / f"deep_dvv_blind_{population}_{leg}_{band_tag}.npz"


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty table: {path}")
    names = fieldnames or list(rows[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def _truth_key(row: dict) -> tuple:
    """Join key for the sealed truth table.

    Must include the group, not just the trial id.  A bare-id join once let a
    duplicate id from a different population overwrite a held-out truth row and
    silently score ~20% of trials against the wrong injected value.
    """
    return (row["population"], row["leg"], row["band"], row["trial_id"])


def _assert_unique_truth(rows: list[dict]) -> dict[tuple, dict]:
    table: dict[tuple, dict] = {}
    for row in rows:
        key = _truth_key(row)
        if key in table:
            raise RuntimeError(f"Duplicate truth key, join would be ambiguous: {key}")
        table[key] = row
    return table


# ---------------------------------------------------------------------------
# Shared loading and preprocessing
# ---------------------------------------------------------------------------
def _leg_channels(n_channels: int, dx: float) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Return per-leg (absolute channel index, leg coordinate) after striding.

    The return leg is reversed so that leg coordinate measures distance from the
    surface end on both legs, matching deep_tube_validation's convention and
    making the two legs directly comparable.
    """
    channels = {}
    outbound = np.arange(0, TURNAROUND_CH, dtype=int)
    inbound = np.arange(TURNAROUND_CH, n_channels, dtype=int)[::-1]
    for name, absolute in (("outbound", outbound), ("return", inbound)):
        coordinate = np.arange(absolute.size, dtype=float) * dx
        keep = (coordinate >= COORD_RANGE_M[0]) & (coordinate <= COORD_RANGE_M[1])
        absolute = absolute[keep][::CHANNEL_STRIDE]
        coordinate = coordinate[keep][::CHANNEL_STRIDE]
        channels[name] = (absolute, coordinate)
    return channels


def _bandpass(section: np.ndarray, fs: float, band: tuple[float, float]) -> np.ndarray:
    sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, np.nan_to_num(section), axis=-1)


def _channel_qc(section: np.ndarray) -> np.ndarray:
    """Boolean keep-mask over channels; ``section`` is (burst, channel, time)."""
    finite_fraction = np.mean(np.isfinite(section), axis=(0, 2))
    dynamic = np.nanmedian(np.nanstd(section, axis=-1), axis=0)
    return (finite_fraction >= 0.99) & np.isfinite(dynamic) & (dynamic > 0)


def _rms_normalize(gather: np.ndarray) -> np.ndarray:
    """Scale each trace to unit RMS.

    A time-invariant per-trace scaling commutes with every time shift in this
    analysis, so it cannot bias a delay.  It stops the shallow high-amplitude
    interval from dominating a beam, the same problem deep_tube_validation had
    to handle for display.
    """
    scale = np.sqrt(np.mean(gather**2, axis=-1, keepdims=True))
    floor = np.median(scale[scale > 0]) * 1e-3 if np.any(scale > 0) else 1.0
    return gather / np.maximum(scale, floor)


def _align(
    filtered: np.ndarray,
    sample_time: np.ndarray,
    coordinate: np.ndarray,
    intercept: float,
    slowness: float,
    q: np.ndarray,
) -> np.ndarray:
    """Sample each channel along ``intercept + slowness * s + q``."""
    aligned = np.empty(filtered.shape[:-1] + (q.size,), dtype=np.float32)
    travel = intercept + slowness * coordinate
    flat = filtered.reshape(-1, filtered.shape[-1])
    out = aligned.reshape(-1, q.size)
    n_channel = coordinate.size
    for index in range(flat.shape[0]):
        out[index] = np.interp(
            travel[index % n_channel] + q, sample_time, flat[index], left=0.0, right=0.0
        )
    return aligned


# ---------------------------------------------------------------------------
# Stage 1: freeze the observable
# ---------------------------------------------------------------------------
def _semblance_profile(
    aligned: np.ndarray, q: np.ndarray, fs: float, half_s: float, intercepts: np.ndarray
) -> np.ndarray:
    """Semblance of a channel-aligned gather for each candidate intercept.

    ``aligned`` is (channel, q).  For a window of ``2*half+1`` samples centred on
    each intercept, semblance is ``sum_t (sum_ch x)^2 / (nch * sum_t sum_ch x^2)``.
    Both sums are sliding, so cumulative sums give every intercept at once.
    """
    half = max(1, int(round(half_s * fs)))
    stacked = np.sum(aligned, axis=0) ** 2
    energy = np.sum(aligned**2, axis=0)
    cum_stacked = np.concatenate(([0.0], np.cumsum(stacked)))
    cum_energy = np.concatenate(([0.0], np.cumsum(energy)))
    centers = np.searchsorted(q, intercepts)
    lo = np.clip(centers - half, 0, q.size)
    hi = np.clip(centers + half + 1, 0, q.size)
    numerator = cum_stacked[hi] - cum_stacked[lo]
    denominator = aligned.shape[0] * (cum_energy[hi] - cum_energy[lo])
    return np.divide(
        numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0
    )


def _fit_trajectory(
    section: np.ndarray, coordinate: np.ndarray, fs: float, config: dict
) -> dict:
    """Grid-search (intercept, slowness) on the discovery-half stack."""
    filtered = _rms_normalize(_bandpass(section, fs, config["band"]))
    sample_time = np.arange(filtered.shape[-1], dtype=float) / fs - PRE_S
    intercepts = np.arange(
        INTERCEPT_GRID_S[0],
        INTERCEPT_GRID_S[1] + 0.5 * INTERCEPT_STEP_S,
        INTERCEPT_STEP_S,
    )
    pad = config["semblance_half_s"] + 2.0 / fs
    q = np.arange(
        int(round(((intercepts[-1] + pad) - (intercepts[0] - pad)) * fs)) + 1, dtype=float
    ) / fs + (intercepts[0] - pad)
    slowness_grid = np.linspace(SLOWNESS_RANGE[0], SLOWNESS_RANGE[1], N_SLOWNESS)

    best = {"semblance": -np.inf}
    profile = np.empty((slowness_grid.size, intercepts.size), dtype=float)
    for islow, slowness in enumerate(slowness_grid):
        aligned = _align(filtered, sample_time, coordinate, 0.0, slowness, q)
        profile[islow] = _semblance_profile(
            aligned, q, fs, config["semblance_half_s"], intercepts
        )
        peak = int(np.argmax(profile[islow]))
        if profile[islow, peak] > best["semblance"]:
            best = {
                "semblance": float(profile[islow, peak]),
                "slowness_s_per_m": float(slowness),
                "intercept_s": float(intercepts[peak]),
            }
    best["velocity_mps"] = 1.0 / best["slowness_s_per_m"]
    best["semblance_grid_max"] = float(np.max(profile))
    return best


def freeze() -> None:
    with np.load(STACKS) as data:
        counts = np.asarray(data["n_common"], dtype=int)
        fs = float(data["fs"])
        dx = float(data["dx_deep"])
        valid = np.flatnonzero(counts > 0)
        discovery_counts = counts.copy()
        discovery_counts[valid[valid % 2 == 0]] = 0  # keep epoch % 2 == 1
        discovery = weighted_stack(data["deep_stacks"], discovery_counts)
        n_channels = discovery.shape[0]

    channels = _leg_channels(n_channels, dx)
    frozen = {
        "generated_by": "deep_dvv_injection_recovery.py --stage freeze",
        "discovery_epochs": int(np.count_nonzero(discovery_counts)),
        "discovery_parity": "epoch % 2 == 1 (inherited from deep_tube_validation.py)",
        "validation_parity": "epoch % 2 == 0",
        "fs_hz": fs,
        "dx_deep_m": dx,
        "turnaround_channel": TURNAROUND_CH,
        "coordinate_range_m": list(COORD_RANGE_M),
        "channel_stride": CHANNEL_STRIDE,
        "slowness_search_range_s_per_m": list(SLOWNESS_RANGE),
        "aperture_m": APERTURE_M,
        "aperture_step_m": STEP_M,
        "min_aperture_correlation": MIN_APERTURE_CORRELATION,
        "min_apertures": MIN_APERTURES,
        "injected_levels": [float(v) for v in INJECTED_DVV],
        "seed": SEED,
        "trajectories": {},
    }
    for leg in LEGS:
        absolute, coordinate = channels[leg]
        section = discovery[absolute]
        for band_tag, config in BAND_CONFIG.items():
            result = _fit_trajectory(section, coordinate, fs, config)
            result["n_channels"] = int(coordinate.size)
            result["lever_arm_s"] = float(
                result["slowness_s_per_m"] * (coordinate[-1] - coordinate[0])
            )
            result["role"] = config["role"]
            frozen["trajectories"][f"{leg}|{band_tag}"] = result
            print(
                f"freeze {leg:9s} {band_tag:6s} "
                f"v={result['velocity_mps']:7.1f} m/s  "
                f"t0={result['intercept_s']:+.3f} s  "
                f"semblance={result['semblance']:.4f}  "
                f"lever arm={result['lever_arm_s']:.3f} s"
            )
    FROZEN.write_text(json.dumps(frozen, indent=2) + "\n")
    print(f"Wrote frozen observable definition to {FROZEN}")


# ---------------------------------------------------------------------------
# Stage 2: inject
# ---------------------------------------------------------------------------
def _extract(
    stacks: np.ndarray,
    epochs: np.ndarray,
    absolute: np.ndarray,
    coordinate: np.ndarray,
    fs: float,
    config: dict,
    trajectory: dict,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return per-burst channel gathers aligned on the frozen trajectory."""
    # np.ix_ selects both axes in one pass.  stacks[epochs][:, absolute] would
    # first materialise a 46 x 3200 x 3500 copy -- 2 GB thrown away immediately.
    raw = np.asarray(stacks[np.ix_(epochs, absolute)], dtype=float)
    keep = _channel_qc(raw)
    if np.sum(keep) < 4 * MIN_APERTURES:
        raise RuntimeError(f"Only {np.sum(keep)} Deep channels pass QC")
    raw = np.nan_to_num(raw[:, keep], copy=False)
    coordinate = coordinate[keep]
    filtered = _rms_normalize(_bandpass(raw, fs, config["band"]))
    sample_time = np.arange(filtered.shape[-1], dtype=float) / fs - PRE_S
    q0, q1 = config["extract_q_s"]
    q = np.arange(int(round((q1 - q0) * fs)), dtype=float) / fs + q0
    aligned = _align(
        filtered,
        sample_time,
        coordinate,
        trajectory["intercept_s"],
        trajectory["slowness_s_per_m"],
        q,
    )
    return aligned, coordinate, q


def _leave_one_out(aligned: np.ndarray, weights: np.ndarray) -> np.ndarray:
    total = float(np.sum(weights))
    weighted_sum = np.tensordot(weights.astype(float), aligned, axes=(0, 0))
    references = np.empty_like(aligned)
    for index, weight in enumerate(weights):
        denominator = total - float(weight)
        if denominator <= 0:
            raise RuntimeError("At least two populated bursts are required")
        references[index] = (weighted_sum - float(weight) * aligned[index]) / denominator
    return references


def _split_half_reference(aligned: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Reference from the opposite time-ordered half of the burst population."""
    n = aligned.shape[0]
    midpoint = n // 2
    first = np.average(aligned[:midpoint], axis=0, weights=weights[:midpoint])
    second = np.average(aligned[midpoint:], axis=0, weights=weights[midpoint:])
    references = np.empty_like(aligned)
    references[:midpoint] = second
    references[midpoint:] = first
    return references


def _shift_gather(gather: np.ndarray, q: np.ndarray, shift: np.ndarray) -> np.ndarray:
    """Delay each channel by ``shift`` seconds (positive = arrives later)."""
    shifted = np.empty_like(gather, dtype=np.float32)
    shift = np.broadcast_to(np.atleast_1d(shift), (gather.shape[0],))
    for index in range(gather.shape[0]):
        shifted[index] = np.interp(
            q - shift[index], q, gather[index], left=0.0, right=0.0
        )
    return shifted


def _inject_population(population: str, epoch_filter) -> None:
    # Seeded per population.  Two default_rng(SEED) instances walk the SAME
    # PCG64 word stream, so drawing trial ids from both produced ~20% duplicate
    # ids and silently cross-joined the two populations' truth rows.
    population_index = 0 if population == "heldout" else 1
    with np.load(STACKS) as data:
        counts = np.asarray(data["n_common"], dtype=int)
        fs = float(data["fs"])
        dx = float(data["dx_deep"])
        valid = np.flatnonzero(counts > 0)
        epochs = valid[epoch_filter(valid)]
        stacks = data["deep_stacks"]
        n_channels = stacks.shape[1]
        channels = _leg_channels(n_channels, dx)
        frozen = json.loads(FROZEN.read_text())
        rng = np.random.default_rng([SEED, population_index])
        truth_rows = []
        for leg in LEGS:
            absolute, coordinate = channels[leg]
            for band_tag, config in BAND_CONFIG.items():
                trajectory = frozen["trajectories"][f"{leg}|{band_tag}"]
                aligned, kept_coordinate, q = _extract(
                    stacks, epochs, absolute, coordinate, fs, config, trajectory
                )
                weights = counts[epochs]
                loo = _leave_one_out(aligned, weights)
                split = _split_half_reference(aligned, weights)
                travel = (
                    trajectory["intercept_s"]
                    + trajectory["slowness_s_per_m"] * kept_coordinate
                )

                specifications = [
                    (i, "dvv", float(v))
                    for i in range(aligned.shape[0])
                    for v in INJECTED_DVV
                ]
                specifications += [
                    (i, "constant_shift", CONSTANT_SHIFT_S)
                    for i in range(aligned.shape[0])
                ]
                specifications += [
                    (i, "source_jitter", float(rng.normal(0.0, JITTER_SIGMA_S)))
                    for i in range(aligned.shape[0])
                ]
                rng.shuffle(specifications)

                trials = np.empty(
                    (len(specifications), aligned.shape[1], aligned.shape[2]),
                    dtype=np.float32,
                )
                target_index = np.empty(len(specifications), dtype=np.int16)
                ids = []
                # The group tag makes the id unique by construction.  It leaks
                # only the group, which the recovery stage already knows from
                # the filename; it says nothing about the injected value.
                group_tag = f"{population[:1]}{leg[:1]}{band_tag.replace('_', '')}"
                for itrial, (index, kind, value) in enumerate(specifications):
                    trial_id = (
                        f"D{group_tag}{rng.integers(0, 2**63, dtype=np.int64):016x}"
                    )
                    ids.append(trial_id)
                    target_index[itrial] = index
                    if kind == "dvv":
                        # eps > 0 speeds the mode up, so the arrival moves earlier
                        # by eps * T0.  Same convention as the Nano injection.
                        trials[itrial] = _shift_gather(
                            aligned[index], q, -value * travel
                        )
                    else:
                        trials[itrial] = _shift_gather(
                            aligned[index], q, np.full(travel.shape, value)
                        )
                    truth_rows.append(
                        {
                            "trial_id": trial_id,
                            "population": population,
                            "leg": leg,
                            "band": band_tag,
                            "target_epoch": int(epochs[index]),
                            "trial_kind": kind,
                            "injected_value": f"{value:.9g}",
                        }
                    )

                np.savez(
                    _blind_path(leg, band_tag, population),
                    trial_id=np.asarray(ids),
                    trial_gather=trials,
                    reference_loo=loo,
                    reference_split=split,
                    target_index=target_index,
                    q_s=q,
                    coordinate_m=kept_coordinate,
                    travel_time_s=travel,
                    fs=fs,
                    intercept_s=trajectory["intercept_s"],
                    slowness_s_per_m=trajectory["slowness_s_per_m"],
                    epochs=epochs,
                    weights=weights,
                )
                print(
                    f"inject {population:9s} {leg:9s} {band_tag:6s} "
                    f"{len(specifications):5d} trials, "
                    f"{kept_coordinate.size} channels"
                )
    return truth_rows


def inject() -> None:
    BLIND_DIR.mkdir(parents=True, exist_ok=True)
    rows = _inject_population("heldout", lambda valid: valid % 2 == 0)
    rows += _inject_population("allbursts", lambda valid: np.ones(valid.size, bool))
    _assert_unique_truth(rows)
    _write_csv(TRUTH, rows)
    print(f"Sealed truth table written separately to {TRUTH}")


# ---------------------------------------------------------------------------
# Stage 3: recover
# ---------------------------------------------------------------------------
def _aperture_layout(coordinate: np.ndarray, width_m: float) -> list[tuple[np.ndarray, float]]:
    layout = []
    start = coordinate[0]
    while start + width_m <= coordinate[-1] + 1e-9:
        members = np.flatnonzero(
            (coordinate >= start) & (coordinate < start + width_m)
        )
        if members.size >= 8:
            layout.append((members, float(coordinate[members].mean())))
        start += STEP_M
    return layout


def _beams(
    gather: np.ndarray,
    q: np.ndarray,
    coordinate: np.ndarray,
    layout: list[tuple[np.ndarray, float]],
    realign: np.ndarray | None,
) -> np.ndarray:
    """Moveout-corrected local beam traces, one per aperture."""
    if realign is not None:
        gather = _shift_gather(gather, q, -realign)
    return np.stack([gather[members].mean(axis=0) for members, _ in layout])


def _estimate(
    target: np.ndarray,
    reference: np.ndarray,
    layout: list[tuple[np.ndarray, float]],
    travel_centers: np.ndarray,
    q: np.ndarray,
    fs: float,
    config: dict,
) -> dict:
    use = (q >= config["recovery_q_s"][0]) & (q <= config["recovery_q_s"][1])
    max_lag = config["max_lag_s"]
    n = len(layout)
    delay = np.full(n, np.nan)
    correlation = np.full(n, np.nan)
    for index in range(n):
        delay[index], correlation[index] = _normalized_correlation_lag(
            reference[index][use], target[index][use], fs, max_lag
        )
    good = (
        np.isfinite(delay)
        & np.isfinite(correlation)
        & (correlation >= MIN_APERTURE_CORRELATION)
        & (np.abs(delay) <= BOUNDARY_FRACTION * max_lag)
    )
    if np.sum(good) < MIN_APERTURES:
        return {
            "estimated_dvv": np.nan,
            "estimated_dvv_se": np.nan,
            "nuisance_intercept_s": np.nan,
            "residual_scale_s": np.nan,
            "median_aperture_correlation": np.nan,
            "n_apertures": int(np.sum(good)),
        }
    weight = np.maximum(correlation[good], 0.0) ** 2
    beta, standard_error, residual_scale = _robust_line(
        travel_centers[good], delay[good], weight
    )
    return {
        "estimated_dvv": -float(beta[1]),
        "estimated_dvv_se": float(standard_error[1]),
        "nuisance_intercept_s": float(beta[0]),
        "residual_scale_s": float(residual_scale),
        "median_aperture_correlation": float(np.median(correlation[good])),
        "n_apertures": int(np.sum(good)),
    }


def _recovery_passes() -> list[dict]:
    """Every recovery configuration applied to the same blinded data."""
    passes = [{"name": "primary", "width": APERTURE_M, "reference": "loo",
               "d_intercept": 0.0, "d_slowness": 0.0}]
    for width in APERTURE_SENSITIVITY_M:
        if width != APERTURE_M:
            passes.append({"name": f"aperture_{int(width)}m", "width": width,
                           "reference": "loo", "d_intercept": 0.0, "d_slowness": 0.0})
    passes.append({"name": "reference_split", "width": APERTURE_M,
                   "reference": "split", "d_intercept": 0.0, "d_slowness": 0.0})
    passes.append({"name": "displaced_intercept", "width": APERTURE_M,
                   "reference": "loo", "d_intercept": DISPLACED_INTERCEPT_S,
                   "d_slowness": 0.0})
    for sign in (+1.0, -1.0):
        passes.append({"name": f"slowness_offset_{'p' if sign > 0 else 'm'}",
                       "width": APERTURE_M, "reference": "loo", "d_intercept": 0.0,
                       "d_slowness": sign * SLOWNESS_OFFSET_S_PER_M})
    return passes


def recover() -> None:
    rows = []
    for population in ("heldout", "allbursts"):
        for leg in LEGS:
            for band_tag, config in BAND_CONFIG.items():
                path = _blind_path(leg, band_tag, population)
                if not path.exists():
                    continue
                with np.load(path) as blind:
                    trial_ids = blind["trial_id"].astype(str)
                    trials = blind["trial_gather"]
                    references = {
                        "loo": blind["reference_loo"],
                        "split": blind["reference_split"],
                    }
                    target_index = blind["target_index"]
                    q = blind["q_s"]
                    coordinate = blind["coordinate_m"]
                    fs = float(blind["fs"])
                    intercept = float(blind["intercept_s"])
                    slowness = float(blind["slowness_s_per_m"])

                    for spec in _recovery_passes():
                        layout = _aperture_layout(coordinate, spec["width"])
                        if len(layout) < MIN_APERTURES:
                            continue
                        assumed_slowness = slowness + spec["d_slowness"]
                        assumed_intercept = intercept + spec["d_intercept"]
                        # Re-alignment needed to move from the frozen trajectory
                        # used at extraction to the trajectory this pass assumes.
                        realign = (
                            spec["d_intercept"] + spec["d_slowness"] * coordinate
                        )
                        realign = None if np.allclose(realign, 0.0) else realign
                        centers = np.asarray([c for _, c in layout])
                        travel_centers = assumed_intercept + assumed_slowness * centers
                        reference_stack = references[spec["reference"]]
                        cache: dict[int, np.ndarray] = {}
                        for itrial, trial_id in enumerate(trial_ids):
                            index = int(target_index[itrial])
                            if index not in cache:
                                cache[index] = _beams(
                                    reference_stack[index], q, coordinate, layout, realign
                                )
                            target = _beams(
                                trials[itrial], q, coordinate, layout, realign
                            )
                            result = _estimate(
                                target, cache[index], layout, travel_centers, q, fs, config
                            )
                            rows.append(
                                {
                                    "trial_id": trial_id,
                                    "population": population,
                                    "leg": leg,
                                    "band": band_tag,
                                    "pass": spec["name"],
                                    "assumed_velocity_mps": f"{1.0 / assumed_slowness:.4f}",
                                    "aperture_m": spec["width"],
                                    "reference": spec["reference"],
                                    **{
                                        key: (
                                            f"{value:.10g}"
                                            if isinstance(value, float)
                                            else value
                                        )
                                        for key, value in result.items()
                                    },
                                }
                            )
                        print(
                            f"recover {population:9s} {leg:9s} {band_tag:6s} "
                            f"{spec['name']:20s} {len(layout):3d} apertures"
                        )
    _write_csv(RECOVERY, rows)
    print(f"Recovery stage estimated {len(rows)} blinded passes and wrote {RECOVERY}")


# ---------------------------------------------------------------------------
# Stage 4: summarize
# ---------------------------------------------------------------------------
def _level_row(truth: float, estimated: np.ndarray, threshold: float) -> dict:
    finite = np.isfinite(estimated)
    values = estimated[finite]
    error = values - truth
    detected = np.abs(values) > threshold
    if truth == 0:
        correct_sign = np.zeros(values.size, dtype=bool)
        correct_detection = detected
    else:
        correct_sign = np.sign(values) == np.sign(truth)
        correct_detection = detected & correct_sign
    median = float(np.median(values)) if values.size else np.nan
    mad = float(1.4826 * np.median(np.abs(values - median))) if values.size else np.nan
    return {
        "injected_dvv": truth,
        "n_trials": int(estimated.size),
        "n_recovered": int(values.size),
        "median_estimated_dvv": median,
        "bias": float(np.mean(error)) if values.size else np.nan,
        "robust_scatter_1p4826_mad": mad,
        "rmse": float(np.sqrt(np.mean(error**2))) if values.size else np.nan,
        "detection_probability": float(np.mean(detected)) if values.size else np.nan,
        "correct_sign_probability": (
            float(np.mean(correct_sign)) if values.size and truth != 0 else np.nan
        ),
        "correct_sign_detection_probability": (
            float(np.mean(correct_detection)) if values.size else np.nan
        ),
    }


def _smallest_reliable(summary_rows: list[dict], key: str) -> float:
    positive = {r["injected_dvv"]: r for r in summary_rows if r["injected_dvv"] > 0}
    negative = {abs(r["injected_dvv"]): r for r in summary_rows if r["injected_dvv"] < 0}
    qualifying = [
        value
        for value in sorted(set(positive) & set(negative))
        if (positive[value][key] or 0) >= 0.95 and (negative[value][key] or 0) >= 0.95
    ]
    return float(min(qualifying)) if qualifying else np.nan


def _analyse_group(injected: np.ndarray, estimated_raw: np.ndarray) -> dict:
    null = estimated_raw[injected == 0]
    finite_null = null[np.isfinite(null)]
    if finite_null.size < 10:
        return {}
    null_center = float(np.median(finite_null))
    estimated = estimated_raw - null_center
    null_calibrated = estimated[injected == 0]
    null_calibrated = null_calibrated[np.isfinite(null_calibrated)]
    threshold = float(np.quantile(np.abs(null_calibrated), 0.95))
    rows = [
        _level_row(float(level), estimated[injected == level], threshold)
        for level in np.unique(injected)
    ]
    return {
        "summary_rows": rows,
        "estimated": estimated,
        "null_calibrated": null_calibrated,
        "null_center": null_center,
        "threshold": threshold,
        "false_positive_rate": float(np.mean(np.abs(null_calibrated) > threshold)),
        "n_null": int(null_calibrated.size),
        "reliable_detection": _smallest_reliable(rows, "correct_sign_detection_probability"),
        "reliable_sign": _smallest_reliable(rows, "correct_sign_probability"),
    }


def summarize() -> None:
    truth = _assert_unique_truth(_read_csv(TRUTH))
    recovery = _read_csv(RECOVERY)
    for row in recovery:
        key = _truth_key(row)
        if key not in truth:
            raise RuntimeError(f"Recovery trial absent from truth table: {key}")
        row.update(
            {
                "trial_kind": truth[key]["trial_kind"],
                "injected_value": float(truth[key]["injected_value"]),
                "target_epoch": int(truth[key]["target_epoch"]),
            }
        )
    # Every dvv level must contain exactly one trial per burst in its population.
    # This is the check that would have caught the bare-id join immediately.
    for population, expected in (("heldout", 23), ("allbursts", 46)):
        counts = {}
        for row in recovery:
            if row["population"] != population or row["trial_kind"] != "dvv":
                continue
            group = (row["leg"], row["band"], row["pass"], row["injected_value"])
            counts[group] = counts.get(group, 0) + 1
        wrong = {group: n for group, n in counts.items() if n != expected}
        if wrong:
            raise RuntimeError(
                f"{population}: {len(wrong)} level groups do not have exactly "
                f"{expected} trials, e.g. {list(wrong.items())[:3]}"
            )

    groups: dict[tuple, dict] = {}
    summary_rows = []
    for population in ("heldout", "allbursts"):
        for leg in LEGS:
            for band_tag in BAND_CONFIG:
                for pass_name in sorted({r["pass"] for r in recovery}):
                    subset = [
                        r
                        for r in recovery
                        if r["population"] == population
                        and r["leg"] == leg
                        and r["band"] == band_tag
                        and r["pass"] == pass_name
                        and r["trial_kind"] == "dvv"
                    ]
                    if not subset:
                        continue
                    injected = np.asarray([r["injected_value"] for r in subset])
                    estimated = np.asarray([float(r["estimated_dvv"]) for r in subset])
                    analysis = _analyse_group(injected, estimated)
                    if not analysis:
                        continue
                    key = (population, leg, band_tag, pass_name)
                    groups[key] = analysis
                    for row in analysis["summary_rows"]:
                        summary_rows.append(
                            {
                                "population": population,
                                "leg": leg,
                                "band": band_tag,
                                "pass": pass_name,
                                **row,
                            }
                        )
    _write_csv(SUMMARY, summary_rows)

    # Timing-jitter controls: eps must stay near zero while the intercept moves.
    control_rows = []
    for kind in ("constant_shift", "source_jitter"):
        for population in ("heldout", "allbursts"):
            for leg in LEGS:
                for band_tag in BAND_CONFIG:
                    subset = [
                        r
                        for r in recovery
                        if r["trial_kind"] == kind
                        and r["population"] == population
                        and r["leg"] == leg
                        and r["band"] == band_tag
                        and r["pass"] == "primary"
                    ]
                    if not subset:
                        continue
                    estimated = np.asarray([float(r["estimated_dvv"]) for r in subset])
                    intercepts = np.asarray(
                        [float(r["nuisance_intercept_s"]) for r in subset]
                    )
                    injected = np.asarray([r["injected_value"] for r in subset])
                    finite = np.isfinite(estimated)
                    key = (population, leg, band_tag, "primary")
                    threshold = groups.get(key, {}).get("threshold", np.nan)
                    control_rows.append(
                        {
                            "control": kind,
                            "population": population,
                            "leg": leg,
                            "band": band_tag,
                            "n_trials": int(estimated.size),
                            "median_injected_shift_s": float(np.median(injected)),
                            "median_recovered_dvv": float(np.median(estimated[finite])),
                            "median_recovered_intercept_s": float(
                                np.median(intercepts[np.isfinite(intercepts)])
                            ),
                            "fraction_beyond_null_threshold": float(
                                np.mean(np.abs(estimated[finite]) > threshold)
                            )
                            if np.isfinite(threshold)
                            else np.nan,
                        }
                    )
    if control_rows:
        _write_csv(CONTROLS, control_rows)

    # Nano/Deep comparison on the primary pass.
    nano = np.load(HERE / "nano_dvv_injection_recovery.npz", allow_pickle=True)
    comparison_rows = [
        {
            "observable": "Nano apparent moveout",
            "band_hz": "30-60",
            "leg": "single",
            "population": "all 46 bursts",
            "lever_arm_s": f"{(440.0 - 80.0) / 2975.0:.4f}",
            "null_threshold_dvv": f"{float(nano['empirical_95pct_threshold_dvv']):.6g}",
            "false_positive_rate": f"{float(nano['empirical_false_positive_rate']):.4f}",
            "reliable_level_detection": f"{float(nano['symmetric_95pct_detection_limit_dvv']):.6g}",
            "reliable_level_sign_only": "not computed by nano script",
        }
    ]
    frozen = json.loads(FROZEN.read_text())
    for population in ("heldout", "allbursts"):
        for leg in LEGS:
            for band_tag, config in BAND_CONFIG.items():
                analysis = groups.get((population, leg, band_tag, "primary"))
                if not analysis:
                    continue
                trajectory = frozen["trajectories"][f"{leg}|{band_tag}"]
                comparison_rows.append(
                    {
                        "observable": f"Deep guided mode ({config['role']})",
                        "band_hz": band_tag.replace("_", "-"),
                        "leg": leg,
                        "population": population,
                        "lever_arm_s": f"{trajectory['lever_arm_s']:.4f}",
                        "null_threshold_dvv": f"{analysis['threshold']:.6g}",
                        "false_positive_rate": f"{analysis['false_positive_rate']:.4f}",
                        "reliable_level_detection": f"{analysis['reliable_detection']:.6g}",
                        "reliable_level_sign_only": f"{analysis['reliable_sign']:.6g}",
                    }
                )
    _write_csv(COMPARISON, comparison_rows)

    _figure(groups, frozen, nano)
    _report(groups, frozen, nano, control_rows)

    np.savez(
        RESULTS,
        **{
            f"{'|'.join(key)}|{name}": value
            for key, analysis in groups.items()
            for name, value in (
                ("estimated", analysis["estimated"]),
                ("null", analysis["null_calibrated"]),
                ("threshold", np.asarray(analysis["threshold"])),
                ("reliable_detection", np.asarray(analysis["reliable_detection"])),
            )
        },
        injected_levels=INJECTED_DVV,
        seed=SEED,
    )
    print(f"Wrote {SUMMARY}, {CONTROLS}, {COMPARISON}, {RESULTS}, {FIGURE}, {REPORT}")


def _figure(groups: dict, frozen: dict, nano) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16.0, 9.0), constrained_layout=True)
    colors = {"outbound": "#2166ac", "return": "#b2182b"}
    levels = np.unique(INJECTED_DVV)

    ax = axes[0, 0]
    for leg in LEGS:
        analysis = groups.get(("heldout", leg, PRIMARY_BAND, "primary"))
        if not analysis:
            continue
        rows = analysis["summary_rows"]
        medians = np.asarray([r["median_estimated_dvv"] for r in rows])
        scatter = np.asarray([r["robust_scatter_1p4826_mad"] for r in rows])
        ax.errorbar(
            [r["injected_dvv"] for r in rows], medians, yerr=scatter, fmt="o",
            color=colors[leg], capsize=2, label=f"Deep {leg}", alpha=0.85,
        )
    limit = 1.05 * float(np.max(np.abs(INJECTED_DVV)))
    ax.plot([-limit, limit], [-limit, limit], color="0.15", ls="--", lw=1, label="one-to-one")
    ax.set(
        xlabel="injected fractional apparent-velocity change",
        ylabel="recovered fractional change",
        title="A  Blind recovery, held-out bursts, 15-30 Hz",
    )
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    for leg in LEGS:
        analysis = groups.get(("heldout", leg, PRIMARY_BAND, "primary"))
        if not analysis:
            continue
        rows = analysis["summary_rows"]
        ax.plot([r["injected_dvv"] for r in rows], [r["bias"] for r in rows],
                "o-", color=colors[leg], label=f"{leg} bias")
        ax.plot([r["injected_dvv"] for r in rows],
                [r["robust_scatter_1p4826_mad"] for r in rows],
                "s--", color=colors[leg], alpha=0.6, label=f"{leg} scatter")
    ax.axhline(0, color="0.25", lw=0.8)
    ax.set(xlabel="injected fractional change", ylabel="fractional change",
           title="B  Bias and inter-burst scatter")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 2]
    for leg in LEGS:
        analysis = groups.get(("heldout", leg, PRIMARY_BAND, "primary"))
        if not analysis:
            continue
        rows = [r for r in analysis["summary_rows"] if r["injected_dvv"] != 0]
        for marker, sign in (("^", 1), ("v", -1)):
            selected = [r for r in rows if np.sign(r["injected_dvv"]) == sign]
            ax.plot(
                [abs(r["injected_dvv"]) for r in selected],
                [r["correct_sign_detection_probability"] for r in selected],
                marker, color=colors[leg], ls="-", alpha=0.8,
                label=f"{leg} {'positive' if sign > 0 else 'negative'}",
            )
        if np.isfinite(analysis["reliable_detection"]):
            ax.axvline(analysis["reliable_detection"], color=colors[leg], ls=":", lw=1.3)
    ax.axhline(0.95, color="0.2", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set(xlabel="absolute injected fractional change",
           ylabel="correct-sign detection probability", ylim=(-0.04, 1.04),
           title="C  Detection reliability")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1, 0]
    for leg in LEGS:
        analysis = groups.get(("heldout", leg, PRIMARY_BAND, "primary"))
        if not analysis:
            continue
        ax.hist(analysis["null_calibrated"], bins=12, alpha=0.55,
                color=colors[leg], label=f"{leg} (n={analysis['n_null']})")
        ax.axvline(analysis["threshold"], color=colors[leg], ls="--", lw=1.2)
        ax.axvline(-analysis["threshold"], color=colors[leg], ls="--", lw=1.2)
    ax.set(xlabel="recovered change for zero injection", ylabel="trials",
           title="D  Zero-injection null and 95% threshold")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    labels, values, bar_colors = [], [], []
    labels.append("Nano threshold")
    values.append(float(nano["empirical_95pct_threshold_dvv"]))
    bar_colors.append("0.45")
    labels.append("Nano reliable")
    values.append(float(nano["symmetric_95pct_detection_limit_dvv"]))
    bar_colors.append("0.7")
    for leg in LEGS:
        analysis = groups.get(("heldout", leg, PRIMARY_BAND, "primary"))
        if not analysis:
            continue
        labels += [f"Deep {leg}\nthreshold", f"Deep {leg}\nreliable"]
        values += [analysis["threshold"], analysis["reliable_detection"]]
        bar_colors += [colors[leg], colors[leg]]
    finite = [v if np.isfinite(v) else np.nan for v in values]
    ax.bar(range(len(finite)), finite, color=bar_colors)
    ax.set_yscale("log")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7, rotation=30, ha="right")
    ax.set(ylabel="fractional apparent-velocity change",
           title="E  Nano vs Deep sensitivity")
    ax.grid(axis="y", alpha=0.25)

    ax = axes[1, 2]
    control_names = ["primary", "reference_split", "aperture_200m", "aperture_800m",
                     "displaced_intercept", "slowness_offset_p", "slowness_offset_m"]
    for leg in LEGS:
        thresholds, present = [], []
        for name in control_names:
            analysis = groups.get(("heldout", leg, PRIMARY_BAND, name))
            if analysis:
                thresholds.append(analysis["threshold"])
                present.append(name)
        if thresholds:
            ax.plot(range(len(thresholds)), thresholds, "o-", color=colors[leg], label=leg)
            ax.set_xticks(range(len(present)))
            ax.set_xticklabels(present, fontsize=7, rotation=30, ha="right")
    control_band = groups.get(("heldout", "outbound", "60_120", "primary"))
    if control_band:
        ax.axhline(control_band["threshold"], color="0.3", ls=":", lw=1.2,
                   label="60-120 Hz control band")
    ax.set_yscale("log")
    ax.set(ylabel="empirical 95% null threshold", title="F  Control passes")
    ax.legend(frameon=False, fontsize=7)

    fig.suptitle(
        "SAFOD AWD Deep guided mode: blind injection--recovery sensitivity\n"
        "recovered quantity is a fractional change in apparent guided-mode speed, "
        "not a formation velocity change",
        fontsize=12,
    )
    fig.savefig(FIGURE, dpi=250)
    plt.close(fig)


def _report(groups: dict, frozen: dict, nano, control_rows: list[dict]) -> None:
    lines = [
        "SAFOD AWD Deep guided-mode blind injection--recovery",
        "=" * 56,
        f"Input: {STACKS.name}",
        f"Held-out bursts: epoch % 2 == 0 ({frozen['validation_parity']}); "
        f"trajectory selected on {frozen['discovery_epochs']} discovery bursts",
        f"Analysis range: {COORD_RANGE_M[0]:.0f}-{COORD_RANGE_M[1]:.0f} m along leg, "
        f"stride {CHANNEL_STRIDE} ({CHANNEL_STRIDE * frozen['dx_deep_m']:.1f} m)",
        f"Apertures: {APERTURE_M:.0f} m wide, {STEP_M:.0f} m step",
        "",
        "Frozen trajectories",
        "-" * 56,
    ]
    for key, trajectory in frozen["trajectories"].items():
        lines.append(
            f"  {key:20s} v={trajectory['velocity_mps']:7.1f} m/s  "
            f"t0={trajectory['intercept_s']:+.3f} s  "
            f"lever arm={trajectory['lever_arm_s']:.3f} s  "
            f"semblance={trajectory['semblance']:.4f}  [{trajectory['role']}]"
        )
    nano_lever = (440.0 - 80.0) / 2975.0
    lines += [
        "",
        "Nano reference (nano_dvv_injection_recovery.npz)",
        "-" * 56,
        f"  lever arm                {nano_lever:.4f} s",
        f"  95% null threshold       {float(nano['empirical_95pct_threshold_dvv']):.6g}",
        f"  false-positive rate      {float(nano['empirical_false_positive_rate']):.4f}",
        f"  reliable level           {float(nano['symmetric_95pct_detection_limit_dvv']):.6g}",
        "",
        "Deep primary pass",
        "-" * 56,
    ]
    for population in ("heldout", "allbursts"):
        for leg in LEGS:
            for band_tag, config in BAND_CONFIG.items():
                analysis = groups.get((population, leg, band_tag, "primary"))
                if not analysis:
                    continue
                trajectory = frozen["trajectories"][f"{leg}|{band_tag}"]
                ratio = float(nano["empirical_95pct_threshold_dvv"]) / analysis["threshold"]
                lines.append(
                    f"  {population:9s} {leg:9s} {band_tag:6s} [{config['role']}]"
                )
                lines.append(
                    f"      lever arm {trajectory['lever_arm_s']:.3f} s "
                    f"({trajectory['lever_arm_s'] / nano_lever:.1f}x Nano); "
                    f"null n={analysis['n_null']}"
                )
                lines.append(
                    f"      95% null threshold  {analysis['threshold']:.6g} "
                    f"({ratio:.1f}x better than Nano)"
                )
                lines.append(
                    f"      false-positive rate {analysis['false_positive_rate']:.4f}"
                )
                lines.append(
                    f"      reliable level      "
                    f"detection {analysis['reliable_detection']:.6g}, "
                    f"sign-only {analysis['reliable_sign']:.6g}"
                )
    lines += ["", "Timing controls (eps must stay inside the null)", "-" * 56]
    for row in control_rows:
        if row["population"] != "heldout" or row["band"] != PRIMARY_BAND:
            continue
        lines.append(
            f"  {row['control']:15s} {row['leg']:9s} "
            f"injected shift {row['median_injected_shift_s'] * 1e3:+.3f} ms -> "
            f"eps {row['median_recovered_dvv']:+.3e}, "
            f"intercept {row['median_recovered_intercept_s'] * 1e3:+.3f} ms, "
            f"beyond null {row['fraction_beyond_null_threshold']:.3f}"
        )
    lines += [
        "",
        "Interpretation limits",
        "-" * 56,
        "  The recovered quantity is a fractional change in the apparent along-fiber",
        "  speed of the selected Deep guided mode. It is not formation Vp or Vs, not",
        "  fault-zone stress, pore pressure, permeability, fracture compliance, or",
        "  tectonic strain, and it carries no depth resolution. Converting it into any",
        "  of those requires guided-wave physics and forward modelling not done here.",
        "  The all-burst population is trajectory-contaminated and is reported only as",
        "  a check that the held-out threshold is not an artefact of 23 null trials.",
        "",
    ]
    report = "\n".join(lines) + "\n"
    REPORT.write_text(report)
    print(report, end="")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("freeze", "inject", "recover", "summarize", "all"),
        default="all",
    )
    args = parser.parse_args()
    if args.stage in ("freeze", "all"):
        freeze()
    if args.stage in ("inject", "all"):
        inject()
    if args.stage in ("recover", "all"):
        recover()
    if args.stage in ("summarize", "all"):
        summarize()


if __name__ == "__main__":
    main()
