"""Blind injection--recovery test for Nano apparent-velocity monitoring.

This script deliberately remains phase-neutral.  It uses the accepted 30--60 Hz
Nano moveout trajectory only as a repeatable observable and estimates fractional
changes in its apparent moveout.  It does not assume that the ridge is direct P
and does not interpret a recovered perturbation as formation velocity change.

The workflow is separated into three stages so recovery never reads the injected
values:

1. ``inject`` writes randomized perturbed gathers and a separate sealed truth CSV.
2. ``recover`` reads only the blinded gather file and estimates differential slope.
3. ``summarize`` joins estimates to truth and measures bias, scatter, false-positive
   rate, and detection probability.

Run ``python nano_dvv_injection_recovery.py --stage all`` for all three stages, or
run the stages as independent jobs for a stricter blind audit.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
BLIND = HERE / "nano_dvv_blind_trials.npz"
TRUTH = HERE / "nano_dvv_blind_truth.csv"
RECOVERY = HERE / "nano_dvv_recovery.csv"
SUMMARY = HERE / "nano_dvv_summary.csv"
RESULTS = HERE / "nano_dvv_injection_recovery.npz"
FIGURE = HERE / "nano_dvv_injection_recovery.png"
REPORT = HERE / "nano_dvv_injection_recovery.txt"

# Accepted phase-neutral trajectory.  The intercept is relative to the first
# retained aperture coordinate, consistent with the moveout/null-test analysis.
BAND_HZ = (30.0, 60.0)
APERTURE_M = (80.0, 440.0)
REFERENCE_SPEED_MPS = 2975.0
REFERENCE_INTERCEPT_S = -0.022
PRE_S = 0.5
CHANNEL_STRIDE = 4
EXTRACT_Q_S = (-0.080, 0.200)
RECOVERY_Q_S = (-0.020, 0.140)
MAX_LAG_S = 0.012
MIN_CHANNEL_CORRELATION = 0.20
MIN_VALID_CHANNELS = 20
SEED = 20260802
INJECTED_DVV = np.asarray(
    [
        -1e-2, -5e-3, -2e-3, -1e-3, -5e-4, -2e-4, -1e-4,
        0.0,
        1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2,
    ],
    dtype=float,
)


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


def _extract_aligned_gathers() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """Return bandpassed burst gathers aligned on the fixed Nano trajectory."""
    with np.load(STACKS) as data:
        fs = float(data["fs"])
        dx = float(data["dx_nano"])
        counts_all = np.asarray(data["n_common"], dtype=int)
        valid_epochs = np.flatnonzero(counts_all > 0)
        c0 = int(np.ceil(APERTURE_M[0] / dx))
        c1 = min(int(np.floor(APERTURE_M[1] / dx)) + 1, data["nano_stacks"].shape[1])
        channels = np.arange(c0, c1, CHANNEL_STRIDE, dtype=int)
        raw = np.asarray(data["nano_stacks"][valid_epochs][:, channels], dtype=float)
        counts = counts_all[valid_epochs]

    finite_fraction = np.mean(np.isfinite(raw), axis=(0, 2))
    dynamic = np.nanmedian(np.nanstd(raw, axis=-1), axis=0)
    keep = (finite_fraction >= 0.99) & np.isfinite(dynamic) & (dynamic > 0)
    if np.sum(keep) < MIN_VALID_CHANNELS:
        raise RuntimeError(f"Only {np.sum(keep)} Nano channels pass finite/dynamic QC")
    raw = np.nan_to_num(raw[:, keep], copy=False)
    channels = channels[keep]

    sos = butter(4, BAND_HZ, btype="bandpass", fs=fs, output="sos")
    filtered = sosfiltfilt(sos, raw, axis=-1)
    absolute_x = channels.astype(float) * dx
    relative_x = absolute_x - absolute_x[0]
    propagation_time = relative_x / REFERENCE_SPEED_MPS
    q = np.arange(
        int(round((EXTRACT_Q_S[1] - EXTRACT_Q_S[0]) * fs)), dtype=float
    ) / fs + EXTRACT_Q_S[0]
    sample_time = np.arange(filtered.shape[-1], dtype=float) / fs - PRE_S
    aligned = np.empty((filtered.shape[0], filtered.shape[1], q.size), dtype=np.float32)
    for iepoch in range(filtered.shape[0]):
        for ichannel, tau in enumerate(propagation_time):
            query = REFERENCE_INTERCEPT_S + tau + q
            aligned[iepoch, ichannel] = np.interp(
                query, sample_time, filtered[iepoch, ichannel], left=0.0, right=0.0
            )
    return aligned, counts, valid_epochs, propagation_time, fs


def _inject_moveout(gather: np.ndarray, q: np.ndarray, tau: np.ndarray, dvv: float) -> np.ndarray:
    """Inject a fractional apparent-speed change by sub-sample channel shifts.

    For ``dvv > 0``, the propagation time shortens by approximately
    ``-dvv * tau``.  Sampling the original trace at ``q + dvv*tau`` places its
    feature at that earlier time.  No injected value is passed to recovery.
    """
    shifted = np.empty_like(gather, dtype=np.float32)
    for ichannel, travel_time in enumerate(tau):
        shifted[ichannel] = np.interp(
            q + dvv * travel_time, q, gather[ichannel], left=0.0, right=0.0
        )
    return shifted


def inject() -> None:
    aligned, counts, epochs, tau, fs = _extract_aligned_gathers()
    q = np.arange(aligned.shape[-1], dtype=float) / fs + EXTRACT_Q_S[0]
    total_weight = float(np.sum(counts))
    weighted_sum = np.tensordot(counts.astype(float), aligned, axes=(0, 0))
    references = np.empty_like(aligned)
    for i, weight in enumerate(counts):
        denominator = total_weight - float(weight)
        if denominator <= 0:
            raise RuntimeError("At least two populated Nano bursts are required")
        references[i] = (weighted_sum - float(weight) * aligned[i]) / denominator

    rng = np.random.default_rng(SEED)
    specifications = [(i, float(dvv)) for i in range(aligned.shape[0]) for dvv in INJECTED_DVV]
    rng.shuffle(specifications)
    ntrial = len(specifications)
    trials = np.empty((ntrial, aligned.shape[1], aligned.shape[2]), dtype=np.float32)
    reference_index = np.empty(ntrial, dtype=np.int16)
    ids = []
    truth_rows = []
    for itrial, (target_index, dvv) in enumerate(specifications):
        trial_id = f"N{rng.integers(0, 2**63, dtype=np.int64):016x}"
        ids.append(trial_id)
        reference_index[itrial] = target_index
        trials[itrial] = _inject_moveout(aligned[target_index], q, tau, dvv)
        truth_rows.append(
            {
                "trial_id": trial_id,
                "target_epoch": int(epochs[target_index]),
                "injected_dvv": f"{dvv:.9g}",
            }
        )

    np.savez_compressed(
        BLIND,
        trial_id=np.asarray(ids),
        trial_gather=trials,
        reference_gather=references,
        reference_index=reference_index,
        q_s=q,
        propagation_time_s=tau,
        fs=fs,
        band_hz=np.asarray(BAND_HZ),
        aperture_m=np.asarray(APERTURE_M),
        reference_speed_mps=REFERENCE_SPEED_MPS,
        reference_intercept_s=REFERENCE_INTERCEPT_S,
        channel_stride=CHANNEL_STRIDE,
        seed=SEED,
    )
    _write_csv(TRUTH, truth_rows)
    print(f"Injection stage wrote {ntrial} randomized trials to {BLIND}")
    print(f"Sealed truth table written separately to {TRUTH}")


def _normalized_correlation_lag(
    reference: np.ndarray, target: np.ndarray, fs: float, max_lag_s: float
) -> tuple[float, float]:
    """Estimate target-minus-reference lag with parabolic sub-sample refinement."""
    max_lag = max(1, int(round(max_lag_s * fs)))
    lags = np.arange(-max_lag, max_lag + 1, dtype=int)
    values = np.full(lags.size, np.nan, dtype=float)
    for ilag, lag in enumerate(lags):
        if lag < 0:
            a, b = target[:lag], reference[-lag:]
        elif lag > 0:
            a, b = target[lag:], reference[:-lag]
        else:
            a, b = target, reference
        a = a - np.mean(a)
        b = b - np.mean(b)
        denominator = np.sqrt(np.dot(a, a) * np.dot(b, b))
        if denominator > 0:
            values[ilag] = np.dot(a, b) / denominator
    if not np.any(np.isfinite(values)):
        return np.nan, np.nan
    peak = int(np.nanargmax(values))
    refined = float(lags[peak])
    if 0 < peak < values.size - 1 and np.all(np.isfinite(values[peak - 1:peak + 2])):
        left, center, right = values[peak - 1:peak + 2]
        curvature = left - 2.0 * center + right
        if curvature < -1e-12:
            refined += 0.5 * (left - right) / curvature
    return refined / fs, float(values[peak])


def _robust_line(x: np.ndarray, y: np.ndarray, base_weight: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Huber-IRLS line fit, returning coefficients, standard errors, and robust scale."""
    design = np.column_stack((np.ones_like(x), x))
    weights = np.asarray(base_weight, dtype=float).copy()
    beta = np.zeros(2, dtype=float)
    for _ in range(20):
        lhs = design.T @ (weights[:, None] * design)
        rhs = design.T @ (weights * y)
        updated = np.linalg.lstsq(lhs, rhs, rcond=None)[0]
        residual = y - design @ updated
        center = np.median(residual)
        scale = 1.4826 * np.median(np.abs(residual - center))
        if not np.isfinite(scale) or scale < 1e-8:
            scale = max(float(np.sqrt(np.average(residual**2, weights=base_weight))), 1e-8)
        huber = np.ones_like(residual)
        large = np.abs(residual) > 1.5 * scale
        huber[large] = 1.5 * scale / np.abs(residual[large])
        new_weights = base_weight * huber
        if np.max(np.abs(updated - beta)) < 1e-10:
            beta, weights = updated, new_weights
            break
        beta, weights = updated, new_weights
    residual = y - design @ beta
    dof = max(1, x.size - 2)
    sigma2 = float(np.sum(weights * residual**2) / dof)
    covariance = sigma2 * np.linalg.pinv(design.T @ (weights[:, None] * design))
    return beta, np.sqrt(np.maximum(np.diag(covariance), 0.0)), float(scale)


def recover() -> None:
    """Recover trial perturbations without opening or accepting the truth table."""
    with np.load(BLIND) as blind:
        trial_ids = blind["trial_id"].astype(str)
        trials = blind["trial_gather"]
        references = blind["reference_gather"]
        reference_index = blind["reference_index"]
        q = blind["q_s"]
        tau = blind["propagation_time_s"]
        fs = float(blind["fs"])
    use_q = (q >= RECOVERY_Q_S[0]) & (q <= RECOVERY_Q_S[1])
    rows = []
    for itrial, trial_id in enumerate(trial_ids):
        reference = references[int(reference_index[itrial])][:, use_q]
        target = trials[itrial][:, use_q]
        lag = np.full(tau.size, np.nan)
        correlation = np.full(tau.size, np.nan)
        for channel in range(tau.size):
            lag[channel], correlation[channel] = _normalized_correlation_lag(
                reference[channel], target[channel], fs, MAX_LAG_S
            )
        good = (
            np.isfinite(lag)
            & np.isfinite(correlation)
            & (correlation >= MIN_CHANNEL_CORRELATION)
        )
        if np.sum(good) >= MIN_VALID_CHANNELS:
            weight = np.maximum(correlation[good], 0.0) ** 2
            beta, standard_error, residual_scale = _robust_line(
                tau[good], lag[good], weight
            )
            estimate = -float(beta[1])
            estimate_se = float(standard_error[1])
            intercept = float(beta[0])
            median_corr = float(np.median(correlation[good]))
        else:
            estimate = estimate_se = intercept = residual_scale = median_corr = np.nan
        rows.append(
            {
                "trial_id": trial_id,
                "estimated_dvv": f"{estimate:.10g}",
                "estimated_dvv_se": f"{estimate_se:.10g}",
                "n_channels": int(np.sum(good)),
                "nuisance_intercept_s": f"{intercept:.10g}",
                "residual_scale_s": f"{residual_scale:.10g}",
                "median_channel_correlation": f"{median_corr:.10g}",
            }
        )
    _write_csv(RECOVERY, rows)
    print(f"Recovery stage estimated {len(rows)} blinded trials and wrote {RECOVERY}")


def _level_row(truth: float, estimated: np.ndarray, threshold: float) -> dict:
    finite = np.isfinite(estimated)
    values = estimated[finite]
    error = values - truth
    detected = np.abs(values) > threshold
    if truth == 0:
        correct_detection = detected
    else:
        correct_detection = detected & (np.sign(values) == np.sign(truth))
    median = float(np.median(values)) if values.size else np.nan
    mad_scatter = float(1.4826 * np.median(np.abs(values - median))) if values.size else np.nan
    return {
        "injected_dvv": truth,
        "n_trials": int(estimated.size),
        "n_recovered": int(values.size),
        "mean_estimated_dvv": float(np.mean(values)) if values.size else np.nan,
        "median_estimated_dvv": median,
        "bias": float(np.mean(error)) if values.size else np.nan,
        "standard_scatter": float(np.std(values, ddof=1)) if values.size > 1 else np.nan,
        "robust_scatter_1p4826_mad": mad_scatter,
        "rmse": float(np.sqrt(np.mean(error**2))) if values.size else np.nan,
        "detection_probability": float(np.mean(detected)) if values.size else np.nan,
        "correct_sign_detection_probability": float(np.mean(correct_detection)) if values.size else np.nan,
    }


def summarize() -> None:
    truth_rows = {row["trial_id"]: row for row in _read_csv(TRUTH)}
    recovery_rows = _read_csv(RECOVERY)
    joined = []
    for row in recovery_rows:
        trial_id = row["trial_id"]
        if trial_id not in truth_rows:
            raise RuntimeError(f"Recovery trial absent from truth table: {trial_id}")
        joined.append(
            (
                trial_id,
                float(truth_rows[trial_id]["injected_dvv"]),
                float(row["estimated_dvv"]),
                float(row["estimated_dvv_se"]),
                int(row["n_channels"]),
            )
        )
    trial_id = np.asarray([row[0] for row in joined])
    injected = np.asarray([row[1] for row in joined])
    estimated_raw = np.asarray([row[2] for row in joined])
    estimated_se = np.asarray([row[3] for row in joined])
    n_channels = np.asarray([row[4] for row in joined])
    null = estimated_raw[injected == 0]
    if np.sum(np.isfinite(null)) < 20:
        raise RuntimeError("Fewer than 20 finite zero-injection trials; null is inadequate")
    null_center = float(np.nanmedian(null))
    estimated = estimated_raw - null_center
    null_calibrated = estimated[injected == 0]
    threshold = float(np.nanquantile(np.abs(null_calibrated), 0.95))
    false_positive_rate = float(np.mean(np.abs(null_calibrated) > threshold))

    summary_rows = []
    levels = np.unique(injected)
    for level in levels:
        summary_rows.append(_level_row(float(level), estimated[injected == level], threshold))
    _write_csv(SUMMARY, summary_rows)

    positive = {row["injected_dvv"]: row for row in summary_rows if row["injected_dvv"] > 0}
    negative = {abs(row["injected_dvv"]): row for row in summary_rows if row["injected_dvv"] < 0}
    common_amplitudes = sorted(set(positive) & set(negative))
    qualifying = [
        value for value in common_amplitudes
        if positive[value]["correct_sign_detection_probability"] >= 0.95
        and negative[value]["correct_sign_detection_probability"] >= 0.95
    ]
    detection_limit = float(min(qualifying)) if qualifying else np.nan

    np.savez(
        RESULTS,
        trial_id=trial_id,
        injected_dvv=injected,
        estimated_dvv_raw=estimated_raw,
        estimated_dvv=estimated,
        estimated_dvv_se=estimated_se,
        n_channels=n_channels,
        null_center_dvv=null_center,
        empirical_95pct_threshold_dvv=threshold,
        empirical_false_positive_rate=false_positive_rate,
        symmetric_95pct_detection_limit_dvv=detection_limit,
        band_hz=np.asarray(BAND_HZ),
        aperture_m=np.asarray(APERTURE_M),
        reference_speed_mps=REFERENCE_SPEED_MPS,
        reference_intercept_s=REFERENCE_INTERCEPT_S,
        injected_levels=INJECTED_DVV,
        seed=SEED,
    )

    medians = np.asarray([row["median_estimated_dvv"] for row in summary_rows])
    robust_scatter = np.asarray([row["robust_scatter_1p4826_mad"] for row in summary_rows])
    bias = np.asarray([row["bias"] for row in summary_rows])
    probability = np.asarray([row["correct_sign_detection_probability"] for row in summary_rows])
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 7.2), constrained_layout=True)
    ax = axes[0, 0]
    ax.scatter(injected, estimated, s=8, color="0.65", alpha=0.35, linewidths=0)
    ax.errorbar(levels, medians, yerr=robust_scatter, fmt="o", color="#2166ac", capsize=2,
                label="median $\\pm$ 1.4826 MAD")
    limit = 1.05 * np.nanmax(np.abs(INJECTED_DVV))
    ax.plot([-limit, limit], [-limit, limit], color="0.15", ls="--", lw=1, label="one-to-one")
    ax.set(xlabel="injected fractional apparent-velocity change",
           ylabel="recovered fractional apparent-velocity change", title="A  Blind recovery")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    ax.plot(levels, bias, "o-", color="#b2182b", label="mean bias")
    ax.plot(levels, robust_scatter, "s-", color="#2166ac", label="robust scatter")
    ax.axhline(0, color="0.25", lw=0.8)
    ax.set(xlabel="injected fractional apparent-velocity change", ylabel="fractional change",
           title="B  Bias and inter-burst scatter")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    colors = np.where(levels < 0, "#4393c3", np.where(levels > 0, "#d6604d", "0.25"))
    ax.scatter(np.abs(levels), probability, c=colors, s=32)
    ax.axhline(0.95, color="0.2", ls="--", lw=1)
    if np.isfinite(detection_limit):
        ax.axvline(detection_limit, color="#7b3294", ls=":", lw=1.3,
                   label=f"two-direction 95% limit = {detection_limit:.1e}")
    ax.set_xscale("symlog", linthresh=8e-5)
    ax.set(xlabel="absolute injected fractional change", ylabel="correct-sign detection probability",
           ylim=(-0.04, 1.04), title="C  Empirical detection probability")
    if np.isfinite(detection_limit):
        ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    ax.hist(null_calibrated, bins=14, color="0.65", edgecolor="white")
    ax.axvline(-threshold, color="#b2182b", ls="--", lw=1.2)
    ax.axvline(threshold, color="#b2182b", ls="--", lw=1.2,
               label=f"empirical 95% threshold = {threshold:.1e}")
    ax.set(xlabel="recovered change for zero injection", ylabel="number of trials",
           title=f"D  Null trials; false-positive rate = {false_positive_rate:.3f}")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Nano 30--60 Hz apparent-moveout blind injection--recovery", fontsize=12)
    fig.savefig(FIGURE, dpi=300)
    plt.close(fig)

    report = (
        "SAFOD AWD Nano phase-neutral apparent-moveout injection--recovery\n"
        f"Input: {STACKS.name}\n"
        f"Band: {BAND_HZ[0]:.0f}-{BAND_HZ[1]:.0f} Hz; aperture: "
        f"{APERTURE_M[0]:.0f}-{APERTURE_M[1]:.0f} m along-fiber coordinate\n"
        f"Fixed reference trajectory: {REFERENCE_SPEED_MPS:.1f} m/s; "
        f"intercept {REFERENCE_INTERCEPT_S:.3f} s relative to p26.cc9.txt UTC_Date\n"
        f"Trials: {injected.size}; zero-injection trials: {null.size}\n"
        f"Null median removed as calibration: {null_center:.6g}\n"
        f"Empirical two-sided 95% null threshold: {threshold:.6g}\n"
        f"Empirical false-positive rate: {false_positive_rate:.6f}\n"
        f"Smallest tested magnitude with >=95% correct-sign detection in both directions: "
        f"{detection_limit if np.isfinite(detection_limit) else 'not reached'}\n"
        "Interpretation: detection applies only to fractional change in the accepted Nano "
        "apparent-moveout observable. It is not yet a formation Vp or direct-P dv/v limit.\n"
    )
    REPORT.write_text(report)
    print(report, end="")
    print(f"Wrote {SUMMARY}, {RESULTS}, {FIGURE}, and {REPORT}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", choices=("inject", "recover", "summarize", "all"), default="all"
    )
    args = parser.parse_args()
    if args.stage in ("inject", "all"):
        inject()
    if args.stage in ("recover", "all"):
        recover()
    if args.stage in ("summarize", "all"):
        summarize()


if __name__ == "__main__":
    main()
