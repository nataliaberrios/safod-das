"""Null tests and visual inspection for the strongest phase-neutral AWD ridge."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
SCAN = HERE / "mode_scan.npz"
BAND = (30.0, 60.0)
APERTURE = (80.0, 440.0)
PRE_S = 0.5
WIN_S = 0.040
N_PERM = 499
SEED = 20260801


def bandpass(x, fs):
    sos = butter(4, BAND, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def weighted_stack(stacks, counts):
    good = counts > 0
    w = counts[good].astype(float)
    return np.tensordot(w, stacks[good], axes=(0, 0)) / w.sum()


def score(section, coordinate, fs, velocity, intercept):
    ns = int(WIN_S * fs)
    idx = int(PRE_S * fs) + np.rint((intercept + coordinate / velocity) * fs).astype(int)
    if idx.min() < 0 or idx.max() + ns > section.shape[1]:
        return np.nan
    gather = np.stack([section[i, j:j + ns] for i, j in enumerate(idx)])
    den = gather.shape[0] * np.sum(gather ** 2)
    return float(np.sum(np.sum(gather, axis=0) ** 2) / den) if den else np.nan


def local_max(section, coordinate, fs, velocities, intercepts):
    best = (-np.inf, np.nan, np.nan)
    for velocity in velocities:
        for intercept in intercepts:
            value = score(section, coordinate, fs, velocity, intercept)
            if np.isfinite(value) and value > best[0]:
                best = (value, velocity, intercept)
    return best


def main():
    data = np.load(STACKS)
    scan = np.load(SCAN)
    fs = float(data["fs"])
    dx = float(data["dx_nano"])
    counts = data["n_common"]
    raw = weighted_stack(data["nano_stacks"], counts)
    filtered = bandpass(raw, fs)

    c0, c1 = int(APERTURE[0] / dx), min(int(APERTURE[1] / dx) + 1, raw.shape[0])
    section = filtered[c0:c1]
    absolute_x = np.arange(c0, c1) * dx
    relative_x = absolute_x - absolute_x[0]

    surface = scan["nano_30_60"]
    iv, it = np.unravel_index(np.argmax(surface), surface.shape)
    v_ref = float(scan["v_grid"][iv])
    t_ref = float(scan["t0_grid"][it])

    # Coarse local look-elsewhere-corrected search, applied identically to data
    # and every depth permutation.
    velocities = np.arange(v_ref - 500, v_ref + 501, 100.0)
    intercepts = np.arange(t_ref - 0.040, t_ref + 0.0401, 0.004)
    stride = 4
    sec_test = section[::stride]
    x_test = relative_x[::stride]
    observed, v_obs, t_obs = local_max(sec_test, x_test, fs, velocities, intercepts)

    rng = np.random.default_rng(SEED)
    permuted = np.empty(N_PERM)
    for i in range(N_PERM):
        permuted[i] = local_max(
            sec_test[rng.permutation(sec_test.shape[0])], x_test, fs,
            velocities, intercepts,
        )[0]
    p_perm = (1 + np.sum(permuted >= observed)) / (N_PERM + 1)

    # Negative-time control: run the same velocity search over intercepts whose
    # predicted windows remain before the GPS source time over the aperture.
    pre_intercepts = np.arange(-0.45, -0.16, 0.004)
    pre_max, pre_v, pre_t = local_max(sec_test, x_test, fs, velocities, pre_intercepts)

    # Inspect unnormalized raw and filtered sections, plus the moveout-aligned
    # filtered gather. No per-channel normalization is used in the first two.
    t = np.arange(raw.shape[1]) / fs - PRE_S
    raw_ap = raw[c0:c1]
    clip_raw = np.percentile(np.abs(raw_ap[:, (t >= -0.1) & (t <= 0.45)]), 99.5)
    clip_bp = np.percentile(np.abs(section[:, (t >= -0.1) & (t <= 0.45)]), 99.5)
    ridge = t_obs + relative_x / v_obs
    shifts = np.rint(ridge * fs).astype(int)
    aligned = np.stack([
        section[i, int(PRE_S * fs) + shift - 40:int(PRE_S * fs) + shift + 120]
        for i, shift in enumerate(shifts)
    ])
    aligned_mean = aligned.mean(axis=0)
    aligned_time = (np.arange(aligned.shape[1]) - 40) / fs

    fig, ax = plt.subplots(1, 3, figsize=(16, 6))
    for axis, image_data, clip, title in (
        (ax[0], raw_ap, clip_raw, "A  Unnormalized 1--100 Hz stack"),
        (ax[1], section, clip_bp, "B  Unnormalized 30--60 Hz stack"),
    ):
        axis.imshow(
            image_data, aspect="auto", origin="upper", cmap="RdBu_r",
            extent=[t[0], t[-1], absolute_x[-1], absolute_x[0]],
            vmin=-clip, vmax=clip,
        )
        axis.plot(ridge, absolute_x, "k--", lw=1.5, label=f"neutral ridge: {v_obs:.0f} m/s")
        axis.axvline(0, color="gold", lw=1, label="GPS source time")
        axis.set_xlim(-0.1, 0.45)
        axis.set_xlabel("time relative to GPS drop (s)")
        axis.set_ylabel("distance along Nano fiber (m)")
        axis.set_title(title)
        axis.legend(fontsize=8)
    ax[2].imshow(
        aligned, aspect="auto", origin="upper", cmap="RdBu_r",
        extent=[aligned_time[0], aligned_time[-1], absolute_x[-1], absolute_x[0]],
        vmin=-clip_bp, vmax=clip_bp,
    )
    ax[2].plot(aligned_time, absolute_x[0] +
               (aligned_mean / (np.max(np.abs(aligned_mean)) + 1e-30)) * 25,
               color="k", lw=1, alpha=0.8, label="mean waveform (scaled)")
    ax[2].axvline(0, color="gold", lw=1)
    ax[2].set_xlabel("time after moveout alignment (s)")
    ax[2].set_ylabel("distance along Nano fiber (m)")
    ax[2].set_title("C  Ridge-aligned 30--60 Hz gather")
    ax[2].legend(fontsize=8)
    fig.suptitle("SAFOD AWD Nano ridge inspection; amplitudes not normalized by channel")
    fig.tight_layout()
    fig.savefig(HERE / "record_section_inspection.png", dpi=180)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].hist(permuted, bins=30, color="0.65", edgecolor="white")
    ax[0].axvline(observed, color="crimson", lw=2, label=f"observed={observed:.3f}")
    ax[0].set_xlabel("maximum local semblance")
    ax[0].set_ylabel("depth permutations")
    ax[0].set_title(f"A  Depth-permutation null, p={p_perm:.4f}")
    ax[0].legend()
    ax[1].plot(scan["v_grid"], surface.max(axis=1), color="navy")
    ax[1].axvline(v_ref, color="crimson", ls="--")
    ax[1].set_xlabel("unsigned apparent speed (m/s)")
    ax[1].set_ylabel("maximum semblance over intercept")
    ax[1].set_title(
        "B  Observed scan profile\n"
        f"pre-source max={pre_max:.3f} at {pre_v:.0f} m/s"
    )
    fig.tight_layout()
    fig.savefig(HERE / "null_tests.png", dpi=180)

    np.savez(
        HERE / "null_tests.npz", observed=observed, v_observed=v_obs,
        t_observed=t_obs, permuted=permuted, p_permutation=p_perm,
        pre_source_max=pre_max, pre_source_v=pre_v, pre_source_t=pre_t,
        band=np.asarray(BAND), aperture=np.asarray(APERTURE), seed=SEED,
    )
    report = (
        "SAFOD AWD phase-neutral null tests\n"
        f"band: {BAND[0]:.0f}-{BAND[1]:.0f} Hz\n"
        f"Nano fiber-coordinate aperture: {APERTURE[0]:.0f}-{APERTURE[1]:.0f} m\n"
        f"observed coarse local-max semblance: {observed:.6f}\n"
        f"observed speed/intercept: {v_obs:.1f} m/s, {t_obs*1e3:.1f} ms\n"
        f"depth permutations: {N_PERM}; seed: {SEED}\n"
        f"permutation p-value (local look-elsewhere corrected): {p_perm:.6f}\n"
        f"largest pre-source control: {pre_max:.6f}\n"
        f"pre-source speed/intercept: {pre_v:.1f} m/s, {pre_t*1e3:.1f} ms\n"
    )
    (HERE / "null_test_report.txt").write_text(report)
    print(report)
    print("Saved record_section_inspection.png, null_tests.png/.npz, and report")


if __name__ == "__main__":
    main()
