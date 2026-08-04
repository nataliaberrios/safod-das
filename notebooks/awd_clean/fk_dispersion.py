"""Signed frequency-slowness and dispersion analysis for SAFOD AWD DAS."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, detrend, sosfiltfilt, windows


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
PRE_S = 0.5
TIME_WINDOW = (-0.08, 0.45)
FREQ_RANGE = (3.0, 100.0)
SLOWNESS = np.linspace(-2.0e-3, 2.0e-3, 801)  # signed s/m; |v| >= 500 m/s
TURNAROUND_CH = 1702  # data-derived provisional Deep hairpin turnaround


def weighted_stack(stacks, counts):
    """Count-weighted mean that does not let isolated NaNs poison a channel."""
    good = counts > 0
    total = np.zeros(stacks.shape[1:], dtype=np.float64)
    weight_sum = np.zeros(stacks.shape[1:], dtype=np.float64)
    for stack, weight in zip(stacks[good], counts[good].astype(float)):
        finite = np.isfinite(stack)
        total[finite] += weight * stack[finite]
        weight_sum[finite] += weight
    return np.divide(
        total, weight_sum, out=np.full_like(total, np.nan), where=weight_sum > 0
    )


def channel_qc(section, coordinate):
    """Remove unusable traces while retaining their true spatial coordinates."""
    finite_fraction = np.mean(np.isfinite(section), axis=1)
    rms = np.sqrt(np.nanmean(section ** 2, axis=1))
    positive = rms[np.isfinite(rms) & (rms > 0)]
    floor = np.nanmedian(positive) * 1e-6 if positive.size else np.inf
    good = (finite_fraction == 1.0) & np.isfinite(rms) & (rms > floor)
    return section[good], coordinate[good], good, finite_fraction, rms


def preprocess(section, fs):
    # Retain a broad band; the frequency axis is resolved by the temporal FFT.
    sos = butter(4, [2.0, 120.0], btype="bandpass", fs=fs, output="sos")
    section = sosfiltfilt(sos, section, axis=-1)
    section = detrend(section, axis=-1, type="linear")
    section *= windows.tukey(section.shape[0], 0.15)[:, None]
    section *= windows.tukey(section.shape[1], 0.15)[None, :]
    return section


def frequency_slowness(section, coordinate, fs):
    """Normalized phase-shift beam power B(f,p), preserving slowness sign."""
    section = preprocess(section.astype(float), fs)
    nfft = int(2 ** np.ceil(np.log2(max(1024, section.shape[1]))))
    spectrum = np.fft.rfft(section, n=nfft, axis=1)
    frequency = np.fft.rfftfreq(nfft, 1 / fs)
    keep = (frequency >= FREQ_RANGE[0]) & (frequency <= FREQ_RANGE[1])
    frequency = frequency[keep]
    spectrum = spectrum[:, keep]

    # Normalize each channel/frequency to unit phase. This estimates phase
    # coherence rather than allowing a few high-amplitude channels to dominate.
    unit = spectrum / (np.abs(spectrum) + 1e-30)
    power = np.empty((frequency.size, SLOWNESS.size), np.float32)
    x = coordinate - coordinate.mean()
    for i, freq in enumerate(frequency):
        steering = np.exp(2j * np.pi * freq * np.outer(SLOWNESS, x))
        beam = steering @ unit[:, i]
        power[i] = (np.abs(beam) / section.shape[0]) ** 2
    return frequency, power


def ridge(frequency, power, minimum_power=0.05):
    index = np.argmax(power, axis=1)
    p = SLOWNESS[index]
    strength = power[np.arange(len(frequency)), index]
    p[strength < minimum_power] = np.nan
    velocity = np.divide(1.0, p, out=np.full_like(p, np.nan), where=np.abs(p) > 1e-8)
    return p, velocity, strength


def segment(section, dx, start_m, stop_m, reverse=False):
    c0 = max(0, int(start_m / dx))
    c1 = min(section.shape[0], int(stop_m / dx) + 1)
    result = section[c0:c1]
    if reverse:
        result = result[::-1]
    coordinate = np.arange(result.shape[0]) * dx
    return result, coordinate, (c0, c1)


def main():
    data = np.load(STACKS)
    counts = data["n_common"]
    fs = float(data["fs"])
    nano = weighted_stack(data["nano_stacks"], counts)
    deep = weighted_stack(data["deep_stacks"], counts)
    dxn, dxd = float(data["dx_nano"]), float(data["dx_deep"])
    i0 = int(round((PRE_S + TIME_WINDOW[0]) * fs))
    i1 = int(round((PRE_S + TIME_WINDOW[1]) * fs))

    nano_match, xn, nano_bounds = segment(nano[:, i0:i1], dxn, 80, 440)
    deep_match, xd, deep_bounds = segment(deep[:, i0:i1], dxd, 80, 440)

    # Analyze Deep legs in a common surface-to-turnaround orientation.  The
    # turnaround channel is provisional; no true-depth label is assigned.
    deep_down = deep[:TURNAROUND_CH, i0:i1]
    x_down = np.arange(deep_down.shape[0]) * dxd
    deep_return = deep[TURNAROUND_CH:, i0:i1][::-1]
    x_return = np.arange(deep_return.shape[0]) * dxd

    definitions = [
        ("Nano 80-440 m", nano_match, xn),
        ("Deep 80-440 m", deep_match, xd),
        ("Deep outbound leg", deep_down, x_down),
        ("Deep return leg (reversed)", deep_return, x_return),
    ]
    products = {}
    summaries = []
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), sharex=True, sharey=True)
    for axis, (name, section, coordinate) in zip(axes.flat, definitions):
        n_input = section.shape[0]
        section, coordinate, good_channels, finite_fraction, channel_rms = channel_qc(
            section, coordinate
        )
        if section.shape[0] < 8:
            raise RuntimeError(f"{name}: only {section.shape[0]} usable channels after QC")
        frequency, power = frequency_slowness(section, coordinate, fs)
        p_ridge, velocity, strength = ridge(frequency, power)
        key = name.lower().replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")
        products[f"{key}_power"] = power
        products[f"{key}_ridge_slowness"] = p_ridge
        products[f"{key}_ridge_velocity"] = velocity
        products[f"{key}_ridge_strength"] = strength
        products[f"{key}_good_channels"] = good_channels
        products[f"{key}_finite_fraction"] = finite_fraction
        products[f"{key}_channel_rms"] = channel_rms
        products[f"{key}_coordinate"] = coordinate

        image = axis.pcolormesh(
            SLOWNESS * 1e3, frequency, power, shading="auto", cmap="magma",
            vmin=0, vmax=max(0.25, np.nanpercentile(power, 99.7)),
        )
        axis.plot(p_ridge * 1e3, frequency, "c.", ms=2, alpha=0.8)
        axis.axvline(0, color="w", lw=0.7, alpha=0.7)
        for speed in (1000, 1500, 3000, 5000):
            for sign in (-1, 1):
                axis.axvline(sign * 1e3 / speed, color="w", ls=":", lw=0.45, alpha=0.35)
        axis.set_title(f"{name}\n{section.shape[0]}/{n_input} channels pass QC")
        axis.set_xlabel("signed slowness (ms/m)")
        axis.set_ylabel("frequency (Hz)")
        fig.colorbar(image, ax=axis, label="normalized phase-beam power")

        band = (frequency >= 20) & (frequency <= 60) & np.isfinite(velocity)
        strong = band & (strength >= np.nanpercentile(strength[band], 60) if band.any() else False)
        median_v = np.nanmedian(velocity[strong]) if np.any(strong) else np.nan
        median_p = np.nanmedian(p_ridge[strong]) if np.any(strong) else np.nan
        summaries.append((name, n_input, section.shape[0], median_p, median_v,
                          np.nanmedian(strength[band]) if band.any() else np.nan))

    fig.suptitle(
        "SAFOD AWD signed frequency-slowness dispersion\n"
        f"stacked common drops, time {TIME_WINDOW[0]:.2f} to {TIME_WINDOW[1]:.2f} s relative to catalog pick"
    )
    fig.tight_layout()
    fig.savefig(HERE / "fk_dispersion.png", dpi=180)
    np.savez(
        HERE / "fk_dispersion.npz", frequency=frequency, slowness=SLOWNESS,
        time_window=np.asarray(TIME_WINDOW), turnaround_ch=TURNAROUND_CH,
        nano_bounds=np.asarray(nano_bounds), deep_bounds=np.asarray(deep_bounds),
        **products,
    )

    lines = [
        "SAFOD AWD signed frequency-slowness summary",
        "Positive/negative signs refer to increasing/decreasing coordinate in each plotted orientation.",
        "Deep return was reversed into a provisional surface-to-turnaround orientation.",
        "20-60 Hz ridge summaries use the stronger 40% of frequency bins and are descriptive, not phase labels.",
        "",
    ]
    for name, n_input, nch, medp, medv, meds in summaries:
        lines.append(
            f"{name}: QC channels={nch}/{n_input}, median signed p={medp*1e3:.5f} ms/m, "
            f"median signed v={medv:.1f} m/s, median ridge power={meds:.4f}"
        )
    report = "\n".join(lines) + "\n"
    (HERE / "fk_dispersion.txt").write_text(report)
    print(report)
    print("Saved fk_dispersion.png/.npz/.txt")


if __name__ == "__main__":
    main()
