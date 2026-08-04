"""Sliding-aperture signed f-p test for depth-varying Deep AWD modes."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, detrend, sosfiltfilt, windows

from fk_dispersion import weighted_stack


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
PRE_S = 0.5
TIME_WINDOW = (-0.08, 2.85)
TURNAROUND_CH = 1702
APERTURE_M = 400.0
STEP_M = 200.0
SLOWNESS = np.linspace(-1.2e-3, 1.2e-3, 481)
BANDS = ((3.0, 15.0), (15.0, 30.0), (30.0, 60.0))


def preprocess(section, fs):
    section = np.nan_to_num(section.astype(float))
    sos = butter(4, [2.0, 80.0], btype="bandpass", fs=fs, output="sos")
    section = sosfiltfilt(sos, section, axis=-1)
    section = detrend(section, axis=-1, type="linear")
    section *= windows.tukey(section.shape[0], 0.15)[:, None]
    section *= windows.tukey(section.shape[1], 0.10)[None, :]
    return section


def band_beam(section, coordinate, fs, band):
    section = preprocess(section, fs)
    nfft = int(2 ** np.ceil(np.log2(section.shape[1])))
    spectrum = np.fft.rfft(section, n=nfft, axis=1)
    frequency = np.fft.rfftfreq(nfft, 1 / fs)
    use = np.flatnonzero((frequency >= band[0]) & (frequency <= band[1]))[::2]
    x = coordinate - coordinate.mean()
    result = np.zeros(SLOWNESS.size)
    for index in use:
        amplitude = np.abs(spectrum[:, index])
        valid = amplitude > np.nanmedian(amplitude) * 1e-6
        if valid.sum() < 8:
            continue
        unit = spectrum[valid, index] / amplitude[valid]
        steering = np.exp(2j * np.pi * frequency[index] * np.outer(SLOWNESS, x[valid]))
        result += (np.abs(steering @ unit) / valid.sum()) ** 2
    return result / max(1, use.size)


def sliding_leg(section, dx, fs):
    aperture = max(16, int(round(APERTURE_M / dx)))
    step = max(8, int(round(STEP_M / dx)))
    starts = np.arange(0, section.shape[0] - aperture + 1, step)
    centers = (starts + 0.5 * (aperture - 1)) * dx
    power = np.empty((len(BANDS), starts.size, SLOWNESS.size), np.float32)
    for j, start in enumerate(starts):
        stop = start + aperture
        x = np.arange(start, stop) * dx
        for ib, band in enumerate(BANDS):
            power[ib, j] = band_beam(section[start:stop], x, fs, band)
    return centers, power


def main():
    with np.load(STACKS) as data:
        fs = float(data["fs"])
        dx = float(data["dx_deep"])
        deep = weighted_stack(data["deep_stacks"], data["n_common"])
    i0 = int(round((PRE_S + TIME_WINDOW[0]) * fs))
    i1 = int(round((PRE_S + TIME_WINDOW[1]) * fs))
    outbound = deep[:TURNAROUND_CH, i0:i1]
    returned = deep[TURNAROUND_CH:, i0:i1][::-1]
    legs = (("Outbound", outbound), ("Return (reversed)", returned))

    products = {}
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey="row")
    for row, (name, section) in enumerate(legs):
        centers, power = sliding_leg(section, dx, fs)
        products[f"{name.split()[0].lower()}_centers_m"] = centers
        products[f"{name.split()[0].lower()}_power"] = power
        for col, band in enumerate(BANDS):
            ax = axes[row, col]
            image = ax.pcolormesh(
                SLOWNESS * 1e3, centers, power[col], shading="auto", cmap="magma"
            )
            ax.axvline(0, color="w", lw=0.6)
            for speed, color in ((1500, "cyan"), (3000, "lime")):
                for sign in (-1, 1):
                    ax.axvline(sign * 1e3 / speed, color=color, ls=":", lw=0.8)
            ax.set_title(f"{name}: {band[0]:.0f}–{band[1]:.0f} Hz")
            ax.set_xlabel("signed slowness (ms/m)")
            ax.set_ylabel("distance from leg start (m)")
            fig.colorbar(image, ax=ax, label="mean normalized beam power")
    fig.suptitle(
        "SAFOD AWD Deep DAS: sliding-aperture signed frequency–slowness\n"
        f"{APERTURE_M:.0f} m apertures every {STEP_M:.0f} m; time {TIME_WINDOW[0]:.2f} to {TIME_WINDOW[1]:.2f} s"
    )
    fig.tight_layout()
    fig.savefig(HERE / "deep_sliding_fk.png", dpi=180)
    np.savez(
        HERE / "deep_sliding_fk.npz", slowness=SLOWNESS,
        bands=np.asarray(BANDS), time_window=np.asarray(TIME_WINDOW),
        aperture_m=APERTURE_M, step_m=STEP_M, turnaround_ch=TURNAROUND_CH,
        **products,
    )
    print("Saved deep_sliding_fk.png and deep_sliding_fk.npz")


if __name__ == "__main__":
    main()
