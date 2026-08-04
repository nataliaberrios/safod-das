"""Signed f-p dispersion of 2005 PGSI axial check-shot records."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import read
from scipy.signal import butter, detrend, sosfiltfilt, windows


HERE = Path(__file__).resolve().parent
PGSI = HERE / "pgsi_reference" / "Check shots"
COORDS = PGSI / "PGSIarray_rec_coords_pos1.txt"
SLOWNESS = np.linspace(-2.0e-3, 2.0e-3, 801)
FREQ_RANGE = (10.0, 120.0)
TIME_RANGE = (-0.10, 0.80)


def beam(section, coordinate, fs):
    section = np.nan_to_num(section.astype(float))
    sos = butter(4, [8.0, 150.0], btype="bandpass", fs=fs, output="sos")
    section = sosfiltfilt(sos, section, axis=-1)
    section = detrend(section, axis=-1)
    section *= windows.tukey(section.shape[0], 0.15)[:, None]
    section *= windows.tukey(section.shape[1], 0.15)[None, :]
    nfft = int(2 ** np.ceil(np.log2(max(4096, section.shape[1]))))
    spectrum = np.fft.rfft(section, n=nfft, axis=1)
    frequency = np.fft.rfftfreq(nfft, 1 / fs)
    keep = (frequency >= FREQ_RANGE[0]) & (frequency <= FREQ_RANGE[1])
    frequency, spectrum = frequency[keep], spectrum[:, keep]
    spectrum /= np.abs(spectrum) + 1e-30
    coordinate = coordinate - coordinate.mean()
    power = np.empty((len(frequency), len(SLOWNESS)), np.float32)
    for i, freq in enumerate(frequency):
        steering = np.exp(2j * np.pi * freq * np.outer(SLOWNESS, coordinate))
        power[i] = (np.abs(steering @ spectrum[:, i]) / len(coordinate)) ** 2
    return frequency, power


def main():
    geometry = np.genfromtxt(COORDS, names=True)
    depth = geometry["REC_DEP"]  # shallow to deep, metres below reference surface
    files = sorted(PGSI.glob("*.dat"))
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), sharex=True, sharey=True)
    products, lines = {}, [
        "2005 PGSI axial check-shot signed frequency-slowness summary",
        "Each record uses its own SEG2 DELAY header; axial channels are reversed into shallow-to-deep order.",
        "",
    ]
    for axis, path in zip(axes.flat, files):
        stream = read(str(path))
        fs = float(stream[0].stats.sampling_rate)
        delay = float(stream[0].stats.seg2.get("DELAY", 0.0))
        # README: deepest level first, axial then two transverse components.
        # One-based axial rows are 1,4,...,238; reverse to shallow-to-deep.
        axial = np.stack([stream[i].data for i in range(0, 240, 3)])[::-1]
        time = np.arange(axial.shape[1]) / fs + delay
        keep = (time >= TIME_RANGE[0]) & (time <= TIME_RANGE[1])
        frequency, power = beam(axial[:, keep], depth, fs)
        ridge_index = np.argmax(power, axis=1)
        ridge_p = SLOWNESS[ridge_index]
        ridge_power = power[np.arange(len(frequency)), ridge_index]

        key = path.stem
        products[f"{key}_power"] = power
        products[f"{key}_ridge_slowness"] = ridge_p
        products[f"{key}_ridge_power"] = ridge_power
        products[f"{key}_delay"] = np.asarray(delay)
        image = axis.pcolormesh(
            SLOWNESS * 1e3, frequency, power, shading="auto", cmap="magma",
            vmin=0, vmax=max(0.25, np.percentile(power, 99.7)),
        )
        axis.plot(ridge_p * 1e3, frequency, "c.", ms=1.5)
        axis.axvline(0, color="w", lw=0.6)
        for speed in (1000, 1500, 3000, 5000):
            for sign in (-1, 1):
                axis.axvline(sign * 1e3 / speed, color="w", ls=":", lw=0.4, alpha=0.35)
        axis.set_title(f"Shot {key}; SEG2 DELAY={delay:+.3f} s")
        axis.set_xlabel("signed slowness (ms/m)")
        axis.set_ylabel("frequency (Hz)")
        fig.colorbar(image, ax=axis, label="normalized phase-beam power")

        band = (frequency >= 20) & (frequency <= 60)
        cutoff = np.percentile(ridge_power[band], 60)
        strong = band & (ridge_power >= cutoff)
        median_p = np.median(ridge_p[strong])
        velocity = 1 / median_p if abs(median_p) >= 1 / 6000 else np.nan
        lines.append(
            f"shot {key}: delay={delay:+.3f} s, median p={median_p*1e3:+.5f} ms/m, "
            f"signed v={velocity:+.1f} m/s, median ridge power={np.median(ridge_power[band]):.4f}"
        )

    fig.suptitle(
        "2005 PGSI check shots: axial-component signed dispersion\n"
        "80 levels, position 1 receiver-depth geometry"
    )
    fig.tight_layout()
    fig.savefig(HERE / "pgsi_checkshot_fk.png", dpi=180)
    np.savez(
        HERE / "pgsi_checkshot_fk.npz", frequency=frequency,
        slowness=SLOWNESS, receiver_depth=depth, **products,
    )
    report = "\n".join(lines) + "\n"
    (HERE / "pgsi_checkshot_fk.txt").write_text(report)
    print(report)
    print("Saved pgsi_checkshot_fk.png/.npz/.txt")


if __name__ == "__main__":
    main()
