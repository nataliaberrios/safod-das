"""Publication-quality frequency-dependent observability for canonical Nano AWD stacks."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import detrend, windows


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
VELOCITY = 2975.0
INTERCEPT_S = -0.022
SIGNAL_RELATIVE = (-0.04, 0.20)
NOISE_WINDOW = (-0.38, -0.14)
DISTANCE_EDGES = np.arange(40.0, 501.0, 40.0)
FREQUENCY_LIMIT = 120.0


def multitaper_psd(traces, fs, nfft=512):
    traces = detrend(traces.astype(float), axis=-1, type="linear")
    tapers = windows.dpss(traces.shape[-1], 2.5, Kmax=4)
    power = 0.0
    for taper in tapers:
        spectrum = np.fft.rfft(traces * taper, n=nfft, axis=-1)
        power += np.abs(spectrum) ** 2
    return power / len(tapers)


def main():
    with np.load(STACKS) as data:
        stacks = data["nano_stacks"]
        counts = data["n_common"]
        fs = float(data["fs"])
        dx = float(data["dx_nano"])
    ns = int(round((SIGNAL_RELATIVE[1] - SIGNAL_RELATIVE[0]) * fs))
    n0 = int(round((0.5 + NOISE_WINDOW[0]) * fs))
    n1 = n0 + ns
    frequency = np.fft.rfftfreq(512, 1 / fs)
    keep_f = frequency <= FREQUENCY_LIMIT
    frequency = frequency[keep_f]
    centers = 0.5 * (DISTANCE_EDGES[:-1] + DISTANCE_EDGES[1:])
    snr = np.full((stacks.shape[0], centers.size, frequency.size), np.nan, np.float32)

    for epoch in np.flatnonzero(counts > 0):
        section = stacks[epoch]
        for iz, (lo, hi) in enumerate(zip(DISTANCE_EDGES[:-1], DISTANCE_EDGES[1:])):
            channels = np.arange(max(0, int(np.ceil(lo / dx))), min(section.shape[0], int(np.floor(hi / dx)) + 1))
            signal = []
            noise = []
            for channel in channels:
                arrival = INTERCEPT_S + channel * dx / VELOCITY
                s0 = int(round((0.5 + arrival + SIGNAL_RELATIVE[0]) * fs))
                s1 = s0 + ns
                if s0 >= 0 and s1 <= section.shape[1] and n0 >= 0 and n1 <= section.shape[1]:
                    signal.append(section[channel, s0:s1])
                    noise.append(section[channel, n0:n1])
            if len(signal) < 8:
                continue
            ps = np.mean(multitaper_psd(np.asarray(signal), fs), axis=0)[keep_f]
            pn = np.mean(multitaper_psd(np.asarray(noise), fs), axis=0)[keep_f]
            snr[epoch, iz] = 10 * np.log10((ps + np.finfo(float).tiny) / (pn + np.finfo(float).tiny))

    median = np.nanmedian(snr, axis=0)
    low, high = np.nanpercentile(snr, [16, 84], axis=0)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.12, 0.9, 4))
    chosen = [1, 4, 7, 10]
    for color, index in zip(colors, chosen):
        axes[0].plot(frequency, median[index], color=color, lw=1.8, label=f"{DISTANCE_EDGES[index]:.0f}–{DISTANCE_EDGES[index+1]:.0f} m")
        axes[0].fill_between(frequency, low[index], high[index], color=color, alpha=0.16, linewidth=0)
    axes[0].axhline(0, color="0.2", ls="--", lw=0.9)
    axes[0].axvspan(30, 60, color="0.5", alpha=0.10)
    axes[0].set(xlim=(2, FREQUENCY_LIMIT), ylim=(-20, 35), xlabel="frequency (Hz)", ylabel="signal-to-noise ratio (dB)", title="A  Burst-median spectral SNR")
    axes[0].legend(frameon=False, title="distance along Nano")
    image = axes[1].pcolormesh(frequency, centers, median, shading="auto", cmap="RdBu_r", vmin=-15, vmax=15)
    axes[1].contour(frequency, centers, median, levels=[0], colors="k", linewidths=1.0)
    axes[1].axvspan(30, 60, color="k", alpha=0.06)
    axes[1].set(xlim=(2, FREQUENCY_LIMIT), xlabel="frequency (Hz)", ylabel="distance along Nano fiber (m)", title="B  Frequency-dependent observability")
    fig.colorbar(image, ax=axes[1], label="median spectral SNR (dB)")
    fig.suptitle("SAFOD AWD Nano: canonical paired-burst frequency observability", fontsize=13)
    fig.savefig(HERE / "nano_frequency_observability.png", dpi=220)
    np.savez(HERE / "nano_frequency_observability.npz", frequency=frequency, distance_centers=centers, distance_edges=DISTANCE_EDGES, snr_db=snr, median_snr_db=median, p16_snr_db=low, p84_snr_db=high, velocity=VELOCITY, intercept_s=INTERCEPT_S, signal_relative=np.asarray(SIGNAL_RELATIVE), noise_window=np.asarray(NOISE_WINDOW))
    band = (frequency >= 30) & (frequency <= 60)
    band_snr = np.nanmedian(median[:, band], axis=1)
    detectable = centers[band_snr > 0]
    limit = np.nanmax(detectable) if detectable.size else np.nan
    report = f"30-60 Hz farthest 40 m bin center with median SNR > 0 dB: {limit:.1f} m\n"
    (HERE / "nano_frequency_observability.txt").write_text(report)
    print(report, end="")
    print("Saved nano_frequency_observability.png/.npz/.txt")


if __name__ == "__main__":
    main()
