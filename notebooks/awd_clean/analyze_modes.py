"""Phase-neutral coherent-moveout scan of canonical AWD epoch stacks."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
OUT = HERE / "mode_scan.npz"
FIG = HERE / "mode_scan.png"

PRE_S = 0.5
BANDS = ((5.0, 15.0), (15.0, 30.0), (30.0, 60.0))
V_GRID = np.arange(500.0, 6000.1, 25.0)
T0_GRID = np.arange(-0.05, 0.301, 0.002)
WINDOW_S = 0.040


def bandpass(data, fs, band):
    sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, data, axis=-1)


def weighted_section(stacks, counts):
    keep = counts > 0
    weights = counts[keep].astype(float)
    return np.tensordot(weights, stacks[keep], axes=(0, 0)) / weights.sum()


def semblance(section, coordinate, fs):
    """Scan unsigned apparent speed and intercept without assigning a phase."""
    nt = section.shape[1]
    ns = max(8, int(WINDOW_S * fs))
    origin = int(PRE_S * fs)
    result = np.zeros((len(V_GRID), len(T0_GRID)), np.float32)
    for iv, velocity in enumerate(V_GRID):
        moveout = np.rint(coordinate / velocity * fs).astype(int)
        for it, intercept in enumerate(T0_GRID):
            indices = origin + int(round(intercept * fs)) + moveout
            if indices.min() < 0 or indices.max() + ns > nt:
                continue
            gather = np.stack([
                section[channel, index:index + ns]
                for channel, index in enumerate(indices)
            ])
            numerator = np.sum(np.sum(gather, axis=0) ** 2)
            denominator = gather.shape[0] * np.sum(gather ** 2)
            result[iv, it] = numerator / denominator if denominator else 0.0
    return result


def distinct_peaks(scan, count=5, separation=250.0):
    profile = scan.max(axis=1)
    selected = []
    for index in np.argsort(profile)[::-1]:
        if all(abs(V_GRID[index] - V_GRID[other]) >= separation for other in selected):
            selected.append(index)
        if len(selected) == count:
            break
    return selected


def main():
    product = np.load(STACKS)
    fs = float(product["fs"])
    counts = product["n_common"]
    nano = weighted_section(product["nano_stacks"], counts)
    deep = weighted_section(product["deep_stacks"], counts)
    dx_nano = float(product["dx_nano"])
    dx_deep = float(product["dx_deep"])

    # Compare matched 80-440 m apertures. Channel coordinate remains distance
    # along fiber; no unverified depth offset enters this scan.
    apertures = {}
    for name, section, dx in (("nano", nano, dx_nano), ("deep", deep, dx_deep)):
        c0, c1 = int(80 / dx), min(int(440 / dx) + 1, section.shape[0])
        coordinate = np.arange(c0, c1) * dx
        apertures[name] = (section[c0:c1], coordinate - coordinate[0], coordinate)

    scans = {}
    fig, axes = plt.subplots(2, len(BANDS), figsize=(16, 9), sharex=True, sharey=True)
    for row, (name, (section, relative, absolute)) in enumerate(apertures.items()):
        for col, band in enumerate(BANDS):
            filtered = bandpass(section, fs, band)
            scan = semblance(filtered, relative, fs)
            scans[f"{name}_{int(band[0])}_{int(band[1])}"] = scan
            peaks = distinct_peaks(scan)
            print(f"{name} {band[0]:.0f}-{band[1]:.0f} Hz")
            for rank, iv in enumerate(peaks, 1):
                it = int(np.argmax(scan[iv]))
                print(
                    f"  {rank}: speed={V_GRID[iv]:.0f} m/s, "
                    f"intercept={T0_GRID[it]*1e3:.1f} ms, "
                    f"semblance={scan[iv, it]:.4f}"
                )
            ax = axes[row, col]
            image = ax.imshow(
                scan.T, origin="lower", aspect="auto", cmap="magma",
                extent=[V_GRID[0], V_GRID[-1], T0_GRID[0] * 1e3, T0_GRID[-1] * 1e3],
            )
            for iv in peaks[:3]:
                it = int(np.argmax(scan[iv]))
                ax.plot(V_GRID[iv], T0_GRID[it] * 1e3, "co", mfc="none")
            ax.set_title(f"{name.capitalize()}, {band[0]:.0f}-{band[1]:.0f} Hz")
            ax.set_xlabel("unsigned apparent speed (m/s)")
            if col == 0:
                ax.set_ylabel("intercept at aperture top (ms)")
            fig.colorbar(image, ax=ax, shrink=0.82, label="semblance")

    fig.suptitle("SAFOD AWD phase-neutral moveout scan, matched 80-440 m fiber coordinates")
    fig.tight_layout()
    fig.savefig(FIG, dpi=180)
    np.savez(
        OUT, v_grid=V_GRID, t0_grid=T0_GRID, bands=np.asarray(BANDS),
        n_common=counts, **scans,
    )
    print(f"Saved {FIG} and {OUT}")


if __name__ == "__main__":
    main()
