"""Measure epoch repeatability of data-selected coherent AWD moveouts."""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
SCANS = HERE / "mode_scan.npz"
OUT_CSV = HERE / "mode_repeatability.csv"
OUT_FIG = HERE / "mode_repeatability.png"
PRE_S = 0.5
WINDOW_S = 0.040


def bp(data, fs, band):
    sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, data, axis=-1)


def score(section, coordinate, fs, velocity, intercept):
    ns = max(8, int(WINDOW_S * fs))
    indices = int(PRE_S * fs) + int(round(intercept * fs)) + np.rint(
        coordinate / velocity * fs
    ).astype(int)
    if indices.min() < 0 or indices.max() + ns > section.shape[1]:
        return np.nan, np.nan
    gather = np.stack([
        section[channel, index:index + ns]
        for channel, index in enumerate(indices)
    ])
    den = gather.shape[0] * np.sum(gather ** 2)
    semblance = np.sum(np.sum(gather, axis=0) ** 2) / den if den else np.nan
    amplitude = np.sqrt(np.mean(gather ** 2))
    return float(semblance), float(amplitude)


def main():
    data = np.load(STACKS)
    scans = np.load(SCANS)
    fs = float(data["fs"])
    counts = data["n_common"]
    bands = scans["bands"]
    vgrid, tgrid = scans["v_grid"], scans["t0_grid"]
    rows = []

    fig, axes = plt.subplots(2, len(bands), figsize=(16, 8), sharex=True)
    for fiber_row, (fiber, key, dx_key) in enumerate((
        ("nano", "nano_stacks", "dx_nano"),
        ("deep", "deep_stacks", "dx_deep"),
    )):
        stack = data[key]
        dx = float(data[dx_key])
        c0, c1 = int(80 / dx), min(int(440 / dx) + 1, stack.shape[1])
        coordinate = np.arange(c0, c1) * dx
        coordinate -= coordinate[0]

        for band_col, band in enumerate(bands):
            scan_key = f"{fiber}_{int(band[0])}_{int(band[1])}"
            global_scan = scans[scan_key]
            iv0, it0 = np.unravel_index(np.argmax(global_scan), global_scan.shape)
            v0, t0 = float(vgrid[iv0]), float(tgrid[it0])
            local_v = vgrid[np.abs(vgrid - v0) <= 300]
            local_t = tgrid[np.abs(tgrid - t0) <= 0.020]

            for epoch in range(stack.shape[0]):
                if counts[epoch] <= 0:
                    continue
                section = bp(stack[epoch, c0:c1], fs, band)
                best = (-np.inf, np.nan, np.nan, np.nan)
                for velocity in local_v:
                    for intercept in local_t:
                        sem, amp = score(section, coordinate, fs, velocity, intercept)
                        if np.isfinite(sem) and sem > best[0]:
                            best = (sem, float(velocity), float(intercept), amp)
                rows.append({
                    "fiber": fiber, "band_low_hz": float(band[0]),
                    "band_high_hz": float(band[1]), "epoch": epoch,
                    "n_common_drops": int(counts[epoch]),
                    "reference_speed_mps": v0,
                    "reference_intercept_s": t0,
                    "epoch_speed_mps": best[1],
                    "epoch_intercept_s": best[2],
                    "epoch_semblance": best[0], "epoch_rms": best[3],
                })

            subset = [r for r in rows if r["fiber"] == fiber and r["band_low_hz"] == band[0]]
            epochs = np.asarray([r["epoch"] for r in subset])
            speeds = np.asarray([r["epoch_speed_mps"] for r in subset])
            sems = np.asarray([r["epoch_semblance"] for r in subset])
            ax = axes[fiber_row, band_col]
            points = ax.scatter(epochs, speeds, c=sems, cmap="viridis", vmin=0, vmax=1)
            ax.axhline(v0, color="k", ls="--", lw=1)
            ax.set_title(
                f"{fiber}, {band[0]:.0f}-{band[1]:.0f} Hz\n"
                f"median={np.nanmedian(speeds):.0f}, MAD="
                f"{np.nanmedian(np.abs(speeds-np.nanmedian(speeds))):.0f} m/s"
            )
            ax.set_xlabel("AWD burst")
            if band_col == 0:
                ax.set_ylabel("best local apparent speed (m/s)")
            fig.colorbar(points, ax=ax, label="semblance")

    with OUT_CSV.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig.suptitle("Epoch stability of phase-neutral AWD moveout peaks")
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=180)
    print(f"Wrote {len(rows)} rows to {OUT_CSV}; saved {OUT_FIG}")


if __name__ == "__main__":
    main()
