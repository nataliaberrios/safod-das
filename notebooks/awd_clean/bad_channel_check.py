"""Are the high-RMS Nano channels an optical artifact or real ground motion?

A channel whose RMS sits far above its neighbours has a short list of possible
causes, and they split cleanly into two groups:

  optical / instrumental          mechanical / real
  ----------------------          -----------------
  fading (low backscatter)        poor or absent coupling, free fiber
  phase-unwrap / 2*pi jumps       casing joint, centraliser, void
  splice, connector, tight bend   local strain concentration, fracture
  damaged fiber                   fluid or tube-wave resonance

Four measurements separate them, and the first is close to decisive.

1. WIDTH AGAINST GAUGE LENGTH. The interrogator integrates strain over a
   16.46 m gauge (13.0 channels at 1.266 m spacing). Any real strain is
   therefore smeared over at least that many channels. An anomaly one or two
   channels wide is narrower than the instrument's own resolution and cannot be
   ground motion.

2. SPECTRAL SHAPE. Fading raises the phase-noise floor at all frequencies, so
   the excess is white to Nyquist. Mechanical energy is band-limited and usually
   peaked.

3. STATIONARITY. Fading follows the Rayleigh speckle pattern, which drifts with
   temperature and laser wavelength over minutes to hours. A splice, joint or
   coupling defect sits at the same channel for the whole survey.

4. NEIGHBOUR COHERENCE. Real strain is spatially coherent over the gauge length,
   so a bad-but-real channel still correlates with its neighbours. Fading noise
   is independent channel to channel.

Reads a few files spread across the survey with all reader processing disabled.

Output
------
figures/awd_2026/plain_look/fig11_bad_channels.png, bad_channels.csv
"""
from __future__ import annotations

import csv
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import welch

SEARCH_DIRS = [
    "/home/groups/edunham/nberrios/safod_das/DAS-utilities/python",
    "/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python",
]
for _p in SEARCH_DIRS:
    if Path(_p).exists() and _p not in sys.path:
        sys.path.insert(0, _p)
from DASutils import readFile_protobuf  # noqa: E402

NANO_DIR = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/Nano")
OUT_DIR = Path(
    "/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026/plain_look"
) / "diagnostic"

DX = 1.26606202          # channel spacing, from the .pb acquisition_stats
GL_M = 16.4588           # gauge length, same source
GL_CH = GL_M / DX        # 13.0 channels
RAW_KW = dict(filter=False, median=False, detrend=False, tapering=False)
N_FILES = 6              # spread across the survey
DUR_S = 60.0             # seconds of each file to use
RMS_FACTOR = 3.0         # flag channels this far above the local median


def local_median(x, half=40):
    """Running median over channels, so the threshold tracks the depth trend."""
    out = np.empty_like(x)
    for i in range(x.size):
        out[i] = np.median(x[max(0, i - half):i + half + 1])
    return out


def anomaly_width(rms, base, idx):
    """Width in channels over which this channel stays above halfway to the peak."""
    peak, floor = rms[idx], base[idx]
    if peak <= floor:
        return 0.0
    half = floor + 0.5 * (peak - floor)
    lo = idx
    while lo > 0 and rms[lo - 1] > half:
        lo -= 1
    hi = idx
    while hi < rms.size - 1 and rms[hi + 1] > half:
        hi += 1
    return float(hi - lo + 1)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    files = sorted(NANO_DIR.glob("*.pb"))
    picks = [files[i] for i in np.linspace(0, len(files) - 1, N_FILES).astype(int)]

    rms_all, spec_all, nb_all, freq = [], [], [], None
    for f in picks:
        print("reading", f.name, flush=True)
        d, info = readFile_protobuf([str(f)], fmin=1.0, fmax=250.0,
                                    desampling=False, **RAW_KW)
        fs = float(info["fs"])
        seg = d[:, :int(DUR_S * fs)].astype(np.float64)
        seg -= seg.mean(axis=1, keepdims=True)
        rms_all.append(np.sqrt(np.mean(seg ** 2, axis=1)))

        freq, p = welch(seg, fs=fs, nperseg=1024, axis=-1)
        spec_all.append(p)

        # correlation with the channel one gauge length away
        k = int(round(GL_CH))
        a, b = seg[:-k], seg[k:]
        num = np.sum(a * b, axis=1)
        den = np.sqrt(np.sum(a * a, axis=1) * np.sum(b * b, axis=1))
        nb = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
        nb_all.append(np.concatenate([nb, np.full(k, np.nan)]))
        del d, seg

    rms_all = np.array(rms_all)
    nb_all = np.array(nb_all)
    z = np.arange(rms_all.shape[1]) * DX
    rms_med = np.median(rms_all, axis=0)
    base = local_median(rms_med)
    flagged = np.where(rms_med > RMS_FACTOR * base)[0]
    print(f"\n{flagged.size} channels above {RMS_FACTOR}x the local median")

    # per-file presence: is the same channel bad every time?
    persist = np.array([
        np.mean([rms_all[i, c] > RMS_FACTOR * local_median(rms_all[i])[c]
                 for i in range(rms_all.shape[0])]) for c in flagged
    ]) if flagged.size else np.array([])

    rows = []
    for j, c in enumerate(flagged):
        w = anomaly_width(rms_med, base, c)
        p = np.median(np.array(spec_all)[:, c, :], axis=0)
        band = freq > 1.0
        lo = (freq > 1) & (freq < 50)
        hi = freq > 200
        flatness = float(np.median(p[hi]) / np.median(p[lo])) if lo.any() else np.nan
        rows.append(dict(
            channel=int(c), z_m=float(z[c]),
            rms_ratio=float(rms_med[c] / base[c]),
            width_ch=w, width_vs_gauge=w / GL_CH,
            persistence=float(persist[j]),
            hi_lo_power=flatness,
            neighbour_corr=float(np.nanmedian(nb_all[:, c])),
            verdict=("instrumental" if w < 0.7 * GL_CH else "wider than gauge; check"),
        ))

    with open(OUT_DIR / "bad_channels.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["channel"])
        wtr.writeheader()
        wtr.writerows(rows)

    print(f"\n{'ch':>5} {'z (m)':>8} {'RMSx':>7} {'width':>7} {'/gauge':>7} "
          f"{'persist':>8} {'hi/lo':>8} {'nbr r':>7}  verdict")
    for r in rows:
        print(f"{r['channel']:5d} {r['z_m']:8.1f} {r['rms_ratio']:7.1f} "
              f"{r['width_ch']:7.1f} {r['width_vs_gauge']:7.2f} "
              f"{r['persistence']:8.2f} {r['hi_lo_power']:8.3f} "
              f"{r['neighbour_corr']:+7.3f}  {r['verdict']}")

    # ------------------------------------------------------------------ figure
    fig, ax = plt.subplots(2, 2, figsize=(15, 9))

    ax[0, 0].semilogy(z, rms_med, lw=0.8, label="median over files")
    ax[0, 0].semilogy(z, base, lw=1.2, color="C1", label="local median (running)")
    ax[0, 0].semilogy(z, RMS_FACTOR * base, lw=1.0, ls="--", color="C3",
                      label=f"{RMS_FACTOR:g}x threshold")
    ax[0, 0].plot(z[flagged], rms_med[flagged], "v", color="C3", ms=6)
    ax[0, 0].set_xlabel("distance along fiber (m)")
    ax[0, 0].set_ylabel("RMS (raw units)")
    ax[0, 0].set_title(f"A  {flagged.size} channels above threshold", fontsize=11)
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(alpha=0.3)

    # B: the decisive test -- width against gauge length
    if rows:
        w = np.array([r["width_ch"] for r in rows])
        rr = np.array([r["rms_ratio"] for r in rows])
        ax[0, 1].scatter(w, rr, s=40, color="C3", zorder=3)
        ax[0, 1].axvline(GL_CH, color="k", lw=1.4,
                         label=f"gauge length = {GL_CH:.1f} ch ({GL_M:g} m)")
        ax[0, 1].axvspan(0, GL_CH, color="C3", alpha=0.10)
        ax[0, 1].text(0.5, 0.03, "narrower than the gauge:\ncannot be ground motion",
                      transform=ax[0, 1].transAxes, ha="center", fontsize=9)
        ax[0, 1].set_yscale("log")
        ax[0, 1].set_xlabel("anomaly width (channels)")
        ax[0, 1].set_ylabel("RMS / local median")
        ax[0, 1].set_title("B  Width against gauge length", fontsize=11)
        ax[0, 1].legend(fontsize=8)
        ax[0, 1].grid(alpha=0.3)

    # C: spectra of flagged channels vs a clean reference
    spec_med = np.median(np.array(spec_all), axis=0)
    clean = np.setdiff1d(np.arange(z.size), flagged)
    ax[1, 0].loglog(freq[1:], np.sqrt(np.median(spec_med[clean], axis=0))[1:],
                    "k", lw=1.6, label="typical channel")
    for c in flagged[:8]:
        ax[1, 0].loglog(freq[1:], np.sqrt(spec_med[c])[1:], lw=1.0, alpha=0.8,
                        label=f"{z[c]:.0f} m")
    ax[1, 0].set_xlabel("frequency (Hz)")
    ax[1, 0].set_ylabel("amplitude spectral density")
    ax[1, 0].set_title("C  Flat to Nyquist = optical; peaked = mechanical", fontsize=11)
    ax[1, 0].legend(fontsize=7, ncol=2)
    ax[1, 0].grid(alpha=0.3, which="both")

    # D: coherence with the channel one gauge length away
    ax[1, 1].plot(z, np.nanmedian(nb_all, axis=0), lw=0.7, color="C0")
    if flagged.size:
        ax[1, 1].plot(z[flagged], np.nanmedian(nb_all, axis=0)[flagged], "v",
                      color="C3", ms=6, label="flagged channels")
        ax[1, 1].legend(fontsize=8)
    ax[1, 1].axhline(0, color="k", lw=0.6)
    ax[1, 1].set_xlabel("distance along fiber (m)")
    ax[1, 1].set_ylabel(f"correlation with channel +{int(round(GL_CH))}")
    ax[1, 1].set_title("D  Real strain stays coherent over a gauge length", fontsize=11)
    ax[1, 1].grid(alpha=0.3)

    fig.suptitle(
        f"11  High-RMS Nano channels: instrumental or real? "
        f"{N_FILES} files across the survey, {DUR_S:g} s each", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig11_bad_channels.png", dpi=140)
    print("\nwrote", OUT_DIR / "fig11_bad_channels.png")


if __name__ == "__main__":
    main()
