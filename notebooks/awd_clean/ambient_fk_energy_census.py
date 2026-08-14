#!/usr/bin/env python3
"""
Is there body-wave energy in the 2024-2025 ambient wavefield at all?

This is the go/no-go that `sanity/README.md` describes as the decision point:
"if there is no body-wave-velocity ridge in the wavefield, no CC trick will
recover a body-wave Green's function". It is a census of the RAW wavefield --
no correlation, no F-K filtering, no velocity selection. It only asks how the
5-20 Hz energy is distributed over apparent velocity and propagation direction.

Method: strain-rate proxy, 5-20 Hz band, 2-D FFT over (channel, time) in 30 s
windows, average power spectrum, then bin |F(f,k)|^2 by apparent velocity
v = f/k and by the sign of f*k (propagation direction along the fiber).

Reference points. Lellouch et al. (2019) report a ~3,200 m/s downgoing P packet
on this fiber. Tube/guided and surface modes at SAFOD sit near 500-1,500 m/s.
If essentially all 5-20 Hz energy is below ~1,500 m/s, the 3,200 m/s target is
absent from the input and the reproduction is data-limited, not method-limited.

Outputs: ambient_fk_energy_census.{npz,png,txt}
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from scipy.signal import butter, sosfiltfilt, detrend
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")
STEM = HERE / "ambient_fk_energy_census"    # per-date suffix appended in main()
FMIN, FMAX = 5.0, 20.0
WIN_S = 30.0
VEL_EDGES = np.array([0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000,
                      5000, 6500, 9000, 1e9])


def corrected_path(p):
    return Path(str(p).replace("/data/SAFODAS1-harddrive-transfer",
                               "/data/SAFOD/SAFODAS1-harddrive-transfer"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="2024-12-20")
    ap.add_argument("--nfiles", type=int, default=12)
    ap.add_argument("--stride", type=int, default=120, help="spread files over the day")
    a = ap.parse_args()
    global STEM
    STEM = HERE / ("ambient_fk_energy_census_%s" % a.date)

    db = pd.read_csv(CSV, sep=r"\s+"); db = db[db.nSamples > 0].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    day = db[db.t.dt.strftime("%Y-%m-%d") == a.date].sort_values("t").reset_index(drop=True)
    picks = day.iloc[::a.stride].head(a.nfiles)

    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    say("Raw 5-20 Hz F-K energy census -- %s" % a.date)
    say("no correlation, no F-K filter, no velocity selection")

    Pacc = None
    nwin = 0
    for row in picks.itertuples(index=False):
        f = corrected_path(row.file)
        if not f.exists():
            continue
        with h5py.File(f, "r") as h:
            g = h["Acquisition/Raw[0]"]
            d = g["RawData"][:].astype(np.float32).T
            fs = float(g.attrs.get("OutputDataRate", 500.0))
            dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 1.0))
        x = np.diff(d, axis=1, prepend=d[:, :1])          # phase -> strain rate
        x = detrend(x, axis=1, type="linear")
        sos = butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos")
        x = sosfiltfilt(sos, x, axis=1)
        nw = int(round(WIN_S * fs))
        for s in range(0, x.shape[1] - nw + 1, nw):
            seg = x[:, s:s + nw]
            seg = seg / (np.std(seg) + 1e-20)
            F = np.fft.fftshift(np.fft.fft2(seg))
            P = np.abs(F) ** 2
            Pacc = P if Pacc is None else Pacc + P
            nwin += 1
        say("  %s : %d windows so far" % (f.name, nwin))

    if Pacc is None:
        raise SystemExit("no usable files")
    P = Pacc / nwin
    nch, nt = P.shape
    k = np.fft.fftshift(np.fft.fftfreq(nch, dx))          # cycles / m
    fr = np.fft.fftshift(np.fft.fftfreq(nt, 1.0 / fs))    # Hz
    KK, FF = np.meshgrid(k, fr, indexing="ij")

    band = (np.abs(FF) >= FMIN) & (np.abs(FF) <= FMAX) & (np.abs(KK) > 1e-9)
    vapp = np.where(band, np.abs(FF) / np.abs(KK), np.nan)
    w = np.where(band, P, 0.0)
    tot = w.sum()

    say("")
    say("--- 5-20 Hz energy by apparent velocity ---")
    say("   velocity band        share of 5-20 Hz energy")
    rows = []
    for lo, hi in zip(VEL_EDGES[:-1], VEL_EDGES[1:]):
        m = band & (vapp >= lo) & (vapp < hi)
        sh = w[m].sum() / tot * 100
        rows.append((lo, hi, sh))
        lbl = "%5.0f - %5.0f m/s" % (lo, hi) if hi < 1e8 else "%5.0f +      m/s" % lo
        say("   %s : %6.2f %%   %s" % (lbl, sh, "#" * int(round(sh / 2))))

    below = sum(s for lo, hi, s in rows if hi <= 1500)
    target = sum(s for lo, hi, s in rows if lo >= 2500 and hi <= 4000)
    say("")
    say("  energy below 1500 m/s          : %.1f %%" % below)
    say("  energy in the 2500-4000 m/s fan : %.1f %%" % target)

    # directionality in the target fan
    fan = band & (vapp >= 2500) & (vapp <= 4000)
    down = fan & (FF * KK < 0)
    up = fan & (FF * KK > 0)
    say("")
    say("--- directionality inside the 2500-4000 m/s fan ---")
    say("  F*K<0 (increasing coordinate / downgoing) : %.1f %%" %
        (w[down].sum() / max(w[fan].sum(), 1e-30) * 100))
    say("  F*K>0 (decreasing coordinate / upgoing)   : %.1f %%" %
        (w[up].sum() / max(w[fan].sum(), 1e-30) * 100))

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    pos = fr >= 0
    im = ax[0].pcolormesh(k, fr[pos], 10 * np.log10(P[:, pos].T + 1e-30),
                          shading="auto", cmap="magma")
    for v, c in ((500, "cyan"), (1500, "lime"), (3200, "white")):
        kk = np.linspace(1e-4, k.max(), 50)
        ax[0].plot(kk, kk * v, "--", color=c, lw=1.2, label="%d m/s" % v)
        ax[0].plot(-kk, kk * v, "--", color=c, lw=1.2)
    ax[0].set_ylim(0, 40); ax[0].set_xlabel("wavenumber (cycles/m)")
    ax[0].set_ylabel("frequency (Hz)"); ax[0].legend(fontsize=7)
    ax[0].set_title("raw F-K power, %s (dB)" % a.date)
    fig.colorbar(im, ax=ax[0])

    lo = [r[0] for r in rows]; sh = [r[2] for r in rows]
    ax[1].barh(range(len(rows)), sh, color="steelblue")
    ax[1].set_yticks(range(len(rows)))
    ax[1].set_yticklabels(["%.0f-%.0f" % (r[0], r[1]) if r[1] < 1e8 else "%.0f+" % r[0]
                           for r in rows], fontsize=8)
    ax[1].axhspan(4.5, 7.5, color="orange", alpha=.25, label="2500-4000 m/s target")
    ax[1].set_xlabel("% of 5-20 Hz energy"); ax[1].legend(fontsize=8)
    ax[1].set_title("where the 5-20 Hz energy lives")
    fig.savefig(str(STEM) + ".png", dpi=200)

    np.savez(str(STEM) + ".npz", P=P, k=k, f=fr, rows=np.array(rows), nwin=nwin)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
