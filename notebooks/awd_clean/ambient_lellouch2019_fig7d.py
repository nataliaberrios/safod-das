#!/usr/bin/env python3
"""Figure 7d: the constant-offset gather, computed from 2024-25 and compared
trace-by-trace against Lellouch's own released 7d correlograms.

WHY THIS IS THE ISOLATION TEST.  Everything tried so far compares our Figure 7c
gather against the PICTURE of his Figure 7c. That is not apples-to-apples. What
is on disk, numerically, is his Figure 7d: 14 correlograms at CONSTANT 50 m
source-receiver separation, labelled by depth, stacked from seven one-day
correlations (github.com/ariellellouch/SAFODDAS, `lellouch_traces/`).

Two facts make 7d the decisive comparison:

  1. Our picker behaves sensibly on HIS 7d traces -- the paper's 3-sample
     parabolic maximum gives 2416 m/s at 50 m rising to 4357 m/s at 700 m, with
     corr(depth, velocity) = +0.948. So the downstream is known-good on
     known-good input.

     CAREFUL, corrected 2026-08-19: r = 0.948 is corr(depth, velocity) for OUR
     picks on HIS traces -- a measure of how monotonic our own picks are. It is
     NOT agreement with his published Figure 9 curve, which has never been
     compared here value-for-value. Earlier wording in this docstring and in two
     status files said this "reproduc[es] his Figure 9"; that is withdrawn as an
     overstatement. The supportable claim is only that on known-good input the
     picker returns a monotonic velocity-depth trend of the expected form, so the
     picking stage is not what fails on 2024-25.
  2. Figure 7d, not 7c, is what the paper's velocity model comes from. Section 4.1:
     "we deduce we are observing a P wave and take the correlation functions shown
     in Figure 7d to estimate its velocity every 50 m."

So computing 7d from 2024-25 and overlaying it on his 7d isolates the difference
to a single, identically-defined observable measured in two epochs.

It is also intrinsically easier than 7c: the arrival sits at t = 50/v, i.e. 12-21
ms, which needs spatial coherence over 50 m rather than across the full 700 m
aperture, and it is a peak-lag measurement rather than a moveout that must be
resolved across offsets.

GEOMETRY.  Unlike 7c the source is not fixed. For each depth z the source is
channel z and the receiver centre is z + 50 m, sliding together down the array.
The paper's neighbour sum still applies, over receivers R-10..R+10 sharing that
source, and for 7d the paper applies travel-time shifts at the average velocity
read off 7c before summing -- "This shift is applied only after an average
velocity has been obtained using simple stacking." Both the unshifted and the
3200 m/s-shifted sums are produced here so the shift is a reported choice, not a
hidden one.

Preprocessing is identical to the validated 7c operator: strain-rate derivative,
optional per-sample median common-mode removal, running-absolute-mean, 30 s
windows at 15 s overlap, and the 5-20 Hz band applied to the stacked correlations.

Output: ambient_transfer/fig7d/ chunks, then ambient_lellouch2019_fig7d.{npz,png,txt}
"""
from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, sosfiltfilt

HERE = Path(__file__).resolve().parent
OUT = HERE / "ambient_transfer" / "fig7d"
STEM = HERE / "ambient_lellouch2019_fig7d"
CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")
LEL = Path("/home/groups/ettore88/nberrios/safod_das_git/lellouch_traces")

FMIN, FMAX = 5.0, 20.0
WIN_S, STEP_S = 30.0, 15.0
MAX_LAG = 0.35
NEIGHBOURS = 10
CH_LO = 23                      # wellhead, G0 registration
OFFSET_M = 50.0                 # the constant source-receiver separation
DEPTHS_M = np.arange(50.0, 700.1, 50.0)   # his released trace labels
V_SHIFT = 3200.0
RAM_S = 0.1


def corrected_path(p):
    return Path(str(p).replace("/data/SAFODAS1-harddrive-transfer",
                               "/data/SAFOD/SAFODAS1-harddrive-transfer"))


def day_rows(date):
    db = pd.read_csv(CSV, sep=r"\s+"); db = db[db.nSamples > 0].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    return db[db.t.dt.strftime("%Y-%m-%d") == date].sort_values("t").reset_index(drop=True)


def run_chunk(a):
    rows = day_rows(a.date)
    sel = rows.iloc[a.start:a.start + a.nfiles + 1]
    if not len(sel):
        raise SystemExit("no files")
    OUT.mkdir(parents=True, exist_ok=True)

    blocks = []; fs = dx = None; used = 0; expect = None; keep = None
    src_ch = rec_ch = None
    for row in sel.itertuples(index=False):
        f = corrected_path(row.file)
        if not f.is_file():
            continue
        with h5py.File(f, "r") as h:
            g = h["Acquisition/Raw[0]"]
            d = g["RawData"][:].astype(np.float32).T
            if fs is None:
                fs = float(g.attrs.get("OutputDataRate", 500.0))
                dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 1.0))
                nch = d.shape[0]
                off = int(round(OFFSET_M / dx))
                src_ch = [CH_LO + int(round(z / dx)) - off for z in DEPTHS_M]
                rec_ch = [s + off for s in src_ch]
                good = [i for i, (s, r) in enumerate(zip(src_ch, rec_ch))
                        if s >= 0 and r + NEIGHBOURS < nch]
                src_ch = [src_ch[i] for i in good]; rec_ch = [rec_ch[i] for i in good]
                keep = sorted({c for s, r in zip(src_ch, rec_ch)
                               for c in list(range(r - NEIGHBOURS, r + NEIGHBOURS + 1)) + [s]})
        t0 = pd.Timestamp(row.startTime)
        if expect is not None and abs((t0 - expect).total_seconds()) > 1.0:
            print("gap before %s, stopping" % f.name, flush=True); break
        expect = t0 + pd.Timedelta(seconds=d.shape[1] / fs)
        blocks.append(d[keep]); used += 1
        if used == a.nfiles + 1:
            break
    if not blocks:
        raise SystemExit("no usable files")
    own = min(used, a.nfiles) * blocks[0].shape[1]

    x = np.concatenate(blocks, axis=1).astype(np.float32); del blocks
    rate = np.empty_like(x); rate[:, 0] = 0.0
    np.subtract(x[:, 1:], x[:, :-1], out=rate[:, 1:]); rate *= np.float32(fs); del x
    if a.common_mode:
        rate -= np.median(rate, axis=0, keepdims=True).astype(np.float32)
    w = uniform_filter1d(np.abs(rate), size=max(3, int(RAM_S * fs) | 1), axis=1, mode="nearest")
    floor = max(np.finfo(np.float32).tiny, np.finfo(np.float32).eps * float(np.nanmedian(w)))
    rate /= np.maximum(w, np.float32(floor)); del w

    idx = {c: i for i, c in enumerate(keep)}
    nw = int(round(WIN_S * fs)); step = int(round(STEP_S * fs))
    nfft = 1 << int(np.ceil(np.log2(2 * nw - 1)))
    nfreq = nfft // 2 + 1
    ndepth = len(src_ch)
    plain = np.zeros((ndepth, nfreq), dtype=np.complex128)   # centre receiver only
    shifted = np.zeros((ndepth, nfreq), dtype=np.complex128) # R+-10, 3200 m/s shifts
    unshift = np.zeros((ndepth, nfreq), dtype=np.complex128) # R+-10, no shifts
    freqs = np.fft.rfftfreq(nfft, 1.0 / fs)
    nwin = 0
    for s0 in range(0, own, step):
        if s0 + nw > rate.shape[1]:
            break
        seg = rate[:, s0:s0 + nw]
        S = np.fft.rfft(seg, n=nfft, axis=1)
        for j, (sc, rc) in enumerate(zip(src_ch, rec_ch)):
            a_ = np.conj(S[idx[sc]])
            plain[j] += a_ * S[idx[rc]]
            for c in range(rc - NEIGHBOURS, rc + NEIGHBOURS + 1):
                cross = a_ * S[idx[c]]
                unshift[j] += cross
                # shift this neighbour back to the centre offset at V_SHIFT
                dt = ((c - sc) * dx - OFFSET_M) / V_SHIFT
                shifted[j] += cross * np.exp(-2j * np.pi * freqs * dt)
        nwin += 1

    stem = "fig7d_%s_start%04d_n%04d%s" % (a.date, a.start, min(used, a.nfiles),
                                           "_cm" if a.common_mode else "")
    np.savez_compressed(OUT / (stem + ".npz"), plain=plain, unshift=unshift,
                        shifted=shifted, nwin=nwin, n_fft=nfft, fs=fs, dx=dx,
                        depths_m=np.asarray(DEPTHS_M[:ndepth]),
                        src_ch=np.asarray(src_ch), rec_ch=np.asarray(rec_ch),
                        nfiles=min(used, a.nfiles), common_mode=a.common_mode)
    print("wrote %s : %d files, %d windows, %d depths"
          % (stem, min(used, a.nfiles), nwin, ndepth), flush=True)


def to_time(spec, nwin, nfft, fs, band=True):
    c = np.fft.irfft(spec / max(nwin, 1), n=nfft, axis=1)
    c = np.fft.fftshift(c, axes=1)
    if band:
        c = sosfiltfilt(butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos"),
                        c, axis=1)
    ml = int(round(MAX_LAG * fs)); mid = nfft // 2
    return np.arange(-ml, ml + 1) / fs, c[:, mid - ml:mid + ml + 1]


def parabolic(y, k, dt):
    if k <= 0 or k >= len(y) - 1:
        return 0.0
    y0, y1, y2 = y[k - 1], y[k], y[k + 1]
    den = y0 - 2 * y1 + y2
    return 0.0 if den == 0 else 0.5 * (y0 - y2) / den * dt


def pick(gather, lags, lo=0.005, hi=0.045):
    """The paper's picker: 3 adjacent samples with the largest value + quadratic."""
    dt = lags[1] - lags[0]
    w = np.where((lags >= lo) & (lags <= hi))[0]
    ts, vs = [], []
    for row in gather:
        k = w[int(np.argmax(row[w]))]
        t = lags[k] + parabolic(row, k, dt)
        ts.append(t); vs.append(OFFSET_M / t if t > 0 else np.nan)
    return np.array(ts), np.array(vs)


def aggregate(a):
    tag = "_cm" if a.common_mode else ""
    fs_ = sorted(glob.glob(str(OUT / ("fig7d_%s_start*_n*%s.npz" % (a.date, tag)))),
                 key=lambda f: int(re.search(r"start(\d+)", f).group(1)))
    fs_ = [f for f in fs_ if ("_cm" in Path(f).stem) == bool(a.common_mode)]
    if not fs_:
        raise SystemExit("no chunks")
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    acc = {}; nwin = 0; nfiles = 0
    for f in fs_:
        d = np.load(f, allow_pickle=True)
        for k in ("plain", "unshift", "shifted"):
            acc[k] = acc.get(k, 0) + d[k]
        nwin += int(d["nwin"]); nfiles += int(d["nfiles"])
        nfft = int(d["n_fft"]); fsr = float(d["fs"]); depths = d["depths_m"]

    say("Figure 7d constant-offset gather -- %s" % a.date)
    say("  %d chunks, %d files, %d windows; source-receiver separation %.0f m"
        % (len(fs_), nfiles, nwin, OFFSET_M))
    say("  common-mode removal: %s" % bool(a.common_mode))
    say("")

    ours = {}
    for k in ("plain", "unshift", "shifted"):
        lags, g = to_time(acc[k], nwin, nfft, fsr)
        ours[k] = g
    say("--- our 2024-25 picks (paper's picker, %.0f-%.0f ms window) ---" % (5, 45))
    say("  %-8s %10s %10s %10s" % ("depth", "centre", "R+-10", "R+-10 shift"))
    picks = {}
    for k in ("plain", "unshift", "shifted"):
        picks[k] = pick(ours[k], lags)
    for i, z in enumerate(depths):
        say("  %6.0f m %10.0f %10.0f %10.0f"
            % (z, picks["plain"][1][i], picks["unshift"][1][i], picks["shifted"][1][i]))

    # Lellouch's released 7d
    lel_z = lel_v = None
    try:
        from obspy import read
        st = read(str(LEL / "CCstackSAFOD.mseed"))
        by = {int(t.stats.station): t.data.astype(float) for t in st}
        lel_z = np.array(sorted(by), float)
        G = np.array([by[int(z)] for z in lel_z])
        n = G.shape[1]; ll = (np.arange(n) - (n - 1) // 2) / 500.0
        _, lel_v = pick(G, ll)
        say("")
        say("--- Lellouch 2017 released 7d, identical picker ---")
        say("  " + "  ".join("%.0f:%0.f" % (z, v) for z, v in zip(lel_z, lel_v)))
        r = np.corrcoef(lel_z, lel_v)[0, 1]
        # Label deliberately does NOT say "his Figure 9 model": r is the
        # monotonicity of OUR picks on HIS traces, not agreement with his
        # published curve, which is not digitised anywhere in this tree.
        say("  fit r = %.3f  (corr(depth, velocity) of OUR picks on HIS "
            "released traces)" % r)
    except Exception as e:
        say("  (released comparison unavailable: %r)" % (e,))

    say("")
    say("--- comparison ---")
    for k, lab in (("shifted", "ours, R+-10 + 3200 m/s shift"),
                   ("unshift", "ours, R+-10 unshifted"),
                   ("plain", "ours, centre receiver only")):
        v = picks[k][1]
        good = np.isfinite(v)
        if good.sum() >= 4:
            rr = np.corrcoef(depths[good], v[good])[0, 1]
            say("  %-30s v(50m)=%5.0f v(700m)=%5.0f  r vs depth = %+.3f"
                % (lab, v[0], v[-1], rr))
    if lel_v is not None:
        say("  %-30s v(50m)=%5.0f v(700m)=%5.0f  r vs depth = %+.3f"
            % ("Lellouch 2017", lel_v[0], lel_v[-1], np.corrcoef(lel_z, lel_v)[0, 1]))
    say("")
    say("  A recovered 7d needs a MONOTONIC velocity-depth trend like his")
    say("  (2416 -> 4357 m/s, corr(depth,v) = +0.948 for OUR picks on HIS traces;")
    say("  that is our picks' monotonicity, NOT agreement with his Figure 9 curve).")
    say("  Scatter or a flat/negative r means the peak lag is not a travel time.")

    fig, ax = plt.subplots(1, 3, figsize=(17, 6), constrained_layout=True)
    for i, (k, lab) in enumerate((("shifted", "ours R+-10 shifted"),
                                  ("unshift", "ours R+-10 unshifted"))):
        g = ours[k]; gn = g / np.maximum(np.abs(g).max(axis=1, keepdims=True), 1e-30)
        for z, row in zip(depths, gn):
            ax[i].plot(lags, -z + row * 22, "k-", lw=0.8)
        ax[i].axvline(OFFSET_M / V_SHIFT, color="crimson", ls="--", lw=1.2,
                      label="50 m / 3200 m/s")
        ax[i].set_xlim(-0.06, 0.09); ax[i].set_xlabel("lag (s)")
        ax[i].set_title(lab); ax[i].legend(fontsize=7)
    ax[0].set_ylabel("depth (m)")
    ax[2].plot(picks["shifted"][1], depths, "o-", color="steelblue", label="ours (shifted)")
    ax[2].plot(picks["unshift"][1], depths, "s--", color="teal", lw=1, label="ours (unshifted)")
    if lel_v is not None:
        ax[2].plot(lel_v, lel_z, "o-", color="crimson", lw=2, label="Lellouch 2017")
    ax[2].invert_yaxis(); ax[2].set_xlim(0, 8000)
    ax[2].set_xlabel("velocity from 50 m / t_peak (m/s)"); ax[2].set_ylabel("depth (m)")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=.3); ax[2].set_title("Figure 9 comparison")
    fig.savefig(str(STEM) + ".png", dpi=190)
    np.savez(str(STEM) + ".npz", lags=lags, depths=depths,
             **{("ours_" + k): v for k, v in ours.items()},
             **{("pick_" + k): picks[k][1] for k in picks},
             lellouch_z=lel_z if lel_z is not None else np.array([]),
             lellouch_v=lel_v if lel_v is not None else np.array([]))
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="2024-12-20")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--nfiles", type=int, default=60)
    ap.add_argument("--common-mode", action="store_true")
    ap.add_argument("--aggregate", action="store_true")
    a = ap.parse_args()
    aggregate(a) if a.aggregate else run_chunk(a)


if __name__ == "__main__":
    main()
