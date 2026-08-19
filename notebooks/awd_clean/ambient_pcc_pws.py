#!/usr/bin/env python3
"""Phase cross-correlation + phase-weighted stacking for the Figure 7c gather.

WHY THIS AND NOT MORE FILTERING.  Every filter tried on this dataset -- pre- and
post-correlation F-K, median common-mode removal, rank-k subspace projection
(global and windowed), and linear Radon tau-p -- fails the same way, because the
correlation gather is dominated by a BROAD zero-lag lobe whose width is
comparable to the whole lag window. A broad lobe stays partially aligned at every
slowness, so it cannot be separated in the slowness domain at all: muting
everything faster than 8000 / 6000 / 5000 / 4200 m/s leaves
corr(apparent velocity, energy) at +0.957 / +0.940 / +0.934 / +0.914. The problem
is upstream of any filter.

PCC and PWS are the literature's standard tools for pulling WEAK body-wave
arrivals out of ambient noise, and they act upstream -- they replace the
correlation and stacking operators rather than filtering their output.

  Schimmel (1999) BSSA 89, 1366   -- phase cross-correlation
  Schimmel & Paulssen (1997) GJI  -- phase-weighted stacking
  Schimmel & Gallart (2007) JGR   -- tf-PWS
  Ventosa, Schimmel & Stutzmann (2017) GJI 211, 30 -- ts-PWS

WHAT THEY CHANGE, and why it bears on our specific failure.

  PCC. Amplitude-unbiased: it correlates the UNIT-MODULUS analytic signals, so
  only instantaneous phase enters. Conventional cross-correlation is quadratic in
  amplitude, so a high-amplitude common component dominates the sum. If the
  broad zero-lag lobe is amplitude-driven, PCC suppresses it by construction.
  Implemented at nu = 2, for which PCC reduces exactly to the cross-correlation
  of unit-modulus analytic signals and is therefore FFT-computable; nu = 1 is not.
  Note RAM normalisation becomes redundant under PCC and is skipped -- amplitude
  is already discarded -- which also removes one processing choice the paper does
  not report.

  PWS. Weights the linear stack by inter-window phase coherence,
  |mean_j exp(i phi_j(tau))|^nu, so a lag at which windows agree in phase is
  enhanced and one where they do not is suppressed, independent of amplitude.
  Ventosa et al.'s practical guidance is followed: phase-stack a modest number of
  already-linearly-stacked groups rather than many raw windows, so the inputs to
  the phase stack have usable SNR. Here each hourly chunk is a group.

CONTROLS. Same receiver-order permutation null as every other branch in this
tree, so the p-value is directly comparable to config 0 (0.1470) and config 3
(0.9220). Same preconditions as ambient_radon_slant_stack.py: no recovery is
declared unless the pedestal is suppressed, the peak lies inside 2500-4000 m/s
and away from a scan edge, the causal side dominates, and p < 0.05.

Chunked like the rest of the tree: hourly products store the summed PCC and the
complex phase sum, both additive, so days aggregate exactly.

Output: ambient_transfer/pcc_pws/ chunks, then {npz,png,txt} beside this file.
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
from scipy.signal import butter, hilbert, sosfiltfilt

HERE = Path(__file__).resolve().parent
OUT = HERE / "ambient_transfer" / "pcc_pws"
STEM = HERE / "ambient_pcc_pws"
CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")

FMIN, FMAX = 5.0, 20.0
WIN_S, STEP_S = 30.0, 15.0
MAX_LAG = 0.35
NEIGHBOURS = 10
CH_LO = 23
OFFSETS_M = np.arange(50.0, 700.1, 50.0)
V_REF, V_LO, V_HI = 3200.0, 2500.0, 4000.0
VELOCITY_GRID = np.arange(1500.0, 6000.1, 25.0)
PWS_NU = 2.0
NULL_COUNT = 10000
SEED = 20260814


def corrected_path(p):
    return Path(str(p).replace("/data/SAFODAS1-harddrive-transfer",
                               "/data/SAFOD/SAFODAS1-harddrive-transfer"))


def day_rows(date):
    db = pd.read_csv(CSV, sep=r"\s+"); db = db[db.nSamples > 0].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    return db[db.t.dt.strftime("%Y-%m-%d") == date].sort_values("t").reset_index(drop=True)


def unit_analytic(x):
    """Unit-modulus analytic signal: discards amplitude, keeps instantaneous phase."""
    a = hilbert(x, axis=1)
    m = np.abs(a)
    return a / np.maximum(m, np.finfo(np.float64).tiny)


def run_chunk(args):
    rows = day_rows(args.date)
    sel = rows.iloc[args.start:args.start + args.nfiles + 1]     # +1 overlap file
    if not len(sel):
        raise SystemExit("no files")
    OUT.mkdir(parents=True, exist_ok=True)

    chunks = []; fs = dx = None; targets = keep = None; used = 0; expect = None
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
                targets = [CH_LO + int(round(o / dx)) for o in OFFSETS_M]
                targets = [t for t in targets if t + NEIGHBOURS < nch]
                keep = sorted({c for t in targets
                               for c in range(t - NEIGHBOURS, t + NEIGHBOURS + 1)} | {CH_LO})
        t0 = pd.Timestamp(row.startTime)
        if expect is not None and abs((t0 - expect).total_seconds()) > 1.0:
            print("gap before %s, stopping" % f.name, flush=True); break
        expect = t0 + pd.Timedelta(seconds=d.shape[1] / fs)
        chunks.append(d[keep]); used += 1
        if used == args.nfiles + 1:
            break
    if not chunks:
        raise SystemExit("no usable files")
    own = min(used, args.nfiles) * chunks[0].shape[1]

    x = np.concatenate(chunks, axis=1).astype(np.float64); del chunks
    x = np.diff(x, axis=1, prepend=x[:, :1])                    # phase -> strain rate
    x -= x.mean(axis=1, keepdims=True)
    sos = butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos")
    x = sosfiltfilt(sos, x, axis=1)
    # NO running-absolute-mean: PCC is amplitude-blind, so RAM is redundant here.
    z = unit_analytic(x); del x

    src_i = keep.index(CH_LO)
    nw = int(round(WIN_S * fs)); step = int(round(STEP_S * fs))
    ml = int(round(MAX_LAG * fs)); lags = np.arange(-ml, ml + 1) / fs
    nfft = 1 << int(np.ceil(np.log2(2 * nw - 1)))

    pcc_sum = np.zeros((len(keep), len(lags)))
    phase_sum = np.zeros((len(keep), len(lags)), dtype=np.complex128)
    nwin = 0
    for s in range(0, own, step):
        if s + nw > z.shape[1]:
            break
        seg = z[:, s:s + nw]
        A = np.fft.fft(seg[src_i], n=nfft)
        B = np.fft.fft(seg, n=nfft, axis=1)
        c = np.fft.ifft(np.conj(A)[None, :] * B, n=nfft, axis=1) / nw
        c = np.concatenate((c[:, -ml:], c[:, :ml + 1]), axis=1)
        real = c.real                       # PCC(nu=2)
        pcc_sum += real
        # instantaneous phase of this window's PCC, for the PWS weight
        ana = hilbert(real, axis=1)
        phase_sum += ana / np.maximum(np.abs(ana), np.finfo(np.float64).tiny)
        nwin += 1

    stem = "pcc_%s_start%04d_n%04d" % (args.date, args.start, min(used, args.nfiles))
    np.savez_compressed(OUT / (stem + ".npz"), lags=lags, pcc_sum=pcc_sum,
                        phase_sum=phase_sum, nwin=nwin, keep=np.array(keep),
                        targets=np.array(targets), dx=dx, fs=fs,
                        nfiles=min(used, args.nfiles), src_channel=CH_LO)
    print("wrote %s : %d files, %d windows" % (stem, min(used, args.nfiles), nwin), flush=True)


def neighbour_stack(corr, keep, targets):
    idx = {c: i for i, c in enumerate(keep)}
    out = np.zeros((len(targets), corr.shape[1])); cnt = []
    for k, R in enumerate(targets):
        mem = [idx[c] for c in range(R - NEIGHBOURS, R + NEIGHBOURS + 1) if c in idx]
        out[k] = corr[mem].sum(axis=0); cnt.append(len(mem))
    return out / np.array(cnt)[:, None]


def moveout_energy(gather, lags, offsets, v, sign=1.0):
    env = np.abs(hilbert(gather, axis=1))
    env /= np.maximum(np.median(env, axis=1, keepdims=True), 1e-30)
    half = max(1, int(0.012 * (len(lags) - 1) / (lags[-1] - lags[0])))
    vals = []
    for row, x in zip(env, offsets):
        k = int(np.argmin(np.abs(lags - sign * x / v)))
        lo, hi = max(0, k - half), min(len(lags), k + half + 1)
        vals.append(row[lo:hi].mean())
    return float(np.median(vals))


def aggregate(args):
    fs_ = sorted(glob.glob(str(OUT / ("pcc_%s_start*_n*.npz" % args.date))),
                 key=lambda f: int(re.search(r"start(\d+)", f).group(1)))
    if not fs_:
        raise SystemExit("no chunks for " + args.date)
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    acc = ph = None; nwin = 0; nfiles = 0
    for f in fs_:
        d = np.load(f, allow_pickle=True)
        if acc is None:
            acc = np.zeros_like(d["pcc_sum"]); ph = np.zeros_like(d["phase_sum"])
            lags = d["lags"]; keep = list(d["keep"]); targets = list(d["targets"])
            dx = float(d["dx"]); src = int(d["src_channel"])
        acc += d["pcc_sum"]; ph += d["phase_sum"]
        nwin += int(d["nwin"]); nfiles += int(d["nfiles"])

    linear = acc / nwin
    coherence = np.abs(ph / nwin)                 # inter-window phase coherence in [0,1]
    pws = linear * coherence ** PWS_NU

    offsets = (np.array(targets) - src) * dx
    g_lin = neighbour_stack(linear, keep, targets)
    g_pws = neighbour_stack(pws, keep, targets)

    say("Phase cross-correlation + phase-weighted stacking -- %s" % args.date)
    say("  %d chunks, %d files, %d windows; PCC nu=2 (unit-modulus analytic), PWS nu=%.0f"
        % (len(fs_), nfiles, nwin, PWS_NU))
    say("  no RAM normalisation (PCC is amplitude-blind); no F-K filter")
    say("  mean inter-window phase coherence: %.4f" % coherence.mean())
    say("")

    rng = np.random.default_rng(SEED)
    results = {}
    for lab, g in (("PCC linear stack", g_lin), ("PCC + PWS", g_pws)):
        sc = np.array([moveout_energy(g, lags, offsets, v) for v in VELOCITY_GRID])
        ac = np.array([moveout_energy(g, lags, offsets, v, sign=-1.0) for v in VELOCITY_GRID])
        k = int(np.argmax(sc))
        nulls = np.array([max(moveout_energy(g[rng.permutation(len(offsets))], lags, offsets, v)
                              for v in VELOCITY_GRID[::4]) for _ in range(NULL_COUNT // 20)])
        p = float((np.sum(nulls >= sc[k]) + 1) / (len(nulls) + 1))
        corr_v = float(np.corrcoef(VELOCITY_GRID, sc)[0, 1])
        c32 = float(np.interp(V_REF, VELOCITY_GRID, sc)); a32 = float(np.interp(V_REF, VELOCITY_GRID, ac))
        results[lab] = dict(sc=sc, ac=ac, peak=sc[k], vpeak=VELOCITY_GRID[k],
                            null95=float(np.percentile(nulls, 95)), p=p,
                            corr_v=corr_v, c32=c32, a32=a32, k=k, nulls=nulls)
        say("--- %s ---" % lab)
        say("  corr(velocity, score) %+.3f   [pedestal indicator; config 0 is +0.976, config 3 -0.381]"
            % corr_v)
        say("  peak %.4f at %.0f m/s | at 3200: causal %.4f acausal %.4f ratio %.2f"
            % (sc[k], VELOCITY_GRID[k], c32, a32, c32 / a32 if a32 else np.nan))
        say("  null95 %.4f | p = %.4f (n=%d)" % (np.percentile(nulls, 95), p, len(nulls)))
        ped = abs(corr_v) < 0.5
        fan = V_LO <= VELOCITY_GRID[k] <= V_HI
        edge = k not in (0, len(sc) - 1)
        dom = (c32 / a32) > 1.0 if a32 else False
        say("  preconditions: pedestal %s | in fan %s | not edge %s | causal %s | p<0.05 %s"
            % (ped, fan, edge, dom, p < 0.05))
        if all((ped, fan, edge, dom, p < 0.05)):
            say("  >>> RECOVERED <<<")
        else:
            say("  not recovered")
        say("")

    fig, ax = plt.subplots(1, 3, figsize=(17, 5.5), constrained_layout=True)
    for a_, g, t in ((ax[0], g_lin, "PCC linear"), (ax[1], g_pws, "PCC + PWS")):
        gn = g / np.abs(g).max(axis=1, keepdims=True)
        for o, row in zip(offsets, gn):
            a_.plot(lags, -o + row * 45, "k-", lw=0.8)
        a_.plot(offsets / V_REF, -offsets, "--", color="crimson", lw=1.5, label="3200 m/s")
        a_.set_xlim(-0.35, 0.35); a_.set_xlabel("lag (s)"); a_.set_title(t)
        a_.legend(fontsize=7)
    ax[0].set_ylabel("depth below wellhead (m)")
    for lab, st in (("PCC linear stack", "-"), ("PCC + PWS", "--")):
        r = results[lab]
        ax[2].plot(VELOCITY_GRID, r["sc"], st, label="%s (p=%.3f)" % (lab, r["p"]))
        ax[2].axhline(r["null95"], color="crimson", ls=":", lw=0.8)
    ax[2].axvline(V_REF, color="steelblue", ls=":"); ax[2].axvspan(V_LO, V_HI, color="orange", alpha=.15)
    ax[2].set_xlabel("trial velocity (m/s)"); ax[2].set_ylabel("moveout energy")
    ax[2].legend(fontsize=7); ax[2].grid(alpha=.3); ax[2].set_title("scan vs null")
    fig.savefig(str(STEM) + ".png", dpi=190)
    np.savez(str(STEM) + ".npz", lags=lags, offsets=offsets, g_lin=g_lin, g_pws=g_pws,
             velocity_grid=VELOCITY_GRID, coherence_mean=coherence.mean(),
             **{("%s_%s" % (k.replace(" ", "_").replace("+", "p"), f)): v[f]
                for k, v in results.items() for f in ("sc", "ac", "p", "corr_v")})
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("wrote %s.{npz,png,txt}" % STEM.name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="2024-12-20")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--nfiles", type=int, default=60)
    ap.add_argument("--aggregate", action="store_true")
    a = ap.parse_args()
    aggregate(a) if a.aggregate else run_chunk(a)


if __name__ == "__main__":
    main()
