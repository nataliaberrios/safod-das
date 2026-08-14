#!/usr/bin/env python3
"""
Reproduce Lellouch et al. (2019) Figure 7c from the 2024-2025 SAFOD archive.

TARGET. Figure 7c is the top-source correlation gather: virtual source fixed at
the top channel, virtual receivers every 50 m down the array, showing "a clear
wave packet with an apparent velocity of about 3,200 m/s". Benchmark trend:

    t(offset) = offset / 3200  ->  50 m: 0.0156 s ,  700 m: 0.2188 s

and the causal side must dominate the acausal side ("the strong signals of these
downgoing waves compared to the upgoing ones suggest that the dominant ambient
field sources originate at the surface").

WHY THIS IS NOT THE EXISTING PIPELINE. `ambient_transfer_test.stack_day` computes
a top-source gather already, but omits three steps that the paper (section 4.1)
or the sibling pipelines in this repository apply, all of them PRE-correlation and
therefore not repairable in the cached `ambient_transfer/transfer_*.npz` products:

  1. NEIGHBOUR STACK.  C_S,R = sum_{Z=R-10}^{R+10} C_S,Z. The paper states this
     "induces a smoothing of the velocity model but is required to extract a
     clear signal". The cached products correlate single pairs -> ~21x less SNR.
     For Figure 7c this is a SIMPLE sum with NO travel-time shifts: the paper
     applies 3,200 m/s shifts only after an average velocity has been obtained
     from simple stacking, and those shifts belong to Figure 7d. Keeping 7c
     shift-free is what makes it an unconditioned observable.

  2. SPECTRAL NORMALISATION.  Equation (6): Im[G] ~ <eps_zs . eps_zR> / |S(w)|^2.
     The cached run had its whitening branch set to None.

  3. COMMON-MODE REMOVAL.  Not in the paper, but the legacy CC pipeline in this
     repository subtracts the per-sample median across channels, and the v19
     audit records that the earlier notebook did median-over-channel removal.
     A zero-delay component common to all channels lands exactly on top of the
     short-offset arrivals. Diagnostic on the cached products: correlations peak
     at t = 0 exactly, C(0) = 0.344 against C(0.015) = 0.196, at every depth.

Each of the three is a switch so their contributions can be separated rather
than assumed. Default is all three ON.

NO F-K FILTER IS APPLIED ANYWHERE IN THIS SCRIPT. The F-K fan was never part of
the paper's ambient workflow, and it failed this project's own pre-filter
channel-scramble gate (see Ambient_FK_QC_workflow.ipynb).

CONTROLS. (a) The receiver-order permutation null: shuffle which correlation
trace is assigned to which offset and rescore, which destroys ordered moveout
while preserving every trace. (b) Causal vs acausal asymmetry, which the paper
predicts on physical grounds. (c) Chunked resumable output so independent
sub-stacks can be compared for convergence.

USAGE
    python ambient_lellouch_fig7c.py --date 2024-12-20 --nfiles 10 --start 0
    python ambient_lellouch_fig7c.py --aggregate --date 2024-12-20
"""
from __future__ import annotations
import argparse, glob, json, os, re
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from scipy.signal import butter, sosfiltfilt, detrend
from scipy.ndimage import uniform_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")
OUT = HERE / "ambient_transfer" / "lellouch_fig7c"

FMIN, FMAX = 5.0, 20.0
RAM_SECONDS = 5.0
WIN_S, OVERLAP_S = 30.0, 15.0      # paper: 30 s windows, 50% overlap
MAX_LAG = 0.35                     # covers 700 m / 2000 m/s = 0.35 s
NEIGHBOURS = 10                    # paper: R-10 .. R+10
OFFSETS_M = np.arange(50.0, 700.1, 50.0)
V_REF = 3200.0                     # paper's average P velocity, for the benchmark only


def corrected_path(p):
    return Path(str(p).replace("/data/SAFODAS1-harddrive-transfer",
                               "/data/SAFOD/SAFODAS1-harddrive-transfer"))


def load_segment(path):
    with h5py.File(path, "r") as h:
        g = h["Acquisition/Raw[0]"]; ds = g["RawData"]
        d = ds[:].astype(np.float32)
        a = h["Acquisition"].attrs
        fs = float(g.attrs.get("OutputDataRate", ds.attrs.get("OutputDataRate", 500.0)))
        dx = float(a.get("SpatialSamplingInterval", 1.0))
    return d.T, fs, dx


def spectral_whiten(x, fs, fmin=FMIN, fmax=FMAX, smooth_hz=1.0):
    """Equation (6): divide by the ambient power spectrum. Implemented as
    per-trace spectral whitening inside the analysis band, with a smoothed
    amplitude spectrum so the operation is a normalisation and not a
    phase-only transform."""
    n = x.shape[1]
    X = np.fft.rfft(x, axis=1)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    amp = np.abs(X)
    w = max(1, int(round(smooth_hz / (f[1] - f[0]))))
    amp_s = uniform_filter1d(amp, size=w, axis=1, mode="nearest")
    band = (f >= fmin) & (f <= fmax)
    floor = np.percentile(amp_s[:, band], 10, axis=1, keepdims=True) * 0.1 + 1e-20
    Xw = np.zeros_like(X)
    Xw[:, band] = X[:, band] / np.maximum(amp_s[:, band], floor)
    return np.fft.irfft(Xw, n=n, axis=1)


def preprocess(x, fs, common_mode=True, whiten=True):
    """Strain -> strain-rate proxy -> optional common mode -> band -> RAM -> optional whitening."""
    x = np.diff(x, axis=1, prepend=x[:, :1])          # strain rate proxy
    if common_mode:
        x = x - np.median(x, axis=0, keepdims=True)    # per-sample median across channels
    x = detrend(x, axis=1, type="linear")
    sos = butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos")
    x = sosfiltfilt(sos, x, axis=1)
    nwin = max(3, int(RAM_SECONDS * fs))
    m = uniform_filter1d(np.abs(x), size=nwin, axis=1, mode="nearest")
    floor = np.percentile(m, 5, axis=1, keepdims=True) * 0.1 + 1e-12
    x = x / np.maximum(m, floor)
    if whiten:
        x = spectral_whiten(x, fs)
    return x


def windowed_correlations(x, fs, src, recv_channels, max_lag=MAX_LAG):
    """Correlate src against every channel in recv_channels, in 30 s windows with
    50% overlap, summing over windows. Returns (lags, [nrecv, nlag]) unnormalised
    sum and the accumulated normalisation, so chunks can be combined correctly."""
    nw = int(round(WIN_S * fs)); step = int(round((WIN_S - OVERLAP_S) * fs))
    ml = int(round(max_lag * fs)); lags = np.arange(-ml, ml + 1) / fs
    nfft = 1 << int(np.ceil(np.log2(2 * nw - 1)))
    acc = np.zeros((len(recv_channels), len(lags)))
    nwin_used = 0
    starts = range(0, x.shape[1] - nw + 1, step)
    for s in starts:
        seg = x[:, s:s + nw]
        a = np.fft.rfft(seg[src], n=nfft)
        B = np.fft.rfft(seg[recv_channels], n=nfft, axis=1)
        C = np.fft.irfft(np.conj(a)[None, :] * B, n=nfft, axis=1)
        # non-negative lags at the head, negative lags wrap to the tail
        C = np.concatenate((C[:, -ml:], C[:, :ml + 1]), axis=1)
        den = np.sqrt(np.sum(seg[src] ** 2) * np.sum(seg[recv_channels] ** 2, axis=1))
        good = den > 0
        C[good] /= den[good][:, None]
        C[~good] = 0.0
        acc += C
        nwin_used += 1
    return lags, acc, nwin_used


def neighbour_stack(corr, recv_channels, targets, n=NEIGHBOURS):
    """C_S,R = sum_{Z=R-n}^{R+n} C_S,Z  (simple sum, no travel-time shifts)."""
    idx = {c: i for i, c in enumerate(recv_channels)}
    out = np.zeros((len(targets), corr.shape[1]))
    counts = []
    for k, R in enumerate(targets):
        members = [idx[c] for c in range(R - n, R + n + 1) if c in idx]
        out[k] = corr[members].sum(axis=0)
        counts.append(len(members))
    return out, np.array(counts)


def moveout_score(gather, lags, offsets, v):
    """Median correlation sampled along a constant apparent velocity."""
    vals = []
    for row, d in zip(gather, offsets):
        t = d / v
        if abs(t) > lags.max():
            vals.append(np.nan); continue
        vals.append(row[int(np.argmin(np.abs(lags - t)))])
    return float(np.nanmedian(vals))


def run_chunk(args):
    db = pd.read_csv(CSV, sep=r"\s+")
    db = db[db.nSamples > 0].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    day = db[db.t.dt.strftime("%Y-%m-%d") == args.date].sort_values("t").reset_index(drop=True)
    rows = day.iloc[args.start:args.start + args.nfiles]
    if not len(rows):
        raise SystemExit("no files for %s at start %d" % (args.date, args.start))

    OUT.mkdir(parents=True, exist_ok=True)
    acc = None; lags = None; nwin_tot = 0; used = []
    recv_channels = targets = None; dx = fs = None

    for i, row in enumerate(rows.itertuples(index=False), 1):
        f = corrected_path(row.file)
        if not f.exists():
            continue
        try:
            x, fs, dx = load_segment(f)
            if targets is None:
                targets = [int(round(o / dx)) for o in OFFSETS_M]
                nch = x.shape[0]
                targets = [t for t in targets if t + NEIGHBOURS < nch]
                keep = sorted({c for t in targets
                               for c in range(t - NEIGHBOURS, t + NEIGHBOURS + 1)
                               if 0 <= c < nch})
                recv_channels = keep
            xp = preprocess(x, fs, common_mode=not args.no_common_mode,
                            whiten=not args.no_whiten)
            sub = np.vstack([xp[0:1], xp[recv_channels]])
            lg, a, nw = windowed_correlations(sub, fs, 0,
                                              np.arange(1, len(recv_channels) + 1))
            if acc is None:
                acc = np.zeros_like(a); lags = lg
            acc += a; nwin_tot += nw; used.append(str(f))
            if i % 5 == 0:
                print("  %d/%d files, %d windows" % (i, len(rows), nwin_tot), flush=True)
        except Exception as e:
            print("  skip", f.name, repr(e), flush=True)

    if not used:
        raise SystemExit("no usable files")

    stem = "fig7c_%s_start%d_n%d%s%s" % (
        args.date, args.start, len(used),
        "" if not args.no_common_mode else "_nocm",
        "" if not args.no_whiten else "_nowh")
    np.savez(OUT / (stem + ".npz"), lags=lags, corr_sum=acc, nwin=nwin_tot,
             recv_channels=np.array(recv_channels), targets=np.array(targets),
             dx=dx, fs=fs, used=np.array(used),
             common_mode=not args.no_common_mode, whiten=not args.no_whiten)
    print("wrote", OUT / (stem + ".npz"))
    print("  %d files, %d windows, gather will be %d offsets"
          % (len(used), nwin_tot, len(targets)))


def aggregate(args):
    pat = str(OUT / ("fig7c_%s_start*_n*.npz" % args.date))
    files = sorted(glob.glob(pat))
    files = [f for f in files
             if ("_nocm" in f) == args.no_common_mode and ("_nowh" in f) == args.no_whiten]
    if not files:
        raise SystemExit("no chunks matching " + pat)
    acc = None; nwin = 0; nfiles = 0
    for f in files:
        d = np.load(f, allow_pickle=True)
        a = d["corr_sum"].astype(float)
        if acc is None:
            acc = np.zeros_like(a); lags = d["lags"].astype(float)
            recv = d["recv_channels"]; targets = d["targets"]; dx = float(d["dx"])
        if a.shape != acc.shape:
            print("  shape mismatch, skipping", os.path.basename(f)); continue
        acc += a; nwin += int(d["nwin"]); nfiles += len(d["used"])

    corr = acc / max(nwin, 1)
    gather, counts = neighbour_stack(corr, list(recv), list(targets))
    gather /= counts[:, None]                       # mean, so panels are comparable
    offsets = np.array(targets) * dx

    ml = np.abs(lags) <= MAX_LAG
    lags = lags[ml]; gather = gather[:, ml]

    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    say("Lellouch Figure 7c reproduction -- %s" % args.date)
    say("  %d chunks, %d files, %d windows; common_mode=%s whiten=%s neighbours=+-%d"
        % (len(files), nfiles, nwin, not args.no_common_mode, not args.no_whiten, NEIGHBOURS))
    say("  NO F-K filter applied.")
    say("")

    # velocity scan
    vv = np.arange(1500, 6001, 25.0)
    sc = np.array([moveout_score(gather, lags, offsets, v) for v in vv])
    kb = int(np.nanargmax(np.abs(sc)))
    say("--- moveout scan ---")
    say("  peak |score| %.4f at %.0f m/s   (paper: ~3200 m/s)" % (abs(sc[kb]), vv[kb]))
    for v in (2500, 3000, 3200, 3500, 4000):
        say("   score at %4d m/s = %+.4f" % (v, moveout_score(gather, lags, offsets, v)))

    # receiver-order permutation null
    rng = np.random.default_rng(20260814)
    null = np.array([max(abs(moveout_score(gather[rng.permutation(len(offsets))],
                                           lags, offsets, v)) for v in vv[::4])
                     for _ in range(2000)])
    obs = abs(sc[kb])
    p = (np.sum(null >= obs) + 1) / (len(null) + 1)
    say("")
    say("--- receiver-order permutation null (2000) ---")
    say("  null median %.4f | 95th pct %.4f | observed %.4f | p = %.4f"
        % (np.median(null), np.percentile(null, 95), obs, p))

    # causal vs acausal (the paper's physical prediction)
    k0 = int(np.argmin(np.abs(lags)))
    cz = np.abs(gather[:, k0 + 1:]).mean(); az = np.abs(gather[:, :k0]).mean()
    say("")
    say("--- causal / acausal asymmetry ---")
    say("  mean|C| causal %.5f | acausal %.5f | ratio %.2f  (paper predicts >1)"
        % (cz, az, cz / az if az > 0 else np.inf))

    # per-offset picks against the benchmark
    say("")
    say("--- per-offset causal pick vs the 3200 m/s benchmark ---")
    picks = []
    for o, row in zip(offsets, gather):
        tb = o / V_REF
        w = np.where((lags > 0.5 * tb) & (lags < 2.0 * tb + 0.02))[0]
        if not len(w):
            picks.append(np.nan); continue
        k = w[int(np.argmax(row[w]))]
        picks.append(lags[k])
        say("   %5.0f m : benchmark %.4f s , picked %.4f s -> v = %6.0f m/s"
            % (o, tb, lags[k], o / lags[k] if lags[k] > 0 else np.nan))
    picks = np.array(picks)
    good = np.isfinite(picks) & (picks > 0)
    if good.sum() >= 4:
        A = np.vstack([offsets[good], np.ones(good.sum())]).T
        sl, ic = np.linalg.lstsq(A, picks[good], rcond=None)[0]
        say("  fitted apparent velocity = %.0f m/s (intercept %.4f s)" % (1 / sl, ic))

    fig, ax = plt.subplots(1, 2, figsize=(13, 7), constrained_layout=True)
    g = gather / np.abs(gather).max(axis=1, keepdims=True)
    for k, (o, row) in enumerate(zip(offsets, g)):
        ax[0].plot(lags, -o + row * 45, "k-", lw=0.8)
    ax[0].plot(offsets / V_REF, -offsets, "--", color="crimson", lw=1.5,
               label="3200 m/s (paper)")
    ax[0].set_xlim(-0.4, 0.7); ax[0].set_xlabel("Time [s]")
    ax[0].set_ylabel("Depth [m]"); ax[0].legend(fontsize=8)
    ax[0].set_title("Figure 7c reproduction, %s\n%d files, no F-K filter" % (args.date, nfiles))

    ax[1].plot(vv, sc, "k-")
    ax[1].axhline(np.percentile(null, 95), color="crimson", ls="--", lw=1,
                  label="perm. null 95th pct")
    ax[1].axhline(-np.percentile(null, 95), color="crimson", ls="--", lw=1)
    ax[1].axvline(V_REF, color="steelblue", ls=":", label="3200 m/s")
    ax[1].set_xlabel("trial apparent velocity (m/s)")
    ax[1].set_ylabel("median correlation along trajectory")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    ax[1].set_title("moveout scan vs receiver-order null")

    tag = "%s%s%s" % (args.date, "_nocm" if args.no_common_mode else "",
                      "_nowh" if args.no_whiten else "")
    fig.savefig(OUT / ("fig7c_%s_summary.png" % tag), dpi=200)
    np.savez(OUT / ("fig7c_%s_aggregate.npz" % tag), lags=lags, gather=gather,
             offsets=offsets, vv=vv, scores=sc, null=null, p=p, nfiles=nfiles, nwin=nwin)
    (OUT / ("fig7c_%s_aggregate.txt" % tag)).write_text("\n".join(log) + "\n")
    say("")
    say("wrote fig7c_%s_{summary.png,aggregate.npz,aggregate.txt}" % tag)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="2024-12-20")
    ap.add_argument("--nfiles", type=int, default=10)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--no-common-mode", action="store_true")
    ap.add_argument("--no-whiten", action="store_true")
    a = ap.parse_args()
    aggregate(a) if a.aggregate else run_chunk(a)


if __name__ == "__main__":
    main()
