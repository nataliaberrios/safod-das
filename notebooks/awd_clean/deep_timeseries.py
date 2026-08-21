#!/usr/bin/env python3
"""Deep-fibre ambient CC across the whole archive: grand stack + time variation.

SCOPE. There are ~947 h of Deep record available:
    ARM A   24.6 h   2026-06-15 -> 06-16, pre-survey (crew on site 06-15)
    ARM B    5.1 h   drop-free gaps inside the AWD survey
    ARM C  917.0 h   2026-05-01 -> 06-08, no personnel on site
Everything so far used at most 5 h of it. This processes blocks spanning the
whole archive, so the arrival can be measured as a FUNCTION OF TIME and the
blocks can also be summed into a grand stack.

DESIGN. Cross-spectra are accumulated per block and written to disk; the block is
released before the next is read, so peak memory is set by the block rather than
the total. Every block is independently scored, and blocks are summed afterwards,
so one pass gives both the time series and the grand stack. Blocks are the natural
unit for a Slurm array: `--task/--ntasks` splits the archive.

FREQUENCY BAND. The default is 5-20 Hz, following Lellouch et al. (2019) and every
other measurement in this project, so results stay comparable. But the AWD
manuscript reports, from the active source on THIS fibre, that the Deep mode sits
"near 1547 m/s, strongest at 15-30 Hz" -- so 5-20 Hz misses the upper half of
where the mode is strongest and carries 5-15 Hz where it is not. `--band` allows
15-30 Hz to be run alongside, which is a declared change of selection region and
is reported as such rather than quietly adopted.

At 1525 m/s (the corrected arrival velocity; 1675 is retracted, see
arrival_velocities.py) 5-20 Hz is 76-305 m wavelength and 15-30 Hz is 51-102 m,
against 2.0419 m channel spacing: both well sampled, worst case 25 samples per
wavelength.

GEOMETRY. Fixed-top virtual source at the wellhead (Deep channel 211, TVD 2 m)
with receivers below it inside the near-vertical section, which is the Lellouch
Figure 7c geometry. Channels 0-210 are surface lead-in and are refused.

Output per task: deep_timeseries_<band>_t<NN>.npz  (block spectra + per-block scores)
Combine with --combine to write deep_timeseries_<band>.{npz,png,txt}
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, hilbert, sosfiltfilt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ambient_lellouch2019_exact_stack as ex
import deep_cc_steps as steps
import arrival_velocities as av
import safod_geometry as geo

STEM = HERE / "deep_timeseries"
SOURCE_CH = 211                 # the wellhead
DEPTH_SPAN = 700.0              # Lellouch's 50-700 m offsets
WINDOW_S, STEP_S = 30.0, 15.0
# Derived, not hand-picked: it must hold DEPTH_SPAN at the slowest velocity
# scored. 0.35 s could not, and truncated everything beyond 534 m at 1525 m/s.
MAX_LAG_S = round(av.required_lag_s(DEPTH_SPAN, min(av.V_DEEP_ARRIVAL, 1000.0)), 2)
V_ARR = av.V_DEEP_ARRIVAL       # retracted from 1675
V_LELL = 3200.0
BANDS = {"5-20": (5.0, 20.0), "15-30": (15.0, 30.0)}
INK, MUTED = "#444444", "#6b6b6b"
C1, C2, C3 = "#0072B2", "#D55E00", "#009E73"


def receivers():
    g = geo.load()
    j = int(np.searchsorted(g["channel"], SOURCE_CH))
    if not g["in_hole"][j]:
        raise SystemExit("source channel %d is surface lead-in" % SOURCE_CH)
    z0 = g["tvd_m"][j]
    m = g["near_vertical"] & (g["tvd_m"] > z0) & (g["tvd_m"] <= z0 + DEPTH_SPAN)
    ch = g["channel"][m].astype(int)
    z = g["tvd_m"][m]
    o = np.argsort(z)
    return ch[o], z[o], z0


def block_spectra(paths, rx_ch, band, target_bytes=4.0e9):
    """Accumulate cross-spectra for one block, streaming over CHANNEL GROUPS.

    Reading a 4 h block whole needed ~40 GB (342 channels x 4 h, doubled by the
    concatenate inside read_records) and OOM-killed 12 of 16 array tasks at
    64 GB. So the block has to be read in pieces -- but WHICH axis is split is a
    correctness question, not a taste one.

    Splitting in TIME was tried first and rejected. It works, but it is only
    approximately equal to reading the block whole: a window straddling a join
    is lost unless a tail is carried across, and the 0.1 s running-abs-mean sees
    a chunk edge rather than real data for 0.05 s either side of every join.
    Measured on synthetic data that came to a 2.0e-3 relative difference, which
    is physically negligible and *therefore the problem*: the only test you can
    write against it is one loose enough to also hide a real bug.

    Splitting in CHANNEL has no such joins. Each receiver's full time series is
    read and processed exactly as before; only the set of receivers held at once
    changes, and the cross-spectra for disjoint receivers are independent. The
    result is bit-identical to reading the block whole, so the accompanying test
    (`test_block_spectra_chunking.py`) can demand equality at 1e-12 and will
    catch anything that is actually wrong. The cost is re-reading the source
    channel and re-walking the files once per group, which is a few passes.

    One deliberate change: the epsilon floor under the running-abs-mean divisor
    is now per CHANNEL rather than per 64-channel block. A floor computed over a
    block would depend on how channels were grouped, which would reintroduce the
    grouping dependence this function exists to avoid. It is a floor of
    eps * median, so it binds only where the divisor is essentially zero.
    """
    probe, fs, _ = steps.read_records(paths[:1], [SOURCE_CH])
    rec_n = int(probe.shape[1]); del probe
    n_win, n_step = int(WINDOW_S * fs), int(STEP_S * fs)
    n_fft = 1 << int(np.ceil(np.log2(2 * n_win - 1)))
    n_ram = ex.odd_ram_samples(0.1, fs)
    total_n = rec_n * len(paths)
    per_group = max(1, int(target_bytes / (total_n * 4 * 3)) - 1)

    acc = np.zeros((len(rx_ch), n_fft // 2 + 1), dtype=np.complex128)
    nw = None
    for g0 in range(0, len(rx_ch), per_group):
        grp = [int(c) for c in rx_ch[g0:g0 + per_group]]
        raw, _, _ = steps.read_records(paths, [SOURCE_CH] + grp)
        x = np.asarray(raw, dtype=np.float32); del raw
        x[:, 1:] = np.diff(x, axis=1) * np.float32(fs); x[:, 0] = 0.0
        w = uniform_filter1d(np.abs(x), size=n_ram, axis=1, mode="nearest")
        floor = (np.float32(np.finfo(np.float32).eps)
                 * np.nanmedian(w, axis=1, keepdims=True).astype(np.float32))
        np.divide(x, np.maximum(w, floor), out=x)
        del w, floor
        nw_g = 0
        for s in range(0, x.shape[1] - n_win + 1, n_step):
            W = np.fft.rfft(x[:, s:s + n_win], n=n_fft, axis=1)
            acc[g0:g0 + len(grp)] += np.conj(W[0]) * W[1:]
            del W
            nw_g += 1
        del x
        if nw is not None and nw_g != nw:
            raise SystemExit("channel groups saw different window counts "
                             "(%d vs %d) -- the records are not aligned"
                             % (nw_g, nw))
        nw = nw_g
    return acc, int(nw or 0), fs, n_fft


def to_gather(acc, nw, fs, n_fft, band):
    full = np.fft.fftshift(np.fft.irfft(acc / max(1, nw), n=n_fft, axis=-1), axes=-1)
    lags = (np.arange(n_fft) - n_fft // 2) / fs
    keep = np.abs(lags) <= MAX_LAG_S
    sos = butter(4, list(band), btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, full, axis=-1)[:, keep], lags[keep]


def score(gather, lags, sep, v):
    """Envelope amplitude along t = sep/v, and the causal/acausal ratio.

    REFUSES rather than clamps. `argmin(|lags - centre|)` returns an EDGE index
    when the target lag lies outside the window, so a separation the window
    cannot reach was silently scored on the edge of the correlation instead of
    failing. With MAX_LAG_S at 0.35 s and separations to 700 m, every separation
    beyond 534 m was scored that way at 1525 m/s -- roughly the deepest quarter
    of the array, and the quarter carrying the most moveout. Because the effect
    grows as v falls, it manufactures a low-velocity trend out of geometry alone.
    """
    env = np.abs(hilbert(gather, axis=1))
    env = env / np.maximum(np.median(env, axis=1, keepdims=True), 1e-30)
    half_s = 0.012
    half = max(1, int(half_s / (lags[1] - lags[0])))
    lag_lo, lag_hi = float(lags[0]), float(lags[-1])
    c, a = [], []
    for row, s in zip(env, sep):
        for sign, out in ((+1.0, c), (-1.0, a)):
            centre = sign * s / v
            if centre - half_s < lag_lo or centre + half_s > lag_hi:
                out.append(np.nan)
                continue
            j = int(np.argmin(np.abs(lags - centre)))
            out.append(row[j - half:j + half + 1].mean())
    if not np.isfinite(c).all() or not np.isfinite(a).all():
        raise SystemExit(
            "lag window +-%.2f s cannot hold separations to %.0f m at %.0f m/s "
            "(needs %.3f s).\n  Widen MAX_LAG_S; do not score a subset of the "
            "offsets, because the survivors are the near ones and they bias the "
            "score upward." % (max(abs(lag_lo), abs(lag_hi)), np.max(sep), v,
                               np.max(sep) / v))
    return float(np.median(c)), float(np.median(a))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="deepC")
    ap.add_argument("--band", choices=sorted(BANDS), default="5-20")
    ap.add_argument("--block-hours", type=float, default=4.0)
    ap.add_argument("--stride-hours", type=float, default=24.0,
                    help="start a block every this many hours; > block-hours "
                         "subsamples the archive for coverage at lower cost")
    ap.add_argument("--task", type=int, default=0)
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--combine", action="store_true")
    a = ap.parse_args()
    band = BANDS[a.band]
    tag = a.band.replace(".", "p")

    if a.combine:
        return combine(a, band, tag)

    rx_ch, rx_z, z0 = receivers()
    sep = rx_z - z0
    rows = ex.deep_rows(a.arm)
    per_block = int(a.block_hours * 60)          # 60 s records
    stride = int(a.stride_hours * 60)
    starts = list(range(0, max(1, len(rows) - per_block + 1), stride))
    mine = starts[a.task::a.ntasks]
    print("arm %s: %.1f h available, %d blocks of %.1f h every %.1f h; this task has %d"
          % (a.arm, len(rows) * 60 / 3600, len(starts), a.block_hours,
             a.stride_hours, len(mine)), flush=True)
    print("  band %s Hz | source ch %d (TVD %.0f m) | %d receivers to %.0f m depth"
          % (a.band, SOURCE_CH, z0, rx_ch.size, rx_z.max()), flush=True)

    out = []
    acc_sum = None; nw_sum = 0; meta = None
    for n, s0 in enumerate(mine):
        blk = rows.iloc[s0:s0 + per_block]
        paths = [ex.corrected_path(r) for r in blk.file]
        try:
            acc, nw, fs, n_fft = block_spectra(paths, rx_ch, band)
        except Exception as exc:
            print("  block %d (%s) failed: %s" % (n, blk.time.iloc[0], exc), flush=True)
            continue
        g, lags = to_gather(acc, nw, fs, n_fft, band)
        c, ac = score(g, lags, sep, V_ARR)
        cl, al = score(g, lags, sep, V_LELL)
        out.append(dict(t=blk.time.iloc[0].value, hours=len(blk) * 60 / 3600.0,
                        nw=nw, c=c, ac=ac, ratio=c / ac if ac else np.nan,
                        c_lell=cl, ratio_lell=cl / al if al else np.nan))
        if acc_sum is None:
            acc_sum = acc.copy(); meta = (fs, n_fft)
        else:
            acc_sum += acc
        nw_sum += nw
        print("  %s | %.1f h | c/a@%.0f = %.2f | c/a@%.0f = %.2f"
              % (blk.time.iloc[0], out[-1]["hours"], V_ARR, out[-1]["ratio"],
                 V_LELL, out[-1]["ratio_lell"]), flush=True)
    if not out:
        raise SystemExit("no blocks processed")
    np.savez_compressed("%s_%s_%s_t%02d.npz" % (STEM, a.arm, tag, a.task),
                        acc_real=acc_sum.real, acc_imag=acc_sum.imag,
                        n_windows=nw_sum, fs=meta[0], n_fft=meta[1],
                        rx_ch=rx_ch, rx_z=rx_z, sep=sep, source_tvd=z0,
                        t=np.array([r["t"] for r in out]),
                        hours=np.array([r["hours"] for r in out]),
                        c=np.array([r["c"] for r in out]),
                        ac=np.array([r["ac"] for r in out]),
                        ratio=np.array([r["ratio"] for r in out]),
                        ratio_lell=np.array([r["ratio_lell"] for r in out]),
                        band=np.array(band))
    print("wrote %s_%s_%s_t%02d.npz (%d blocks, %.1f h)"
          % (STEM.name, a.arm, tag, a.task, len(out),
             sum(r["hours"] for r in out)), flush=True)


def combine(a, band, tag):
    import pandas as pd
    files = sorted(glob.glob("%s_%s_%s_t*.npz" % (STEM, a.arm, tag)))
    if not files:
        raise SystemExit("no task products for arm %s band %s" % (a.arm, tag))
    acc = None; nw = 0; t = []; hrs = []; rat = []; ratl = []; cs = []
    for f in files:
        d = np.load(f, allow_pickle=True)
        A = d["acc_real"] + 1j * d["acc_imag"]
        acc = A.copy() if acc is None else acc + A
        nw += int(d["n_windows"])
        t += list(d["t"]); hrs += list(d["hours"])
        rat += list(d["ratio"]); ratl += list(d["ratio_lell"]); cs += list(d["c"])
        fs = float(d["fs"]); n_fft = int(d["n_fft"]); sep = d["sep"]; rx_z = d["rx_z"]
    order = np.argsort(t)
    t = np.array(t)[order]; hrs = np.array(hrs)[order]
    rat = np.array(rat)[order]; ratl = np.array(ratl)[order]; cs = np.array(cs)[order]
    total = float(hrs.sum())
    g, lags = to_gather(acc, nw, fs, n_fft, band)
    c, ac_ = score(g, lags, sep, V_ARR)
    cl, al = score(g, lags, sep, V_LELL)

    log = ["Deep-fibre ambient CC, arm %s, band %s Hz" % (a.arm, a.band),
           "  GRAND STACK: %.1f h in %d blocks, %d windows" % (total, hrs.size, nw),
           "  source channel %d (TVD %.0f m), receivers to %.0f m depth"
           % (SOURCE_CH, float(np.min(rx_z) - sep[0]), float(rx_z.max())),
           "",
           "  causal/acausal at %.0f m/s : %.3f  (causal %.4f)" % (V_ARR, c / ac_ if ac_ else np.nan, c),
           "  causal/acausal at %.0f m/s : %.3f  (Lellouch velocity)" % (V_LELL, cl / al if al else np.nan),
           "",
           "  per-block c/a at %.0f m/s : median %.3f, IQR %.3f-%.3f, max %.3f"
           % (V_ARR, np.nanmedian(rat), np.nanpercentile(rat, 25),
              np.nanpercentile(rat, 75), np.nanmax(rat))]
    tt = pd.to_datetime(t).tz_localize("UTC").tz_convert("America/Los_Angeles")
    log += ["  span %s to %s local" % (tt[0].strftime("%Y-%m-%d %H:%M"),
                                       tt[-1].strftime("%Y-%m-%d %H:%M"))]
    print("\n".join(log), flush=True)

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6), constrained_layout=True)
    fig.suptitle("Deep fibre, arm %s, %.0f h stacked, %s Hz   |   %s to %s  (local)"
                 % (a.arm, total, a.band, tt[0].strftime("%a %d %b %Y"),
                    tt[-1].strftime("%a %d %b %Y")), fontsize=11)
    fig.text(0.995, 0.002, "UTC: %s to %s"
             % (pd.to_datetime(t[0]).strftime("%Y-%m-%d %H:%M"),
                pd.to_datetime(t[-1]).strftime("%Y-%m-%d %H:%M")),
             ha="right", va="bottom", fontsize=6.5, color="#9a9a9a")
    ax[0].plot(tt, rat, "o-", color=C1, ms=3, lw=0.8, label="%.0f m/s" % V_ARR)
    ax[0].plot(tt, ratl, "s-", color=C2, ms=3, lw=0.8, alpha=.7,
               label="%.0f m/s (Lellouch)" % V_LELL)
    ax[0].axhline(1.0, color=INK, ls="--", lw=1.0)
    ax[0].set(xlabel="date and time, LOCAL (America/Los_Angeles)",
              ylabel="causal / acausal",
              title="(a) Directionality over time")
    ax[0].legend(); ax[0].tick_params(axis="x", rotation=30, labelsize=7)
    ax[0].grid(alpha=.3)
    ax[1].hist(rat[np.isfinite(rat)], bins=30, color=C1, alpha=.85)
    ax[1].axvline(1.0, color=INK, ls="--", lw=1.2, label="balanced")
    ax[1].set(xlabel="causal/acausal at %.0f m/s" % V_ARR, ylabel="blocks",
              title="(b) Distribution over %d blocks" % hrs.size)
    ax[1].legend(); ax[1].grid(alpha=.3)
    tn = g / np.maximum(np.abs(g).max(axis=1, keepdims=True), 1e-30)
    lim = float(np.percentile(np.abs(tn), 97.0))
    ax[2].imshow(tn, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                 interpolation="nearest",
                 extent=[lags[0], lags[-1], rx_z[-1], rx_z[0]])
    ax[2].plot(sep / V_ARR, rx_z, "-", color=C3, lw=1.5, label="%.0f m/s" % V_ARR)
    _XLIM = av.plot_lag_limit(lags, float(np.max(sep)), V_ARR)
    ax[2].set(xlim=(-_XLIM, _XLIM), xlabel="correlation lag (s)",
              ylabel="true vertical depth (m)",
              title="(c) Grand stack, %.0f h" % total)
    ax[2].legend(); ax[2].grid(False)
    fig.savefig("%s_%s_%s.png" % (STEM, a.arm, tag), dpi=300, bbox_inches="tight")
    np.savez_compressed("%s_%s_%s.npz" % (STEM, a.arm, tag), gather=g, lags=lags,
                        rx_z=rx_z, sep=sep, t=t, hours=hrs, ratio=rat,
                        ratio_lell=ratl, total_hours=total, n_windows=nw,
                        band=np.array(band))
    Path("%s_%s_%s.txt" % (STEM, a.arm, tag)).write_text("\n".join(log) + "\n")
    print("wrote %s_%s_%s.{npz,png,txt}" % (STEM.name, a.arm, tag), flush=True)


if __name__ == "__main__":
    main()
