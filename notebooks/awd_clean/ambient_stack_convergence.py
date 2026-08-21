#!/usr/bin/env python3
"""Does the moveout score grow like sqrt(N) as the stack lengthens?

THE QUESTION THIS ANSWERS, and it is the right question to ask of an
ambient-noise result.

Ambient-noise interferometry converges toward the Green's function as independent
noise realisations are averaged, and the signal-to-noise of the reconstruction
grows as sqrt(T). That is why people stack for months, and it is the reason a
year-long archive ought to beat Lellouch et al. (2019)'s single day. This script
tests whether it does, by stacking 1, 2, 4, ... 24 hourly chunks of the same day
and scoring each partial stack identically.

WHAT THE THREE POSSIBLE SHAPES MEAN

  score grows ~ sqrt(N)    an arrival is present and stacking is working. Extend
                           the stack.
  score flat               the statistic is limited by something that averages
                           exactly like the signal does -- a REPEATABLE
                           contaminant. Stacking cannot help; the contaminant
                           must be removed instead.
  score flat AND the       there is nothing to accumulate. sqrt(N) times zero is
  null-relative            zero, so more data yields a more precise measurement
  detectability flat       of absence, not a detection.

The discriminating quantity is therefore NOT the raw score but the score divided
by its own null threshold ("detectability"). A raw score can rise simply because
a coherent contaminant accumulates; detectability rises only if the observation
outruns what permuted receiver order produces from the same data.

BOTH BRANCHES ARE RUN. The baseline (no common-mode removal) is dominated by the
static pattern, so its raw score should climb while its detectability does not --
the signature of accumulating junk. The common-mode-removed branch has the
pattern suppressed (pedestal diagnostic -0.381 against +0.976), so if an arrival
is present at all, this is the branch in which sqrt(N) growth should appear.

A fitted exponent is reported: score ~ N^b. b ~ 0.5 is stacking working, b ~ 0 is
a flat statistic. The fit is on log-log, and the exponent for DETECTABILITY is the
one to read.

Reads the existing hourly chunk products only -- no raw data, no reprocessing.

Output: ambient_stack_convergence.{npz,png,txt}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, hilbert, sosfiltfilt

import ambient_lellouch2019_exact_stack as X

import safod_geometry as geo
HERE = Path(__file__).resolve().parent
STEM = HERE / "ambient_stack_convergence"
CHUNKS = HERE / "ambient_transfer" / "lellouch2019_exact_stack"
DATE = "2024-12-20"
NULL_COUNT = 2000
SEED = 20260819
V_REF = 3200.0
GATE_S = 0.012


def chunk_files(suffix: str) -> list[Path]:
    pat = ("chunk_%s_src23_ram0p1_cross_correlation_ordered_r0%s_start*.npz"
           % (DATE, suffix))
    return sorted(CHUNKS.glob(pat))


def combine(paths):
    """Sum the stored cross-spectra, exactly as the aggregate step does."""
    central = neighbor = None
    windows = 0
    meta = None
    for p in paths:
        with np.load(p, allow_pickle=False) as d:
            if central is None:
                central = np.zeros_like(d["central_cross_spectrum_sum"])
                neighbor = np.zeros_like(d["neighbor_cross_spectrum_sum"])
                meta = dict(n_fft=int(d["n_fft"]), fs=float(d["fs_hz"]),
                            dx=float(d["dx_m"]),
                            offsets=d["offsets_m"].copy())
            central += d["central_cross_spectrum_sum"]
            neighbor += d["neighbor_cross_spectrum_sum"]
            windows += int(d["n_windows"])
    return central, neighbor, windows, meta


def correlation(spectrum_sum, windows, meta):
    """Average, inverse-transform, band-pass the FULL correlation, then crop."""
    avg = spectrum_sum / max(1, windows)
    full = np.fft.irfft(avg, n=meta["n_fft"], axis=-1)
    full = np.fft.fftshift(full, axes=-1)
    fs = meta["fs"]
    sos = butter(4, list(X.OUTPUT_BAND_HZ), btype="bandpass", fs=fs, output="sos")
    full = sosfiltfilt(sos, full, axis=-1)
    lags = (np.arange(meta["n_fft"]) - meta["n_fft"] // 2) / fs
    keep = np.abs(lags) <= X.MAX_LAG_SECONDS
    return full[:, keep], lags[keep]


def score_curve(gather, lags, offsets, grid, sign=1.0):
    env = np.abs(hilbert(gather, axis=1))
    env = env / np.maximum(np.median(env, axis=1, keepdims=True), 1e-30)
    half = max(1, int(GATE_S / (lags[1] - lags[0])))
    out = np.empty(grid.size)
    for i, v in enumerate(grid):
        vals = []
        for row, x in zip(env, offsets):
            k = int(np.argmin(np.abs(lags - sign * x / v)))
            lo, hi = max(0, k - half), min(len(lags), k + half + 1)
            vals.append(row[lo:hi].mean())
        out[i] = np.median(vals)
    return out


def evaluate(gather, lags, offsets, rng):
    grid = X.VELOCITY_GRID_M_S
    cs = score_curve(gather, lags, offsets, grid)
    obs = float(cs.max())
    coarse = grid[::4]
    nulls = np.empty(NULL_COUNT)
    for i in range(NULL_COUNT):
        perm = rng.permutation(len(offsets))
        nulls[i] = score_curve(gather[perm], lags, offsets, coarse).max()
    n95 = float(np.percentile(nulls, 95))
    return dict(obs=obs, at=float(grid[int(np.argmax(cs))]),
                null95=n95, detect=obs / n95 if n95 else np.nan,
                p=float((np.sum(nulls >= obs) + 1) / (NULL_COUNT + 1)),
                pedestal=float(np.corrcoef(grid, cs)[0, 1]))


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    rng = np.random.default_rng(SEED)

    say("Stack convergence: does the score grow like sqrt(N)?")
    say("  %s, hourly chunks, both branches" % DATE)
    say("  The quantity to read is DETECTABILITY = score / its own null 95th.")
    say("  A raw score can rise merely because a repeatable contaminant")
    say("  accumulates; detectability rises only if the observation outruns what")
    say("  permuted receiver order produces from the same data.")
    say("")

    results = {}
    for suffix, label in (("", "baseline (no common-mode removal)"),
                          ("_cm", "common-mode removed")):
        files = chunk_files(suffix)
        if len(files) < 4:
            say("  %-34s only %d chunks, skipped" % (label, len(files)))
            continue
        say("=== %s: %d hourly chunks ===" % (label, len(files)))
        counts, rows = [], []
        n = 1
        while n <= len(files):
            central, neighbor, windows, meta = combine(files[:n])
            gather, lags = correlation(neighbor, windows, meta)
            r = evaluate(gather, lags, meta["offsets"], rng)
            counts.append(n); rows.append(r)
            say("  %2d chunk(s) %6d windows | score %7.4f at %4.0f m/s | "
                "null95 %7.4f | detect %.3f | p %.4f | pedestal %+.3f"
                % (n, windows, r["obs"], r["at"], r["null95"], r["detect"],
                   r["p"], r["pedestal"]))
            n *= 2
        if counts[-1] != len(files):
            central, neighbor, windows, meta = combine(files)
            gather, lags = correlation(neighbor, windows, meta)
            r = evaluate(gather, lags, meta["offsets"], rng)
            counts.append(len(files)); rows.append(r)
            say("  %2d chunk(s) %6d windows | score %7.4f at %4.0f m/s | "
                "null95 %7.4f | detect %.3f | p %.4f | pedestal %+.3f"
                % (len(files), windows, r["obs"], r["at"], r["null95"],
                   r["detect"], r["p"], r["pedestal"]))
        c = np.array(counts, float)
        raw = np.array([r["obs"] for r in rows])
        det = np.array([r["detect"] for r in rows])
        b_raw = float(np.polyfit(np.log(c), np.log(raw), 1)[0])
        ok = np.isfinite(det) & (det > 0)
        b_det = float(np.polyfit(np.log(c[ok]), np.log(det[ok]), 1)[0]) if ok.sum() > 2 else np.nan
        say("  fitted exponent, raw score      : N^%+.3f" % b_raw)
        say("  fitted exponent, DETECTABILITY  : N^%+.3f   (0.5 = stacking works,"
            " 0 = flat)" % b_det)
        say("")
        results[label] = dict(counts=c, raw=raw, det=det,
                              p=np.array([r["p"] for r in rows]),
                              ped=np.array([r["pedestal"] for r in rows]),
                              b_raw=b_raw, b_det=b_det)

    if not results:
        raise SystemExit("no chunk products found under %s" % CHUNKS)

    say("=== reading ===")
    for label, r in results.items():
        grows = r["b_det"] > 0.25
        say("  %-34s detectability exponent %+.3f -> %s"
            % (label, r["b_det"],
               "GROWING: stacking is working, extend it" if grows
               else "FLAT: more data does not help"))
    say("")
    best = max(results.values(), key=lambda r: r["b_det"] if np.isfinite(r["b_det"]) else -9)
    if np.isfinite(best["b_det"]) and best["b_det"] > 0.25:
        say("  At least one branch shows sqrt(N)-like growth. The arrival may be")
        say("  present but under-stacked, and a longer stack is justified.")
    else:
        say("  NEITHER branch shows sqrt(N) growth in detectability. This is the")
        say("  quantitative answer to 'more data should win': it should, and it")
        say("  does not. Averaging suppresses INCOHERENT noise as 1/sqrt(N), so")
        say("  detectability grows as sqrt(N) ONLY if there is a coherent arrival")
        say("  to accumulate. A flat detectability curve means either the limiting")
        say("  contaminant averages exactly as the signal would -- i.e. it repeats")
        say("  window to window -- or there is no arrival to accumulate. The")
        say("  common-mode-removed branch discriminates: with the repeatable")
        say("  pattern suppressed (pedestal %+.3f at full stack), a flat curve"
            % results.get("common-mode removed", {"ped": np.array([np.nan])})["ped"][-1])
        say("  leaves absence as the explanation.")
        say("")
        say("  Consistent with Behm (2016), where 30 s sufficed under adequate")
        say("  illumination, and with the four-day stack (p = 0.9184) sitting")
        say("  further from significance than the best single day (p = 0.1345).")
    say("")
    say("  LIMITS: one day, hourly granularity, so this measures within-day")
    say("  convergence; the four-day result covers the across-day case. The null")
    say("  is receiver-order permutation, which does not contain any pre-")
    say("  correlation operator.")

    # PI convention: hours stacked and the LOCAL date range belong in the title,
    # not just in the log, so a reader can check the run against site activity.
    _hrs = max(max(r["counts"]) for r in results.values()) if results else 0.0
    _t0 = pd.Timestamp(DATE + " 00:00", tz="UTC")
    _title, _foot = geo.figure_label(_t0, _t0 + pd.Timedelta(hours=float(_hrs)),
                                     float(_hrs), fibre="Nano (2024-25 archive)",
                                     extra="stack convergence, hourly chunks")
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6), constrained_layout=True)
    cols = {"baseline (no common-mode removal)": "#D55E00",
            "common-mode removed": "#0072B2"}
    for label, r in results.items():
        c = cols.get(label, "k")
        ax[0].loglog(r["counts"], r["raw"], "o-", color=c, ms=5, label=label)
        ax[1].semilogx(r["counts"], r["det"], "o-", color=c, ms=5,
                       label="%s (N^%+.2f)" % (label, r["b_det"]))
        ax[2].semilogx(r["counts"], r["p"], "o-", color=c, ms=5, label=label)
    ref = np.array([1, 24], float)
    first = list(results.values())[0]
    ax[0].loglog(ref, first["raw"][0] * np.sqrt(ref), "k--", lw=1.2,
                 label=r"$\sqrt{N}$ reference")
    ax[0].set(xlabel="hourly chunks stacked", ylabel="raw moveout score",
              title="(a) Raw score\n(can rise from accumulating junk)")
    ax[0].legend(fontsize=7); ax[0].grid(alpha=.3, which="both")
    ax[1].axhline(1.0, color="k", ls="--", lw=1.2, label="detection threshold")
    ax[1].set(xlabel="hourly chunks stacked",
              ylabel="detectability = score / null 95th",
              title="(b) THE ONE TO READ\nflat = more data does not help")
    ax[1].legend(fontsize=7); ax[1].grid(alpha=.3, which="both")
    ax[2].axhline(0.05, color="k", ls="--", lw=1.2, label="p = 0.05")
    ax[2].set(xlabel="hourly chunks stacked", ylabel="p value", yscale="log",
              title="(c) Significance vs stack length")
    ax[2].legend(fontsize=7); ax[2].grid(alpha=.3, which="both")
    fig.suptitle(_title, fontsize=10.5)
    fig.text(0.995, 0.002, _foot, ha="right", va="bottom",
             fontsize=6.5, color="#9a9a9a")
    fig.savefig(str(STEM) + ".png", dpi=190, bbox_inches="tight")

    np.savez(str(STEM) + ".npz",
             **{("%s_%s" % (k.split()[0], f)): v[f]
                for k, v in results.items()
                for f in ("counts", "raw", "det", "p", "ped")})
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
