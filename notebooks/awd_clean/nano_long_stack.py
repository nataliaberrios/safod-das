#!/usr/bin/env python3
"""Long Nano stack: the FIRST test of the cemented fibre from a downhole source.

WHY, REWRITTEN 2026-08-20. Everything this file previously rested on is void.
The 3.0 h "p = 0.0474 at 2950 m/s" hint AND the 25 h "N^+0.023" null that
superseded it were both computed with the virtual source at channel 10, which is
in the AIR -- nano_find_wellhead.py puts fibre entry at channel 73. An air source
correlated against downhole channels returns only what they share (instrumental
common mode at zero lag) and cannot show moveout whether or not moveout exists,
so it manufactures exactly the null it appeared to measure. Neither result was
evidence. THE CEMENTED FIBRE HAS NEVER BEEN TESTED FROM A DOWNHOLE SOURCE; this
run is the first.

The Deep comparison in the old text is also retracted: that arrival is at
1350-1550 m/s, not 1675 (see arrival_velocities.py). The point it made still
stands -- a guided mode in borehole fluid is a far easier target than a body wave
through rock, so a null here is weaker evidence than a null there.

There are ~25 h of Nano record before the first weight drop. This stacks all of
it and asks the one question that separates a weak real arrival from a fluctuation:

    DOES DETECTABILITY AT 3200 m/s GROW LIKE sqrt(N)?

A coherent arrival accumulating against incoherent noise grows as sqrt(N) in
detectability -- the score divided by the 95th percentile of its OWN null rebuilt
at the same stack length. A fluctuation does not. This is the same diagnostic
that showed the 2024-25 main-hole fibre contains nothing (N^+0.042 where an
arrival requires N^+0.50), applied here to a pre-specified velocity.

3200 m/s IS PRE-SPECIFIED, and is not a number chosen after seeing this data: it
is Lellouch et al. (2019)'s own value, and independently it is what this fibre's
measured dispersion extrapolates to inside 5-20 Hz (2950 m/s at 35-80 Hz and
3300 at 15-25 Hz, slowness trend 0.483 us/m/Hz -> 3102-3238 m/s in band). The
earlier a-priori of 2950 m/s was a 30-60 Hz number applied to a 5-20 Hz analysis,
which was simply the wrong target. The scan-wide per-velocity result is reported
alongside but is the weaker statement.

MEMORY. 25 h of 732 channels at 250 Hz is ~66 GB as one float64 array, and the
pipeline needs several copies. Cross-spectra are therefore ACCUMULATED batch by
batch and the raw data for each batch is released before the next is read, so
peak memory is set by the batch, not by the total stack. This also makes the
convergence curve free: the partial sums after each batch are the shorter stacks.

Output: nano_long_stack.{npz,png,txt}
"""
from __future__ import annotations

import argparse
import sys
from datetime import timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, hilbert, sosfiltfilt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ambient_lellouch2019_exact_stack as ex
import nano_ambient_cc as nano
import arrival_velocities as av

STEM = HERE / "nano_long_stack"
LA = ZoneInfo("America/Los_Angeles")

# A-PRIORI VELOCITY, corrected 2026-08-20. This was 2950 m/s, but
# nano_mode_identification.txt measures 2950 at 35-80 Hz and 3300 at 15-25 Hz,
# a slowness trend of 0.483 us/m/Hz that extrapolates to 3102-3238 m/s inside
# the 5-20 Hz analysis band. 2950 was therefore the WRONG a-priori target for
# this band. The in-band expectation coincides with Lellouch's 3200 m/s.
V_ACTIVE = av.V_NANO_INBAND
NANO_WELLHEAD = 73         # nano_find_wellhead.py; 0-72 are surface/air
V_LELL = 3200.0
OFFSETS_M = np.arange(50.0, 700.1, 50.0)
NEIGHBOUR = 10
WINDOW_S, STEP_S = 30.0, 15.0
RAM_S = 0.1
BAND = (5.0, 20.0)
# Derived: must hold the largest offset at the slowest velocity SCORED. The old
# 0.35 s could not reach 700 m below 1934 m/s, so most of V_GRID was scored on a
# shrinking subset of near offsets -- an upward bias that grows as v falls.
MAX_LAG_S = None                # set below, once V_GRID exists
V_GRID = np.arange(300.0, 6000.1, 25.0)
MAX_LAG_S = round(av.required_lag_s(OFFSETS_M.max(), V_GRID.min()), 2)
GATE_S = 0.012
NULLS = 400
SEED = 20260820


def moveout_curve(gather, lags, offsets, sign=1.0):
    env = np.abs(hilbert(gather, axis=1))
    env = env / np.maximum(np.median(env, axis=1, keepdims=True), 1e-30)
    half = max(1, int(GATE_S / (lags[1] - lags[0])))
    lo_lag, hi_lag = float(lags[0]), float(lags[-1])
    out = np.empty(V_GRID.size)
    for i, v in enumerate(V_GRID):
        centres = sign * np.asarray(offsets, dtype=float) / v
        # score a velocity only if EVERY offset's gate fits: argmin would other-
        # wise clamp to the window edge, and the offsets that survive are the
        # near ones, which sit closest to the zero-lag lobe and read high
        if (centres - GATE_S < lo_lag).any() or (centres + GATE_S > hi_lag).any():
            out[i] = np.nan
            continue
        vals = [row[max(0, k - half):k + half + 1].mean()
                for row, k in zip(env, (np.abs(lags[None, :] - centres[:, None])
                                        ).argmin(axis=1))]
        out[i] = np.median(vals)
    return out


def spectra_to_gather(acc, n_used, n_fft, fs):
    avg = acc / max(1, n_used)
    full = np.fft.fftshift(np.fft.irfft(avg, n=n_fft, axis=-1), axes=-1)
    lags = (np.arange(n_fft) - n_fft // 2) / fs
    sos = butter(4, list(BAND), btype="bandpass", fs=fs, output="sos")
    keep = np.abs(lags) <= MAX_LAG_S
    return sosfiltfilt(sos, full, axis=-1)[:, keep], lags[keep]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfiles", type=int, default=300, help="5-min records (300 = 25 h)")
    ap.add_argument("--batch", type=int, default=6, help="records read at once")
    ap.add_argument("--source", type=int, default=NANO_WELLHEAD,
                    help="virtual source; must be at or below the wellhead")
    ap.add_argument("--allow-air-source", action="store_true",
                    help=argparse.SUPPRESS)     # only to reproduce the void runs
    a = ap.parse_args()

    # A 25 h run finished on 2026-08-20 with --source 10 and was reported as a
    # null result for the cemented fibre. Channel 10 is in the AIR: the Nano
    # fibre enters the hole at channel 73 (nano_find_wellhead.py, from a step in
    # RMS, a step in neighbour coherence, and the HF ratio). Correlating an air
    # channel against downhole channels returns only what they share, which is
    # instrumental common mode at zero lag, and CANNOT show moveout whether or
    # not moveout exists. So that run was not evidence either way, and neither
    # was the p = 0.0474 it superseded. Refusing rather than warning, because
    # the run costs 1.5 node hours and reads as a result.
    if a.source < NANO_WELLHEAD and not a.allow_air_source:
        raise SystemExit(
            "REFUSING to run: channel %d is in the AIR.\n"
            "  The Nano fibre enters the hole at channel %d; a virtual source\n"
            "  above that cannot produce moveout and the output would look\n"
            "  like a null result rather than an invalid one." % (a.source, NANO_WELLHEAD))

    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    rng = np.random.default_rng(SEED)

    files = sorted(nano.NANO_DIR.glob("*.pb"))
    pre = [(p, t) for p, t in ((p, nano.nano_time(p.name)) for p in files)
           if t is not None and t < nano.FIRST_DROP]
    use = pre[: a.nfiles]
    hours_total = len(use) * 300 / 3600.0
    l0 = use[0][1].astimezone(LA); l1 = use[-1][1].astimezone(LA)
    say("Nano long stack (CEMENTED fibre), source channel %d" % a.source)
    say("  %.2f h available before the first weight drop" % hours_total)
    say("  %s to %s LOCAL" % (l0.strftime("%Y-%m-%d %H:%M"), l1.strftime("%Y-%m-%d %H:%M")))
    say("  %s to %s UTC" % (use[0][1].strftime("%Y-%m-%d %H:%M"),
                            use[-1][1].strftime("%Y-%m-%d %H:%M")))
    say("  a-priori test velocity %.0f m/s (Lellouch, and this fibre's in-band" % V_ACTIVE)
    say("    dispersion extrapolation; 2950 m/s is a 30-60 Hz value)")
    say("")

    acc = None
    n_used = 0
    n_fft = fs = dx = None
    centres = offs = rows_nb = None
    conv = []          # (hours, detectability at V_ACTIVE, p at V_ACTIVE)

    for i in range(0, len(use), a.batch):
        chunk = [str(p) for p, _ in use[i:i + a.batch]]
        arr, info = nano.readFile_protobuf(chunk, fmin=1.0, fmax=100.0,
                                           desampling=True, **nano.RAW_KW)
        x = np.asarray(arr, dtype=np.float32); del arr
        if fs is None:
            fs = float(info["fs"]); dx = float(info.get("dx", nano.DX_NANO))
            centres = a.source + np.rint(OFFSETS_M / dx).astype(int)
            keep = centres + NEIGHBOUR < x.shape[0]
            centres, offs = centres[keep], OFFSETS_M[keep]
            rows_nb = [np.arange(c - NEIGHBOUR, c + NEIGHBOUR + 1) for c in centres]
            n_win = int(WINDOW_S * fs)
            n_fft = 1 << int(np.ceil(np.log2(2 * n_win - 1)))
            acc = np.zeros((len(rows_nb), n_fft // 2 + 1), dtype=np.complex128)
            say("  fs %.0f Hz, dx %.4f m, %d receivers %.0f-%.0f m"
                % (fs, dx, len(rows_nb), offs[0], offs[-1]))
        n_win, n_step = int(WINDOW_S * fs), int(STEP_S * fs)

        # differentiate + RAM in place, then correlate this batch
        x[:, 1:] = np.diff(x, axis=1) * np.float32(fs)
        x[:, 0] = 0.0
        n_ram = ex.odd_ram_samples(RAM_S, fs)
        for j in range(0, x.shape[0], 64):
            blk = x[j:j + 64]
            w = uniform_filter1d(np.abs(blk), size=n_ram, axis=1, mode="nearest")
            sc = float(np.nanmedian(w)) or 1.0
            np.divide(blk, np.maximum(w, np.float32(np.finfo(np.float32).eps * sc)),
                      out=blk)
            del w
        for s in range(0, x.shape[1] - n_win + 1, n_step):
            S = np.conj(np.fft.rfft(x[a.source, s:s + n_win], n=n_fft))
            for r, rr in enumerate(rows_nb):
                acc[r] += S * np.fft.rfft(x[rr, s:s + n_win].sum(axis=0), n=n_fft)
            n_used += 1
        del x

        hours = (i + len(chunk)) * 300 / 3600.0
        # convergence point: score / its own null 95th, at the a-priori velocity
        gather, lags = spectra_to_gather(acc, n_used, n_fft, fs)
        cs = moveout_curve(gather, lags, offs)
        k = int(np.argmin(np.abs(V_GRID - V_ACTIVE)))
        null_at = np.empty(NULLS)
        for m in range(NULLS):
            g = gather[rng.permutation(len(offs))]
            env = np.abs(hilbert(g, axis=1))
            env = env / np.maximum(np.median(env, axis=1, keepdims=True), 1e-30)
            half = max(1, int(GATE_S / (lags[1] - lags[0])))
            vals = []
            for row, xo in zip(env, offs):
                centre = xo / V_ACTIVE
                if centre - GATE_S < lags[0] or centre + GATE_S > lags[-1]:
                    raise SystemExit(
                        "lag window +-%.2f s cannot hold %.0f m at %.0f m/s"
                        % (MAX_LAG_S, xo, V_ACTIVE))
                kk = int(np.argmin(np.abs(lags - centre)))
                lo, hi = max(0, kk - half), min(len(lags), kk + half + 1)
                vals.append(row[lo:hi].mean())
            null_at[m] = np.median(vals)
        thr = float(np.percentile(null_at, 95.0))
        pv = float((np.sum(null_at >= cs[k]) + 1) / (NULLS + 1))
        conv.append((hours, cs[k] / thr if thr else np.nan, pv, n_used))
        say("  %5.2f h | %5d windows | score@%.0f %.4f | detect %.3f | p %.4f"
            % (hours, n_used, V_ACTIVE, cs[k], cs[k] / thr if thr else np.nan, pv))

    # final full-scan statistics
    gather, lags = spectra_to_gather(acc, n_used, n_fft, fs)
    cs = moveout_curve(gather, lags, offs)
    ac = moveout_curve(gather, lags, offs, sign=-1.0)
    null = np.empty((NULLS, V_GRID.size))
    for m in range(NULLS):
        null[m] = moveout_curve(gather[rng.permutation(len(offs))], lags, offs)
    thresh = np.percentile(null, 95.0, axis=0)
    pv_all = (np.sum(null >= cs[None, :], axis=0) + 1.0) / (NULLS + 1.0)
    k = int(np.argmax(cs))
    ka = int(np.argmin(np.abs(V_GRID - V_ACTIVE)))
    clears = cs > thresh

    conv = np.array([(c[0], c[1], c[2]) for c in conv], dtype=float)
    b = (np.polyfit(np.log(conv[:, 0]), np.log(np.maximum(conv[:, 1], 1e-9)), 1)[0]
         if len(conv) > 2 else np.nan)

    say("")
    say("=== full stack: %.2f h, %d windows ===" % (hours_total, n_used))
    say("  scan peak %.4f at %.0f m/s (p = %.4f)" % (cs[k], V_GRID[k], pv_all[k]))
    say("  A-PRIORI test at %.0f m/s: score %.4f, p = %.4f"
        % (V_ACTIVE, cs[ka], pv_all[ka]))
    say("  at %.0f m/s (Lellouch): p = %.4f"
        % (V_LELL, float(np.interp(V_LELL, V_GRID, pv_all))))
    say("  velocities clearing their own null: %d of %d" % (int(clears.sum()), clears.size))
    if clears.any():
        say("    band %.0f-%.0f m/s" % (V_GRID[clears].min(), V_GRID[clears].max()))
    say("  causal/acausal at %.0f m/s: %.3f"
        % (V_ACTIVE, cs[ka] / ac[ka] if ac[ka] else np.nan))
    say("")
    say("=== convergence at the a-priori velocity ===")
    say("  detectability exponent: N^%+.3f   (0.50 = a coherent arrival, 0 = flat)" % b)
    if b > 0.25:
        say("  GROWING. Consistent with a real arrival accumulating; the 3 h result")
        say("  was under-stacked, not absent.")
    elif b < 0.10:
        say("  FLAT. More hours do not help, which is what an absent arrival looks")
        say("  like -- the same behaviour as the 2024-25 main-hole fibre.")
    else:
        say("  INTERMEDIATE (%+.3f). Not decisive either way." % b)

    loc = "%s %s-%s local" % (l0.strftime("%Y-%m-%d"), l0.strftime("%H:%M"),
                              l1.strftime("%H:%M")) if l0.date() == l1.date() else \
          "%s to %s local" % (l0.strftime("%Y-%m-%d %H:%M"), l1.strftime("%Y-%m-%d %H:%M"))
    utc = "%s to %s UTC" % (use[0][1].strftime("%Y-%m-%d %H:%M"),
                            use[-1][1].strftime("%Y-%m-%d %H:%M"))

    INK, C1, C2, C3 = "#444444", "#0072B2", "#D55E00", "#009E73"
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6), constrained_layout=True)
    fig.suptitle("Nano fibre (cemented), %.1f h stacked   |   %s   (%s)"
                 % (hours_total, loc, utc), fontsize=10.5, y=1.02)
    g = gather / np.maximum(np.abs(gather).max(axis=1, keepdims=True), 1e-30)
    for o, row in zip(offs, g):
        ax[0].plot(lags, -o + row * 40, "-", color=INK, lw=0.7)
    ax[0].plot(offs / V_ACTIVE, -offs, "-", color=C2, lw=1.5,
               label="%.0f m/s (active source)" % V_ACTIVE)
    ax[0].plot(offs / V_LELL, -offs, "--", color=C1, lw=1.3,
               label="%.0f m/s (Lellouch)" % V_LELL)
    ax[0].set(xlim=(-0.35, 0.35), xlabel="correlation lag (s)",
              ylabel="offset from virtual source (m)", title="(a) Gather")
    ax[0].legend(fontsize=7)
    ax[1].plot(V_GRID / 1e3, cs, color=INK, lw=1.6, label="causal")
    ax[1].plot(V_GRID / 1e3, ac, color="#8a8a8a", lw=1.0, label="acausal")
    ax[1].plot(V_GRID / 1e3, thresh, "--", color=C2, lw=1.2, label="null 95th")
    ax[1].fill_between(V_GRID / 1e3, cs, thresh, where=clears, color=C1, alpha=.22)
    ax[1].axvline(V_ACTIVE / 1e3, color=C3, ls=":", lw=1.5)
    ax[1].set(xlabel="trial velocity (km/s)", ylabel="moveout score",
              title="(b) Scan, p=%.4f at %.0f m/s" % (pv_all[ka], V_ACTIVE))
    ax[1].legend(fontsize=7)
    ax[2].semilogx(conv[:, 0], conv[:, 1], "o-", color=C1, ms=5,
                   label="detectability at %.0f m/s" % V_ACTIVE)
    ax[2].axhline(1.0, color=C2, ls="--", lw=1.2, label="detection threshold")
    ax[2].set(xlabel="hours stacked", ylabel="score / own null 95th",
              title="(c) Convergence: N^%+.3f" % b)
    ax[2].legend(fontsize=7)
    fig.savefig(str(STEM) + ".png", dpi=300, bbox_inches="tight")

    np.savez_compressed(str(STEM) + ".npz", gather=gather, lags=lags, offsets=offs,
                        v_grid=V_GRID, causal=cs, acausal=ac, thresh=thresh,
                        p_per_velocity=pv_all, conv=conv, exponent=b,
                        hours=hours_total, n_windows=n_used, fs=fs, dx=dx,
                        source_channel=a.source, t_first=str(use[0][1]),
                        t_last=str(use[-1][1]))
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
