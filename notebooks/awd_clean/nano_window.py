#!/usr/bin/env python3
"""Nano over one chosen LOCAL-time window, so crew-active hours can be isolated.

WHY. The 25 h Nano stack is flat, which is surprising: the same window on the
wireline Deep fibre gives causal/acausal 1.1-4.0 in 12 of 12 blocks, and the
CEMENTED fibre is the one that ought to be the better sensor for a downgoing body
wave, because it is tied to the formation rather than hanging in fluid.

A 25 h grand stack averages over everything. On 2026-06-15 a crew arrived and
spent the afternoon and evening drilling, digging and deploying nodes -- and that
surface work is the plausible source for a downgoing wavefield. Averaging it
together with a quiet night can only dilute it.

So this stacks ONE local-time window, and is meant to be run as a PAIR: an
active window and a quiet one of the SAME LENGTH on the same fibre with the same
processing. Matched duration matters -- comparing 5 h against 20 h would confound
illumination with stack depth, which is the error this project has made before.

NOTE ON COVERAGE. Nano recording begins 2026-06-15 15:48 local, so there is no
morning of the 15th to look at. A "09:00-21:00" request yields 15:48-20:58.

Output: nano_window_<tag>.{npz,png,txt}
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.signal import hilbert

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ambient_lellouch2019_exact_stack as ex
import arrival_velocities as av
import nano_ambient_cc as nano
import nano_long_stack as nls

LA = ZoneInfo("America/Los_Angeles")
V_BODY = av.V_NANO_INBAND        # 3200 m/s, the in-band expectation
V_DEEPMODE = 1400.0              # what the Deep recovers, tested here too
NULLS = 400
INK, C1, C2, C3 = "#444444", "#0072B2", "#D55E00", "#009E73"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True, help="local, e.g. 2026-06-15T15:00")
    ap.add_argument("--end", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--source", type=int, default=nls.NANO_WELLHEAD)
    ap.add_argument("--max-records", type=int, default=0,
                    help="truncate to this many, so a pair can be length-matched")
    a = ap.parse_args()
    if a.source < nls.NANO_WELLHEAD:
        raise SystemExit("channel %d is in the AIR; the hole starts at %d"
                         % (a.source, nls.NANO_WELLHEAD))
    stem = HERE / ("nano_window_%s" % a.tag)

    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    t0 = datetime.fromisoformat(a.start).replace(tzinfo=LA)
    t1 = datetime.fromisoformat(a.end).replace(tzinfo=LA)
    files = sorted(nano.NANO_DIR.glob("*.pb"))
    use = [(p, t) for p, t in ((p, nano.nano_time(p.name)) for p in files)
           if t is not None and t0 <= t.astimezone(LA) < t1 and t < nano.FIRST_DROP]
    if not use:
        raise SystemExit("no Nano records in %s .. %s local" % (a.start, a.end))
    if a.max_records:
        use = use[: a.max_records]
    hours = len(use) * 300 / 3600.0

    say("Nano window '%s', source channel %d (the wellhead)" % (a.tag, a.source))
    say("  requested %s to %s local" % (a.start, a.end))
    say("  got %d records = %.2f h, %s to %s LOCAL"
        % (len(use), hours,
           use[0][1].astimezone(LA).strftime("%a %d %b %H:%M"),
           use[-1][1].astimezone(LA).strftime("%a %d %b %H:%M")))
    say("")

    acc = None; nw = 0; fs = dx = n_fft = None; centres = offs = None
    for i in range(0, len(use), 6):
        arr, info = nano.readFile_protobuf([str(p) for p, _ in use[i:i + 6]],
                                           fmin=1.0, fmax=100.0, desampling=True,
                                           **nano.RAW_KW)
        x = np.asarray(arr, dtype=np.float32); del arr
        if fs is None:
            fs = float(info["fs"]); dx = float(info.get("dx", nano.DX_NANO))
            centres = a.source + np.rint(nls.OFFSETS_M / dx).astype(int)
            keep = centres + nls.NEIGHBOUR < x.shape[0]
            centres, offs = centres[keep], nls.OFFSETS_M[keep]
            n_win = int(nls.WINDOW_S * fs)
            n_fft = 1 << int(np.ceil(np.log2(2 * n_win - 1)))
            acc = np.zeros((len(centres), n_fft // 2 + 1), dtype=np.complex128)
            say("  fs %.0f Hz, dx %.4f m, %d receivers %.0f-%.0f m"
                % (fs, dx, len(centres), offs[0], offs[-1]))
        n_win, n_step = int(nls.WINDOW_S * fs), int(nls.STEP_S * fs)
        x[:, 1:] = np.diff(x, axis=1) * np.float32(fs); x[:, 0] = 0.0
        n_ram = ex.odd_ram_samples(nls.RAM_S, fs)
        for k in range(0, x.shape[0], 64):
            b = x[k:k + 64]
            w = uniform_filter1d(np.abs(b), size=n_ram, axis=1, mode="nearest")
            np.divide(b, np.maximum(w, np.float32(np.finfo(np.float32).eps)), out=b)
            del w
        for s in range(0, x.shape[1] - n_win + 1, n_step):
            S = np.conj(np.fft.rfft(x[a.source, s:s + n_win], n=n_fft))
            for j, c in enumerate(centres):
                seg = x[c - nls.NEIGHBOUR:c + nls.NEIGHBOUR + 1, s:s + n_win].sum(axis=0)
                acc[j] += S * np.fft.rfft(seg, n=n_fft)
            nw += 1
        del x

    gather, lags = nls.spectra_to_gather(acc, nw, n_fft, fs)
    cs = nls.moveout_curve(gather, lags, offs)
    say("  %d windows stacked" % nw)
    say("")

    rng = np.random.default_rng(nls.SEED)
    null = np.empty((NULLS, nls.V_GRID.size))
    for m in range(NULLS):
        null[m] = nls.moveout_curve(gather[rng.permutation(len(offs))], lags, offs)
    thr = np.nanpercentile(null, 95.0, axis=0)
    pv = (np.nansum(null >= cs[None, :], axis=0) + 1.0) / (NULLS + 1.0)
    clears = np.flatnonzero(np.isfinite(cs) & (cs > thr))

    env = np.abs(hilbert(gather, axis=1))
    env = env / np.maximum(np.median(env, axis=1, keepdims=True), 1e-30)
    half = max(1, int(nls.GATE_S / (lags[1] - lags[0])))
    def ca(v):
        c, ac = [], []
        for row, x_ in zip(env, offs):
            for sign, out in ((+1.0, c), (-1.0, ac)):
                k = int(np.argmin(np.abs(lags - sign * x_ / v)))
                out.append(row[k - half:k + half + 1].mean())
        return float(np.median(c)) / max(float(np.median(ac)), 1e-30)

    k3 = int(np.argmin(np.abs(nls.V_GRID - V_BODY)))
    say("=== %.2f h, %d windows ===" % (hours, nw))
    say("  a-priori %0.f m/s : score %.4f, per-velocity p = %.4f, c/a = %.3f"
        % (V_BODY, cs[k3], pv[k3], ca(V_BODY)))
    say("  Deep mode %.0f m/s : c/a = %.3f" % (V_DEEPMODE, ca(V_DEEPMODE)))
    if np.isfinite(cs).any():
        kp = int(np.nanargmax(cs))
        say("  scan peak %.4f at %.0f m/s (p = %.4f)" % (cs[kp], nls.V_GRID[kp], pv[kp]))
    say("  velocities clearing their own null: %d of %d (chance ~%.1f)"
        % (clears.size, np.isfinite(cs).sum(), 0.05 * np.isfinite(cs).sum()))
    if clears.size:
        say("    %.0f-%.0f m/s" % (nls.V_GRID[clears].min(), nls.V_GRID[clears].max()))

    fig, ax = plt.subplots(1, 2, figsize=(12.6, 5.2), constrained_layout=True)
    fig.suptitle("Nano fibre, %.2f h, source ch %d   |   %s to %s  (local)   |   %s"
                 % (hours, a.source,
                    use[0][1].astimezone(LA).strftime("%a %d %b %H:%M"),
                    use[-1][1].astimezone(LA).strftime("%a %d %b %H:%M"), a.tag),
                 fontsize=10.5)
    fig.text(0.995, 0.002, "UTC: %s to %s"
             % (use[0][1].strftime("%Y-%m-%d %H:%M"), use[-1][1].strftime("%Y-%m-%d %H:%M")),
             ha="right", va="bottom", fontsize=6.5, color="#9a9a9a")
    tn = gather / np.maximum(np.abs(gather).max(axis=1, keepdims=True), 1e-30)
    lim = float(np.percentile(np.abs(tn), 97.0))
    show = min(float(np.abs(lags).max()), 1.3 * offs.max() / V_DEEPMODE)
    ax[0].imshow(tn, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                 interpolation="nearest", extent=[lags[0], lags[-1], offs[-1], offs[0]])
    for v, c, ls in ((V_BODY, C1, "--"), (V_DEEPMODE, C2, "-")):
        ax[0].plot(offs / v, offs, ls, color=c, lw=1.5, label="%.0f m/s" % v)
    ax[0].set(xlim=(-show, show), xlabel="correlation lag (s)",
              ylabel="offset from virtual source (m)", title="(a) Gather")
    ax[0].legend(fontsize=8); ax[0].grid(False)
    ax[1].plot(nls.V_GRID, cs, "-", color=C1, lw=1.4, label="observed")
    ax[1].plot(nls.V_GRID, thr, "-", color=INK, lw=1.0, alpha=.7, label="null 95th")
    ax[1].axvline(V_BODY, color=C2, ls="--", lw=1.2, label="%.0f m/s" % V_BODY)
    ax[1].set(xlabel="trial velocity (m/s)", ylabel="moveout score",
              title="(b) Scan vs its own per-velocity null")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    fig.savefig(str(stem) + ".png", dpi=200, bbox_inches="tight")
    np.savez(str(stem) + ".npz", gather=gather, lags=lags, offsets=offs,
             v_grid=nls.V_GRID, score=cs, thresh=thr, p=pv, hours=hours,
             n_windows=nw, source_channel=a.source,
             t_first=str(use[0][1]), t_last=str(use[-1][1]))
    Path(str(stem) + ".txt").write_text("\n".join(log) + "\n")
    say("wrote %s.{npz,png,txt}" % stem.name)


if __name__ == "__main__":
    main()
