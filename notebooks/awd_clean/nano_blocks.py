#!/usr/bin/env python3
"""Nano, block by block, on the days people were on site. Is it washed out?

THE PUZZLE. On 2026-06-15/16 a crew was at SAFOD drilling, digging and deploying
nodes. Over essentially the same window (Nano 15:48-16:43 local, Deep
16:09-14:09), the wireline Deep fibre recovers a coherent arrival in 12 of 12
two-hour blocks, causal/acausal 1.1-4.0. The CEMENTED Nano fibre -- the one that
ought to be the better sensor for a downgoing body wave, because it is tied to
the formation -- gives a flat null over the same 25 h: p = 0.10 at 3200 m/s,
detectability N^+0.042, 5 of 229 velocities clearing against 11.5 by chance.

That is surprising and worth attacking before it is believed.

HYPOTHESIS TESTED HERE. The 25 h Nano number is a single grand stack. The Deep
time series varies by a factor of 3.5 between blocks and is strongly diurnal, so
if the Nano signal is concentrated in a few hours of peak surface activity, a
flat 25 h average would bury it. Averaging destroys a transient.

So: the same 2 h blocking the Deep got, and per block

    c/a at 3200 m/s   the a-priori body-wave velocity (in-band value)
    c/a at 1400 m/s   the velocity the Deep recovers, in case the Nano carries
                      the same mode more weakly rather than a body wave
    scan peak         where the moveout score actually maximises, unconstrained
    band RMS          a plain activity proxy -- DOES THIS FIBRE EVEN SEE THE
                      CREW? If the noise level never rises while people are
                      working, the null is about coupling or depth, not about
                      interferometry, and that is a different problem.

The RMS column is the one that decides which explanation is live, so it is
computed even though it is not an interferometry statistic.

Output: nano_blocks.{npz,png,txt}
"""
from __future__ import annotations

import argparse
import sys
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
import arrival_velocities as av
import nano_ambient_cc as nano
import nano_long_stack as nls

STEM = HERE / "nano_blocks"
LA = ZoneInfo("America/Los_Angeles")
V_BODY = av.V_NANO_INBAND        # 3200 m/s
V_DEEPMODE = 1400.0              # what the Deep recovers, tested here too
INK = "#444444"
C1, C2, C3 = "#0072B2", "#D55E00", "#009E73"


def ca_at(gather, lags, offs, v):
    """Causal and acausal envelope amplitude along t = +-x/v."""
    env = np.abs(hilbert(gather, axis=1))
    env = env / np.maximum(np.median(env, axis=1, keepdims=True), 1e-30)
    half = max(1, int(nls.GATE_S / (lags[1] - lags[0])))
    c, a = [], []
    for row, x in zip(env, offs):
        for sign, out in ((+1.0, c), (-1.0, a)):
            centre = sign * x / v
            if centre - nls.GATE_S < lags[0] or centre + nls.GATE_S > lags[-1]:
                continue
            k = int(np.argmin(np.abs(lags - centre)))
            out.append(row[k - half:k + half + 1].mean())
    if not c or not a:
        return np.nan, np.nan
    return float(np.median(c)), float(np.median(a))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block-hours", type=float, default=2.0)
    ap.add_argument("--source", type=int, default=nls.NANO_WELLHEAD)
    a = ap.parse_args()
    if a.source < nls.NANO_WELLHEAD:
        raise SystemExit("channel %d is in the AIR; the hole starts at %d"
                         % (a.source, nls.NANO_WELLHEAD))

    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    files = sorted(nano.NANO_DIR.glob("*.pb"))
    pre = [(p, t) for p, t in ((p, nano.nano_time(p.name)) for p in files)
           if t is not None and t < nano.FIRST_DROP]
    per = max(1, int(round(a.block_hours * 3600 / 300)))     # 300 s records
    blocks = [pre[i:i + per] for i in range(0, len(pre), per)]
    blocks = [b for b in blocks if len(b) == per]

    say("Nano, block by block, source channel %d (the wellhead)" % a.source)
    say("  %d blocks of %.1f h, %s to %s LOCAL"
        % (len(blocks), a.block_hours,
           pre[0][1].astimezone(LA).strftime("%Y-%m-%d %H:%M"),
           pre[-1][1].astimezone(LA).strftime("%Y-%m-%d %H:%M")))
    say("  crew on site 2026-06-15 (nodal deployment: drilling, digging)")
    say("")
    say("  %-18s %8s %8s %9s %10s" % ("block start (local)", "c/a 3200",
                                      "c/a 1400", "peak m/s", "band RMS"))

    rows = []
    for blk in blocks:
        paths = [str(p) for p, _ in blk]
        arr, info = nano.readFile_protobuf(paths, fmin=1.0, fmax=100.0,
                                           desampling=True, **nano.RAW_KW)
        x = np.asarray(arr, dtype=np.float32); del arr
        fs = float(info["fs"]); dx = float(info.get("dx", nano.DX_NANO))
        x[:, 1:] = np.diff(x, axis=1) * np.float32(fs); x[:, 0] = 0.0

        centres = a.source + np.rint(nls.OFFSETS_M / dx).astype(int)
        keep = centres + nls.NEIGHBOUR < x.shape[0]
        centres, offs = centres[keep], nls.OFFSETS_M[keep]

        sos = butter(4, list(nls.BAND), btype="bandpass", fs=fs, output="sos")
        rms = float(np.sqrt(np.mean(sosfiltfilt(sos, x[a.source]) ** 2)))

        n_ram = ex.odd_ram_samples(nls.RAM_S, fs)
        for i in range(0, x.shape[0], 64):
            b = x[i:i + 64]
            w = uniform_filter1d(np.abs(b), size=n_ram, axis=1, mode="nearest")
            np.divide(b, np.maximum(w, np.float32(np.finfo(np.float32).eps)), out=b)
            del w

        n_win, n_step = int(nls.WINDOW_S * fs), int(nls.STEP_S * fs)
        n_fft = 1 << int(np.ceil(np.log2(2 * n_win - 1)))
        acc = np.zeros((len(centres), n_fft // 2 + 1), dtype=np.complex128)
        nw = 0
        for s in range(0, x.shape[1] - n_win + 1, n_step):
            S = np.conj(np.fft.rfft(x[a.source, s:s + n_win], n=n_fft))
            for j, c in enumerate(centres):
                seg = x[c - nls.NEIGHBOUR:c + nls.NEIGHBOUR + 1, s:s + n_win].sum(axis=0)
                acc[j] += S * np.fft.rfft(seg, n=n_fft)
            nw += 1
        del x

        gather, lags = nls.spectra_to_gather(acc, nw, n_fft, fs)
        cs = nls.moveout_curve(gather, lags, offs)
        pk = float(nls.V_GRID[np.nanargmax(cs)]) if np.isfinite(cs).any() else np.nan
        c3, a3 = ca_at(gather, lags, offs, V_BODY)
        c1, a1 = ca_at(gather, lags, offs, V_DEEPMODE)
        r3 = c3 / a3 if a3 else np.nan
        r1 = c1 / a1 if a1 else np.nan
        t_loc = blk[0][1].astimezone(LA)
        say("  %-18s %8.3f %8.3f %9.0f %10.3e"
            % (t_loc.strftime("%a %d %b %H:%M"), r3, r1, pk, rms))
        rows.append((blk[0][1], r3, r1, pk, rms, nw))

    t = np.array([str(r[0]) for r in rows])
    r3 = np.array([r[1] for r in rows]); r1 = np.array([r[2] for r in rows])
    pk = np.array([r[3] for r in rows]); rms = np.array([r[4] for r in rows])

    say("")
    say("=== reading ===")
    say("  c/a at 3200 m/s : median %.3f, max %.3f, %d of %d blocks above 1.5"
        % (np.nanmedian(r3), np.nanmax(r3), int(np.nansum(r3 > 1.5)), r3.size))
    say("  c/a at 1400 m/s : median %.3f, max %.3f, %d of %d blocks above 1.5"
        % (np.nanmedian(r1), np.nanmax(r1), int(np.nansum(r1 > 1.5)), r1.size))
    say("  band RMS varies %.2fx across the window (max/min)"
        % (np.nanmax(rms) / max(np.nanmin(rms), 1e-30)))
    say("  Deep, same days, same blocking: median 2.694, 11 of 12 above 1.5")
    say("")
    if np.nansum(r3 > 1.5) == 0 and np.nansum(r1 > 1.5) == 0:
        say("  NOT a washout. No individual block shows a one-sided arrival at")
        say("  either velocity, so the flat 25 h stack is not hiding a transient.")
        if np.nanmax(rms) / max(np.nanmin(rms), 1e-30) > 3.0:
            say("  The fibre DOES see the surface activity (RMS varies strongly),")
            say("  so it is receiving the noise and not converting it into a")
            say("  recoverable arrival -- a wavefield or coupling question, not a")
            say("  question about the stack.")
        else:
            say("  The band RMS barely moves, so this fibre may not be registering")
            say("  the surface work at all. Check coupling and depth before")
            say("  concluding anything about interferometry.")
    else:
        say("  SOME blocks are one-sided. The 25 h grand stack averaged over a")
        say("  time-varying signal; report the blocks, not the stack.")

    fig, ax = plt.subplots(3, 1, figsize=(11.0, 8.4), sharex=True,
                           constrained_layout=True)
    import pandas as pd
    tl = pd.to_datetime(t, utc=True).tz_convert(LA)
    fig.suptitle("Nano fibre, %.1f h blocks, source ch %d   |   %s to %s  (local)"
                 "   |   crew on site 15 Jun"
                 % (a.block_hours, a.source, tl[0].strftime("%a %d %b %H:%M"),
                    tl[-1].strftime("%a %d %b %H:%M")), fontsize=10.5)
    ax[0].plot(tl, r3, "o-", color=C1, label="c/a at %.0f m/s" % V_BODY)
    ax[0].plot(tl, r1, "s-", color=C2, label="c/a at %.0f m/s (Deep mode)" % V_DEEPMODE)
    ax[0].axhline(1.0, color=INK, ls="--", lw=1.1, label="balanced")
    ax[0].axhline(2.694, color=C3, ls=":", lw=1.4, label="Deep median, same days")
    ax[0].set(ylabel="causal / acausal"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    ax[1].plot(tl, pk, "o-", color=INK)
    ax[1].axhline(V_BODY, color=C1, ls="--", lw=1.1)
    ax[1].set(ylabel="scan peak (m/s)"); ax[1].grid(alpha=.3)
    ax[2].semilogy(tl, rms, "o-", color=C3)
    ax[2].set(ylabel="band RMS (activity)", xlabel="local time"); ax[2].grid(alpha=.3)
    fig.autofmt_xdate()
    fig.savefig(str(STEM) + ".png", dpi=200, bbox_inches="tight")
    np.savez(str(STEM) + ".npz", t=t, ca_3200=r3, ca_1400=r1, peak=pk, rms=rms,
             source_channel=a.source, block_hours=a.block_hours)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
