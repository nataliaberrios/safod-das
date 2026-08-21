#!/usr/bin/env python3
"""The whole Deep fibre, one gather: wellhead source, every channel to the bottom.

Everything so far used 700 m of offset, which is Lellouch's Figure 7c aperture
inherited from the reproduction. The Deep outbound limb runs channel 211 to 1700
-- 3051 m along the hole -- so that was 23% of the fibre. This uses all of it.

Distance is measured ALONG THE HOLE (MD from the geometry spreadsheet), because
that is the distance the correlation lag actually corresponds to for anything
travelling down the borehole, and because it is what the channel axis is. Below
channel 949 the hole deviates to ~55 degrees, so MD and true vertical depth
diverge by 503 m (20%) by the bottom -- 0.36 s at 1400 m/s, which is 30x the
0.012 s picking gate. TVD is carried alongside so the difference is visible
rather than assumed.

No velocity is imposed and nothing is scanned. Read the moveout off the picture.

Output: deep_full_limb.{npz,png,txt}
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, hilbert, sosfiltfilt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ambient_lellouch2019_exact_stack as ex
import arrival_velocities as av
import deep_timeseries as dts
import safod_geometry as geo

STEM = HERE / "deep_full_limb"
SOURCE_CH = 211
BAND = (5.0, 20.0)
INK = "#444444"
C1, C2, C3 = "#0072B2", "#D55E00", "#009E73"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfiles", type=int, default=300, help="60 s records (300 = 5 h)")
    ap.add_argument("--step", type=int, default=4, help="use every Nth channel")
    ap.add_argument("--max-lag", type=float, default=3.5)
    a = ap.parse_args()

    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    g = geo.load()
    rev = int(g["reversal_channel"])
    rx = np.arange(SOURCE_CH + a.step, rev + 1, a.step)
    j = np.searchsorted(g["channel"], rx)
    md = g["md_m"][j] - g["md_m"][int(np.searchsorted(g["channel"], SOURCE_CH))]
    tvd = g["tvd_m"][j]
    inc = g["inclination_deg"][j]

    say("Deep fibre, whole outbound limb from the wellhead")
    say("  source ch %d (MD 2 m); receivers ch %d-%d, every %d"
        % (SOURCE_CH, rx[0], rx[-1], a.step))
    say("  %d traces, %.0f m along the hole, %.0f m vertically"
        % (rx.size, md.max(), tvd.max() - tvd.min()))
    say("  inclination %.1f deg at the top to %.1f deg at the bottom"
        % (inc.min(), inc.max()))
    say("")

    rows = ex.deep_rows("deepA").iloc[: a.nfiles]
    # the manifest column is `file`, and some rows predate a data move, so every
    # path goes through corrected_path -- same as deep_timeseries does
    paths = [ex.corrected_path(r) for r in rows.file]
    hours = len(paths) * 60 / 3600.0
    say("  %.2f h stacked (%d x 60 s records)" % (hours, len(paths)))

    dts.SOURCE_CH = SOURCE_CH
    acc, nw, fs, n_fft = dts.block_spectra(paths, rx, BAND)
    say("  %d windows at %.0f Hz" % (nw, fs))

    full = np.fft.fftshift(np.fft.irfft(acc / max(1, nw), n=n_fft, axis=-1), axes=-1)
    lags = (np.arange(n_fft) - n_fft // 2) / fs
    keep = np.abs(lags) <= a.max_lag
    sos = butter(4, list(BAND), btype="bandpass", fs=fs, output="sos")
    gather = sosfiltfilt(sos, full, axis=-1)[:, keep]
    lags = lags[keep]
    say("  lag window +-%.2f s (needs %.2f s for %.0f m at 1400 m/s)"
        % (a.max_lag, md.max() / 1400.0, md.max()))
    say("")

    # pick the causal envelope peak per trace, away from the zero-lag lobe
    env = np.abs(hilbert(gather, axis=1))
    m = lags > 0.05
    pk = lags[m][env[:, m].argmax(axis=1)]

    # fit over the near-vertical part only, then see whether it EXTRAPOLATES
    # into the deviated part -- that is the whole question, and it needs no scan
    nv = inc < 5.0
    fit = nv & (md > 100)
    A = np.polyfit(md[fit], pk[fit], 1)
    v_md = 1.0 / A[0]
    r = float(np.corrcoef(md[fit], pk[fit])[0, 1])
    say("=== fit on the NEAR-VERTICAL part only (%d traces, %.0f-%.0f m) ==="
        % (fit.sum(), md[fit].min(), md[fit].max()))
    say("  velocity along the hole : %.0f m/s" % v_md)
    say("  intercept               : %+.3f s" % A[1])
    say("  r                       : %+.4f" % r)
    say("")

    dev = inc > 20.0
    if dev.any():
        pred_md = np.polyval(A, md[dev])
        # if the wave went straight down through rock instead, lag would track TVD
        tvd_from_src = tvd[dev] - tvd[0]
        pred_tvd = A[1] + tvd_from_src / v_md
        r_md = float(np.sqrt(np.mean((pk[dev] - pred_md) ** 2)))
        r_tvd = float(np.sqrt(np.mean((pk[dev] - pred_tvd) ** 2)))
        say("=== extrapolated into the DEVIATED part (%d traces, inc %.0f-%.0f deg) ==="
            % (dev.sum(), inc[dev].min(), inc[dev].max()))
        say("  the two predictions separate by up to %.2f s there"
            % float(np.max(np.abs(pred_md - pred_tvd))))
        say("  residual if it follows the HOLE      : %.3f s" % r_md)
        say("  residual if it goes straight DOWN    : %.3f s" % r_tvd)
        if min(r_md, r_tvd) > 0.15:
            say("  NEITHER fits. The picks in the deviated section are probably not")
            say("  the same arrival -- do not read a path from this.")
        elif r_md < r_tvd:
            say("  -> tracks the BOREHOLE, i.e. it is guided along the hole.")
        else:
            say("  -> tracks VERTICAL DEPTH, i.e. it is not following the hole.")
        say("")

    fig, ax = plt.subplots(1, 2, figsize=(13.5, 6.4), constrained_layout=True)
    t0, t1 = rows.time.iloc[0], rows.time.iloc[-1]
    title, foot = geo.figure_label(t0, t1, hours, fibre="Deep",
                                   extra="whole outbound limb, source ch %d" % SOURCE_CH)
    fig.suptitle(title, fontsize=10.5)
    fig.text(0.995, 0.002, foot, ha="right", va="bottom", fontsize=6.5, color="#9a9a9a")

    tn = gather / np.maximum(np.abs(gather).max(axis=1, keepdims=True), 1e-30)
    lim = float(np.percentile(np.abs(tn), 97.0))
    ax[0].imshow(tn, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                 interpolation="nearest",
                 extent=[lags[0], lags[-1], md[-1], md[0]])
    ax[0].plot(np.polyval(A, md), md, "-", color=C3, lw=1.6,
               label="%.0f m/s fitted above 5 deg" % v_md)
    ax[0].axhline(md[np.argmax(inc > 5.0)], color=INK, ls=":", lw=1.2,
                  label="hole starts to deviate")
    ax[0].set(xlim=(-a.max_lag, a.max_lag), xlabel="correlation lag (s)",
              ylabel="distance along the hole from the wellhead (m)",
              title="(a) Whole limb, %d traces" % rx.size)
    ax[0].legend(fontsize=8, loc="lower left"); ax[0].grid(False)

    ax[1].plot(pk, md, "o", ms=2.5, color=C1, label="picked envelope peak")
    ax[1].plot(np.polyval(A, md), md, "-", color=C3, lw=1.6,
               label="follows the hole (%.0f m/s)" % v_md)
    ax[1].plot(A[1] + (tvd - tvd[0]) / v_md, md, "--", color=C2, lw=1.6,
               label="straight down at the same speed")
    ax[1].axhline(md[np.argmax(inc > 5.0)], color=INK, ls=":", lw=1.2)
    ax[1].set(xlabel="lag (s)", ylabel="distance along the hole (m)",
              title="(b) Picks vs the two paths")
    ax[1].invert_yaxis(); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)

    fig.savefig(str(STEM) + ".png", dpi=300, bbox_inches="tight")
    np.savez_compressed(str(STEM) + ".npz", gather=gather, lags=lags, md=md,
                        tvd=tvd, inclination=inc, channels=rx, picks=pk,
                        hours=hours, n_windows=nw, fs=fs, fit=A,
                        source_channel=SOURCE_CH, t_first=str(t0), t_last=str(t1))
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
