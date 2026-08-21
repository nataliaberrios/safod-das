#!/usr/bin/env python3
"""Dense virtual-source record section on the Deep fibre -- see the moveout.

WHY THIS EXISTS. Every Deep-fibre figure so far has been a QC diagnostic: score
against trial velocity, p-value against velocity. Those answer "is it
significant", not "what does it look like", and they are the wrong figure to
judge a seismic arrival by.

The figure a seismologist wants is the GATHER: correlation lag on x, distance
along the fibre on y, with the wavefront visible as a coherent alignment. That is
what Lellouch et al. (2019) Figure 7c is.

Two things were wrong with the earlier plots:
  - they plotted the SCAN, not the section;
  - the published geometry uses only 14 receivers at 50 m spacing, which is far
    too sparse to see a wavefront even when one is present.

This correlates the virtual source against EVERY channel over the offset range,
giving a dense section (hundreds of traces) in which moveout is either visible or
it is not. No velocity assumption enters the picture -- the moveout line is drawn
on top afterwards, so the reader can judge the alignment independently.

Output: deep_record_section.{npz,png}
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ambient_lellouch2019_exact_stack as ex
import deep_cc_steps as steps

STEM = HERE / "deep_record_section"   # --tag appends a suffix


def window_label(t_first, t_last):
    """UTC and LOCAL span for a figure title.

    Local time is included deliberately: whether personnel were at SAFOD is a
    local-clock question, and these figures are cross-checked against site
    activity (e.g. the nodal deployment of 2026-06-15, the AWD weight drops of
    2026-06-16/17). A UTC-only label makes that check error-prone, since UTC-7
    pushes a California afternoon onto the following UTC date.
    """
    import pandas as pd
    a = pd.Timestamp(t_first); b = pd.Timestamp(t_last)
    if a.tzinfo is None:
        a = a.tz_localize("UTC")
    if b.tzinfo is None:
        b = b.tz_localize("UTC")
    la = a.tz_convert("America/Los_Angeles"); lb = b.tz_convert("America/Los_Angeles")
    same_day = la.date() == lb.date()
    loc = ("%s %s-%s local" % (la.strftime("%Y-%m-%d"), la.strftime("%H:%M"),
                               lb.strftime("%H:%M"))
           if same_day else
           "%s to %s local" % (la.strftime("%Y-%m-%d %H:%M"),
                               lb.strftime("%Y-%m-%d %H:%M")))
    utc = "%s to %s UTC" % (a.strftime("%Y-%m-%d %H:%M"), b.strftime("%Y-%m-%d %H:%M"))
    return loc, utc
SRC = 400
MAX_OFFSET_M = 750.0
WINDOW_S, STEP_S = 30.0, 15.0
V_MARK = 1675.0          # the recovered arrival
V_LELL = 3200.0          # what Lellouch reports on the cemented main hole


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfiles", type=int, default=120)
    ap.add_argument("--arm", default="deepA")
    ap.add_argument("--source", type=int, default=SRC)
    ap.add_argument("--common-mode", action="store_true",
                    help="subtract the across-channel MEDIAN at each time sample. "
                         "This is the simple fix for the zero-lag stripe. It is NOT "
                         "part of Lellouch et al. (2019) section 4.1 and is a "
                         "labelled sensitivity branch in this project.")
    ap.add_argument("--neighbour", type=int, default=0,
                    help="half-width of the published R+-10 receiver-neighbour sum, "
                         "in channels. 0 (default) uses single channels, which is a "
                         "DEVIATION from Lellouch made to keep the section dense; "
                         "10 reproduces the published sum.")
    ap.add_argument("--tag", default="")
    ap.add_argument("--fibre", choices=("deep", "nano"), default="deep",
                    help="nano reads the CEMENTED protobuf fibre instead")
    ap.add_argument("--v-mark", type=float, default=None,
                    help="velocity to draw; defaults to 1675 (deep) / 2950 (nano)")
    a = ap.parse_args()

    global V_MARK
    if a.v_mark:
        V_MARK = a.v_mark
    elif a.fibre == "nano":
        V_MARK = 2950.0        # the mode the AWD active source measures on Nano

    # HARD GUARD, on both fibres. A virtual source outside the borehole
    # correlates air against rock: it returns instrumental common mode at zero
    # lag and cannot show moveout, so it produces a null that says nothing about
    # the wavefield. Deep channel 98 did exactly that, and every Nano product
    # written before 2026-08-20 used channel 10 -- 63 channels above the entry
    # measured by nano_find_wellhead.py. deep_section_depth.py already refuses;
    # this script did not, and deep_record_section_nano{,cm}.npz are the result.
    if a.fibre == "nano":
        import nano_ambient_cc as nano
        nano.require_downhole(a.source)
    else:
        import safod_geometry as geo
        _g = geo.load()
        _j = int(np.clip(np.searchsorted(_g["channel"], a.source), 0,
                         _g["channel"].size - 1))
        if not _g["in_hole"][_j]:
            raise SystemExit(
                "REFUSING to run: Deep channel %d is SURFACE LEAD-IN (not in "
                "the hole).\n  The first downhole channel is %d (TVD %.1f m)."
                % (a.source, int(_g["channel"][_g["in_hole"]][0]),
                   float(_g["tvd_m"][_g["in_hole"]][0])))

    if a.fibre == "nano":
        files = sorted(nano.NANO_DIR.glob("*.pb"))
        pre = [(p_, t) for p_, t in ((p_, nano.nano_time(p_.name)) for p_ in files)
               if t is not None and t < nano.FIRST_DROP]
        use = pre[: a.nfiles]
        print("Nano (CEMENTED), %d pre-drop records = %.2f h"
              % (len(use), len(use) * 300 / 3600), flush=True)
        print("  %s -> %s UTC" % (use[0][1], use[-1][1]), flush=True)
        batches, fs, dx = [], None, None
        for i in range(0, len(use), 6):
            arr, info = nano.readFile_protobuf([str(q) for q, _ in use[i:i + 6]],
                                               fmin=1.0, fmax=100.0,
                                               desampling=True, **nano.RAW_KW)
            if fs is None:
                fs = float(info["fs"]); dx = float(info.get("dx", nano.DX_NANO))
            batches.append(np.asarray(arr, dtype=np.float32)); del arr
        allch = np.concatenate(batches, axis=1); del batches
        n_out = int(MAX_OFFSET_M / dx)
        hi = min(a.source + n_out + 1, allch.shape[0])
        raw = allch[a.source:hi]; del allch
        channels = list(range(a.source, hi))
        t_first, t_last = str(use[0][1]), str(use[-1][1])
        print("  %d receiver channels, %.3f m spacing, out to %.0f m, fs=%.0f Hz"
              % (len(channels), dx, (hi - a.source) * dx, fs), flush=True)
    else:
        rows = ex.deep_rows(a.arm).iloc[: a.nfiles]
        paths = [ex.corrected_path(r) for r in rows.file]
        print("arm %s, %d records = %.2f h" % (a.arm, len(paths), len(paths)*60/3600),
              flush=True)
        print("  %s -> %s UTC" % (rows.time.iloc[0], rows.time.iloc[-1]), flush=True)
        with __import__("h5py").File(paths[0], "r") as h:
            dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 2.0419))
        n_out = int(MAX_OFFSET_M / dx)
        channels = list(range(a.source, a.source + n_out + 1))
        print("  %d receiver channels, %.1f m spacing, out to %.0f m"
              % (len(channels), dx, n_out * dx), flush=True)
        raw, fs, _ = steps.read_records(paths, channels)
        t_first, t_last = str(rows.time.iloc[0]), str(rows.time.iloc[-1])

    # MEMORY. A dense section is 369 channels x nfiles x 60 s x 1000 Hz. At
    # float64 that is ~21 GB per array for 120 records, and holding raw + rate +
    # normalised at once OOM-killed a 64 GB job. Work in float32 and overwrite in
    # place so only ONE full-size array is ever live.
    x = np.asarray(raw, dtype=np.float32); del raw
    hours = x.shape[1] / fs / 3600.0
    print("  %.2f h stacked  (array %s float32 = %.1f GB)"
          % (hours, x.shape, x.nbytes / 2**30), flush=True)
    # differentiate in place
    x[:, 1:] = np.diff(x, axis=1) * np.float32(fs)
    x[:, 0] = 0.0
    # running-absolute-mean, in place, one channel block at a time
    from scipy.ndimage import uniform_filter1d
    n_ram = ex.odd_ram_samples(0.1, fs)
    BLK = 64
    for i in range(0, x.shape[0], BLK):
        blk = x[i:i + BLK]
        w = uniform_filter1d(np.abs(blk), size=n_ram, axis=1, mode="nearest")
        scale = float(np.nanmedian(w)) or 1.0
        np.divide(blk, np.maximum(w, np.float32(np.finfo(np.float32).eps * scale)),
                  out=blk)
        del w
    if a.common_mode:
        # The simple fix: remove what is identical on every channel at each
        # instant. One line, and it is what kills the zero-lag stripe.
        med = np.median(x, axis=0)
        x -= med.astype(np.float32)[None, :]
        print("  common mode removed (across-channel median per sample)", flush=True)
    normed = x

    n_win, n_step = int(WINDOW_S * fs), int(STEP_S * fs)
    nb = int(a.neighbour)
    if nb > 0:
        idx = [i for i in range(1 + nb, len(channels) - nb)]
        rx = [np.arange(i - nb, i + nb + 1) for i in idx]
        offsets = np.array(idx, dtype=float) * dx
        print("  R+-%d neighbour sum (published), %d receivers" % (nb, len(rx)), flush=True)
    else:
        idx = list(range(1, len(channels)))
        rx = [np.array([i]) for i in idx]
        offsets = np.array(idx, dtype=float) * dx
        print("  single channels (dense, deviates from the published R+-10)", flush=True)
    gather, lags, nw = steps.correlate(normed, 0, rx, fs, n_win, n_step)
    print("  %.2f h in %d windows, gather %s" % (hours, nw, gather.shape), flush=True)

    stem = Path(str(STEM) + (("_" + a.tag) if a.tag else ""))
    np.savez_compressed(str(stem) + ".npz", gather=gather, lags=lags,
                        offsets=offsets, fs=fs, dx=dx, n_windows=nw,
                        source_channel=a.source, arm=a.arm, hours=hours,
                        t_first=t_first, t_last=t_last, fibre=a.fibre)

    # ---------------- the figure ----------------
    INK, MUTED = "#444444", "#6b6b6b"
    C_MARK, C_ALT = "#D55E00", "#0072B2"
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "font.size": 9.5,
        "axes.titlesize": 10, "axes.labelsize": 9.5, "axes.edgecolor": "#b0b0b0",
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "axes.labelcolor": INK, "legend.frameon": False, "legend.fontsize": 8,
    })

    # trace-normalise so distant traces are visible at all
    tn = gather / np.maximum(np.abs(gather).max(axis=1, keepdims=True), 1e-30)
    env = np.abs(hilbert(gather, axis=1))
    env = env / np.maximum(env.max(axis=1, keepdims=True), 1e-30)

    loc_lab, utc_lab = window_label(t_first, t_last)
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 6.0), constrained_layout=True,
                           sharey=True)
    fig.suptitle("%s fibre, %.1f h stacked   |   %s   (%s)"
                 % (a.fibre.capitalize(), hours, loc_lab, utc_lab),
                 fontsize=10.5, y=1.005)

    # (a) amplitude, diverging map with neutral zero, clipped so the arrival shows
    lim = float(np.percentile(np.abs(tn), 97.0))
    ax[0].imshow(tn, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                 interpolation="nearest",
                 extent=[lags[0], lags[-1], offsets[-1], offsets[0]])
    ax[0].set(xlim=(-0.35, 0.35), xlabel="correlation lag (s)",
              ylabel="distance from virtual source along fibre (m)",
              title="(a) Record section, %.1f h stacked, %d traces"
                    % (hours, len(offsets)))

    # (b) same, with the moveout lines drawn ON TOP so alignment is judged by eye
    ax[1].imshow(tn, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                 interpolation="nearest",
                 extent=[lags[0], lags[-1], offsets[-1], offsets[0]])
    for v, c, ls, lab in ((V_MARK, C_MARK, "-", "%.0f m/s (recovered)" % V_MARK),
                          (V_LELL, C_ALT, "--", "%.0f m/s (Lellouch)" % V_LELL)):
        ax[1].plot(offsets / v, offsets, ls, color=c, lw=1.5, alpha=0.9, label=lab)
        ax[1].plot(-offsets / v, offsets, ls, color=c, lw=1.0, alpha=0.5)
    ax[1].set(xlim=(-0.35, 0.35), xlabel="correlation lag (s)",
              title="(b) With moveout curves overlaid")
    ax[1].legend(loc="lower right")

    # (c) envelope: coherent energy is easier to see than polarity
    ax[2].imshow(env, aspect="auto", cmap="magma", vmin=0,
                 vmax=float(np.percentile(env, 99.0)), interpolation="nearest",
                 extent=[lags[0], lags[-1], offsets[-1], offsets[0]])
    ax[2].plot(offsets / V_MARK, offsets, "-", color="w", lw=1.3, alpha=0.85,
               label="%.0f m/s" % V_MARK)
    ax[2].set(xlim=(-0.35, 0.35), xlabel="correlation lag (s)",
              title="(c) Envelope")
    ax[2].legend(loc="lower right")
    fig.savefig(str(stem) + ".png", bbox_inches="tight")
    print("wrote %s.{npz,png}" % stem.name, flush=True)

    # a wiggle version, decimated in trace count so individual traces are legible
    step = max(1, len(offsets) // 60)
    fig2, ax2 = plt.subplots(figsize=(7.6, 8.4), constrained_layout=True)
    scale = 2.2 * step * dx
    for o, row in zip(offsets[::step], tn[::step]):
        y = o - row * scale
        ax2.plot(lags, y, "-", color=INK, lw=0.6)
        ax2.fill_between(lags, o, y, where=(y < o), color=INK, alpha=0.55, lw=0)
    ax2.plot(offsets / V_MARK, offsets, "-", color=C_MARK, lw=1.6,
             label="%.0f m/s" % V_MARK)
    ax2.plot(offsets / V_LELL, offsets, "--", color=C_ALT, lw=1.4,
             label="%.0f m/s (Lellouch)" % V_LELL)
    ax2.set(xlim=(-0.25, 0.25), ylim=(offsets[-1], 0),
            xlabel="correlation lag (s)",
            ylabel="distance from virtual source along fibre (m)",
            title="%s fibre virtual-source gather, %.1f h stacked\n%s\n%s"
                  % (a.fibre.capitalize(), hours, loc_lab, utc_lab))
    ax2.legend(loc="lower right")
    fig2.savefig(str(stem) + "_wiggle.png", bbox_inches="tight")
    print("wrote %s_wiggle.png" % stem.name, flush=True)


if __name__ == "__main__":
    main()
