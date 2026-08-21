#!/usr/bin/env python3
"""Deep-fibre virtual-source gather plotted against TRUE VERTICAL DEPTH.

Supersedes `deep_record_section.py`, which plotted distance ALONG FIBRE because
no depth registration existed. `SAFOD_Phase2_GeoReferenced_Channels.xlsx` now
supplies it (read through `safod_geometry.py`), and that changes what can be
claimed as well as what is drawn:

  - the array is restricted to the NEAR-VERTICAL section, channels 211-949,
    inclination < 5 deg, TVD 2-1513 m. Outside it the hole reaches 55 deg and an
    along-fibre apparent velocity is not a vertical velocity;
  - channels 0-210 are SURFACE LEAD-IN and are excluded. Source channel 98 in the
    earlier scan was surface fibre, which is why it alone returned nothing;
  - the vertical axis is TVD in metres, so the moveout line is a real velocity.

RETURN-LIMB CONTROL. The fibre reverses at channel 1700 and retraces the same
depths, so every outbound channel has a return-limb partner at the same TVD on a
physically different piece of fibre with the same interrogator. `--limb return`
runs the identical measurement there. An arrival appearing at the same DEPTH on
both limbs cannot be a processing artefact; one appearing on only one limb is a
warning.

Figure labels follow the PI's convention via `safod_geometry.figure_label`:
duration stacked in HOURS (not file counts, which do not compare across fibres
with different record lengths), and the date range in BOTH local time and UTC
(because whether personnel were on site is a local-clock question).

Output: deep_section_depth[_<tag>].{npz,png}
"""
from __future__ import annotations

import argparse
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
import safod_geometry as geo
import arrival_velocities as av

STEM = HERE / "deep_section_depth"
WINDOW_S, STEP_S = 30.0, 15.0
BAND = (5.0, 20.0)
# Derived from the aperture, not hand-picked. At 0.35 s a 700 m separation is
# unreachable below 1934 m/s, so the arrival simply was not in the window -- and
# the plot's xlim was hard-set to the same 0.35 s, so it was not on the figure
# either. Bound below by the slowest velocity any curve is drawn at.
MAX_LAG_S = None                # set per run in main(), from the actual span
V_MARK = av.V_DEEP_ARRIVAL       # retracted from 1675
V_LELL = 3200.0

INK, MUTED, NULLC = "#444444", "#6b6b6b", "#8a8a8a"
C1, C2, C3 = "#0072B2", "#D55E00", "#009E73"   # validated colourblind-safe


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfiles", type=int, default=96)
    ap.add_argument("--arm", default="deepA")
    ap.add_argument("--source", type=int, default=400)
    ap.add_argument("--limb", choices=("outbound", "return"), default="outbound")
    ap.add_argument("--max-depth-span", type=float, default=760.0,
                    help="metres of TVD below the source to include")
    ap.add_argument("--common-mode", action="store_true")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    global V_MARK, MAX_LAG_S
    v_src = av.deep_velocity_for_source(a.source)
    if v_src is not None:
        V_MARK = v_src          # this source's own measured velocity
    # the window must reach the deepest receiver at the SLOWEST curve drawn
    MAX_LAG_S = round(av.required_lag_s(a.max_depth_span, min(V_MARK, V_LELL)), 2)

    g = geo.load()
    src = a.source
    if a.limb == "return":
        src = geo.return_limb_partner(a.source)
        print("return limb: outbound ch %d -> ch %d (same TVD)" % (a.source, src), flush=True)
    print(geo.describe(src), flush=True)
    # HARD GUARD. Do not correlate from a channel that is in the air. An earlier
    # scan used Deep channel 98, which is surface lead-in, and it duly returned
    # nothing -- a null that said something about the cable, not the wavefield.
    _g = geo.load()
    _j = int(np.searchsorted(_g["channel"], src))
    if not _g["in_hole"][_j]:
        raise SystemExit(
            "REFUSING to run: channel %d is SURFACE LEAD-IN (not in the hole).\n"
            "  The Deep fibre's first downhole channel is %d (TVD %.1f m)."
            % (src, int(_g["channel"][_g["in_hole"]][0]),
               _g["tvd_m"][_g["in_hole"]][0]))

    # receivers: same limb, increasing depth, inside the near-vertical section
    ch_all = g["channel"].astype(int)
    j_src = int(np.searchsorted(ch_all, src))
    tvd_src = g["tvd_m"][j_src]
    if a.limb == "outbound":
        mask = g["near_vertical"] & (g["tvd_m"] > tvd_src) & \
               (g["tvd_m"] <= tvd_src + a.max_depth_span)
    else:
        # on the return limb depth DECREASES with channel index
        mask = g["return_limb"] & (g["inclination_deg"] < geo.NEAR_VERTICAL_MAX_INC) & \
               (g["tvd_m"] > tvd_src) & (g["tvd_m"] <= tvd_src + a.max_depth_span)
    rx_ch = ch_all[mask]
    rx_tvd = g["tvd_m"][mask]
    order = np.argsort(rx_tvd)
    rx_ch, rx_tvd = rx_ch[order], rx_tvd[order]
    if rx_ch.size < 20:
        raise SystemExit("only %d receivers in the near-vertical span; widen "
                         "--max-depth-span or pick a shallower source" % rx_ch.size)
    print("  %d receivers, TVD %.0f-%.0f m (source at %.0f m)"
          % (rx_ch.size, rx_tvd.min(), rx_tvd.max(), tvd_src), flush=True)

    rows = ex.deep_rows(a.arm).iloc[: a.nfiles]
    paths = [ex.corrected_path(r) for r in rows.file]
    hours = len(paths) * 60 / 3600.0
    print("  %.2f h stacked (%d x 60 s records)" % (hours, len(paths)), flush=True)

    channels = [src] + list(rx_ch)
    raw, fs, dx = steps.read_records(paths, channels)
    x = np.asarray(raw, dtype=np.float32); del raw
    x[:, 1:] = np.diff(x, axis=1) * np.float32(fs); x[:, 0] = 0.0
    n_ram = ex.odd_ram_samples(0.1, fs)
    for i in range(0, x.shape[0], 64):
        blk = x[i:i + 64]
        w = uniform_filter1d(np.abs(blk), size=n_ram, axis=1, mode="nearest")
        sc = float(np.nanmedian(w)) or 1.0
        np.divide(blk, np.maximum(w, np.float32(np.finfo(np.float32).eps * sc)), out=blk)
        del w
    if a.common_mode:
        x -= np.median(x, axis=0).astype(np.float32)[None, :]
        print("  common mode removed", flush=True)

    n_win, n_step = int(WINDOW_S * fs), int(STEP_S * fs)
    n_fft = 1 << int(np.ceil(np.log2(2 * n_win - 1)))
    acc = np.zeros((len(rx_ch), n_fft // 2 + 1), dtype=np.complex128)
    nw = 0
    for s in range(0, x.shape[1] - n_win + 1, n_step):
        S = np.conj(np.fft.rfft(x[0, s:s + n_win], n=n_fft))
        for r in range(1, x.shape[0]):
            acc[r - 1] += S * np.fft.rfft(x[r, s:s + n_win], n=n_fft)
        nw += 1
    del x
    full = np.fft.fftshift(np.fft.irfft(acc / max(1, nw), n=n_fft, axis=-1), axes=-1)
    lags = (np.arange(n_fft) - n_fft // 2) / fs
    keep = np.abs(lags) <= MAX_LAG_S
    gather = sosfiltfilt(butter(4, list(BAND), btype="bandpass", fs=fs, output="sos"),
                         full, axis=-1)[:, keep]
    lags = lags[keep]
    sep = rx_tvd - tvd_src           # vertical separation, the moveout variable
    print("  %.2f h in %d windows" % (hours, nw), flush=True)

    stem = Path(str(STEM) + (("_" + a.tag) if a.tag else ""))
    np.savez_compressed(str(stem) + ".npz", gather=gather, lags=lags,
                        tvd=rx_tvd, separation=sep, source_channel=src,
                        source_tvd=tvd_src, hours=hours, n_windows=nw,
                        limb=a.limb, fs=fs,
                        t_first=str(rows.time.iloc[0]), t_last=str(rows.time.iloc[-1]))

    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "font.size": 9.5,
        "axes.titlesize": 10, "axes.labelsize": 9.5, "axes.edgecolor": "#b0b0b0",
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "axes.labelcolor": INK, "legend.frameon": False, "legend.fontsize": 8,
    })
    # Spell the configuration out. The filename tags (src400, cm, top211, pubcm)
    # are shorthand nobody should have to decode from a figure.
    bits = ["virtual source: channel %d at %.0f m depth" % (src, tvd_src),
            "%s limb" % a.limb,
            "receivers %.0f-%.0f m depth (near-vertical section)"
            % (rx_tvd.min(), rx_tvd.max()),
            "common mode removed" if a.common_mode else "no common-mode removal"]
    label, utc_note = geo.figure_label(rows.time.iloc[0], rows.time.iloc[-1], hours, "Deep",
                             extra="; ".join(bits))
    tn = gather / np.maximum(np.abs(gather).max(axis=1, keepdims=True), 1e-30)
    env = np.abs(hilbert(gather, axis=1))
    env = env / np.maximum(env.max(axis=1, keepdims=True), 1e-30)

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 6.0), constrained_layout=True, sharey=True)
    fig.suptitle(label, fontsize=11)
    fig.text(0.995, 0.002, utc_note, ha="right", va="bottom",
             fontsize=6.5, color="#9a9a9a")
    lim = float(np.percentile(np.abs(tn), 97.0))
    ext = [lags[0], lags[-1], rx_tvd[-1], rx_tvd[0]]
    for k, (dat, cmap, vmin, vmax, ttl) in enumerate((
            (tn, "RdBu_r", -lim, lim, "(a) Gather vs true vertical depth"),
            (tn, "RdBu_r", -lim, lim, "(b) With moveout curves"),
            (env, "magma", 0.0, float(np.percentile(env, 99.0)), "(c) Envelope"))):
        ax[k].imshow(dat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                     interpolation="nearest", extent=ext)
        ax[k].set(xlim=(-MAX_LAG_S, MAX_LAG_S),
                  xlabel="correlation lag (s)", title=ttl)
        ax[k].grid(False)
    ax[0].set_ylabel("true vertical depth (m)")
    for v, c, ls, lab in ((V_MARK, C2, "-", "%.0f m/s" % V_MARK),
                          (V_LELL, C1, "--", "%.0f m/s (Lellouch)" % V_LELL)):
        ax[1].plot(sep / v, rx_tvd, ls, color=c, lw=1.6, label=lab)
        ax[1].plot(-sep / v, rx_tvd, ls, color=c, lw=1.0, alpha=0.5)
    ax[1].legend(loc="lower right")
    ax[2].plot(sep / V_MARK, rx_tvd, "-", color="w", lw=1.3, alpha=0.85)
    fig.savefig(str(stem) + ".png", bbox_inches="tight")
    print("wrote %s.{npz,png}" % stem.name, flush=True)

    step = max(1, rx_ch.size // 55)
    fig2, ax2 = plt.subplots(figsize=(7.8, 8.8), constrained_layout=True)
    scale = 2.2 * step * float(np.median(np.diff(rx_tvd)))
    for z, row in zip(rx_tvd[::step], tn[::step]):
        y = z - row * scale
        ax2.plot(lags, y, "-", color=INK, lw=0.6)
        ax2.fill_between(lags, z, y, where=(y < z), color=INK, alpha=0.55, lw=0)
    ax2.plot(sep / V_MARK, rx_tvd, "-", color=C2, lw=1.7, label="%.0f m/s" % V_MARK)
    ax2.plot(sep / V_LELL, rx_tvd, "--", color=C1, lw=1.4,
             label="%.0f m/s (Lellouch)" % V_LELL)
    ax2.set(xlim=(-0.25, 0.25), ylim=(rx_tvd[-1], rx_tvd[0]),
            xlabel="correlation lag (s)", ylabel="true vertical depth (m)",
            title=label)
    fig2.text(0.995, 0.002, utc_note, ha="right", va="bottom", fontsize=6.5,
              color="#9a9a9a")
    ax2.legend(loc="lower right")
    fig2.savefig(str(stem) + "_wiggle.png", bbox_inches="tight")
    print("wrote %s_wiggle.png" % stem.name, flush=True)


if __name__ == "__main__":
    main()
