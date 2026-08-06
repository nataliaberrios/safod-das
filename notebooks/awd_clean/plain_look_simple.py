"""The same plain looks as `plain_look.py`, redrawn to be readable on first sight.

`plain_look.py` produces diagnostic figures: six panels at a time, dense, with no
annotation. This produces the same information as four presentation figures --
one result each, annotated, with the measured quantity in the title.

Register: written for a scientist who is not an active-source seismologist.
Acquisition terms are used, not avoided, but glossed on first use, and every
claim carries the number behind it. The point is to remove the jargon barrier,
not the content.

Nothing is recomputed differently. Same data, same burst, same bands. Only the
drawing changes. The detailed versions stay where they are.

Design rules used here, since they are the whole point:
  * one result per figure, quantified in the title;
  * annotation on the image identifying the arrival and the aperture limit;
  * every axis and colour bar in physical units;
  * a self-contained caption, so the figure survives being pasted into a poster
    with no surrounding text.

Outputs
-------
figures/awd_2026/plain_look/simple01..simple04*.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch  # noqa: F401  (kept for caption box styling)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plain_look import (  # noqa: E402
    DX_NANO, NANO_DIR, OUT_DIR, PRE_S, RAW_KW,
    bandpass, ch_at, nano_time, norm_rows, pick_burst, read_manifest,
    section_around,
)
from DASutils import readFile_protobuf  # noqa: E402

BAND = (20.0, 50.0)
V_DIRECT = 2975.0          # measured elsewhere in this repo; drawn only as a guide
CAPTION_KW = dict(fontsize=11, va="top", ha="left", wrap=True,
                  bbox=dict(boxstyle="round,pad=0.6", fc="#f4f4f4", ec="#999999"))


def caption(fig, text, height=0.26):
    """Put the caption in its own axes at the bottom, so it always has room.

    fig.text with va='top' grows downward off the canvas and silently truncates,
    which is how the first version lost four of its six lines."""
    fig.subplots_adjust(bottom=height + 0.06)
    ax = fig.add_axes([0.02, 0.01, 0.96, height])
    ax.axis("off")
    ax.text(0.0, 1.0, text, transform=ax.transAxes, **CAPTION_KW)


def section_image(ax, sec, fs, dx, band_label, zmax=None, cmap="seismic"):
    """One record section on a single absolute scale, with units on everything.

    Deliberately NOT normalised per trace: scaling every row to its own peak
    lifts pure-noise channels to full colour and the whole panel reads as
    static. On one shared scale the arrival is bright and the noise stays
    faint, which is both easier to read and more honest."""
    t = np.arange(sec.shape[1]) / fs - PRE_S
    z = np.arange(sec.shape[0]) * dx
    v = np.percentile(np.abs(sec), 99.0)
    im = ax.pcolormesh(t, z, sec, cmap=cmap, vmin=-v, vmax=v, shading="auto")
    ax.invert_yaxis()
    if zmax:
        ax.set_ylim(zmax, 0)
    ax.set_xlabel("time since drop (s)", fontsize=11)
    ax.set_ylabel("distance along fiber (m)", fontsize=11)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(f"strain rate, {band_label}\n"
                 r"($\mu\varepsilon$ s$^{-1}$, common scale)", fontsize=9)
    return im


# --------------------------------------------------------------------------

def simple01(sec, fs, tag):
    """What the experiment records, with the arrival pointed at."""
    fig, ax = plt.subplots(1, 2, figsize=(16, 11),
                           gridspec_kw=dict(width_ratios=[1.5, 1]))

    filt = bandpass(sec, fs, BAND)
    section_image(ax[0], filt, fs, DX_NANO, "20-50 Hz", zmax=600)
    zz = np.linspace(0, 600, 50)
    ax[0].plot(zz / V_DIRECT, zz, "k--", lw=1.6, alpha=0.8)
    ax[0].annotate(f"direct arrival, {V_DIRECT:.0f} m s$^{{-1}}$ apparent",
                   xy=(0.10, 300), xytext=(0.55, 190), fontsize=12,
                   arrowprops=dict(arrowstyle="->", lw=2, color="k"),
                   bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="k", alpha=0.9))
    ax[0].annotate("aperture limit: arrival indistinguishable\nfrom noise below ~450 m",
                   xy=(0.6, 520), xytext=(1.05, 470), fontsize=12,
                   arrowprops=dict(arrowstyle="->", lw=2, color="0.3"),
                   bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.3", alpha=0.9))
    ax[0].set_xlim(-0.3, 2.0)
    ax[0].set_title("Record section: one row per channel, 1.27 m spacing", fontsize=13)

    # the same thing as ordinary wiggly lines, which is what it really is
    depths = (75, 150, 225, 300, 375, 450, 525)
    t = np.arange(filt.shape[1]) / fs - PRE_S
    win = (t >= -0.1) & (t <= 0.6)
    rows = [filt[ch_at(zq, DX_NANO, filt.shape[0])][win] for zq in depths]
    scale = 55.0 / max(np.max(np.abs(r)) for r in rows)   # one scale for all
    for zq, tr in zip(depths, rows):
        ax[1].plot(t[win], scale * tr + zq, "k", lw=1.1)
        ax[1].text(-0.095, zq - 22, f"{zq} m", fontsize=10, color="0.35")
    ax[1].plot(zz / V_DIRECT, zz, "r--", lw=2.0, alpha=0.85,
               label=f"{V_DIRECT:.0f} m s$^{{-1}}$ moveout")
    ax[1].invert_yaxis()
    ax[1].set_xlim(-0.1, 0.6)
    ax[1].set_ylim(600, 0)
    ax[1].set_xlabel("time since drop (s)", fontsize=11)
    ax[1].set_ylabel("distance along fiber (m)", fontsize=11)
    ax[1].legend(fontsize=10, loc="lower right")
    ax[1].set_title("Seven of those channels as individual traces", fontsize=13)

    fig.suptitle("1.  Surface weight drop recorded on borehole fiber: direct "
                 f"arrival at {V_DIRECT:.0f} m s$^{{-1}}$ apparent velocity, "
                 "usable to ~450 m", fontsize=15)
    caption(fig,
            "A \"record section\" is the wavefield u(z, t) as an image: one row per "
            "fiber channel, time across,\namplitude in colour. The right panel is the "
            "same data as seven individual traces, which is what\neach row actually "
            "is.\n\n"
            "Arrival time increases linearly with depth, t(z) = z/v, so the arrival "
            "appears as a straight line whose\nslope is 1/v. Fitting that slope gives "
            f"{V_DIRECT:.0f} m/s. \"Apparent\" because DAS is sensitive only to "
            "strain along\nthe fiber, so a wave crossing it obliquely measures faster "
            "than it travels.\n\n"
            f"{tag}")
    fig.savefig(OUT_DIR / "simple01_what_we_record.png", dpi=140,
                bbox_inches="tight")
    plt.close(fig)


def simple02(sec, fs, tag):
    """Which frequencies carry the signal -- four panels, each with a verdict."""
    picks = [(None, "broadband, unfiltered", "noise-dominated"),
             ((1.0, 5.0), "1-5 Hz", "coherent, but zero moveout: not a body wave"),
             ((20.0, 50.0), "20-50 Hz", "signal band: clear moveout to ~450 m"),
             ((100.0, 250.0), "100-250 Hz", "below the noise floor")]
    fig, axes = plt.subplots(1, 4, figsize=(19, 10.5), sharey=True)
    for a, (band, label, verdict) in zip(axes, picks):
        filt = bandpass(sec, fs, band) if band else sec
        t = np.arange(filt.shape[1]) / fs - PRE_S
        z = np.arange(filt.shape[0]) * DX_NANO
        # each panel on its own absolute scale, so a band with no signal looks
        # empty rather than being stretched up to look like signal
        v = np.percentile(np.abs(filt), 99.5)
        a.pcolormesh(t, z, filt, cmap="seismic", vmin=-v, vmax=v, shading="auto")
        a.set_xlim(-0.2, 1.5)
        a.set_ylim(600, 0)
        a.set_xlabel("time (s)", fontsize=11)
        good = verdict.startswith("signal band")
        a.set_title(f"{label}\n{verdict}", fontsize=12,
                    color=("darkgreen" if good else "0.35"),
                    fontweight=("bold" if good else "normal"))
        for sp in a.spines.values():
            sp.set_edgecolor("darkgreen" if good else "0.8")
            sp.set_linewidth(3 if good else 1)
    axes[0].set_ylabel("distance along fiber (m)", fontsize=11)
    sm = plt.cm.ScalarMappable(cmap="seismic", norm=plt.Normalize(-1, 1))
    cb = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.01)
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels(["$-$", "0", "$+$"])
    cb.set_label(r"strain rate, each panel scaled to its own 99.5th percentile",
                 fontsize=10)

    fig.suptitle("2.  Band decomposition: coherent moveout is confined to "
                 "20-50 Hz", fontsize=15)
    caption(fig,
            "One drop, bandpass-filtered into four ranges (4th-order Butterworth, "
            "zero-phase). Each panel is\nscaled independently, so a band carrying no "
            "signal reads as empty rather than being stretched to\nlook like signal.\n\n"
            "Broadband is noise-dominated. The 1-5 Hz energy is coherent across the "
            "array but arrives with zero\nmoveout and starts before t=0, so it is "
            "ambient rather than source-generated. Above 100 Hz there is\nnothing above "
            "the noise floor. The 20-50 Hz band is the only one showing a "
            "moveout-consistent\narrival, which is why the rest of the analysis uses "
            "it -- an empirical choice, not an inherited one.\n\n"
            f"{tag}")
    fig.savefig(OUT_DIR / "simple02_where_the_signal_is.png", dpi=140,
                bbox_inches="tight")
    plt.close(fig)


def simple03(nano_all, fs, t_file, drops, tag):
    """How repeatable the source is, and how deep it stays usable."""
    depths = [150.0, 350.0, 550.0]
    verdicts = None   # computed below from the data
    # Filter a LONG window and display a short one. Filtering the display window
    # directly puts the biggest excursions on the figure at its two edges, and
    # those are filter ringing, not data.
    show_lo, show_hi = -0.05, 0.35
    fig, axes = plt.subplots(1, 3, figsize=(17, 10), sharey=True)
    per_depth = []
    for zq in depths:
        c = ch_at(zq, DX_NANO, nano_all.shape[0])
        kept = []
        for d in drops:
            s_ = section_around(nano_all[c:c + 1], fs, t_file, d["utc_time"],
                                pre=0.5, post=1.0)
            if s_ is None:
                continue
            tr = bandpass(s_, fs, BAND)[0]
            t = np.arange(tr.size) / fs - 0.5
            win = (t >= show_lo) & (t <= show_hi)
            kept.append((t[win], tr[win]))
        per_depth.append(kept)

    # one shared vertical scale, so the amplitude really does fall with depth
    vmax = max(np.max(np.abs(tr)) for kept in per_depth for _, tr in kept)

    # quantify rather than assert: median pairwise normalised cross-correlation
    # between drops, and the coefficient of variation of their peak amplitudes
    verdicts, stats = [], []
    for kept in per_depth:
        M = np.asarray([tr for _, tr in kept])
        M = M - M.mean(axis=1, keepdims=True)
        nrm = np.linalg.norm(M, axis=1)
        C = (M @ M.T) / np.outer(nrm, nrm)
        iu = np.triu_indices(len(M), k=1)
        cc = float(np.median(C[iu]))
        pk = np.max(np.abs(M), axis=1)
        cv = float(pk.std() / pk.mean())
        stats.append((cc, cv))
        verdicts.append(f"CC = {cc:.2f}    CV = {cv:.0%}")

    for a, zq, verdict, kept in zip(axes, depths, verdicts, per_depth):
        for t, tr in kept:
            a.plot(t, tr, lw=0.9, alpha=0.6)
        good = stats[depths.index(zq)][0] > 0.9
        a.set_title(f"{zq:.0f} m along fiber\n{verdict}", fontsize=12, pad=12,
                    color=("darkgreen" if good else "0.35"),
                    fontweight=("bold" if good else "normal"))
        a.set_xlabel("time since drop (s)", fontsize=11)
        a.grid(alpha=0.3)
        a.set_xlim(show_lo, show_hi)
        a.set_ylim(-1.1 * vmax, 1.1 * vmax)
        if good:
            for sp in a.spines.values():
                sp.set_edgecolor("darkgreen")
                sp.set_linewidth(2.5)
    axes[0].set_ylabel(r"strain rate ($\mu\varepsilon$ s$^{-1}$)", fontsize=11)
    fig.suptitle(f"3.  Source repeatability against depth: {len(drops)} drops from one "
                 f"burst, CC {stats[0][0]:.2f} at 150 m falling to "
                 f"{stats[-1][0]:.2f} at 550 m", fontsize=15)
    caption(fig,
            f"All {len(drops)} drops of one burst overlaid, at three fiber positions, "
            "20-50 Hz. Perfect source and\npath repeatability would put the traces "
            "exactly on top of each other, so the spread is the combined\n"
            "source-plus-medium variability.\n\n"
            "CC is the median pairwise normalised cross-correlation between drops; "
            "CV is the coefficient of\nvariation of their peak amplitudes. All three panels share one vertical scale, so the\namplitude "
            "decay with depth is real rather than an artefact of autoscaling. Traces "
            "are filtered over a\nlonger window than is displayed, so the plotted "
            "range is free of filter edge transients.\n\n"
            f"{tag}")
    fig.savefig(OUT_DIR / "simple03_how_repeatable.png", dpi=140,
                bbox_inches="tight")
    plt.close(fig)


def simple04(sec, fs, tag):
    """Does the standard cleanup step change the answer? No."""
    filt = bandpass(sec, fs, BAND)
    med = np.median(filt, axis=0)
    after = filt - med[None, :]
    removed_pct = 100 * np.sqrt(np.mean(med ** 2)) / np.sqrt(np.mean(filt ** 2))

    fig, ax = plt.subplots(1, 3, figsize=(17, 10.5), sharey=True)
    for a, data, title in [(ax[0], filt, "before common-mode removal"),
                           (ax[1], after, "after common-mode removal"),
                           (ax[2], filt - after, "the removed component")]:
        t = np.arange(data.shape[1]) / fs - PRE_S
        z = np.arange(data.shape[0]) * DX_NANO
        v = np.percentile(np.abs(filt), 99.0)
        im = a.pcolormesh(t, z, data, cmap="seismic", vmin=-v, vmax=v, shading="auto")
        a.set_ylim(600, 0)
        a.set_xlim(-0.3, 2.0)
        a.set_xlabel("time since drop (s)", fontsize=11)
        a.set_title(title, fontsize=13)
    ax[0].set_ylabel("distance along fiber (m)", fontsize=11)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label(r"strain rate ($\mu\varepsilon$ s$^{-1}$), common scale", fontsize=10)

    fig.suptitle("4.  Common-mode removal accounts for "
                 f"{removed_pct:.1f}% of section RMS at 20-50 Hz", fontsize=15)
    caption(fig,
            "The DAS reader subtracts the across-channel median at every time sample by "
            "default, on the\nassumption that anything simultaneous on all channels is "
            "interrogator noise. It is applied silently,\nso it is worth knowing what "
            "it costs.\n\n"
            "The concern is that it also removes any genuine arrival with zero moveout. "
            f"Here the removed\ncomponent is {removed_pct:.2f}% of the section RMS and "
            "the first two panels are visually identical, so no\nresult in the "
            "20-50 Hz band depends on this choice either way. It matters more at low "
            "frequency,\nwhere the common mode carries most of the coherent energy.\n\n"
            f"{tag}")
    fig.savefig(OUT_DIR / "simple04_does_cleanup_matter.png", dpi=140,
                bbox_inches="tight")
    plt.close(fig)
    return removed_pct


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_manifest()
    burst_id, drops = pick_burst(rows)
    drops.sort(key=lambda r: r["utc_time"])
    name = drops[len(drops) // 2]["nano_file"]
    drops = [d for d in drops if d["nano_file"] == name]
    t_file = nano_time(name)

    print(f"burst {burst_id}, {len(drops)} drops, file {name}", flush=True)
    nano, info = readFile_protobuf([str(NANO_DIR / name)], fmin=1.0, fmax=250.0,
                                   desampling=False, **RAW_KW)
    fs = float(info["fs"])
    print(f"  {nano.shape} at {fs} Hz", flush=True)

    t_drop = drops[len(drops) // 2]["utc_time"]
    sec = section_around(nano, fs, t_file, t_drop)
    tag = (f"Data: burst {burst_id}, {len(drops)} drops, "
           f"{t_drop:%Y-%m-%d %H:%M} UTC, cemented (Nano) fiber.")

    simple01(sec, fs, tag)
    simple02(sec, fs, tag)
    simple03(nano, fs, t_file, drops, tag)
    pct = simple04(sec, fs, tag)
    print(f"cleanup removes {pct:.2f}% of the section RMS")
    print("wrote", OUT_DIR)


if __name__ == "__main__":
    main()
