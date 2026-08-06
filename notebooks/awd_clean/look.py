"""One figure, one thing. Single-panel looks at the AWD data.

`plain_look.py` is the dense diagnostic set (six panels a figure).
`plain_look_simple.py` was an attempt at a presentation set, but it still put
three or four unrelated subplots on one canvas, which is the thing that made it
hard to read.

This is the version to actually look at: **each file is one panel showing one
thing.** No subplot grids. Put them side by side yourself if you want to compare;
that way you choose the pairing rather than inheriting mine.

Language is at the same register as the last pass -- real terms, real numbers,
glossed on first use -- because that part was right. Only the composition
changed.

Outputs (each a single panel)
-----------------------------
look01_section_20_50Hz.png     the record section, annotated
look02_section_unfiltered.png  the same drop with no filter, for contrast
look03_traces.png              seven channels as ordinary traces
look04_overlay_150m.png        20 drops overlaid, shallow
look05_overlay_350m.png        20 drops overlaid, mid
look06_overlay_550m.png        20 drops overlaid, deep
look07_common_mode.png         what the automatic cleanup removes
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plain_look import (  # noqa: E402
    DX_NANO, NANO_DIR, OUT_DIR, PRE_S, RAW_KW,
    bandpass, ch_at, nano_time, pick_burst, read_manifest, section_around,
)
from DASutils import readFile_protobuf  # noqa: E402

BAND = (20.0, 50.0)
V_DIRECT = 2975.0
FIGSIZE = (11, 8.5)          # one panel, sized for a poster block
NOTE = dict(fontsize=10.5, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.5", fc="#f5f5f5", ec="#aaaaaa"))


def despike_channels(sec, factor=3.0, half=40):
    """Replace channels whose RMS is far above their neighbours' by interpolation.

    fig11_bad_channels.py identified these as optical artifacts narrower than the
    16.5 m gauge length, so they cannot be ground motion. Left in, they draw
    horizontal stripes across the whole record and are the first thing the eye
    lands on. Replaced by the mean of the nearest clean channels either side."""
    rms = np.sqrt(np.mean(np.asarray(sec, float) ** 2, axis=1))
    base = np.array([np.median(rms[max(0, i - half):i + half + 1])
                     for i in range(rms.size)])
    bad = rms > factor * base
    out = np.array(sec, dtype=float, copy=True)
    good = np.flatnonzero(~bad)
    for i in np.flatnonzero(bad):
        nb = good[np.argsort(np.abs(good - i))[:2]]
        out[i] = out[nb].mean(axis=0)
    return out, int(bad.sum())


def note(fig, what, why):
    """Two short lines under the panel: what it is, then why it matters.

    The 'why' line is the one that was missing. A reader should not have to
    reconstruct the purpose of a figure from its axes."""
    fig.subplots_adjust(bottom=0.24)
    fig.text(0.02, 0.175, f"{what}\n\nWHY IT MATTERS:  {why}", **NOTE)


def draw_section(sec, fs, band, title, fname, what, why, annotate=False):
    sec, n_bad = despike_channels(sec)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    t = np.arange(sec.shape[1]) / fs - PRE_S
    z = np.arange(sec.shape[0]) * DX_NANO
    v = np.percentile(np.abs(sec), 99.7)   # only the strongest 0.3% saturates,
    im = ax.pcolormesh(t, z, sec, cmap="seismic",   # so ambient noise stays pale
                       vmin=-v, vmax=v, shading="auto")
    ax.set_ylim(600, 0)
    ax.set_xlim(-0.1, 0.8)                 # the arrival lives here; 0.8-2.0 s was
                                           # half the panel and all of it noise
    ax.set_xlabel("time since drop (s)", fontsize=12)
    ax.set_ylabel("distance along fiber (m)", fontsize=12)
    ax.set_title(title, fontsize=14, pad=12)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(r"strain rate ($\mu\varepsilon$ s$^{-1}$)", fontsize=11)
    if annotate:
        # stop the guide line where the arrival is still traceable, so it does
        # not run through the region the annotation calls empty
        zz = np.linspace(0, 450, 50)
        ax.plot(zz / V_DIRECT, zz, "k--", lw=1.8, alpha=0.85)
        ax.annotate(f"direct arrival, {V_DIRECT:.0f} m s$^{{-1}}$",
                    xy=(0.105, 300), xytext=(0.33, 200), fontsize=13,
                    arrowprops=dict(arrowstyle="->", lw=2, color="k"),
                    bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="k"))
        ax.annotate("no arrival distinguishable\nbelow ~450 m",
                    xy=(0.22, 510), xytext=(0.36, 540), fontsize=12,
                    ha="left", annotation_clip=True,
                    arrowprops=dict(arrowstyle="->", lw=1.8, color="0.3"),
                    bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.3"))
    note(fig, what + f"\n{n_bad} instrumental bad channels interpolated over.", why)
    fig.savefig(OUT_DIR / fname, dpi=150)
    plt.close(fig)


def draw_traces(sec, fs, fname):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    filt = bandpass(sec, fs, BAND)
    depths = (75, 150, 225, 300, 375, 450, 525)
    t = np.arange(filt.shape[1]) / fs - PRE_S
    win = (t >= -0.1) & (t <= 0.6)
    rows = [filt[ch_at(zq, DX_NANO, filt.shape[0])][win] for zq in depths]
    scale = 55.0 / max(np.max(np.abs(r)) for r in rows)
    for zq, tr in zip(depths, rows):
        ax.plot(t[win], scale * tr + zq, "k", lw=1.3)
        ax.text(-0.093, zq - 20, f"{zq} m", fontsize=11, color="0.35")
    zz = np.linspace(0, 600, 50)
    ax.plot(zz / V_DIRECT, zz, "r--", lw=2.0, alpha=0.85,
            label=f"{V_DIRECT:.0f} m s$^{{-1}}$ moveout")
    ax.set_ylim(600, 0)
    ax.set_xlim(-0.1, 0.6)
    ax.set_xlabel("time since drop (s)", fontsize=12)
    ax.set_ylabel("distance along fiber (m)", fontsize=12)
    ax.set_title("The same data as seven ordinary seismograms, 20-50 Hz",
                 fontsize=14, pad=12)
    ax.legend(fontsize=11, loc="lower right")
    note(fig,
         "Each trace is one row of the record section, drawn as an ordinary\n"
         "seismogram. The arrival moves to later time with depth.",
         "the colourful image is not exotic -- it is these same traces stacked up.\n"
         "The slope of the red line is 1/velocity, which is how the 2975 m/s is measured.")
    fig.savefig(OUT_DIR / fname, dpi=150)
    plt.close(fig)


def draw_overlay(nano_all, fs, t_file, drops, zq, fname, vmax):
    """All drops of one burst at one depth. One depth per figure."""
    c = ch_at(zq, DX_NANO, nano_all.shape[0])
    kept = []
    for d in drops:
        s_ = section_around(nano_all[c:c + 1], fs, t_file, d["utc_time"],
                            pre=0.5, post=1.0)
        if s_ is None:
            continue
        tr = bandpass(s_, fs, BAND)[0]
        t = np.arange(tr.size) / fs - 0.5
        w = (t >= -0.05) & (t <= 0.35)
        kept.append((t[w], tr[w]))

    M = np.asarray([tr for _, tr in kept])
    M = M - M.mean(axis=1, keepdims=True)
    nrm = np.linalg.norm(M, axis=1)
    C = (M @ M.T) / np.outer(nrm, nrm)
    iu = np.triu_indices(len(M), k=1)
    cc = float(np.median(C[iu]))
    pk = np.max(np.abs(M), axis=1)
    cv = float(pk.std() / pk.mean())

    fig, ax = plt.subplots(figsize=FIGSIZE)
    for t, tr in kept:
        ax.plot(t, tr, lw=1.0, alpha=0.65)
    ax.set_xlim(-0.05, 0.35)
    ax.set_ylim(-1.1 * vmax, 1.1 * vmax)
    ax.grid(alpha=0.3)
    ax.set_xlabel("time since drop (s)", fontsize=12)
    ax.set_ylabel(r"strain rate ($\mu\varepsilon$ s$^{-1}$)", fontsize=12)
    ax.set_title(f"How repeatable is the source at {zq:.0f} m? "
                 f"{len(kept)} drops overlaid\n"
                 f"median pairwise CC = {cc:.2f},   peak-amplitude CV = {cv:.0%}",
                 fontsize=14, pad=12)
    note(fig,
         f"All {len(kept)} drops of one burst at this one depth, 20-50 Hz. All three\n"
         "overlay figures share a vertical scale, so the decay with depth is real.",
         "if the source repeated perfectly the lines would coincide. How much they\n"
         "scatter sets the smallest velocity change this experiment could ever detect.")
    fig.savefig(OUT_DIR / fname, dpi=150)
    plt.close(fig)
    return cc, cv


def draw_common_mode(sec, fs, fname):
    filt = bandpass(sec, fs, BAND)
    med = np.median(filt, axis=0)
    pct = 100 * np.sqrt(np.mean(med ** 2)) / np.sqrt(np.mean(filt ** 2))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    t = np.arange(filt.shape[1]) / fs - PRE_S
    z = np.arange(filt.shape[0]) * DX_NANO
    v = np.percentile(np.abs(filt), 99.0)
    im = ax.pcolormesh(t, z, np.tile(med, (filt.shape[0], 1)),
                       cmap="seismic", vmin=-v, vmax=v, shading="auto")
    ax.set_ylim(600, 0)
    ax.set_xlim(-0.3, 2.0)
    ax.set_xlabel("time since drop (s)", fontsize=12)
    ax.set_ylabel("distance along fiber (m)", fontsize=12)
    ax.set_title("Does the automatic cleanup step matter? No:\n"
                 f"it removes {pct:.2f}% of the signal", fontsize=14, pad=12)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(r"strain rate ($\mu\varepsilon$ s$^{-1}$)", fontsize=11)
    note(fig,
         "The reading software subtracts the across-channel median at every time\n"
         "sample, automatically. This is that removed part, on the data's own scale.",
         "it is a hidden processing step nobody chose deliberately. Drawn here it is\n"
         "almost invisible, so no 20-50 Hz result depends on it either way.")
    fig.savefig(OUT_DIR / fname, dpi=150)
    plt.close(fig)
    return pct


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_manifest()
    burst_id, drops = pick_burst(rows)
    drops.sort(key=lambda r: r["utc_time"])
    name = drops[len(drops) // 2]["nano_file"]
    drops = [d for d in drops if d["nano_file"] == name]
    t_file = nano_time(name)
    print(f"burst {burst_id}, {len(drops)} drops, {name}", flush=True)

    nano, info = readFile_protobuf([str(NANO_DIR / name)], fmin=1.0, fmax=250.0,
                                   desampling=False, **RAW_KW)
    fs = float(info["fs"])
    t_drop = drops[len(drops) // 2]["utc_time"]
    sec = section_around(nano, fs, t_file, t_drop)

    draw_section(bandpass(sec, fs, BAND), fs, BAND,
                 "One weight drop, recorded all the way down the fiber\n"
                 f"(20-50 Hz, burst {burst_id})",
                 "look01_section_20_50Hz.png",
                 "One row per fiber channel, 1.27 m apart. Colour is strain rate.\n"
                 "The arrival leans over because it reaches deeper channels later.",
                 "this is the raw material for everything else in the project. It shows\n"
                 "the experiment works, and that it stops working below about 450 m.",
                 annotate=True)

    draw_section(sec, fs, None,
                 "The same drop unfiltered: the arrival survives only to ~100 m\n"
                 f"(burst {burst_id})",
                 "look02_section_unfiltered.png",
                 "Identical data, no bandpass. Compare with the previous figure, which\n"
                 "is this same drop filtered to 20-50 Hz.",
                 "filtering is what buys the aperture. Unfiltered the arrival is clear\n"
                 "in the top ~100 m and lost below it; at 20-50 Hz it is traceable to\n"
                 "~450 m. Same data, four times the usable range.")

    draw_traces(sec, fs, "look03_traces.png")

    # one shared vertical scale across the three overlay figures
    vmax = 0.0
    for zq in (150.0, 350.0, 550.0):
        c = ch_at(zq, DX_NANO, nano.shape[0])
        for d in drops:
            s_ = section_around(nano[c:c + 1], fs, t_file, d["utc_time"],
                                pre=0.5, post=1.0)
            if s_ is None:
                continue
            tr = bandpass(s_, fs, BAND)[0]
            t = np.arange(tr.size) / fs - 0.5
            w = (t >= -0.05) & (t <= 0.35)
            vmax = max(vmax, float(np.max(np.abs(tr[w]))))

    for zq, fname in [(150.0, "look04_overlay_150m.png"),
                      (350.0, "look05_overlay_350m.png"),
                      (550.0, "look06_overlay_550m.png")]:
        cc, cv = draw_overlay(nano, fs, t_file, drops, zq, fname, vmax)
        print(f"  {zq:.0f} m: CC {cc:.2f}, CV {cv:.0%}")

    pct = draw_common_mode(sec, fs, "look07_common_mode.png")
    print(f"  common mode: {pct:.2f}% of section RMS")
    print("wrote", OUT_DIR)


if __name__ == "__main__":
    main()
