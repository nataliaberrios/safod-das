"""Plain, unprocessed looks at the raw SAFOD AWD data.

Everything else in this directory starts from a bandpass and a stack. This does
not. It reads a small number of raw Nano and Deep files with every reader option
turned off and makes the boring plots -- record sections, spectra, spectrograms,
per-channel statistics, drops overlaid on each other -- so the data can be looked
at before any statistic is computed on it.

What the reader does before you see anything
--------------------------------------------
``DASutils.readFile_protobuf`` and ``readFile_HDF`` both default to
``detrend=True``, ``tapering=True`` (Tukey, alpha=0.4), ``filter=True`` (the
``fmin``/``fmax`` you pass) and ``median=True``. ``median=True`` calls
``preprocess_medfilt``, which subtracts the across-channel median at every time
sample -- common-mode removal, usually applied to suppress laser and
interrogator noise shared by all channels.

Whether that is the right choice here is a question this script asks rather than
answers. Common-mode removal also deletes any real arrival that reaches every
channel at nearly the same time, so it is a physics assumption, not a cleanup
step. Figure 7 draws the whole reader chain one step at a time and figure 9
compares the section, the traces and the spectra with and without it. The wider
problem is that all of this is invisible: the steps happen inside one function
call, the order matters, and several AWD scripts then apply a *second* bandpass
to data the reader has already filtered.

Reader settings used here: ``filter=False, median=False, detrend=False,
tapering=False, desampling=False``, so figures 1-6 show what is actually on
disk. ``fmin``/``fmax`` are positional in the DASutils signature so they are
still passed, but with ``filter=False`` they are inert.

Outputs
-------
figures/awd_2026/plain_look/fig01..fig08*.png and plain_look.npz
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt, welch, spectrogram

SEARCH_DIRS = [
    "/home/groups/edunham/nberrios/safod_das/DAS-utilities/python",
    "/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python",
]
for _p in SEARCH_DIRS:
    if Path(_p).exists() and _p not in sys.path:
        sys.path.insert(0, _p)
from DASutils import readFile_protobuf, readFile_HDF  # noqa: E402

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "awd_manifest.csv"
NANO_DIR = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/Nano")
DEEP_DIR = Path(
    "/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/"
    "01_--_recording_2026-06-15T230629Z_--_active_source"
)
OUT_DIR = Path(
    "/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026/plain_look"
)

PRE_S, POST_S = 0.5, 3.0
DX_NANO, DX_DEEP = 1.26606202, 2.0419

# Deliberately spans the whole recorded band rather than the inherited working
# bands. Nyquist is 500 Hz at 1 kHz sampling; the top band is clipped to it.
BANDS = [None, (1.0, 5.0), (5.0, 20.0), (20.0, 50.0), (50.0, 100.0), (100.0, 250.0)]
# Distances along fiber to pull single traces from, in metres.
TRACE_DEPTHS_M = [50.0, 150.0, 250.0, 350.0, 450.0, 550.0]

RAW_KW = dict(filter=False, median=False, detrend=False, tapering=False)


# ----------------------------------------------------------------- utilities

def nano_time(name: str) -> datetime | None:
    parts = os.path.basename(name).split("_")
    try:
        return datetime.strptime(
            f"{parts[1]}T{parts[2].replace('.', ':')}Z", "%Y-%m-%dT%H:%M:%SZ"
        ).replace(tzinfo=timezone.utc)
    except (IndexError, ValueError):
        return None


def deep_time(name: str) -> datetime | None:
    m = re.search(r"_(\d{4}-\d{2}-\d{2}T\d{6}Z)\.h5$", os.path.basename(name))
    if not m:
        return None
    return datetime.strptime(m.group(1), "%Y-%m-%dT%H%M%SZ").replace(tzinfo=timezone.utc)


def bandpass(x, fs, band):
    if band is None:
        return x
    hi = min(band[1], 0.45 * fs)
    sos = butter(4, [band[0], hi], btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def read_manifest():
    """Return the manifest rows as a list of dicts, times parsed."""
    import csv

    rows = []
    with open(MANIFEST) as fh:
        for r in csv.DictReader(fh):
            r["utc_time"] = datetime.fromisoformat(r["utc_time"])
            r["burst_id"] = int(r["burst_id"])
            r["nano_available"] = r["nano_available"] == "1"
            r["deep_available"] = r["deep_available"] == "1"
            rows.append(r)
    return rows


def pick_burst(rows):
    """The burst with the most drops available on both fibers from one file each.

    Chosen mechanically, not by how the data look."""
    from collections import Counter

    ok = [r for r in rows if r["nano_available"] and r["deep_available"]]
    counts = Counter(r["burst_id"] for r in ok)
    for burst_id, _ in counts.most_common():
        drops = [r for r in ok if r["burst_id"] == burst_id]
        if len({d["nano_file"] for d in drops}) == 1:
            return burst_id, drops
    burst_id = counts.most_common(1)[0][0]
    return burst_id, [r for r in ok if r["burst_id"] == burst_id]


def section_around(data, fs, t_file, t_drop, pre=PRE_S, post=POST_S):
    i = int((t_drop - t_file).total_seconds() * fs)
    i0, i1 = i - int(pre * fs), i + int(post * fs)
    if i0 < 0 or i1 > data.shape[1]:
        return None
    return data[:, i0:i1]


def norm_rows(sec):
    """Per-trace peak normalisation, for display only."""
    peak = np.max(np.abs(sec), axis=1, keepdims=True)
    return sec / np.where(peak > 0, peak, 1.0)


def image(ax, sec, fs, dx, clip=99.0, pre=PRE_S):
    t = np.arange(sec.shape[1]) / fs - pre
    z = np.arange(sec.shape[0]) * dx
    v = np.percentile(np.abs(sec), clip)
    im = ax.pcolormesh(t, z, sec, cmap="seismic", vmin=-v, vmax=v, shading="auto")
    ax.invert_yaxis()
    ax.set_xlabel("time after drop (s)")
    ax.set_ylabel("distance along fiber (m)")
    return im


def ch_at(depth_m, dx, n_ch):
    return min(int(round(depth_m / dx)), n_ch - 1)


# ----------------------------------------------------------------- figures

def fig01_raw_sections(nano_raw, deep_raw, fs_n, fs_d, tag):
    """The record section with nothing done to it, and the same trace-normalised."""
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    for row, (sec, fs, dx, name) in enumerate(
        [(nano_raw, fs_n, DX_NANO, "Nano (cemented)"),
         (deep_raw, fs_d, DX_DEEP, "Deep (wireline)")]
    ):
        if sec is None:
            continue
        im = image(ax[row, 0], sec, fs, dx)
        ax[row, 0].set_title(f"{name}: raw amplitude, no filter, no median removal")
        plt.colorbar(im, ax=ax[row, 0], label="raw units")
        im = image(ax[row, 1], norm_rows(sec), fs, dx)
        ax[row, 1].set_title(f"{name}: same section, each trace scaled to its own peak")
        plt.colorbar(im, ax=ax[row, 1], label="normalised")
    fig.suptitle(f"1  One weight drop, unprocessed -- {tag}", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig01_raw_record_section.png", dpi=140)
    plt.close(fig)


def fig02_bands(sec, fs, dx, name, tag):
    """The same drop seen through a ladder of frequency bands."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for a, band in zip(axes.ravel(), BANDS):
        filt = bandpass(sec, fs, band)
        im = image(a, norm_rows(filt), fs, dx)
        label = "unfiltered" if band is None else f"{band[0]:g}-{band[1]:g} Hz"
        a.set_title(label)
        plt.colorbar(im, ax=a)
    fig.suptitle(
        f"2  {name}: one drop, trace-normalised, band by band -- {tag}", fontsize=13
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig02_frequency_bands.png", dpi=140)
    plt.close(fig)


def fig03_psd(nano_raw, deep_raw, fs_n, fs_d, tag):
    """Welch spectra to Nyquist, signal window and pre-drop noise window."""
    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    store = {}
    for row, (sec, fs, dx, name) in enumerate(
        [(nano_raw, fs_n, DX_NANO, "Nano"), (deep_raw, fs_d, DX_DEEP, "Deep")]
    ):
        if sec is None:
            continue
        i0 = int(PRE_S * fs)
        nper = min(1024, i0, sec.shape[1] - i0)
        for depth in TRACE_DEPTHS_M:
            c = ch_at(depth, dx, sec.shape[0])
            f, p = welch(sec[c, i0:], fs=fs, nperseg=nper)
            ax[row, 0].loglog(f[1:], np.sqrt(p[1:]), lw=1.2, label=f"{depth:.0f} m")
            store[f"{name}_sig_{depth:.0f}"] = p
        f, p = welch(sec[:, :i0], fs=fs, nperseg=nper, axis=-1)
        ax[row, 0].loglog(f[1:], np.sqrt(np.mean(p, axis=0))[1:], "k--", lw=1.3,
                          label="pre-drop noise, all ch")
        store[f"{name}_noise_mean"] = np.mean(p, axis=0)
        store[f"{name}_freq"] = f
        for line in (60.0, 120.0, 180.0):
            if line < fs / 2:
                ax[row, 0].axvline(line, color="0.7", lw=0.8, ls=":", zorder=0)
        ax[row, 0].set_xlabel("frequency (Hz)")
        ax[row, 0].set_ylabel("amplitude spectral density")
        ax[row, 0].set_title(f"{name}: spectra to Nyquist (dotted = 60 Hz harmonics)")
        ax[row, 0].legend(fontsize=8)
        ax[row, 0].grid(alpha=0.3, which="both")

        # noise PSD against channel -- where the bad channels and lines show up
        f, pall = welch(sec[:, :i0], fs=fs, nperseg=nper, axis=-1)
        z = np.arange(sec.shape[0]) * dx
        db = 10 * np.log10(pall[:, 1:] + 1e-30)
        im = ax[row, 1].pcolormesh(f[1:], z, db, cmap="magma", shading="auto",
                                   vmin=np.percentile(db, 2), vmax=np.percentile(db, 98))
        ax[row, 1].set_xscale("log")
        ax[row, 1].invert_yaxis()
        ax[row, 1].set_xlabel("frequency (Hz)")
        ax[row, 1].set_ylabel("distance along fiber (m)")
        ax[row, 1].set_title(f"{name}: pre-drop noise PSD by channel")
        plt.colorbar(im, ax=ax[row, 1], label="dB")
        store[f"{name}_psd_by_channel"] = pall
    fig.suptitle(f"3  Noise and signal spectra, full band -- {tag}", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig03_noise_psd.png", dpi=140)
    plt.close(fig)
    return store


def fig04_spectrogram(nano_file_raw, fs_n, t_file, drops, tag):
    """Whole-file spectrograms at a few depths, with the drop times marked."""
    depths = [50.0, 250.0, 550.0]
    fig, axes = plt.subplots(len(depths), 1, figsize=(14, 10), sharex=True)
    for a, depth in zip(axes, depths):
        c = ch_at(depth, DX_NANO, nano_file_raw.shape[0])
        f, t, S = spectrogram(nano_file_raw[c], fs=fs_n, nperseg=2048, noverlap=1024)
        db = 10 * np.log10(S + 1e-30)
        im = a.pcolormesh(t, f, db, cmap="magma", shading="auto",
                          vmin=np.percentile(db, 5), vmax=np.percentile(db, 99.5))
        for d in drops:
            a.axvline((d["utc_time"] - t_file).total_seconds(), color="c", lw=0.7,
                      alpha=0.8)
        a.set_ylabel("frequency (Hz)")
        a.set_title(f"Nano at {depth:.0f} m (cyan = manifest drop times)")
        plt.colorbar(im, ax=a, label="dB")
    axes[-1].set_xlabel("seconds into file")
    fig.suptitle(f"4  Whole raw file, unfiltered spectrogram -- {tag}", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig04_spectrogram.png", dpi=140)
    plt.close(fig)


def fig05_channel_qc(nano_file_raw, fs_n, tag):
    """Per-channel statistics: which channels are dead, clipped, or odd."""
    x = nano_file_raw
    z = np.arange(x.shape[0]) * DX_NANO
    rms = np.sqrt(np.mean(x.astype(np.float64) ** 2, axis=1))
    dc = np.mean(x, axis=1)
    # lag-1 autocorrelation: near zero means the trace is white, i.e. no signal
    xc = x - dc[:, None]
    denom = np.sum(xc * xc, axis=1)
    lag1 = np.sum(xc[:, 1:] * xc[:, :-1], axis=1) / np.where(denom > 0, denom, 1.0)
    amax = np.max(np.abs(x), axis=1)
    # fraction of samples within 0.1% of that channel's own extreme value
    near = np.mean(np.abs(x) >= 0.999 * amax[:, None], axis=1)
    dead = rms <= 0

    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    ax[0, 0].semilogy(z, np.where(rms > 0, rms, np.nan), lw=0.8)
    ax[0, 0].set_xlabel("distance along fiber (m)")
    ax[0, 0].set_ylabel("RMS (raw units)")
    ax[0, 0].set_title(f"RMS per channel ({int(dead.sum())} channels exactly zero)")
    ax[0, 0].grid(alpha=0.3)

    ax[0, 1].plot(z, dc, lw=0.8)
    ax[0, 1].set_xlabel("distance along fiber (m)")
    ax[0, 1].set_ylabel("mean (raw units)")
    ax[0, 1].set_title("DC offset per channel")
    ax[0, 1].grid(alpha=0.3)

    ax[1, 0].plot(z, lag1, lw=0.8)
    ax[1, 0].axhline(0, color="k", lw=0.6)
    ax[1, 0].set_xlabel("distance along fiber (m)")
    ax[1, 0].set_ylabel("lag-1 autocorrelation")
    ax[1, 0].set_title("Whiteness: values near 0 are noise-only channels")
    ax[1, 0].grid(alpha=0.3)

    ax[1, 1].semilogy(z, np.where(near > 0, near, np.nan), lw=0.8)
    ax[1, 1].set_xlabel("distance along fiber (m)")
    ax[1, 1].set_ylabel("fraction of samples at the channel extreme")
    ax[1, 1].set_title("Flat-topping / clipping indicator")
    ax[1, 1].grid(alpha=0.3)

    fig.suptitle(f"5  Nano per-channel QC over one whole raw file -- {tag}", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig05_channel_qc.png", dpi=140)
    plt.close(fig)
    return dict(z=z, rms=rms, dc=dc, lag1=lag1, near_extreme=near)


def fig06_repeatability(nano_file_raw, fs_n, t_file, drops, tag):
    """Every drop in the burst, overlaid. No stacking, no metric."""
    depths = [150.0, 350.0, 550.0]
    fig, ax = plt.subplots(len(depths), 2, figsize=(14, 10))
    peaks, times = [], []
    for row, depth in enumerate(depths):
        c = ch_at(depth, DX_NANO, nano_file_raw.shape[0])
        for band, col in [(None, 0), ((20.0, 50.0), 1)]:
            a = ax[row, col]
            n_used = 0
            for d in drops:
                sec = section_around(nano_file_raw[c:c + 1], fs_n, t_file,
                                     d["utc_time"], pre=0.1, post=0.4)
                if sec is None:
                    continue
                tr = bandpass(sec, fs_n, band)[0]
                t = np.arange(tr.size) / fs_n - 0.1
                a.plot(t, tr, lw=0.6, alpha=0.55)
                n_used += 1
                if col == 1:
                    peaks.append(np.max(np.abs(tr)))
                    times.append(d["utc_time"])
            a.set_xlabel("time after drop (s)")
            a.set_ylabel("raw units")
            label = "unfiltered" if band is None else "20-50 Hz"
            a.set_title(f"{depth:.0f} m, {label} -- {n_used} drops overlaid")
            a.grid(alpha=0.3)
    fig.suptitle(
        f"6  Drop-to-drop repeatability, every trace drawn -- {tag}", fontsize=13
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig06_drop_repeatability.png", dpi=140)
    plt.close(fig)
    return np.asarray(peaks), times


def fig07_processing_ladder(nano_raw, nano_default, fs_n, tag, band=(20.0, 50.0)):
    """The reader's default chain, one step at a time, in the order it happens.

    Reproduces detrend -> taper -> bandpass -> median removal by hand so each
    step can be seen on its own. The last panel is the reader's own default
    output, as a check that the hand-built ladder lands in the same place."""
    from scipy.signal import detrend as sp_detrend
    from scipy.signal.windows import tukey

    t = np.arange(nano_raw.shape[1]) / fs_n - PRE_S
    steps, rms = [], []

    x = nano_raw.astype(np.float64)
    steps.append(("1. as stored on disk", x.copy()))
    x = sp_detrend(x, axis=-1, type="linear")
    steps.append(("2. + detrend (linear)", x.copy()))
    x = x * tukey(x.shape[1], alpha=0.4)[None, :]
    steps.append((r"3. + Tukey taper ($\alpha$=0.4)", x.copy()))
    x = bandpass(x, fs_n, band)
    steps.append((f"4. + bandpass {band[0]:g}-{band[1]:g} Hz", x.copy()))
    med = np.median(x, axis=0)
    x = x - med[None, :]
    steps.append(("5. + median removal (common mode)", x.copy()))
    for _, s in steps:
        rms.append(float(np.sqrt(np.mean(s ** 2))))

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for a, (title, s) in zip(axes.ravel(), steps):
        im = image(a, s, fs_n, DX_NANO)
        a.set_title(title, fontsize=10)
        plt.colorbar(im, ax=a)

    a = axes.ravel()[5]
    a.plot(t, med, lw=0.9, label="common mode removed at step 5")
    a.set_xlabel("time after drop (s)")
    a.set_ylabel(r"$\mu\varepsilon$/s")
    a.set_title("the across-channel median trace itself", fontsize=10)
    a.grid(alpha=0.3)
    a.legend(fontsize=8)

    note = "  ".join(f"{i + 1}:{v:.3g}" for i, v in enumerate(rms))
    fig.suptitle(
        f"7  What each preprocessing step does, in order -- {tag}\n"
        f"section RMS after each step -- {note}", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig07_processing_ladder.png", dpi=140)
    plt.close(fig)

    # separate check: does the hand-built ladder match the reader's own output?
    if nano_default is not None:
        fig, ax = plt.subplots(1, 3, figsize=(16, 5))
        im = image(ax[0], steps[-1][1], fs_n, DX_NANO)
        ax[0].set_title("ladder rebuilt by hand (step 5)")
        plt.colorbar(im, ax=ax[0])
        im = image(ax[1], nano_default, fs_n, DX_NANO)
        ax[1].set_title("readFile_protobuf with its defaults")
        plt.colorbar(im, ax=ax[1])
        ax[2].plot(t, steps[-1][1][ch_at(250.0, DX_NANO, nano_raw.shape[0])],
                   lw=0.9, label="hand-built")
        ax[2].plot(t, nano_default[ch_at(250.0, DX_NANO, nano_default.shape[0])],
                   lw=0.9, ls="--", label="reader default")
        ax[2].set_xlabel("time after drop (s)")
        ax[2].set_title("one trace at 250 m; bands differ, so these need not overlie")
        ax[2].legend(fontsize=8)
        ax[2].grid(alpha=0.3)
        fig.suptitle(f"7b  Hand-built ladder vs the reader -- {tag}", fontsize=12)
        fig.tight_layout()
        fig.savefig(OUT_DIR / "fig07b_ladder_check.png", dpi=140)
        plt.close(fig)
    return rms


def fig09_common_mode(nano_raw, fs_n, tag, band=(20.0, 50.0)):
    """Common-mode removal, before and after, so the choice can be judged.

    Subtracting the across-channel median suppresses whatever every channel
    shares. That is laser and interrogator noise, and it is also any real
    arrival with little or no moveout across the array. Which of those dominates
    is an empirical question, so this shows the section, single traces, spectra
    and per-channel energy both ways, plus how much of each channel's energy the
    step actually takes."""
    x = bandpass(nano_raw.astype(np.float64), fs_n, band)
    med = np.median(x, axis=0)
    y = x - med[None, :]
    t = np.arange(x.shape[1]) / fs_n - PRE_S
    z = np.arange(x.shape[0]) * DX_NANO
    i0 = int(PRE_S * fs_n)
    nper = min(512, x.shape[1] - i0)

    fig, ax = plt.subplots(2, 3, figsize=(17, 9))

    v = np.percentile(np.abs(x), 99.0)
    for a, s, title in [(ax[0, 0], x, "before: common mode kept"),
                        (ax[0, 1], y, "after: across-channel median removed"),
                        (ax[0, 2], x - y, "difference: exactly what was removed")]:
        im = a.pcolormesh(t, z, s, cmap="seismic", vmin=-v, vmax=v, shading="auto")
        a.invert_yaxis()
        a.set_xlabel("time after drop (s)")
        a.set_ylabel("distance along fiber (m)")
        a.set_title(title, fontsize=10)
        plt.colorbar(im, ax=a)

    # single traces, same vertical scale, so the change is readable
    for depth, off in zip([150.0, 350.0, 550.0], [0, 1, 2]):
        c = ch_at(depth, DX_NANO, x.shape[0])
        sc = 2.5 * np.std(x[c])
        ax[1, 0].plot(t, x[c] / sc + off, lw=0.7, color="C0")
        ax[1, 0].plot(t, y[c] / sc + off, lw=0.7, color="C3", alpha=0.8)
        ax[1, 0].text(t[0], off + 0.35, f"{depth:.0f} m", fontsize=8)
    ax[1, 0].plot([], [], color="C0", label="before")
    ax[1, 0].plot([], [], color="C3", label="after")
    ax[1, 0].set_xlabel("time after drop (s)")
    ax[1, 0].set_yticks([])
    ax[1, 0].set_title(f"traces, {band[0]:g}-{band[1]:g} Hz, common scale", fontsize=10)
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].grid(alpha=0.3)

    for depth in [150.0, 350.0, 550.0]:
        c = ch_at(depth, DX_NANO, x.shape[0])
        for arr, ls in [(x, "-"), (y, "--")]:
            f, p = welch(arr[c, i0:], fs=fs_n, nperseg=nper)
            ax[1, 1].loglog(f[1:], np.sqrt(p[1:]), ls, lw=1.0,
                            label=f"{depth:.0f} m {'before' if ls == '-' else 'after'}")
    ax[1, 1].set_xlabel("frequency (Hz)")
    ax[1, 1].set_ylabel("amplitude spectral density")
    ax[1, 1].set_title("spectra either side of the step", fontsize=10)
    ax[1, 1].legend(fontsize=7, ncol=2)
    ax[1, 1].grid(alpha=0.3, which="both")

    e_before = np.sqrt(np.mean(x ** 2, axis=1))
    e_after = np.sqrt(np.mean(y ** 2, axis=1))
    frac = 1.0 - np.divide(e_after, e_before, out=np.zeros_like(e_after),
                           where=e_before > 0)
    ax[1, 2].plot(z, 100 * frac, lw=0.8)
    ax[1, 2].axhline(0, color="k", lw=0.6)
    ax[1, 2].set_xlabel("distance along fiber (m)")
    ax[1, 2].set_ylabel("% of channel RMS removed")
    ax[1, 2].set_title("how much each channel loses", fontsize=10)
    ax[1, 2].grid(alpha=0.3)

    fig.suptitle(
        f"9  Common-mode removal, before and after -- {tag}\n"
        f"median removed carries {np.sqrt(np.mean(med ** 2)) / np.sqrt(np.mean(x ** 2)):.1%} "
        f"of the section RMS; per-channel loss {100 * np.median(frac):.1f}% median, "
        f"{100 * frac.max():.1f}% worst", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig09_common_mode_before_after.png", dpi=140)
    plt.close(fig)
    return dict(cm_median_trace=med, cm_frac_removed=frac, cm_z=z)


def fig08_units(deep_raw, fs_d, tag):
    """Deep is stored as strain, Nano as strain rate. What that looks like.

    OptaSense stores strain; Sintela stores strain rate. readFile_HDF only
    differentiates when it is passed diff=True *and* an explicit
    system='OptaSense'. awd_spectra.py, wireline_tube_look.py and
    tube_wave_gate.py reconcile this with np.gradient; the paired-stack and
    record-section paths do not. A factor of omega is a 90 degree phase
    rotation plus a slope across the band, so it changes waveform shape and
    cross-fiber amplitude, but not the position of an F-K ridge or a moveout
    velocity."""
    if deep_raw is None:
        return
    rate = np.gradient(deep_raw, 1.0 / fs_d, axis=-1)
    c = ch_at(250.0, DX_DEEP, deep_raw.shape[0])
    t = np.arange(deep_raw.shape[1]) / fs_d - PRE_S
    i0 = int(PRE_S * fs_d)
    nper = min(1024, deep_raw.shape[1] - i0)

    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    im = image(ax[0, 0], norm_rows(deep_raw), fs_d, DX_DEEP)
    ax[0, 0].set_title("Deep as read: strain")
    plt.colorbar(im, ax=ax[0, 0])
    im = image(ax[0, 1], norm_rows(rate), fs_d, DX_DEEP)
    ax[0, 1].set_title("Deep differentiated: strain rate, matching Nano")
    plt.colorbar(im, ax=ax[0, 1])

    ax[1, 0].plot(t, deep_raw[c] / np.max(np.abs(deep_raw[c])), lw=0.9, label="strain")
    ax[1, 0].plot(t, rate[c] / np.max(np.abs(rate[c])), lw=0.9, ls="--",
                  label="strain rate")
    ax[1, 0].set_xlim(-0.1, 0.6)
    ax[1, 0].set_xlabel("time after drop (s)")
    ax[1, 0].set_ylabel("normalised")
    ax[1, 0].set_title("one Deep trace at 250 m, both conventions")
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].grid(alpha=0.3)

    for arr, lab in [(deep_raw, "strain"), (rate, "strain rate")]:
        f, p = welch(arr[c, i0:], fs=fs_d, nperseg=nper)
        ax[1, 1].loglog(f[1:], np.sqrt(p[1:]), lw=1.2, label=lab)
    ax[1, 1].set_xlabel("frequency (Hz)")
    ax[1, 1].set_ylabel("amplitude spectral density")
    ax[1, 1].set_title(r"the difference is a factor of $\omega$: a slope across the band")
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].grid(alpha=0.3, which="both")

    fig.suptitle(f"8  Deep stores strain, Nano stores strain rate -- {tag}", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig08_units.png", dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------- driver

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--burst", type=int, default=None,
                    help="manifest burst_id; default picks the fullest single-file burst")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_manifest()
    if args.burst is None:
        burst_id, drops = pick_burst(rows)
    else:
        burst_id = args.burst
        drops = [r for r in rows
                 if r["burst_id"] == burst_id and r["nano_available"]]
    drops.sort(key=lambda r: r["utc_time"])
    nano_name = drops[0]["nano_file"]
    tag = f"burst {burst_id}, {len(drops)} drops, {drops[0]['utc_time']:%Y-%m-%d %H:%M} UTC"
    print(f"burst {burst_id}: {len(drops)} drops, nano file {nano_name}")

    nano_path = NANO_DIR / nano_name
    t_nano = nano_time(nano_name)
    print("reading Nano raw ...", flush=True)
    nano_all, info_n = readFile_protobuf([str(nano_path)], fmin=1.0, fmax=250.0,
                                         desampling=False, **RAW_KW)
    fs_n, dx_n = float(info_n["fs"]), float(info_n["dx"])
    print(f"  nano {nano_all.shape}  fs={fs_n}  dx={dx_n}  "
          f"min={nano_all.min():.4g} max={nano_all.max():.4g}", flush=True)

    print("reading Nano with reader defaults for comparison ...", flush=True)
    try:
        nano_def_all, _ = readFile_protobuf([str(nano_path)], fmin=1.0, fmax=250.0,
                                            desampling=False)
    except Exception as exc:                                    # noqa: BLE001
        print("  default read failed:", exc)
        nano_def_all = None

    t_drop = drops[len(drops) // 2]["utc_time"]
    nano_sec = section_around(nano_all, fs_n, t_nano, t_drop)
    nano_def_sec = (section_around(nano_def_all, fs_n, t_nano, t_drop)
                    if nano_def_all is not None else None)

    deep_sec, fs_d = None, None
    deep_name = drops[len(drops) // 2].get("deep_file", "")
    if deep_name:
        try:
            print("reading Deep raw ...", flush=True)
            deep_all, info_d = readFile_HDF([str(DEEP_DIR / deep_name)], fmin=1.0,
                                            fmax=250.0, desampling=False,
                                            verbose=False, **RAW_KW)
            fs_d = float(info_d["fs"])
            print(f"  deep {deep_all.shape}  fs={fs_d}", flush=True)
            deep_sec = section_around(deep_all, fs_d, deep_time(deep_name), t_drop)
            del deep_all
        except Exception as exc:                                # noqa: BLE001
            print("  deep read failed:", exc)

    print("figures ...", flush=True)
    fig01_raw_sections(nano_sec, deep_sec, fs_n, fs_d, tag)
    fig02_bands(nano_sec, fs_n, DX_NANO, "Nano", tag)
    psd = fig03_psd(nano_sec, deep_sec, fs_n, fs_d, tag)
    fig04_spectrogram(nano_all, fs_n, t_nano, drops, tag)
    qc = fig05_channel_qc(nano_all, fs_n, tag)
    peaks, ptimes = fig06_repeatability(nano_all, fs_n, t_nano, drops, tag)
    rms_ladder = fig07_processing_ladder(nano_sec, nano_def_sec, fs_n, tag)
    fig08_units(deep_sec, fs_d, tag)
    cm = fig09_common_mode(nano_sec, fs_n, tag)

    np.savez(OUT_DIR / "plain_look.npz",
             burst_id=burst_id, tag=tag, fs_nano=fs_n, dx_nano=dx_n,
             nano_section=nano_sec.astype(np.float32),
             deep_section=(deep_sec.astype(np.float32) if deep_sec is not None
                           else np.zeros((0, 0), np.float32)),
             drop_peaks=peaks,
             drop_times=np.array([str(t) for t in ptimes]),
             ladder_rms=np.asarray(rms_ladder),
             **{f"qc_{k}": v for k, v in qc.items()},
             **{f"psd_{k}": v for k, v in psd.items()},
             **cm)

    # the numbers worth reading in the log rather than off a plot
    print("\n--- plain numbers ---")
    print(f"nano raw range      {nano_all.min():.6g} .. {nano_all.max():.6g}")
    print(f"channels exactly 0  {int((qc['rms'] <= 0).sum())} / {qc['rms'].size}")
    lo = np.percentile(qc["rms"][qc["rms"] > 0], 1)
    print(f"1st-pct RMS         {lo:.4g}  (median {np.median(qc['rms']):.4g})")
    print(f"|lag1| < 0.05 on    {int((np.abs(qc['lag1']) < 0.05).sum())} channels")
    if peaks.size:
        print(f"drop peak amplitude {peaks.min():.4g} .. {peaks.max():.4g}  "
              f"(median {np.median(peaks):.4g}, CV {peaks.std() / peaks.mean():.2%})")
    print("\nwrote", OUT_DIR)


if __name__ == "__main__":
    main()
