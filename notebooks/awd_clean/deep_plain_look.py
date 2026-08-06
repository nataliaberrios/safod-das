"""Plain looks at the raw Deep (wireline) fiber, and a direct test of the taper ramp.

The Nano counterpart is `plain_look.py`. This is the Deep half, which was the
gap: figures 2 and 4-7 and 9 and 11 there are cemented-fiber only, so every
statement about bad channels, band content, repeatability and common mode so far
applies to Nano and not to this fiber.

Deep needs its own script rather than a loop over both, for three reasons.

* **The fiber hairpins at channel 1702** (`deep_target_scan.py:38`). Channel
  index is not monotonic in depth: channels 0-1702 go down the hole and
  1702-3200 come back up. Plotting `channel * dx` straight through, as
  `plain_look.py` figures 1 and 3 do, mislabels everything past the turnaround.
  Here the two legs are drawn separately and then folded onto a common
  coordinate.
* **Files are 62 s, not 300 s.** The read-time Tukey taper therefore bites much
  harder: 3.5 s of cut window is 5.6% of the file, so a drop inside the ramp
  gets a strongly time-varying gain rather than a scale factor.
* **Deep is the fiber with the open claims on it** -- the slow modes, the
  time-gated branch search, the conditional registration -- and they live at
  late times, which is exactly what a decaying ramp suppresses.

Figure 6 is the point of the script: two Deep drops, one with a flat window and
one with a steep ramp, drawn against each other. If the late-time energy tracks
the taper rather than the Earth, it shows up there.

Reader settings: `filter=False, median=False, detrend=False, tapering=False`.

Outputs
-------
figures/awd_2026/plain_look/deep_fig01..deep_fig06*.png
"""
from __future__ import annotations

import csv
import re
import sys
from datetime import datetime, timezone
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
from DASutils import readFile_HDF  # noqa: E402

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "awd_manifest.csv"
DEEP_DIR = Path(
    "/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/"
    "01_--_recording_2026-06-15T230629Z_--_active_source"
)
OUT_DIR = Path(
    "/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026/plain_look"
)

DX = 2.0419                 # Deep channel spacing
TURNAROUND_CH = 1702        # deep_target_scan.py:38
PRE_S, POST_S = 0.5, 3.0
DEEP_DURATION_S = 62.0
ALPHA = 0.4
FS_NOMINAL = 1000.0
RAW_KW = dict(filter=False, median=False, detrend=False, tapering=False)

# Deep working bands in this repo are 3-15, 5-20, 15-30 and 20-50 Hz.
BANDS = [None, (1.0, 5.0), (3.0, 15.0), (5.0, 20.0), (20.0, 50.0), (50.0, 100.0)]


def deep_time(name: str):
    m = re.search(r"_(\d{4}-\d{2}-\d{2}T\d{6}Z)\.h5$", name)
    return (datetime.strptime(m.group(1), "%Y-%m-%dT%H%M%SZ").replace(tzinfo=timezone.utc)
            if m else None)


def bandpass(x, fs, band):
    if band is None:
        return x
    sos = butter(4, [band[0], min(band[1], 0.45 * fs)], btype="bandpass",
                 fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def taper_weight(offset_s, dur=DEEP_DURATION_S, fs=FS_NOMINAL, alpha=ALPHA):
    from scipy.signal.windows import tukey
    n = int(round(dur * fs))
    w = tukey(n, alpha=alpha)
    return w[int(np.clip(offset_s * fs, 0, n - 1))]


def load_manifest():
    rows = []
    with open(MANIFEST) as fh:
        for r in csv.DictReader(fh):
            if r["deep_available"] != "1":
                continue
            off = float(r["deep_offset_s"])
            w0 = taper_weight(max(off - PRE_S, 0.0))
            w1 = taper_weight(min(off + POST_S, DEEP_DURATION_S))
            rows.append(dict(
                burst=int(r["burst_id"]),
                t=datetime.fromisoformat(r["utc_time"]),
                file=r["deep_file"], off=off,
                w=taper_weight(off),
                ramp=(w1 / w0 if w0 > 1e-6 else np.inf),
            ))
    return rows


def read_deep(name):
    d, info = readFile_HDF([str(DEEP_DIR / name)], fmin=1.0, fmax=250.0,
                           desampling=False, verbose=False, **RAW_KW)
    return d, float(info["fs"])


def cut(data, fs, t_file, t_drop, pre=PRE_S, post=POST_S):
    i = int((t_drop - t_file).total_seconds() * fs)
    i0, i1 = i - int(pre * fs), i + int(post * fs)
    if i0 < 0 or i1 > data.shape[1]:
        return None
    return data[:, i0:i1]


def norm_rows(sec):
    p = np.max(np.abs(sec), axis=1, keepdims=True)
    return sec / np.where(p > 0, p, 1.0)


def legs(sec):
    """Split at the hairpin. Return leg is reversed so both run away from surface."""
    out = sec[:TURNAROUND_CH]
    ret = sec[TURNAROUND_CH:][::-1]
    return out, ret


def image(ax, sec, fs, coord, clip=99.0):
    t = np.arange(sec.shape[1]) / fs - PRE_S
    v = np.percentile(np.abs(sec), clip)
    im = ax.pcolormesh(t, coord, sec, cmap="seismic", vmin=-v, vmax=v, shading="auto")
    ax.invert_yaxis()
    ax.set_xlabel("time after drop (s)")
    ax.set_ylabel("distance along leg (m)")
    return im


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_manifest()

    # a clean reference drop: flat window, from the fullest burst available
    clean = [r for r in rows if abs(r["ramp"] - 1.0) < 0.02 and r["w"] > 0.99]
    ref = sorted(clean, key=lambda r: r["t"])[len(clean) // 2]
    print(f"reference drop: burst {ref['burst']} at {ref['t']}, "
          f"offset {ref['off']:.1f} s, weight {ref['w']:.3f}, ramp {ref['ramp']:.3f}")

    d_all, fs = read_deep(ref["file"])
    print(f"  deep {d_all.shape} fs={fs}", flush=True)
    sec = cut(d_all, fs, deep_time(ref["file"]), ref["t"])
    tag = (f"burst {ref['burst']}, {ref['t']:%Y-%m-%d %H:%M} UTC, "
           f"taper weight {ref['w']:.2f}, window ramp {ref['ramp']:.2f}")

    out, ret = legs(sec)
    z_out = np.arange(out.shape[0]) * DX
    z_ret = np.arange(ret.shape[0]) * DX

    # ---------------------------------------------------------------- fig 1
    fig, ax = plt.subplots(1, 3, figsize=(17, 6))
    im = image(ax[0], norm_rows(out), fs, z_out)
    ax[0].set_title(f"outbound leg, ch 0-{TURNAROUND_CH} (down the hole)")
    plt.colorbar(im, ax=ax[0])
    im = image(ax[1], norm_rows(ret), fs, z_ret)
    ax[1].set_title(f"return leg, ch {TURNAROUND_CH}-{sec.shape[0]} (reversed)")
    plt.colorbar(im, ax=ax[1])
    n = min(out.shape[0], ret.shape[0])
    im = image(ax[2], norm_rows(out[:n] + ret[:n]), fs, z_out[:n])
    ax[2].set_title("legs folded and summed")
    plt.colorbar(im, ax=ax[2])
    fig.suptitle(f"Deep 1  Raw, unfiltered, hairpin split at ch {TURNAROUND_CH} -- {tag}",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "deep_fig01_raw_legs.png", dpi=140)
    plt.close(fig)

    # ---------------------------------------------------------------- fig 2
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for a, band in zip(axes.ravel(), BANDS):
        im = image(a, norm_rows(bandpass(out, fs, band)), fs, z_out)
        a.set_title("unfiltered" if band is None else f"{band[0]:g}-{band[1]:g} Hz")
        plt.colorbar(im, ax=a)
    fig.suptitle(f"Deep 2  Outbound leg band by band -- {tag}", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "deep_fig02_bands.png", dpi=140)
    plt.close(fig)

    # ---------------------------------------------------------------- fig 3
    x = d_all.astype(np.float64)
    rms = np.sqrt(np.mean(x ** 2, axis=1))
    dc = np.mean(x, axis=1)
    xc = x - dc[:, None]
    den = np.sum(xc * xc, axis=1)
    lag1 = np.sum(xc[:, 1:] * xc[:, :-1], axis=1) / np.where(den > 0, den, 1.0)
    ch = np.arange(x.shape[0])

    fig, ax = plt.subplots(2, 2, figsize=(15, 9))
    ax[0, 0].semilogy(ch, np.where(rms > 0, rms, np.nan), lw=0.6)
    ax[0, 0].axvline(TURNAROUND_CH, color="C3", lw=1.2, label="hairpin")
    ax[0, 0].set_xlabel("channel")
    ax[0, 0].set_ylabel("RMS")
    ax[0, 0].set_title(f"RMS per channel ({int((rms <= 0).sum())} exactly zero)")
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(alpha=0.3)

    ax[0, 1].plot(ch, dc, lw=0.6)
    ax[0, 1].axvline(TURNAROUND_CH, color="C3", lw=1.2)
    ax[0, 1].set_xlabel("channel")
    ax[0, 1].set_ylabel("mean")
    ax[0, 1].set_title("DC offset per channel")
    ax[0, 1].grid(alpha=0.3)

    ax[1, 0].plot(ch, lag1, lw=0.6)
    ax[1, 0].axvline(TURNAROUND_CH, color="C3", lw=1.2)
    ax[1, 0].axhline(0, color="k", lw=0.6)
    ax[1, 0].set_xlabel("channel")
    ax[1, 0].set_ylabel("lag-1 autocorrelation")
    ax[1, 0].set_title("Whiteness: near 0 means noise-only")
    ax[1, 0].grid(alpha=0.3)

    # do the two legs see the same thing at the same depth?
    n = min(TURNAROUND_CH, x.shape[0] - TURNAROUND_CH)
    ax[1, 1].semilogy(z_out[:n], rms[:n], lw=0.7, label="outbound")
    ax[1, 1].semilogy(z_out[:n], rms[TURNAROUND_CH:][::-1][:n], lw=0.7, label="return")
    ax[1, 1].set_xlabel("distance along leg (m)")
    ax[1, 1].set_ylabel("RMS")
    ax[1, 1].set_title("Same fiber, both legs: should agree if depth-mapped")
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].grid(alpha=0.3)

    fig.suptitle(f"Deep 3  Per-channel QC over one whole 62 s file -- {tag}", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "deep_fig03_channel_qc.png", dpi=140)
    plt.close(fig)

    # ---------------------------------------------------------------- fig 4
    fig, ax = plt.subplots(2, 2, figsize=(15, 9))
    i0 = int(PRE_S * fs)
    nper = min(1024, sec.shape[1] - i0)
    for zq in [100.0, 500.0, 1500.0, 3000.0]:
        c = min(int(zq / DX), out.shape[0] - 1)
        f, p = welch(out[c, i0:], fs=fs, nperseg=nper)
        ax[0, 0].loglog(f[1:], np.sqrt(p[1:]), lw=1.1, label=f"{zq:.0f} m")
    for line in (60.0, 120.0, 180.0):
        ax[0, 0].axvline(line, color="0.7", lw=0.8, ls=":", zorder=0)
    ax[0, 0].set_xlabel("frequency (Hz)")
    ax[0, 0].set_ylabel("amplitude spectral density")
    ax[0, 0].set_title("Signal-window spectra (dotted = 60 Hz harmonics)")
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(alpha=0.3, which="both")

    f, pall = welch(out[:, :i0], fs=fs, nperseg=nper, axis=-1)
    db = 10 * np.log10(pall[:, 1:] + 1e-30)
    im = ax[0, 1].pcolormesh(f[1:], z_out, db, cmap="magma", shading="auto",
                             vmin=np.percentile(db, 2), vmax=np.percentile(db, 98))
    ax[0, 1].set_xscale("log")
    ax[0, 1].invert_yaxis()
    ax[0, 1].set_xlabel("frequency (Hz)")
    ax[0, 1].set_ylabel("distance along leg (m)")
    ax[0, 1].set_title("Pre-drop noise PSD by channel, outbound leg")
    plt.colorbar(im, ax=ax[0, 1], label="dB")

    c = min(int(500.0 / DX), x.shape[0] - 1)
    ff, tt, S = spectrogram(x[c], fs=fs, nperseg=2048, noverlap=1024)
    sdb = 10 * np.log10(S + 1e-30)
    im = ax[1, 0].pcolormesh(tt, ff, sdb, cmap="magma", shading="auto",
                             vmin=np.percentile(sdb, 5), vmax=np.percentile(sdb, 99.5))
    ax[1, 0].axvline(ref["off"], color="c", lw=1.0)
    for frac in (0.2, 0.8):
        ax[1, 0].axvline(frac * DEEP_DURATION_S, color="w", ls="--", lw=1.0)
    ax[1, 0].set_xlabel("seconds into the 62 s file")
    ax[1, 0].set_ylabel("frequency (Hz)")
    ax[1, 0].set_title("Whole file at 500 m; cyan = drop, white = taper flat zone")
    plt.colorbar(im, ax=ax[1, 0], label="dB")

    # common mode, before and after, on this fiber
    y = bandpass(out, fs, (5.0, 20.0))
    med = np.median(y, axis=0)
    e0 = np.sqrt(np.mean(y ** 2, axis=1))
    e1 = np.sqrt(np.mean((y - med[None, :]) ** 2, axis=1))
    frac = 1.0 - np.divide(e1, e0, out=np.zeros_like(e1), where=e0 > 0)
    ax[1, 1].plot(z_out, 100 * frac, lw=0.7)
    ax[1, 1].axhline(0, color="k", lw=0.6)
    ax[1, 1].set_xlabel("distance along leg (m)")
    ax[1, 1].set_ylabel("% of channel RMS removed")
    ax[1, 1].set_title(f"Common-mode removal at 5-20 Hz: median "
                       f"{100 * np.median(frac):.1f}%, worst {100 * frac.max():.1f}%")
    ax[1, 1].grid(alpha=0.3)

    fig.suptitle(f"Deep 4  Spectra, spectrogram, common mode -- {tag}", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "deep_fig04_spectra_commonmode.png", dpi=140)
    plt.close(fig)

    # ---------------------------------------------------------------- fig 5
    burst_drops = sorted([r for r in rows if r["burst"] == ref["burst"]],
                         key=lambda r: r["t"])
    same_file = [r for r in burst_drops if r["file"] == ref["file"]]
    fig, ax = plt.subplots(3, 2, figsize=(14, 10))
    for row, zq in enumerate([100.0, 500.0, 1500.0]):
        c = min(int(zq / DX), out.shape[0] - 1)
        for col, band in enumerate([None, (5.0, 20.0)]):
            a = ax[row, col]
            used = 0
            for r in same_file:
                s = cut(d_all, fs, deep_time(r["file"]), r["t"], pre=0.1, post=0.6)
                if s is None:
                    continue
                tr = bandpass(s[c:c + 1], fs, band)[0]
                a.plot(np.arange(tr.size) / fs - 0.1, tr, lw=0.6, alpha=0.55)
                used += 1
            a.set_xlabel("time after drop (s)")
            a.set_ylabel("microstrain")
            a.set_title(f"{zq:.0f} m, "
                        f"{'unfiltered' if band is None else '5-20 Hz'} -- {used} drops")
            a.grid(alpha=0.3)
    fig.suptitle(f"Deep 5  Drop-to-drop repeatability, every trace drawn -- {tag}",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "deep_fig05_repeatability.png", dpi=140)
    plt.close(fig)
    del d_all

    # ---------------------------------------------------------------- fig 6
    # the flagged risk, tested directly: flat window vs steep ramp
    ramped = [r for r in rows if np.isfinite(r["ramp"]) and r["ramp"] < 0.5]
    if not ramped:
        print("no strongly ramped Deep drop found; skipping figure 6")
        return
    bad = sorted(ramped, key=lambda r: r["ramp"])[len(ramped) // 2]
    print(f"ramped drop: burst {bad['burst']} at {bad['t']}, offset {bad['off']:.1f} s, "
          f"weight {bad['w']:.3f}, ramp {bad['ramp']:.3f}")

    d_bad, fs2 = read_deep(bad["file"])
    sec_bad = cut(d_bad, fs2, deep_time(bad["file"]), bad["t"])
    out_bad, _ = legs(sec_bad)

    fig, ax = plt.subplots(2, 2, figsize=(15, 9))
    for col, (s, r, name) in enumerate([
        (out, ref, "flat window"), (out_bad, bad, "ramped window")
    ]):
        im = image(ax[0, col], norm_rows(bandpass(s, fs, (5.0, 20.0))), fs, z_out)
        ax[0, col].set_title(f"{name}: burst {r['burst']}, weight {r['w']:.2f}, "
                             f"ramp {r['ramp']:.2f}", fontsize=10)
        plt.colorbar(im, ax=ax[0, col])

    # energy against time: what the taper does to the late record
    for s, r, lab in [(out, ref, "flat"), (out_bad, bad, "ramped")]:
        y = bandpass(s, fs, (5.0, 20.0))
        env = np.sqrt(np.mean(y ** 2, axis=0))
        t = np.arange(env.size) / fs - PRE_S
        ax[1, 0].semilogy(t, env / env[t < -0.1].mean(), lw=1.0,
                          label=f"{lab} (ramp {r['ramp']:.2f})")
    ax[1, 0].set_xlabel("time after drop (s)")
    ax[1, 0].set_ylabel("array RMS / pre-drop level")
    ax[1, 0].set_title("Array energy vs time: a ramp tilts this curve")
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].grid(alpha=0.3, which="both")

    # the taper the reader would have applied, drawn over the same window
    from scipy.signal.windows import tukey
    n = int(round(DEEP_DURATION_S * FS_NOMINAL))
    w = tukey(n, alpha=ALPHA)
    for r, lab in [(ref, "flat"), (bad, "ramped")]:
        i = int(r["off"] * FS_NOMINAL)
        lo, hi = i - int(PRE_S * FS_NOMINAL), i + int(POST_S * FS_NOMINAL)
        seg = w[max(lo, 0):min(hi, n)]
        ax[1, 1].plot(np.linspace(-PRE_S, POST_S, seg.size), seg / seg[0], lw=1.4,
                      label=f"{lab} (x{r['ramp']:.2f} across window)")
    ax[1, 1].axhline(1.0, color="k", lw=0.8, ls="--")
    ax[1, 1].set_xlabel("time after drop (s)")
    ax[1, 1].set_ylabel("reader gain, normalised to window start")
    ax[1, 1].set_title("The gain the default reader applies across the cut window")
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].grid(alpha=0.3)

    fig.suptitle("Deep 6  Does the taper ramp change what the record looks like?",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "deep_fig06_taper_ramp_test.png", dpi=140)
    plt.close(fig)
    print("\nwrote", OUT_DIR)


if __name__ == "__main__":
    main()
