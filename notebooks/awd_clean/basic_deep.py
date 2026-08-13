"""The same standard first-look set as `basic.py`, for the Deep (wireline) fiber.

`basic.py` is Nano only. This is the other half.

Three things differ from the Nano version and none of them are cosmetic.

**Units.** OptaSense (Deep) records strain; Sintela (Nano) records strain rate.
Everything here is differentiated so the two fibers are the same quantity and
the numbers can be compared at all.

**The fiber hairpins at channel 1702** (`deep_target_scan.py:38`): it goes down
the hole and comes back up, so channel index is not monotonic in position. Only
the outbound leg is used, which is the same choice `multiscale_record_sections.py`
makes.

**Band.** The Deep working band in this project is 3-15 Hz, not the Nano 20-50 Hz.
The band figure shows why.

Outputs, in figures/awd_2026/plain_look/basic_deep/
  dbasic01  shot gather, one drop
  dbasic02  shot gather, all drops stacked
  dbasic03  traces down the leg
  dbasic04  amplitude spectrum, signal vs pre-drop noise
  dbasic05  f-k spectrum
  dbasic06  RMS amplitude against distance
  dbasic07  SNR against distance
  dbasic08  relative traveltime by inter-channel cross-correlation
  dbasic09  one drop bandpassed into six bands
  dbasic10  spectrogram of a whole 60 s raw file
  dbasic11  amplitude spectrum against distance
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import spectrogram as _spec

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plain_look import PLAIN_LOOK_DIR  # noqa: E402
from deep_plain_look import (  # noqa: E402
    DX, TURNAROUND_CH, PRE_S, bandpass, cut, deep_time, load_manifest,
    read_deep, to_rate,
)

OUT_DIR = PLAIN_LOOK_DIR / "basic_deep"
HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"

BAND = (3.0, 15.0)          # the Deep working band in this project
V_REF = 1450.0              # the Deep slow-mode candidate, drawn as a guide only
FIG = (9.5, 7.5)


def despike(sec, factor=3.0, half=40):
    rms = np.sqrt(np.mean(np.asarray(sec, float) ** 2, axis=1))
    base = np.array([np.median(rms[max(0, i - half):i + half + 1])
                     for i in range(rms.size)])
    bad = rms > factor * base
    out = np.array(sec, float, copy=True)
    good = np.flatnonzero(~bad)
    for i in np.flatnonzero(bad):
        out[i] = out[good[np.argsort(np.abs(good - i))[:2]]].mean(axis=0)
    return out


def gather(sec, fs, z, title, fname, clip=99.7):
    fig, ax = plt.subplots(figsize=FIG)
    t = np.arange(sec.shape[1]) / fs - PRE_S
    v = np.percentile(np.abs(sec), clip)
    im = ax.pcolormesh(t, z, sec, cmap="seismic", vmin=-v, vmax=v, shading="auto")
    ax.invert_yaxis()
    ax.set_xlim(-0.2, 2.5)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("distance along outbound leg (m)")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, label=r"strain rate ($\mu\varepsilon$ s$^{-1}$)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / fname, dpi=150)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    d = np.load(STACKS)
    fs = float(d["fs"])
    counts = d["n_common"]
    g = counts > 0
    w = counts[g].astype(float)
    ds = d["deep_stacks"]
    start = int(d["deep_start_ch"])
    stack_strain = np.tensordot(w, ds[g], axes=(0, 0)) / w.sum()
    n_drops = int(w.sum())
    print(f"deep stacks {ds.shape}, start channel {start}, {n_drops} drops",
          flush=True)

    # outbound leg only, in strain rate
    n_out = max(0, min(stack_strain.shape[0], TURNAROUND_CH - start))
    if n_out < 50:
        raise SystemExit(f"only {n_out} outbound channels saved; nothing to plot")
    S = despike(bandpass(to_rate(stack_strain[:n_out], fs), fs, BAND))
    z = (start + np.arange(n_out)) * DX
    i0 = int(PRE_S * fs)
    print(f"  outbound leg: {n_out} channels, {z[0]:.0f}-{z[-1]:.0f} m", flush=True)

    # ---- 01 one drop -----------------------------------------------------
    rows = load_manifest()
    ok = [r for r in rows
          if r["off"] - PRE_S >= 0 and r["off"] + 3.0 <= 60.0 and r["w"] > 0.99]
    ref = sorted(ok, key=lambda r: r["t"])[len(ok) // 2]
    raw_all, fs2 = read_deep(ref["file"])
    one_strain = cut(raw_all, fs2, deep_time(ref["file"]), ref["t"])
    one = to_rate(one_strain[start:start + n_out], fs2)
    gather(despike(bandpass(one, fs2, BAND)), fs2, z,
           f"Deep shot gather, one drop (burst {ref['burst']}, "
           f"{BAND[0]:g}-{BAND[1]:g} Hz)", "dbasic01_gather_one_drop.png")

    # ---- 02 stack --------------------------------------------------------
    gather(S, fs, z, f"Deep shot gather, {n_drops} drops stacked "
                     f"({BAND[0]:g}-{BAND[1]:g} Hz)", "dbasic02_gather_stack.png")

    # ---- 03 traces -------------------------------------------------------
    fig, ax = plt.subplots(figsize=FIG)
    t = np.arange(S.shape[1]) / fs - PRE_S
    win = (t >= -0.2) & (t <= 2.5)
    want = np.linspace(z[0], z[-1], 12)
    idx = [int(np.argmin(np.abs(z - q))) for q in want]
    span = (z[-1] - z[0]) / 14
    scale = span / max(np.max(np.abs(S[i][win])) for i in idx)
    for i in idx:
        ax.plot(t[win], scale * S[i][win] + z[i], "k", lw=0.9)
    ax.invert_yaxis()
    ax.set_xlim(-0.2, 2.5)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("distance along outbound leg (m)")
    ax.set_title(f"Deep stacked traces ({BAND[0]:g}-{BAND[1]:g} Hz)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "dbasic03_traces.png", dpi=150)
    plt.close(fig)

    # ---- 04 spectrum -----------------------------------------------------
    raw = despike(to_rate(stack_strain[:n_out], fs))
    sig, noi = raw[:, i0:], raw[:, :i0]
    n = min(sig.shape[1], noi.shape[1])
    f = np.fft.rfftfreq(n, 1 / fs)
    A_s = np.mean(np.abs(np.fft.rfft(sig[:, :n], axis=1)), axis=0)
    A_n = np.mean(np.abs(np.fft.rfft(noi[:, :n], axis=1)), axis=0)
    fig, ax = plt.subplots(figsize=FIG)
    ax.loglog(f[1:], A_s[1:], lw=1.4, label="signal window")
    ax.loglog(f[1:], A_n[1:], lw=1.4, label="pre-drop noise")
    ax.axvspan(*BAND, color="0.85", zorder=0, label=f"{BAND[0]:g}-{BAND[1]:g} Hz")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("mean amplitude spectrum")
    ax.set_title("Deep amplitude spectrum, averaged over the outbound leg")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "dbasic04_spectrum.png", dpi=150)
    plt.close(fig)

    # ---- 05 f-k ----------------------------------------------------------
    seg = raw[:, i0:i0 + int(2.0 * fs)]
    F = np.fft.fftshift(np.fft.fft2(seg), axes=0)
    fq = np.fft.rfftfreq(seg.shape[1], 1 / fs)
    k = np.fft.fftshift(np.fft.fftfreq(seg.shape[0], DX))
    P = 20 * np.log10(np.abs(F[:, :fq.size]) + 1e-12)
    fig, ax = plt.subplots(figsize=FIG)
    im = ax.pcolormesh(fq, k, P, cmap="magma", shading="auto",
                       vmin=np.percentile(P, 60), vmax=np.percentile(P, 99.9))
    for vv, ls in [(V_REF, "-"), (2975.0, "--")]:
        ax.plot(fq, fq / vv, "c" + ls, lw=1.2, label=f"{vv:.0f} m/s")
        ax.plot(fq, -fq / vv, "c" + ls, lw=1.2)
    ax.set_xlim(0, 60)
    ax.set_ylim(-0.05, 0.05)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel(r"wavenumber (m$^{-1}$)")
    ax.set_title("Deep f-k spectrum, outbound leg, 2 s window")
    ax.legend(loc="upper right", fontsize=9)
    plt.colorbar(im, ax=ax, label="dB")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "dbasic05_fk.png", dpi=150)
    plt.close(fig)

    # ---- 06/07 amplitude and SNR ------------------------------------------
    noise_rms = np.sqrt(np.mean(S[:, :i0] ** 2, axis=1))
    sig_rms = np.sqrt(np.mean(S[:, i0:i0 + int(1.5 * fs)] ** 2, axis=1))
    snr_db = 20 * np.log10(np.maximum(sig_rms, 1e-30) / np.maximum(noise_rms, 1e-30))

    fig, ax = plt.subplots(figsize=FIG)
    ax.semilogy(z, sig_rms, lw=1.4, label="signal window (0-1.5 s)")
    ax.semilogy(z, noise_rms, lw=1.4, label="pre-drop noise")
    ax.set_xlabel("distance along outbound leg (m)")
    ax.set_ylabel(r"RMS strain rate ($\mu\varepsilon$ s$^{-1}$)")
    ax.set_title(f"Deep amplitude against distance, {n_drops}-drop stack")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "dbasic06_amplitude_vs_distance.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=FIG)
    ax.plot(z, snr_db, lw=1.4)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(10, color="0.5", lw=0.8, ls="--")
    ax.set_xlabel("distance along outbound leg (m)")
    ax.set_ylabel("SNR (dB)")
    ax.set_title(f"Deep SNR against distance, {n_drops}-drop stack")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "dbasic07_snr_vs_distance.png", dpi=150)
    plt.close(fig)

    # ---- 08 moveout by inter-channel cross-correlation ---------------------
    step = 4
    v_min = 800.0
    half = max(3, int(np.ceil(step * DX / v_min * fs)))
    s0 = int(np.argmax(snr_db))
    ref_c, tt = s0, {s0: 0.0}
    for c in range(s0 + step, n_out, step):
        cc = np.correlate(S[c] - S[c].mean(), S[ref_c] - S[ref_c].mean(), mode="same")
        mid = cc.size // 2
        w2 = cc[mid:mid + half]
        if w2.size < 3 or not np.any(w2):
            break
        tt[c] = tt[ref_c] + int(np.argmax(w2)) / fs
        ref_c = c
    cs = np.array(sorted(tt))
    zz, ts = z[cs], np.array([tt[c] for c in cs])
    keep = snr_db[cs] > 6.0
    fig, ax = plt.subplots(figsize=FIG)
    ax.plot(ts[keep] * 1e3, zz[keep], "k.", ms=5, label="accumulated lag")
    if keep.sum() > 10:
        A = np.vstack([ts[keep], np.ones(int(keep.sum()))]).T
        vel, icpt = np.linalg.lstsq(A, zz[keep], rcond=None)[0]
        tl = np.array([ts[keep].min(), ts[keep].max()])
        ax.plot(tl * 1e3, vel * tl + icpt, "r-", lw=1.8, label=f"fit: {vel:.0f} m/s")
        print(f"  deep moveout fit: {vel:.0f} m/s over {int(keep.sum())} points "
              f"({zz[keep].min():.0f}-{zz[keep].max():.0f} m)")
    ax.invert_yaxis()
    ax.set_xlabel("relative traveltime (ms)")
    ax.set_ylabel("distance along outbound leg (m)")
    ax.set_title("Deep relative traveltime by inter-channel cross-correlation")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "dbasic08_moveout.png", dpi=150)
    plt.close(fig)

    # ---- 09 bands ----------------------------------------------------------
    bands = [None, (1.0, 5.0), (3.0, 15.0), (5.0, 20.0), (20.0, 50.0), (50.0, 100.0)]
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True)
    t1 = np.arange(one.shape[1]) / fs2 - PRE_S
    for a, bd in zip(axes.ravel(), bands):
        x = despike(bandpass(one, fs2, bd) if bd else one)
        v = np.percentile(np.abs(x), 99.7)
        im = a.pcolormesh(t1, z, x, cmap="seismic", vmin=-v, vmax=v, shading="auto")
        a.set_title("unfiltered" if bd is None else f"{bd[0]:g}-{bd[1]:g} Hz")
        plt.colorbar(im, ax=a, label=r"$\mu\varepsilon$ s$^{-1}$")
    axes[0, 0].invert_yaxis()
    axes[0, 0].set_xlim(-0.2, 2.5)
    for a in axes[1]:
        a.set_xlabel("time (s)")
    for a in axes[:, 0]:
        a.set_ylabel("distance along leg (m)")
    fig.suptitle(f"Deep, one drop bandpassed into six bands (burst {ref['burst']})")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "dbasic09_bands.png", dpi=150)
    plt.close(fig)

    # ---- 10 spectrogram ----------------------------------------------------
    c_spec = start + n_out // 4
    ff, ts_, Sxx = _spec(to_rate(raw_all[c_spec], fs2), fs=fs2,
                         nperseg=2048, noverlap=1024)
    db = 10 * np.log10(Sxx + 1e-30)
    fig, ax = plt.subplots(figsize=(12, 6.5))
    im = ax.pcolormesh(ts_, ff, db, cmap="magma", shading="auto",
                       vmin=np.percentile(db, 5), vmax=np.percentile(db, 99.5))
    ax.axvline(ref["off"], color="c", lw=1.0)
    ax.set_ylim(0, 200)
    ax.set_xlabel("seconds into the 60 s raw file")
    ax.set_ylabel("frequency (Hz)")
    ax.set_title(f"Deep spectrogram at {c_spec * DX:.0f} m (cyan = the drop)")
    plt.colorbar(im, ax=ax, label="dB")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "dbasic10_spectrogram.png", dpi=150)
    plt.close(fig)
    del raw_all

    # ---- 11 spectrum vs distance -------------------------------------------
    A = np.abs(np.fft.rfft(sig, axis=1))
    fq2 = np.fft.rfftfreq(sig.shape[1], 1 / fs)
    m = (fq2 >= 1) & (fq2 <= 120)
    AdB = 20 * np.log10(A[:, m] + 1e-30)
    fig, ax = plt.subplots(figsize=(10, 7.5))
    im = ax.pcolormesh(fq2[m], z, AdB, cmap="magma", shading="auto",
                       vmin=np.percentile(AdB, 5), vmax=np.percentile(AdB, 99.5))
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("distance along outbound leg (m)")
    ax.set_title(f"Deep amplitude spectrum against distance, {n_drops}-drop stack")
    plt.colorbar(im, ax=ax, label="dB")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "dbasic11_spectrum_vs_distance.png", dpi=150)
    plt.close(fig)

    print(f"  SNR range {snr_db.min():.1f} to {snr_db.max():.1f} dB")
    print("wrote", OUT_DIR)


if __name__ == "__main__":
    main()
