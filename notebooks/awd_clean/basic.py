"""The standard first-look figures for a borehole active-source DAS dataset.

No annotations, no caption boxes, no arrows. Just the plots anyone handed this
kind of data would make first, one per file, properly labelled.

  basic01  shot gather, one drop
  basic02  shot gather, all 859 drops stacked
  basic03  a few traces from the stack
  basic04  amplitude spectrum, signal window vs pre-drop noise
  basic05  f-k spectrum of the stack
  basic06  RMS amplitude against distance, signal and noise
  basic07  SNR against distance
  basic08  relative traveltime by inter-channel cross-correlation
  basic09  the same drop bandpassed into six bands
  basic10  spectrogram of one whole raw file at one channel
  basic11  amplitude spectrum against distance
  basic12  drops per burst across the survey

Sources: `canonical_epoch_stacks_paired_deep_all.npz` for anything stacked, and
one raw Nano file for the single drop.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plain_look import (  # noqa: E402
    DX_NANO, NANO_DIR, PLAIN_LOOK_DIR, PRE_S, RAW_KW,
    bandpass, nano_time, pick_burst, read_manifest, section_around,
)
from DASutils import readFile_protobuf  # noqa: E402

OUT_DIR = PLAIN_LOOK_DIR / "basic"

HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
BAND = (20.0, 50.0)
ZMAX = 600.0
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


def gather(sec, fs, title, fname, clip=99.7):
    fig, ax = plt.subplots(figsize=FIG)
    t = np.arange(sec.shape[1]) / fs - PRE_S
    z = np.arange(sec.shape[0]) * DX_NANO
    v = np.percentile(np.abs(sec), clip)
    im = ax.pcolormesh(t, z, sec, cmap="seismic", vmin=-v, vmax=v, shading="auto")
    ax.set_ylim(ZMAX, 0)
    ax.set_xlim(-0.1, 0.8)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("distance along fiber (m)")
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
    good = counts > 0
    w = counts[good].astype(float)
    stack = np.tensordot(w, d["nano_stacks"][good], axes=(0, 0)) / w.sum()
    n_drops = int(w.sum())
    nch = min(int(ZMAX / DX_NANO), stack.shape[0])
    z = np.arange(nch) * DX_NANO
    S = despike(bandpass(stack[:nch], fs, BAND))
    i0 = int(PRE_S * fs)
    print(f"stack {stack.shape}, {n_drops} drops, fs={fs}", flush=True)

    # ---- 01 single drop -----------------------------------------------
    rows = read_manifest()
    burst_id, drops = pick_burst(rows)
    drops.sort(key=lambda r: r["utc_time"])
    name = drops[len(drops) // 2]["nano_file"]
    nano, info = readFile_protobuf([str(NANO_DIR / name)], fmin=1.0, fmax=250.0,
                                   desampling=False, **RAW_KW)
    one = section_around(nano, fs, nano_time(name),
                         drops[len(drops) // 2]["utc_time"])[:nch]
    gather(despike(bandpass(one, fs, BAND)), fs,
           f"Shot gather, one drop (burst {burst_id}, {BAND[0]:g}-{BAND[1]:g} Hz)",
           "basic01_gather_one_drop.png")

    # ---- 02 stacked gather --------------------------------------------
    gather(S, fs, f"Shot gather, {n_drops} drops stacked "
                  f"({BAND[0]:g}-{BAND[1]:g} Hz)", "basic02_gather_stack.png")

    # ---- 03 traces ------------------------------------------------------
    fig, ax = plt.subplots(figsize=FIG)
    t = np.arange(S.shape[1]) / fs - PRE_S
    win = (t >= -0.1) & (t <= 0.5)
    depths = np.arange(50, 601, 50)
    idx = [min(int(zq / DX_NANO), nch - 1) for zq in depths]
    scale = 42.0 / max(np.max(np.abs(S[i][win])) for i in idx)
    for zq, i in zip(depths, idx):
        ax.plot(t[win], scale * S[i][win] + zq, "k", lw=1.0)
    ax.set_ylim(650, 0)
    ax.set_xlim(-0.1, 0.5)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("distance along fiber (m)")
    ax.set_title(f"Stacked traces every 50 m ({BAND[0]:g}-{BAND[1]:g} Hz)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "basic03_traces.png", dpi=150)
    plt.close(fig)

    # ---- 04 amplitude spectrum -----------------------------------------
    raw = despike(stack[:nch])
    sig, noi = raw[:, i0:], raw[:, :i0]
    n = min(sig.shape[1], noi.shape[1])
    f = np.fft.rfftfreq(n, 1 / fs)
    A_s = np.mean(np.abs(np.fft.rfft(sig[:, :n], axis=1)), axis=0)
    A_n = np.mean(np.abs(np.fft.rfft(noi[:, :n], axis=1)), axis=0)
    fig, ax = plt.subplots(figsize=FIG)
    ax.loglog(f[1:], A_s[1:], lw=1.4, label="signal window")
    ax.loglog(f[1:], A_n[1:], lw=1.4, label="pre-drop noise")
    ax.axvspan(BAND[0], BAND[1], color="0.85", zorder=0, label="20-50 Hz")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("mean amplitude spectrum")
    ax.set_title("Amplitude spectrum, averaged over channels 0-600 m")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "basic04_spectrum.png", dpi=150)
    plt.close(fig)

    # ---- 05 f-k ---------------------------------------------------------
    seg = raw[:, i0:i0 + int(1.0 * fs)]
    F = np.fft.fftshift(np.fft.fft2(seg), axes=0)
    freq = np.fft.rfftfreq(seg.shape[1], 1 / fs)
    k = np.fft.fftshift(np.fft.fftfreq(seg.shape[0], DX_NANO))
    P = 20 * np.log10(np.abs(F[:, :freq.size]) + 1e-12)
    fig, ax = plt.subplots(figsize=FIG)
    im = ax.pcolormesh(freq, k, P, cmap="magma", shading="auto",
                       vmin=np.percentile(P, 60), vmax=np.percentile(P, 99.9))
    for v_line, ls in [(2975.0, "-"), (1500.0, "--")]:
        ax.plot(freq, freq / v_line, "c" + ls, lw=1.2,
                label=f"{v_line:.0f} m/s")
        ax.plot(freq, -freq / v_line, "c" + ls, lw=1.2)
    ax.set_xlim(0, 120)
    ax.set_ylim(-0.12, 0.12)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel(r"wavenumber (m$^{-1}$)")
    ax.set_title("f-k spectrum of the stacked gather, 1 s window")
    ax.legend(loc="upper right", fontsize=9)
    plt.colorbar(im, ax=ax, label="dB")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "basic05_fk.png", dpi=150)
    plt.close(fig)

    # ---- 06/07 amplitude and SNR vs distance ----------------------------
    env = np.abs(hilbert(S, axis=1))
    noise_rms = np.sqrt(np.mean(S[:, :i0] ** 2, axis=1))
    sig_rms = np.sqrt(np.mean(S[:, i0:i0 + int(0.25 * fs)] ** 2, axis=1))
    snr_db = 20 * np.log10(np.maximum(sig_rms, 1e-30) /
                           np.maximum(noise_rms, 1e-30))

    fig, ax = plt.subplots(figsize=FIG)
    ax.semilogy(z, sig_rms, lw=1.4, label="signal window (0-0.25 s)")
    ax.semilogy(z, noise_rms, lw=1.4, label="pre-drop noise")
    ax.set_xlabel("distance along fiber (m)")
    ax.set_ylabel(r"RMS strain rate ($\mu\varepsilon$ s$^{-1}$)")
    ax.set_title(f"Amplitude against distance, {n_drops}-drop stack")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "basic06_amplitude_vs_distance.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=FIG)
    ax.plot(z, snr_db, lw=1.4)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(10, color="0.5", lw=0.8, ls="--")
    ax.set_xlabel("distance along fiber (m)")
    ax.set_ylabel("SNR (dB)")
    ax.set_title(f"Signal-to-noise ratio against distance, {n_drops}-drop stack")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "basic07_snr_vs_distance.png", dpi=150)
    plt.close(fig)

    # ---- 08 relative traveltime by channel-to-channel cross-correlation ----
    # A threshold pick on the envelope does not work here: the pre-drop level is
    # already exceeded at t=0 on many channels, so picks pin to zero and the
    # picker latches onto different cycles at different depths. Instead track the
    # lag between neighbouring channels by cross-correlation and accumulate it.
    # That is differential timing -- immune to origin-time error, and the same
    # argument dvv_core.py makes for interval velocity.
    step = 4                                   # channels per correlation step
    # The search window must be bounded by physics, not left wide open. Over
    # step*dx = 5.1 m, an apparent velocity of 1000 m/s is 5.1 ms; anything
    # slower is not a body wave here. A 120 ms window let the tracker jump whole
    # cycles, and because the lags accumulate those skips propagated.
    v_min = 1000.0
    half = max(3, int(np.ceil(step * DX_NANO / v_min * fs)))
    start = int(np.argmax(snr_db[:nch]))
    ref_c = start
    tt = {start: 0.0}
    lag_ok = True
    for c in range(start + step, nch, step):
        a = S[ref_c]
        b = S[c]
        cc = np.correlate(b - b.mean(), a - a.mean(), mode="same")
        mid = cc.size // 2
        seg = cc[mid:mid + half]               # causal: deeper arrives later
        if seg.size < 3 or not np.any(seg):
            lag_ok = False
            break
        shift = int(np.argmax(seg)) / fs
        tt[c] = tt[ref_c] + shift
        ref_c = c
    cs = np.array(sorted(tt))
    zz = cs * DX_NANO
    ts = np.array([tt[c] for c in cs])
    keep = snr_db[cs] > 6.0

    fig, ax = plt.subplots(figsize=FIG)
    ax.plot(ts[keep] * 1e3, zz[keep], "k.", ms=5,
            label="accumulated inter-channel lag")
    vel = np.nan
    if keep.sum() > 10:
        A = np.vstack([ts[keep], np.ones(int(keep.sum()))]).T
        vel, icpt = np.linalg.lstsq(A, zz[keep], rcond=None)[0]
        tl = np.array([ts[keep].min(), ts[keep].max()])
        ax.plot(tl * 1e3, vel * tl + icpt, "r-", lw=1.8,
                label=f"fit: {vel:.0f} m/s")
        print(f"  moveout fit: {vel:.0f} m/s over {int(keep.sum())} points "
              f"({zz[keep].min():.0f}-{zz[keep].max():.0f} m)")
    ax.invert_yaxis()
    ax.set_xlabel("relative traveltime (ms)")
    ax.set_ylabel("distance along fiber (m)")
    ax.set_title("Relative traveltime by inter-channel cross-correlation")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "basic08_moveout.png", dpi=150)
    plt.close(fig)


    # ---- 09 band comparison ---------------------------------------------
    bands = [None, (1.0, 5.0), (5.0, 20.0), (20.0, 50.0), (50.0, 100.0),
             (100.0, 250.0)]
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True)
    tt1 = np.arange(one.shape[1]) / fs - PRE_S
    zz1 = np.arange(one.shape[0]) * DX_NANO
    for a, bd in zip(axes.ravel(), bands):
        x = despike(bandpass(one, fs, bd) if bd else one)
        v = np.percentile(np.abs(x), 99.7)
        im = a.pcolormesh(tt1, zz1, x, cmap="seismic", vmin=-v, vmax=v,
                          shading="auto")
        a.set_title("unfiltered" if bd is None else f"{bd[0]:g}-{bd[1]:g} Hz")
        plt.colorbar(im, ax=a, label=r"$\mu\varepsilon$ s$^{-1}$")
    axes[0, 0].set_ylim(ZMAX, 0)
    axes[0, 0].set_xlim(-0.1, 0.8)
    for a in axes[1]:
        a.set_xlabel("time (s)")
    for a in axes[:, 0]:
        a.set_ylabel("distance along fiber (m)")
    fig.suptitle(f"One drop bandpassed into six bands (burst {burst_id})")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "basic09_bands.png", dpi=150)
    plt.close(fig)

    # ---- 10 spectrogram --------------------------------------------------
    from scipy.signal import spectrogram as _spec
    c_spec = min(int(150.0 / DX_NANO), nano.shape[0] - 1)
    ff, ts_, Sxx = _spec(nano[c_spec].astype(float), fs=fs, nperseg=2048,
                         noverlap=1024)
    db = 10 * np.log10(Sxx + 1e-30)
    fig, ax = plt.subplots(figsize=(12, 6.5))
    im = ax.pcolormesh(ts_, ff, db, cmap="magma", shading="auto",
                       vmin=np.percentile(db, 5), vmax=np.percentile(db, 99.5))
    for dr in drops:
        ax.axvline((dr["utc_time"] - nano_time(name)).total_seconds(),
                   color="c", lw=0.6, alpha=0.75)
    ax.set_ylim(0, 200)
    ax.set_xlabel("seconds into the 5-minute raw file")
    ax.set_ylabel("frequency (Hz)")
    ax.set_title(f"Spectrogram at {c_spec * DX_NANO:.0f} m "
                 f"(cyan = drop times from the manifest)")
    plt.colorbar(im, ax=ax, label="dB")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "basic10_spectrogram.png", dpi=150)
    plt.close(fig)
    del nano

    # ---- 11 spectrum against distance ------------------------------------
    nfft = sig.shape[1]
    A = np.abs(np.fft.rfft(sig, axis=1))
    fq = np.fft.rfftfreq(nfft, 1 / fs)
    m = (fq >= 1) & (fq <= 200)
    AdB = 20 * np.log10(A[:, m] + 1e-30)
    fig, ax = plt.subplots(figsize=(10, 7.5))
    im = ax.pcolormesh(fq[m], z, AdB, cmap="magma", shading="auto",
                       vmin=np.percentile(AdB, 5), vmax=np.percentile(AdB, 99.5))
    ax.set_ylim(ZMAX, 0)
    ax.set_xscale("log")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("distance along fiber (m)")
    ax.set_title(f"Signal amplitude spectrum against distance, "
                 f"{n_drops}-drop stack")
    plt.colorbar(im, ax=ax, label="dB")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "basic11_spectrum_vs_distance.png", dpi=150)
    plt.close(fig)

    # ---- 12 survey fold ---------------------------------------------------
    from collections import Counter
    per = Counter(r["burst_id"] for r in rows if r["nano_available"])
    bids = sorted(per)
    t0 = min(r["utc_time"] for r in rows)
    hrs = []
    for b in bids:
        ts = [r["utc_time"] for r in rows if r["burst_id"] == b]
        hrs.append((min(ts) - t0).total_seconds() / 3600.0)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(hrs, [per[b] for b in bids], width=0.25, color="0.35")
    ax.set_xlabel("hours into the survey")
    ax.set_ylabel("drops in burst")
    ax.set_title(f"Survey timeline: {sum(per.values())} drops in "
                 f"{len(bids)} bursts")
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "basic12_survey_fold.png", dpi=150)
    plt.close(fig)

    j = np.argmin(np.abs(z - 460))
    below = np.flatnonzero(snr_db < 10)
    where = f"{z[below[0]]:.0f} m" if below.size else "not within 0-600 m"
    print(f"  SNR at 460 m: {snr_db[j]:.1f} dB;  first drops below 10 dB at {where}")
    print("wrote", OUT_DIR)


if __name__ == "__main__":
    main()
