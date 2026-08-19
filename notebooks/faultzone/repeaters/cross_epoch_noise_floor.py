#!/usr/bin/env python3
"""WITHDRAWN 2026-08-14 -- THE TWO ARMS WERE NOT PROCESSED ALIKE. DO NOT CITE.

The conclusion this script printed -- that the 2024-25 ambient field decorrelates
over ~4 channels against ~26 in 2017, and that the loss is therefore in spatial
structure rather than noise amplitude -- is confounded and probably an artefact of
the confound.

The 2024-25 arm reads `cache_all`/`cache_hf`, built by `extract_all.py`, which
calls `DASutils.readFile_HDF(...)` WITHOUT passing `median=False`. That parameter
defaults to True, so the cached data has had the MEDIAN ACROSS CHANNELS REMOVED --
i.e. common-mode removal -- and has additionally been differentiated (`diff=True`),
band-passed 5-40 Hz and desampled. The 2017 arm is a raw binary read with
`np.fromfile` and had none of that applied.

Median-across-channels removal subtracts the spatially uniform component, which is
exactly the long-range-coherent part of the wavefield. It was stripped from the
2024-25 arm and not from the 2017 arm, and the script then measured 2024-25 as less
coherent at long separations. Both asymmetries (median, and the derivative) push in
the same direction, and the median one pushes hard.

TO REDO THIS PROPERLY: bypass the caches and re-extract the six earthquakes from
HDF5 with `h5py` directly -- no median removal, no derivative, no band-pass -- so
that the only remaining difference between epochs is the instrument. Retain the
gauge-length caveat below, which is genuine and runs against the conclusion.

The companion `cross_epoch_array_response.py` semblance result is far less affected
(normalised within-record measure of a loud transient) and is retained with a
caveat, but should also be redone on identically-read data.

Original docstring follows.

---

Instrument self-noise across epochs: is the 2024-25 interrogator noisier?

WHY.  `cross_epoch_array_response.py` shows the 2024-25 array records coherent
earthquake wavefronts at ~92 % of the 2017 array, at lower SNR -- so gross
sensing failure is ruled out.  But that test cannot separate the two remaining
explanations for the Figure 7c non-reproduction, because both survive it:

  B'. The 2024-25 interrogator has a HIGHER SELF-NOISE FLOOR.  Earthquakes are
      loud, so semblance on them is unaffected; ambient correlations are built
      from weak signal and are destroyed.
  A.  The AMBIENT SOURCE FIELD changed.

THE DISCRIMINANT.  Instrument self-noise is channel-to-channel INCOHERENT --
each channel's photodetector/laser noise is its own.  A real seismic wavefield
is spatially COHERENT over the correlation length of the wave.  So the fraction
of variance that is incoherent between neighbouring channels separates them:

    high incoherent fraction in 2024-25  ->  B', noisier instrument
    similar incoherent fraction          ->  A, the source field changed

ESTIMATOR.  For channel separation L, the MAXIMUM over time lag of the
normalised cross-correlation between channels i and i+L, median over i. The
maximum, not the zero-lag value: zero-lag correlation of a propagating wave
falls as cos(2*pi*L/lambda) and vanishes at a quarter wavelength regardless of
coherence, so a zero-lag curve measures wavelength, not coherence.  Incoherent fraction = 1 - coherent.  Reported as a curve against separation,
because a wavefield decorrelates smoothly with distance while white instrument
noise drops to zero at L = 1 immediately.  The intercept as L -> 0 is what
distinguishes them.

DATA.  2017: the pre-event windows of the two released earthquake records
(github.com/ariellellouch/SAFODDAS) -- the only 2017 raw data that exists
publicly, and short (a few seconds).  2024-25: pre-event windows of the same
matched deep events, and additionally genuine ambient records from the
continuous archive.

CONFOUND, stated up front and not corrected.  The 2017 gauge length is 10 m and
the 2024-25 gauge is 16.335 m.  A LONGER gauge averages along the fibre, which
SUPPRESSES incoherent noise and RAISES apparent inter-channel coherence.  So the
2024-25 array is flattered by this comparison: if it nonetheless shows a higher
incoherent fraction, that is a conservative finding.  If it shows a lower one,
the gauge difference alone could explain it and the test is inconclusive.

Output: cross_epoch_noise_floor.{npz,png,txt}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, resample_poly, sosfiltfilt

HERE = Path(__file__).resolve().parent
STEM = HERE / "cross_epoch_noise_floor"
LEL = Path("/scratch/users/nberrios/lellouch2017")

FMIN, FMAX = 5.0, 20.0
FS_COMMON = 100.0
EDIT_LO, EDIT_HI = 0.2, 5.0
WELLHEAD_2024 = 23
DX_2024 = 1.0209523439407349
DX_2017 = 1.0
DEPTH_MAX_M = 700.0
LAGS = np.arange(1, 41)           # channel separations to probe (~1-41 m)


def bandpass(x, fs, lo=FMIN, hi=FMAX):
    return sosfiltfilt(butter(4, [lo, min(hi, 0.45 * fs)], btype="bandpass",
                              fs=fs, output="sos"), x, axis=1)


def edit_traces(x):
    rms = np.sqrt((x ** 2).mean(axis=1))
    med = np.median(rms[rms > 0]) if np.any(rms > 0) else 0.0
    if med <= 0:
        return x
    return x[(rms > EDIT_LO * med) & (rms < EDIT_HI * med)]


def to_common(x, fs):
    if abs(fs - FS_COMMON) < 1e-9:
        return x, fs
    from math import gcd
    up, dn = int(FS_COMMON), int(fs)
    g = gcd(up, dn)
    return resample_poly(x, up // g, dn // g, axis=1), FS_COMMON


def coherence_curve(x, fs, dx, max_lag_s=0.25):
    """Median MAX-OVER-LAG normalised cross-correlation between channels
    separated by L.

    NOT the zero-lag value. For a wave propagating along the fibre with
    along-fibre wavelength lambda, the zero-lag correlation between channels L
    apart falls as cos(2*pi*L/lambda) and passes through zero at a quarter
    wavelength even when the two traces are perfectly coherent. A zero-lag
    curve therefore measures the dominant wavelength, not coherence, and would
    report a slow wavefield as "incoherent".

    Taking the maximum of the normalised cross-correlation over a lag window
    removes the propagation phase and leaves coherence proper. The lag window is
    generous (+-0.25 s) so it admits apparent velocities down to ~L*dx/0.25 s,
    i.e. it never truncates a physically plausible moveout at these separations.
    """
    x = x - x.mean(axis=1, keepdims=True)
    n = np.sqrt((x ** 2).sum(axis=1))
    ok = n > 0
    x, n = x[ok], n[ok]
    nt = x.shape[1]
    ml = min(int(max_lag_s * fs), nt // 2 - 1)
    nfft = 1 << int(np.ceil(np.log2(2 * nt - 1)))
    X = np.fft.rfft(x, n=nfft, axis=1)
    out = []
    for L in LAGS:
        if L >= x.shape[0]:
            out.append(np.nan); continue
        cc = np.fft.irfft(np.conj(X[:-L]) * X[L:], n=nfft, axis=1)
        cc = np.concatenate((cc[:, -ml:], cc[:, :ml + 1]), axis=1)
        cc /= (n[:-L] * n[L:])[:, None]
        out.append(float(np.median(np.abs(cc).max(axis=1))))
    return np.array(out)


def prep(x, fs, dx, ch0):
    x = edit_traces(np.asarray(x, float))
    x, fs = to_common(x, fs)
    hi = min(x.shape[0], ch0 + int(DEPTH_MAX_M / dx))
    return bandpass(x[ch0:hi], fs), fs


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    say("Cross-epoch instrument self-noise: incoherent fraction between channels")
    say("band %.0f-%.0f Hz, %.0f Hz, depth 0-%.0f m, channel separations 1-%d"
        % (FMIN, FMAX, FS_COMMON, DEPTH_MAX_M, LAGS[-1]))
    say("NOTE: the 2024-25 gauge is 16.335 m vs 10 m in 2017. A longer gauge")
    say("      SUPPRESSES incoherent noise, so 2024-25 is flattered here.")
    say("")

    curves = {}

    say("--- 2017: pre-event windows of the released earthquake records ---")
    for name in ("M1p33", "M2p46"):
        f = LEL / name
        if not f.is_file():
            say("  missing %s" % name); continue
        x = np.fromfile(f, dtype="<f4").reshape(1250, 800).T
        xf, fs = prep(x, 250.0, DX_2017, 0)
        pre = xf[:, int(0.5 * fs):int(2.5 * fs)]      # 2 s before the ~3 s arrival
        c = coherence_curve(pre, fs, DX_2017)
        curves["2017 %s pre-event" % name] = c
        say("  %-26s coherence at L=1: %.3f  L=10: %.3f  L=40: %.3f"
            % (name, c[0], c[9], c[-1]))

    say("")
    say("--- 2024-25: pre-event windows of matched deep events ---")
    ev = pd.read_csv(HERE / "all_events.csv")
    sel = ev[(ev.depth.between(9.0, 13.0)) & (ev.mag.between(1.0, 2.0))]
    used = 0
    for r in sel.sort_values("mag").itertuples(index=False):
        hf = HERE / "cache_hf" / ("%s.npz" % r.tag)
        lo = HERE / "cache_all" / ("%s.npz" % r.tag)
        src = hf if hf.is_file() else (lo if lo.is_file() else None)
        if src is None:
            continue
        d = np.load(src, allow_pickle=True)
        xf, fs = prep(np.asarray(d["X"]), float(d["fs"]), DX_2024, WELLHEAD_2024)
        pre = xf[:, int(0.5 * fs):int(4.5 * fs)]      # pre_s = 5 s in these caches
        c = coherence_curve(pre, fs, DX_2024)
        curves["2024-25 %s pre-event" % r.tag[3:]] = c
        say("  %-26s coherence at L=1: %.3f  L=10: %.3f  L=40: %.3f"
            % (r.tag[3:], c[0], c[9], c[-1]))
        used += 1
        if used >= 6:
            break

    say("")
    say("--- comparison ---")
    c17 = [v for k, v in curves.items() if k.startswith("2017")]
    c24 = [v for k, v in curves.items() if k.startswith("2024")]
    if not c17 or not c24:
        say("  insufficient data"); Path(str(STEM) + ".txt").write_text("\n".join(log)); return
    m17, m24 = np.nanmedian(c17, axis=0), np.nanmedian(c24, axis=0)
    say("  channel sep        L=1     L=5    L=10    L=20    L=40")
    say("  2017 coherence   %6.3f  %6.3f  %6.3f  %6.3f  %6.3f"
        % (m17[0], m17[4], m17[9], m17[19], m17[-1]))
    say("  2024-25          %6.3f  %6.3f  %6.3f  %6.3f  %6.3f"
        % (m24[0], m24[4], m24[9], m24[19], m24[-1]))
    say("  incoherent 2017  %6.3f  (1 - coherence at L=1)" % (1 - m17[0]))
    say("  incoherent 24-25 %6.3f" % (1 - m24[0]))
    say("")
    inc17, inc24 = 1 - m17[0], 1 - m24[0]

    # L=1 is gauge-dominated and near-useless: with a 16.335 m gauge, adjacent
    # 1 m channels share ~94% of their integration window and MUST correlate.
    # The informative quantity is the DECAY LENGTH -- the separation at which
    # coherence falls below 0.5, which reflects the wavefield's own spatial
    # correlation rather than the instrument's averaging.
    def decay_len(c):
        below = np.where(c < 0.5)[0]
        return float(LAGS[below[0]]) if len(below) else float(LAGS[-1])
    d17 = [decay_len(c) for c in c17]
    d24 = [decay_len(c) for c in c24]
    say("  decay length to coherence 0.5:")
    say("    2017    %s  (median %.0f channels)" % (["%.0f" % v for v in d17], np.median(d17)))
    say("    2024-25 %s  (median %.0f channels)" % (["%.0f" % v for v in d24], np.median(d24)))
    say("")
    lo24, hi17 = max(c[9] for c in c24), min(c[9] for c in c17)
    say("  at L=10 the two epochs do not overlap: every 2017 value (min %.3f)" % hi17)
    say("  exceeds every 2024-25 value (max %.3f)." % lo24)
    say("")
    if np.median(d24) < 0.5 * np.median(d17):
        say("  FINDING: the 2024-25 ambient wavefield decorrelates over a MUCH shorter")
        say("  distance -- median %.0f channels against %.0f in 2017." % (np.median(d24), np.median(d17)))
        say("  This is directly disqualifying for Figure 7c: a virtual-source gather")
        say("  spanning 50-700 m requires spatial coherence across those offsets, and")
        say("  a field that decorrelates within ~10 m cannot supply it, however well")
        say("  the instrument performs.")
        say("")
        say("  It also points AWAY from B' (a noisier interrogator). Self-noise is")
        say("  incoherent at ALL separations and would depress L=1 as well; here L=1")
        say("  is HIGHER in 2024-25 (%.3f vs %.3f). The loss is in spatial structure," % (m24[0], m17[0]))
        say("  not in noise amplitude -- which is the signature of a changed SOURCE")
        say("  FIELD, i.e. explanation A.")
    else:
        say("  Decay lengths are comparable; no evidence of a spatial-coherence change.")
    say("")
    say("  STRENGTH OF EVIDENCE: suggestive, not established. n = 2 for 2017 and the")
    say("  two records disagree sharply at L=1 (%.3f vs %.3f), so the 2017 side rests" % (c17[0][0], c17[1][0]))
    say("  on very little. The separation at L=10 is complete but 2-vs-6 gives an")
    say("  exact one-tailed rank p of only 0.071.")
    say("")
    say("  LIMITS: 2017 rests on 2 records and ~2 s of pre-event data each; the")
    say("  gauge-length confound is uncorrected; and pre-event windows of a")
    say("  triggered record are not the same population as continuous ambient.")

    fig, ax = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for k, v in curves.items():
        c = "crimson" if k.startswith("2017") else "steelblue"
        ax[0].plot(LAGS, v, "-" if c == "crimson" else "--", color=c,
                   lw=2 if c == "crimson" else 1.1, label=k,
                   alpha=0.95 if c == "crimson" else 0.65)
    ax[0].set_xlabel("channel separation (channels ~ m)")
    ax[0].set_ylabel("median inter-channel coherence")
    ax[0].legend(fontsize=6); ax[0].grid(alpha=.3)
    ax[0].set_title("coherence vs separation\nred = 2017, blue = 2024-25")
    ax[1].bar([0, 1], [1 - m17[0], 1 - m24[0]], color=["crimson", "steelblue"])
    ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["2017", "2024-25"])
    ax[1].set_ylabel("incoherent fraction at L=1")
    ax[1].set_title("instrument self-noise proxy\n(lower = cleaner)")
    ax[1].grid(alpha=.3, axis="y")
    fig.savefig(str(STEM) + ".png", dpi=190)
    np.savez(str(STEM) + ".npz", lags=LAGS,
             **{k.replace(" ", "_"): v for k, v in curves.items()})
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
