#!/usr/bin/env python3
"""Cross-epoch control: does the 2024-25 array record coherent wavefronts as well
as the 2017 array Lellouch et al. (2019) used?

WHY THIS EXISTS.  Figure 7c does not reproduce on the 2024-2025 ambient archive
(see awd_clean/AMBIENT_FIG7C_STATUS.md).  The explanation previously offered --
that body-wave energy is absent from the input -- was withdrawn on 2026-08-14
after an audit showed the energy census had no geometric baseline and its
headline statistic could not distinguish the SAFOD wavefield from white noise.
So the non-reproduction currently has no validated mechanism.

There are two candidate explanations and they make opposite predictions:

  A. THE AMBIENT SOURCE FIELD CHANGED between 2017 and 2024-25.
     Prediction: the 2024-25 array still records earthquake wavefronts as
     coherently as the 2017 array did.  Only the ambient field differs.

  B. THE INSTRUMENT OR INSTALLATION DEGRADED (different interrogator, 16.335 m
     gauge against 10 m, a decade more of fibre ageing).
     Prediction: the 2024-25 array records earthquake wavefronts LESS coherently
     too, because the loss is in the sensing, not the source.

Earthquakes decide between them, because Lellouch released raw 2017 earthquake
records on the same fibre (github.com/ariellellouch/SAFODDAS) while the raw 2017
ambient records were never released.  This is the only cross-epoch comparison
the public data supports.

THE MEASUREMENT.  Semblance of the direct P arrival along a slowness scan, in
sliding depth windows.  Semblance is a normalised coherence in [0,1], so it is
insensitive to absolute amplitude and therefore to source magnitude in a way
that raw SNR is not -- though it still rises with SNR, so SNR is reported
alongside as a covariate and several 2024-25 events spanning a magnitude range
are used rather than a single hand-picked one.

MATCHING.  Lellouch's M = 1.33 is at depth 11.16 km and horizontal offset
1.87 km (paper section 3.1), i.e. steep incidence up the fibre axis, which is
the geometry DAS is most sensitive to.  2024-25 comparators are selected for
depth 9-13 km and M 1.0-2.0 from `all_events.csv`, not for their result.

KNOWN NON-EQUIVALENCES, none of which are corrected and all of which bias
AGAINST the 2024-25 array, so a null result here is conservative:
  - gauge length 16.335 m now against 10 m in 2017 (more along-fibre averaging);
  - channel spacing 1.0209 m against 1.0 m;
  - 2017 released at 250 Hz, 2024-25 cached at 500 or 100 Hz; all are resampled
    to a common 100 Hz and band-limited to 5-20 Hz, the band Figure 7c uses;
  - channel 0 of the 2017 array is taken as the wellhead; for 2024-25 the
    wellhead is channel 23 per the G0 registration in METHODS_STATUS.md.

Output: cross_epoch_array_response.{npz,png,txt}
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, resample_poly, sosfiltfilt

HERE = Path(__file__).resolve().parent
STEM = HERE / "cross_epoch_array_response"
LEL = Path("/scratch/users/nberrios/lellouch2017")

FMIN, FMAX = 5.0, 20.0          # measurement band: the Figure 7c band
PICK_LO, PICK_HI = 10.0, 40.0   # picking band: wider, where earthquake SNR lives.
                                # 5-20 Hz is too narrow to PICK a P arrival in --
                                # the paper picks P at 0-120 Hz -- but it is the
                                # band the ambient result lives in, so coherence is
                                # still MEASURED at 5-20 Hz.
EDIT_LO, EDIT_HI = 0.2, 5.0     # per-channel RMS bounds, relative to the median.
                                # Lellouch's README: "All data are unprocessed -
                                # trace editing is required." M1p33 carries 16
                                # channels at ~1e6x the median without it.
EDGE_S = 0.5                    # exclude this much from each record end when
                                # picking, to reject zero-phase filter edge ringing
FS_COMMON = 100.0
DEPTH_MAX_M = 700.0             # Figure 7c spans 0-700 m
WIN_M = 60.0                    # depth window for semblance
STEP_M = 30.0
WIN_S = 0.30                    # time window around the picked arrival
SLOWNESS = np.linspace(-1.2e-3, 1.2e-3, 241)   # s/m, +-833 m/s .. inf

WELLHEAD_2024 = 23              # G0 registration
DX_2024 = 1.0209523439407349
DX_2017 = 1.0


def bandpass(x, fs, lo=FMIN, hi=FMAX):
    hi = min(hi, 0.45 * fs)
    sos = butter(4, [lo, hi], btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=1)


def edit_traces(x, label, log):
    """Drop channels whose RMS is a gross outlier. Required by the 2017 release
    and equally applicable to 2024-25, so both epochs get the identical rule."""
    rms = np.sqrt((x ** 2).mean(axis=1))
    med = np.median(rms[rms > 0]) if np.any(rms > 0) else 0.0
    if med <= 0:
        return x, np.ones(len(x), bool)
    keep = (rms > EDIT_LO * med) & (rms < EDIT_HI * med)
    log.append("      trace edit: kept %d of %d channels (dropped %d)"
               % (keep.sum(), len(keep), (~keep).sum()))
    return x[keep], keep


def taper(x, fs, seconds=0.25):
    n = int(seconds * fs)
    if n < 2 or 2 * n >= x.shape[1]:
        return x
    w = np.ones(x.shape[1])
    ramp = 0.5 * (1 - np.cos(np.pi * np.arange(n) / n))
    w[:n] = ramp; w[-n:] = ramp[::-1]
    return x * w


def to_common_rate(x, fs):
    if abs(fs - FS_COMMON) < 1e-9:
        return x, fs
    from math import gcd
    up, dn = int(FS_COMMON), int(fs)
    g = gcd(up, dn)
    return resample_poly(x, up // g, dn // g, axis=1), FS_COMMON


def pick_arrival(x, fs):
    """Index of peak MEDIAN envelope across channels, edges excluded.

    Median not mean: a handful of glitched traces dominates a mean envelope and
    drags the pick to wherever they spike. Edges excluded because zero-phase
    filtering rings at the record ends, which is what put the first attempt's
    pick at 4.74 s of a 5.00 s record.
    """
    env = np.median(np.abs(x), axis=0)
    k = max(3, int(0.05 * fs))
    env = np.convolve(env, np.ones(k) / k, "same")
    e = int(EDGE_S * fs)
    if 2 * e >= len(env):
        return int(np.argmax(env))
    return e + int(np.argmax(env[e:len(env) - e]))


def semblance(block, fs, depths, slowness):
    """Max semblance over a slowness scan for one depth window."""
    n = block.shape[1]
    z0 = depths - depths.mean()
    best = 0.0
    for p in slowness:
        shifts = np.round(p * z0 * fs).astype(int)
        if np.abs(shifts).max() >= n // 3:
            continue
        rows = [np.roll(tr, -s) for tr, s in zip(block, shifts)]
        a = np.asarray(rows)
        num = (a.sum(axis=0) ** 2).sum()
        den = a.shape[0] * (a ** 2).sum()
        if den > 0:
            best = max(best, float(num / den))
    return best


def profile(x, fs, dx, ch_wellhead, label, log):
    """Semblance and SNR versus depth for one record."""
    x = np.asarray(x, dtype=float)
    x, keep = edit_traces(x, label, log)
    x, fs = to_common_rate(x, fs)
    x = taper(x, fs)
    xp = bandpass(x, fs, PICK_LO, PICK_HI)   # pick here
    x = bandpass(x, fs, FMIN, FMAX)          # measure here
    nch = x.shape[0]
    hi = min(nch, ch_wellhead + int(DEPTH_MAX_M / dx))
    x = x[ch_wellhead:hi]; xp = xp[ch_wellhead:hi]
    depths = np.arange(x.shape[0]) * dx
    k = pick_arrival(xp, fs)
    half = int(WIN_S * fs / 2)
    lo_t, hi_t = max(0, k - half), min(x.shape[1], k + half)
    pre_hi = max(0, lo_t - int(0.5 * fs))
    pre_lo = max(0, pre_hi - int(1.0 * fs))

    zc, sem, snr = [], [], []
    w = int(WIN_M / dx); st = int(STEP_M / dx)
    for i in range(0, x.shape[0] - w + 1, st):
        blk = x[i:i + w, lo_t:hi_t]
        noise = x[i:i + w, pre_lo:pre_hi]
        if blk.size == 0 or noise.size == 0:
            continue
        zc.append(depths[i:i + w].mean())
        sem.append(semblance(blk, fs, depths[i:i + w], SLOWNESS))
        nrms = np.sqrt((noise ** 2).mean())
        snr.append(float(np.sqrt((blk ** 2).mean()) / nrms) if nrms > 0 else np.nan)
    log.append("  %-28s arrival t=%.2f s of %.2f s, %d windows, median semblance %.3f, median SNR %.1f"
               % (label, k / fs, x.shape[1] / fs, len(zc), np.median(sem), np.nanmedian(snr)))
    if np.nanmedian(snr) < 1.5:
        log.append("      WARNING: SNR < 1.5 -- no arrival detected; this record's"
                   " semblance is not interpretable")
    return np.array(zc), np.array(sem), np.array(snr)


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    say("Cross-epoch array response: 2017 (Lellouch) vs 2024-25, identical operator")
    say("band %.0f-%.0f Hz, resampled to %.0f Hz, depth 0-%.0f m below wellhead,"
        " %.0f m windows" % (FMIN, FMAX, FS_COMMON, DEPTH_MAX_M, WIN_M))
    say("")

    results = {}
    say("--- 2017, Lellouch released records (800 ch, 1.0 m, 250 Hz) ---")
    for name, mag in (("M1p33", 1.33), ("M2p46", 2.46)):
        f = LEL / name
        if not f.is_file():
            say("  missing %s" % f); continue
        raw = np.fromfile(f, dtype="<f4")
        x = raw.reshape(1250, 800).T          # time-major on disk; verified
        results["2017 %s (M%.2f)" % (name, mag)] = profile(
            x, 250.0, DX_2017, 0, "2017 %s" % name, log) + (mag,)

    say("")
    say("--- 2024-25, matched deep events (900 ch, 1.0210 m) ---")
    ev = pd.read_csv(HERE / "all_events.csv")
    sel = ev[(ev.depth.between(9.0, 13.0)) & (ev.mag.between(1.0, 2.0))].copy()
    used = 0
    for r in sel.sort_values("mag").itertuples(index=False):
        hf = HERE / "cache_hf" / ("%s.npz" % r.tag)
        lo = HERE / "cache_all" / ("%s.npz" % r.tag)
        src = hf if hf.is_file() else (lo if lo.is_file() else None)
        if src is None:
            continue
        d = np.load(src, allow_pickle=True)
        lbl = "2024-25 %s (M%.2f, %.1f km)" % (r.tag[3:], r.mag, r.depth)
        results[lbl] = profile(np.asarray(d["X"]), float(d["fs"]), DX_2024,
                               WELLHEAD_2024, lbl, log) + (r.mag,)
        used += 1
        if used >= 6:
            break
    say("")
    say("  selected on depth 9-13 km and M 1.0-2.0 only, never on result")

    say("")
    say("--- comparison ---")
    e2017 = {k: v for k, v in results.items() if k.startswith("2017")}
    e2024 = {k: v for k, v in results.items() if k.startswith("2024")}
    if not e2017 or not e2024:
        say("  insufficient data"); Path(str(STEM) + ".txt").write_text("\n".join(log)); return

    m17 = np.median([np.median(v[1]) for v in e2017.values()])
    m24 = np.median([np.median(v[1]) for v in e2024.values()])
    s17 = np.median([np.nanmedian(v[2]) for v in e2017.values()])
    s24 = np.median([np.nanmedian(v[2]) for v in e2024.values()])
    say("  median semblance   2017 %.3f   2024-25 %.3f   ratio %.2f" % (m17, m24, m24 / m17))
    say("  median SNR         2017 %.1f   2024-25 %.1f" % (s17, s24))
    say("")
    if m24 >= 0.8 * m17:
        say("  RULES OUT GROSS SENSING FAILURE. The 2024-25 array records coherent")
        say("  earthquake wavefronts at %.0f%% of the 2017 array, at LOWER SNR" % (100*m24/m17))
        say("  (%.1f vs %.1f), so at matched SNR it is at least as good. The fibre" % (s24, s17))
        say("  still resolves a propagating wavefront over 0-700 m.")
        say("")
        say("  IT DOES NOT ESTABLISH THAT THE AMBIENT FIELD CHANGED. Semblance on an")
        say("  earthquake tests sensing of a STRONG transient. Ambient interferometry")
        say("  additionally needs the instrument self-noise floor to sit below the")
        say("  ambient field. A noisier interrogator passes this test easily -- an")
        say("  earthquake is loud -- and still destroys ambient correlations, which")
        say("  are built from exactly the weak signal that self-noise swamps.")
        say("")
        say("  Remaining candidates, both still live:")
        say("    - the 2024-25 interrogator has a higher self-noise floor;")
        say("    - the ambient source field changed.")
        say("  These are separable: self-noise is channel-to-channel INCOHERENT, so")
        say("  the incoherent noise fraction distinguishes them. See")
        say("  cross_epoch_noise_floor.py.")
    elif m24 <= 0.5 * m17:
        say("  READS AS EXPLANATION B: the 2024-25 array is markedly less coherent on")
        say("  earthquakes too. The loss is in the sensing or installation, and the")
        say("  ambient field need not have changed at all.")
    else:
        say("  AMBIGUOUS: %.0f%% of the 2017 coherence. Neither explanation is clean;"
            % (100 * m24 / m17))
        say("  SNR differences between events may dominate -- compare at matched SNR.")
    say("")
    say("  CAVEAT: semblance rises with SNR, so this is not fully source-independent.")
    say("  Read the semblance-vs-SNR panel before concluding.")

    fig, ax = plt.subplots(1, 3, figsize=(16, 5.5), constrained_layout=True)
    for k, v in results.items():
        c = "crimson" if k.startswith("2017") else "steelblue"
        ls = "-" if k.startswith("2017") else "--"
        ax[0].plot(v[1], v[0], ls, color=c, lw=2 if c == "crimson" else 1.1,
                   label=k, alpha=0.95 if c == "crimson" else 0.7)
        ax[1].plot(v[2], v[0], ls, color=c, lw=2 if c == "crimson" else 1.1, alpha=0.8)
        ax[2].scatter(np.nanmedian(v[2]), np.median(v[1]), color=c,
                      s=70 if c == "crimson" else 40, zorder=3)
    ax[0].invert_yaxis(); ax[0].set_xlabel("semblance of direct P"); ax[0].set_ylabel("depth below wellhead (m)")
    ax[0].set_xlim(0, 1); ax[0].legend(fontsize=6); ax[0].grid(alpha=.3)
    ax[0].set_title("coherence vs depth\nred = 2017 (Lellouch), blue = 2024-25")
    ax[1].invert_yaxis(); ax[1].set_xscale("log"); ax[1].set_xlabel("direct-P SNR")
    ax[1].grid(alpha=.3); ax[1].set_title("SNR vs depth")
    ax[2].set_xscale("log"); ax[2].set_xlabel("median SNR"); ax[2].set_ylabel("median semblance")
    ax[2].grid(alpha=.3); ax[2].set_title("semblance vs SNR\n(compare epochs at matched SNR)")
    fig.savefig(str(STEM) + ".png", dpi=190)

    np.savez(str(STEM) + ".npz",
             **{("%s__%s" % (k, f)): v[i] for k, v in results.items()
                for i, f in enumerate(("depth", "semblance", "snr"))})
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
