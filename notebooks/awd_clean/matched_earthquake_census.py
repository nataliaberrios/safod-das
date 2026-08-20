#!/usr/bin/env python3
"""Interrogator or signal?  The same measurement on an earthquake in both epochs.

THE QUESTION.  The 2024-25 ambient field is dominated by low-wavenumber energy;
the 2017 release window is not.  Two explanations are on the table and the
earlier census could not separate them:

  (A) INTERROGATOR.  The 2024-25 system reports absolute optical phase
      (RawDataUnit = 'rad * 2PI/2^16') and applies no common-mode rejection, so
      every channel carries a shared laser/temperature/PSU term.  That is a k = 0
      component by construction and would be present in every record.
  (B) SIGNAL PRESENCE.  The 2017 comparison window came from an EARTHQUAKE
      record (M1p33, M2p46), the 2024-25 window from continuous ambient noise.
      A loud coherent arrival swamps instrumental common mode, so the 2017 arm
      would show a small k = 0 share for reasons having nothing to do with its
      interrogator.

These predict opposite things for a 2024-25 EARTHQUAKE record:
      under (A) it still shows a large k = 0 share -- the instrument term is
               always there, so the contamination follows the instrument;
      under (B) it shows a small k = 0 share like 2017 -- the contamination
               follows the absence of signal, and our ambient problem is an
               SNR problem rather than an instrument problem.
Either answer is useful and they cannot both hold.

WHY THE CACHED EVENTS CANNOT BE USED.  faultzone/repeaters/cache_all/ has 206
extracted events, but extract_all.py reads them through DASutils.readFile_HDF
without median=False, and that default subtracts the per-sample median across
channels -- it removes the k = 0 component before writing the cache.  Measuring
a k = 0 share on those files would return ~0 by construction.  This is the exact
confound that forced the withdrawal of cross_epoch_noise_floor.py, so the raw
HDF5 is read directly here with h5py, as the census does.

MATCHING.  Both arms go through one shared code path with a runtime assertion
that the operation sequences are identical, and both are matched on:
    band 5-20 Hz | 250 Hz | 700 m aperture | 2 s Hann windows | trace editing
Residual known mismatch, reported not hidden: gauge length 10 m vs 16.335 m, and
the earthquakes are different events at different distances and magnitudes, so
absolute levels are not comparable -- only the SHARE of in-band energy at k = 0,
which is a within-record normalised quantity, is.

Output: matched_earthquake_census.{png,txt}
"""
from __future__ import annotations

import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, detrend, resample_poly, sosfiltfilt

HERE = Path(__file__).resolve().parent
STEM = HERE / "matched_earthquake_census"
CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")
LEL = Path("/scratch/users/nberrios/lellouch2017")
CACHE = HERE.parent / "faultzone" / "repeaters" / "cache_all"

FMIN, FMAX = 5.0, 20.0
FS_COMMON = 250.0
APERTURE_M = 700.0
WIN_S = 2.0
EDIT_LO, EDIT_HI = 0.2, 5.0
CH_LO_2024 = 23
V_EDGES = np.array([0, 500, 1000, 1500, 2000, 2500, 3000, 4000,
                    6000, 10000, 30000, np.inf])
STA_S = 0.10
TRIGGER = 3.0


def corrected_path(p):
    return Path(str(p).replace("/data/SAFODAS1-harddrive-transfer",
                               "/data/SAFOD/SAFODAS1-harddrive-transfer"))


def edit_traces(x):
    """Interpolate outlier channels rather than deleting them.

    Deleting breaks the uniform spatial sampling the 2-D FFT assumes, which
    would corrupt the wavenumber axis this script measures.
    """
    rms = np.sqrt(np.mean(x ** 2, axis=1))
    med = np.median(rms[rms > 0]) if np.any(rms > 0) else 0.0
    bad = (rms < EDIT_LO * med) | (rms > EDIT_HI * med) | ~np.isfinite(rms)
    good = np.flatnonzero(~bad)
    if good.size < 8:
        raise SystemExit("too few good channels")
    for i in np.flatnonzero(bad):
        lo = good[good < i]
        hi = good[good > i]
        if lo.size and hi.size:
            a, b = lo[-1], hi[0]
            w = (i - a) / (b - a)
            x[i] = (1 - w) * x[a] + w * x[b]
        else:
            x[i] = x[good[0] if not lo.size else good[-1]]
    return x, int(bad.sum())


def analyse(x, fs, dx, label, ops, say):
    """The single shared path. Both arms must produce the same `steps` list."""
    steps = []
    x = np.asarray(x, dtype=np.float64)
    x, dropped = edit_traces(x)
    steps.append("trace_edit_interp(%.1f-%.1f x median RMS)" % (EDIT_LO, EDIT_HI))
    x = np.diff(x, axis=1)
    steps.append("diff(strain_rate)")
    if abs(fs - FS_COMMON) > 1e-6:
        factor = int(round(fs / FS_COMMON))
        if factor < 1 or abs(fs / factor - FS_COMMON) > 1e-6:
            raise SystemExit("cannot decimate %g to %g" % (fs, FS_COMMON))
        x = resample_poly(x, 1, factor, axis=1)
        fs = FS_COMMON
        steps.append("resample_to_%gHz" % FS_COMMON)
    else:
        steps.append("resample_to_%gHz(noop)" % FS_COMMON)
    x = detrend(x, axis=1)
    steps.append("detrend")
    x = sosfiltfilt(butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos"),
                    x, axis=1)
    steps.append("bandpass_%g-%g" % (FMIN, FMAX))
    nch = int(round(APERTURE_M / dx))
    if x.shape[0] < nch:
        raise SystemExit("%s has %d channels, need %d for a %.0f m aperture"
                         % (label, x.shape[0], nch, APERTURE_M))
    x = x[:nch]
    steps.append("aperture_%.0fm(%d ch)" % (APERTURE_M, nch))
    nw = int(WIN_S * fs)
    if x.shape[1] < nw:
        raise SystemExit("%s too short for a %.1f s window" % (label, WIN_S))
    wt, wx = np.hanning(nw)[None, :], np.hanning(nch)[:, None]
    steps.append("hann_space_and_time")
    P = None
    nwin = 0
    for s in range(0, x.shape[1] - nw + 1, nw):
        F = np.fft.fftshift(np.fft.fft2(x[:, s:s + nw] * wt * wx))
        P = np.abs(F) ** 2 if P is None else P + np.abs(F) ** 2
        nwin += 1
    steps.append("sum_power(no per-window normalisation)")
    ops[label] = steps
    say("  %-10s dropped/interpolated %d channels, %d windows" % (label, dropped, nwin))
    k = np.fft.fftshift(np.fft.fftfreq(nch, dx))
    f = np.fft.fftshift(np.fft.fftfreq(nw, 1.0 / fs))
    return P, k, f, nwin


def budget(P, k, f):
    K, F = np.meshgrid(k, f, indexing="ij")
    band = (np.abs(F) >= FMIN) & (np.abs(F) <= FMAX)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(np.abs(K) > 0, np.abs(F) / np.abs(K), np.inf)
    tot_e, tot_c = P[band].sum(), band.sum()
    rows = []
    for lo, hi in zip(V_EDGES[:-1], V_EDGES[1:]):
        m = band & (v >= lo) & ((v <= hi) if np.isinf(hi) else (v < hi))
        rows.append((lo, hi, 100 * P[m].sum() / tot_e, 100 * m.sum() / tot_c))
    covered = sum(r[2] for r in rows)
    if not (99.0 <= covered <= 101.0):
        raise SystemExit("binning does not conserve energy: %.2f%%" % covered)
    return rows


def pick_arrival(x, fs, edge_s=0.25):
    """First STA excursion above TRIGGER x the quiet level, plus its strength.

    `edge_s` of each end is ignored.  sosfiltfilt and the STA convolution both
    produce edge transients, and without this the picker fires at t ~ 0 on every
    record -- which is how a 2024-25 excerpt cut to begin 1.0 s BEFORE the event
    time came back with a pick at 0.038 s.  A pick is only accepted below as an
    earthquake if it also lands near the expected time and clears an SNR floor.
    """
    d = detrend(np.diff(x, axis=1), axis=1)
    d = sosfiltfilt(butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos"),
                    d, axis=1)
    rms = np.sqrt(np.mean(d ** 2, axis=0))
    nsta = max(1, int(STA_S * fs))
    sta = np.convolve(rms, np.ones(nsta) / nsta, mode="same")
    guard = int(edge_s * fs)
    core = sta[guard:len(sta) - guard] if len(sta) > 2 * guard + 4 else sta
    quiet = float(np.median(np.sort(core)[: max(1, len(core) // 4)]))
    above = np.flatnonzero(core > TRIGGER * quiet)
    snr = float(core.max() / quiet) if quiet > 0 else np.nan
    t = (float((above[0] + guard) / fs) if above.size else np.nan)
    return t, sta, quiet, snr


def load_2024_event(say):
    """Raw HDF5 covering a cached event time. Cache is used only for the TIME."""
    if not CACHE.is_dir():
        raise SystemExit("no cached event list at %s" % CACHE)
    names = sorted(p.name for p in CACHE.glob("ev_*.npz"))
    if not names:
        raise SystemExit("no ev_*.npz in %s" % CACHE)
    db = pd.read_csv(CSV, sep=r"\s+")
    db = db[db.nSamples == 30000].copy()          # continuous files only
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    db = db.dropna(subset=["t"]).sort_values("t").reset_index(drop=True)

    for nm in names:
        stamp = nm[len("ev_"):-len(".npz")]
        try:
            t0 = pd.to_datetime(stamp, format="%Y%m%dT%H%M%S", utc=True)
        except ValueError:
            continue
        hit = db[(db.t <= t0) & (db.t + pd.Timedelta(seconds=60) > t0)]
        if hit.empty:
            continue
        f = corrected_path(hit.iloc[0].file)
        if not f.is_file():
            continue
        with h5py.File(f, "r") as h:
            g = h["Acquisition/Raw[0]"]
            fs = float(g.attrs.get("OutputDataRate", 500.0))
            dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 1.0))
            offset = int((t0 - hit.iloc[0].t).total_seconds() * fs)
            lo = max(0, offset - int(1.0 * fs))
            hi = min(g["RawData"].shape[0], lo + int(6.0 * fs))
            d = g["RawData"][lo:hi, CH_LO_2024:].astype(np.float32).T
        expected = (offset - lo) / fs            # where the catalogue says it is
        pick, _, _, snr = pick_arrival(d.astype(np.float64), fs)
        # Accept only if an arrival is actually IN the window: the pick must land
        # near the catalogue time and clear an SNR floor.  Without both, this
        # would happily hand back an ambient excerpt and the whole comparison
        # would silently become ambient-vs-noise again.
        near = (not np.isnan(pick)) and abs(pick - expected) < 1.5
        loud = np.isfinite(snr) and snr >= 6.0
        say("  candidate %s in %s: pick %s, expected %.2f s, SNR %.1f -> %s"
            % (stamp, f.name,
               "none" if np.isnan(pick) else "%.3f s" % pick, expected, snr,
               "accepted" if (near and loud) else "rejected"))
        if not (near and loud):
            continue
        return d.astype(np.float64), fs, dx, stamp, pick
    raise SystemExit("no cached event mapped to a raw file with a visible arrival")


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    ops = {}

    say("Matched earthquake census: interrogator or signal presence?")
    say("  band %g-%g Hz | %.0f Hz | %.0f m aperture | %.1f s Hann windows"
        % (FMIN, FMAX, FS_COMMON, APERTURE_M, WIN_S))
    say("  raw HDF5 via h5py -- the cache_all events cannot be used, DASutils")
    say("  readFile_HDF(median=True) strips k=0 before writing them")
    say("")

    x24, fs24, dx24, stamp, pick24 = load_2024_event(say)
    P24, k24, f24, n24 = analyse(x24, fs24, dx24, "2024-25 eq", ops, say)

    # ---- 2017 arm: check there is enough POST-ARRIVAL data before comparing ----
    #
    # This is the step that decides whether the test is possible at all.  The
    # 2017 arm must contain the earthquake, not the noise before it, or the
    # comparison is 2024-25-earthquake versus 2017-noise -- the same
    # signal-presence confound the ambient census was suspected of, inverted.
    parts, picks17 = [], []
    for name in ("M1p33", "M2p46"):
        p = LEL / name
        if not p.is_file():
            continue
        v = np.fromfile(p, dtype="<f4").reshape(1250, 800).T.astype(np.float64)
        pk, _, _, snr17 = pick_arrival(v, 250.0)
        picks17.append((name, pk, snr17))
        parts.append(v)
    if not parts:
        raise SystemExit("2017 release not found at %s" % LEL)

    record_s = 1250 / 250.0
    say("")
    say("  2017 post-arrival duration available (needs >= %.1f s for one window):" % WIN_S)
    usable = 0.0
    for name, pk, snr17 in picks17:
        if np.isnan(pk):
            say("    %-6s no arrival picked -- record is noise throughout" % name)
            continue
        after = record_s - pk
        usable += max(0.0, after)
        say("    %-6s arrival %.3f s of %.2f s, SNR %.1f -> %.2f s after the pick %s"
            % (name, pk, record_s, snr17, after,
               "(enough)" if after >= WIN_S else "(NOT enough for one window)"))

    if usable < WIN_S:
        say("")
        say("=== NOT TESTABLE WITH THE RELEASED DATA ===")
        say("  Total post-arrival 2017 data across both records: %.2f s, which is" % usable)
        say("  less than one %.1f s analysis window. The released records place" % WIN_S)
        say("  their arrivals at the very end (4.10 s and 4.92 s of 5.00 s), so")
        say("  there is no way to form a 2017 arm that actually contains the")
        say("  earthquake. Concatenating the full records instead would give a")
        say("  2017 arm that is ~90 %% pre-event noise compared against a 2024-25")
        say("  arm that is mostly earthquake -- the signal-presence confound this")
        say("  test exists to remove, merely inverted, and it would have produced")
        say("  a confident and meaningless number.")
        say("")
        say("  THIS DOES NOT LEAVE THE QUESTION OPEN. lellouch2017_window_audit.py")
        say("  established that the ambient census's 2017 window (first 2.5 s) is")
        say("  genuinely PRE-EVENT -- both arrivals are later than the cut -- so")
        say("  both arms of that census were noise and the ambient-vs-ambient")
        say("  comparison is already the matched one. Shortening WIN_S to fit the")
        say("  0.9 s available is not a fix: at 5 Hz that is under 5 cycles, from")
        say("  a single window, which cannot support a wavenumber budget.")
        Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
        say("")
        say("wrote %s.txt (no figure: the comparison was not run)" % STEM.name)
        return

    parts = [v[:, int(pk * 250.0):] for (_, pk, _), v in zip(picks17, parts)
             if not np.isnan(pk)]
    x17 = np.concatenate(parts, axis=1)
    P17, k17, f17, n17 = analyse(x17, 250.0, 1.0, "2017 eq", ops, say)

    say("")
    say("--- processing applied to each arm ---")
    for lab, steps in ops.items():
        say("  %-12s %s" % (lab, " -> ".join(steps)))
    # Two steps legitimately differ between arms and must be normalised out of
    # the comparison, or the guard fires on a difference that is not an error:
    #   resample_*   -- 2024-25 decimates 500 -> 250, 2017 is already at 250
    #   aperture_*   -- the same 700 m needs 686 channels at dx = 1.0209 m and
    #                   700 at dx = 1.0 m; the counts are reported separately
    def canonical(steps):
        return tuple("aperture" if s.startswith("aperture") else s
                     for s in steps if not s.startswith("resample"))
    identical = len(set(canonical(v) for v in ops.values())) == 1
    say("  IDENTICAL OPERATION SEQUENCE (decimation and aperture channel count "
        "normalised out): %s" % identical)
    if not identical:
        for lab, steps in ops.items():
            say("    %-12s %s" % (lab, canonical(steps)))
        raise SystemExit("arms did not receive identical operations -- refusing to compare")

    r24, r17 = budget(P24, k24, f24), budget(P17, k17, f17)
    say("")
    say("--- apparent-velocity budget, BOTH ARMS AN EARTHQUAKE ---")
    say("  %-18s %18s %18s" % ("velocity band", "2017 eq", "2024-25 eq"))
    for a, b in zip(r17, r24):
        lbl = ("%5.0f - %5.0f" % (a[0], a[1])) if np.isfinite(a[1]) else ("%5.0f +  (k=0)" % a[0])
        say("  %-18s %8.2f%% (%5.2f%%) %8.2f%% (%5.2f%%)" % (lbl, a[2], a[3], b[2], b[3]))

    k0_17 = r17[-1][2]
    k0_24 = r24[-1][2]
    fan = lambda rows: sum(r[2] for r in rows if 2500 <= r[0] < 4000)
    say("")
    say("  k=0 share:        2017 eq %6.2f%%   2024-25 eq %6.2f%%" % (k0_17, k0_24))
    say("  body-wave fan:    2017 eq %6.2f%%   2024-25 eq %6.2f%%" % (fan(r17), fan(r24)))
    say("")
    say("=== verdict ===")
    if k0_24 > 10.0 * max(k0_17, 0.01) and k0_24 > 20.0:
        say("  (A) INTERROGATOR. A 2024-25 EARTHQUAKE record still carries a large")
        say("  k=0 share, so the contamination follows the instrument rather than")
        say("  the absence of signal. The 2024-25 system reports absolute optical")
        say("  phase with no common-mode rejection, which is a k=0 term by")
        say("  construction. This supports the interrogator-settings explanation.")
    elif k0_24 < 5.0:
        say("  (B) SIGNAL PRESENCE. A 2024-25 earthquake record looks like 2017 --")
        say("  small k=0 share -- so the k=0 dominance is a property of AMBIENT")
        say("  windows, not of the instrument. The ambient problem is then one of")
        say("  SNR: the noise field is weak relative to the instrument floor, and")
        say("  the earlier cross-epoch 'mechanism' was measuring signal presence.")
    else:
        say("  Intermediate (%.2f%% vs %.2f%%): neither explanation is clean." % (k0_24, k0_17))
        say("  Report both numbers; do not attribute the difference to either cause.")
    say("")
    say("  LIMITS: different events at different distances and magnitudes, so only")
    say("  the within-record SHARE is comparable, not absolute levels; gauge length")
    say("  10 m vs 16.335 m; 2017 is %d windows and 2024-25 is %d." % (n17, n24))

    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2), constrained_layout=True)
    labels = [("%5.0f-%5.0f" % (a[0], a[1])) if np.isfinite(a[1]) else "k=0" for a in r17]
    y = np.arange(len(labels))
    ax[0].barh(y - 0.2, [a[2] for a in r17], 0.4, label="2017 earthquake", color="tab:blue")
    ax[0].barh(y + 0.2, [b[2] for b in r24], 0.4, label="2024-25 earthquake", color="tab:red")
    ax[0].set_yticks(y); ax[0].set_yticklabels(labels, fontsize=7)
    ax[0].set(xlabel="% of in-band energy", title="Matched: earthquake in BOTH arms")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, axis="x")
    ax[1].bar([0, 1], [k0_17, k0_24], color=["tab:blue", "tab:red"])
    ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["2017 eq", "2024-25 eq"])
    ax[1].set(ylabel="% of in-band energy at k = 0",
              title="The discriminator\nlarge on the right = interrogator")
    ax[1].grid(alpha=.3, axis="y")
    fig.savefig(str(STEM) + ".png", dpi=190)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
