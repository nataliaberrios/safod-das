#!/usr/bin/env python3
"""Is the interrogator to blame for the static fixed-k pattern?

AMBIENT_LOWK_MECHANISM.md established WHAT the contaminant is -- a static
spatial pattern at fixed wavenumber holding ~39 % of the 5-20 Hz energy -- and
left WHERE IT COMES FROM open, on the grounds that the matched
earthquake-versus-earthquake test is impossible with Lellouch's release (only
0.90 s of post-arrival 2017 data exists).  That was giving up too early: the
question can be attacked with 2024-25 data alone, because an instrumental term
and an earth/site term make different predictions that do not need 2017 at all.

    an INSTRUMENTAL term (laser phase, PSU, interrogator electronics, a fixed
      optical response per channel) has a FIXED spatial fingerprint.  Lasers and
      power supplies do not change with the weather or the season.
    an EARTH or SITE term varies with noise conditions, which change hour to
      hour and month to month.

THREE TESTS, all on 2024-25 raw HDF5 read with h5py (never through
DASutils.readFile_HDF, whose median=True default subtracts the per-sample median
across channels and would delete the very thing being measured -- the confound
that forced the withdrawal of cross_epoch_noise_floor.py).

T1  STATIONARITY OF THE SPATIAL PATTERN -- the strong test.
    For each of several widely separated days, take the leading left singular
    vector u1 of the channel-by-time matrix: the dominant spatial pattern.
    Correlate u1 between days, using |corr| because the sign of a singular
    vector is arbitrary.  A pattern that is near-identical across months is
    instrumental.  A pattern that decorrelates between days is not.
    Control: the same statistic between u1 of one day and u2..u5 of another,
    which bounds what "unrelated spatial patterns" scores on this data.

T2  PRESENCE IN THE SURFACE LEAD-IN.
    Channel 23 is the wellhead (gate G0); channels 0-22 are the surface lead-in.
    Reported, but interpreted with care: the lead-in is at the surface and WILL
    record genuine surface ground motion, so its energy is not by itself proof
    of an instrumental origin.  What is diagnostic is whether u1 has comparable
    amplitude there AND in the deep cemented section with a consistent shape,
    since a body wave at depth has no reason to appear in a surface lead-in in
    phase with deep channels.

T3  SPECTRAL CHARACTER OF THE COMMON COMPONENT.
    The projection of the data onto u1 is a time series.  A seismic ambient
    field is broadband; instrumental artefacts commonly carry narrowband lines
    (mains harmonics, interrogator internals).  Spectral flatness and the
    strongest line-to-median ratio are reported.  A strong narrow line is
    positive evidence of an instrumental origin; its absence is not evidence
    against, since a drifting laser term is also broadband.

NO SINGLE TEST IS CONCLUSIVE and the script does not pretend otherwise: it
prints per-test evidence and a combined reading that says which way each test
points, with explicit non-conclusions where a test cannot decide.  Four false
automated verdicts have been printed in this thread already.

Output: interrogator_blame_test.{png,txt}
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, detrend, resample_poly, sosfiltfilt, welch

HERE = Path(__file__).resolve().parent
STEM = HERE / "interrogator_blame_test"
CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")

FMIN, FMAX = 5.0, 20.0
FS_COMMON = 250.0
WELLHEAD_CH = 23              # G0
N_CH = 700                    # from channel 0, so the lead-in is included
SECONDS = 30.0
N_DAYS = 6
EDIT_LO, EDIT_HI = 0.2, 5.0
SEED = 20260819


def corrected_path(p):
    return Path(str(p).replace("/data/SAFODAS1-harddrive-transfer",
                               "/data/SAFOD/SAFODAS1-harddrive-transfer"))


def edit_traces(x):
    """Interpolate outlier channels; never delete (deleting breaks uniform dx)."""
    rms = np.sqrt(np.mean(x ** 2, axis=1))
    med = np.median(rms[rms > 0]) if np.any(rms > 0) else 0.0
    bad = (rms < EDIT_LO * med) | (rms > EDIT_HI * med) | ~np.isfinite(rms)
    good = np.flatnonzero(~bad)
    if good.size < 8:
        return x, -1
    for i in np.flatnonzero(bad):
        lo, hi = good[good < i], good[good > i]
        if lo.size and hi.size:
            a, b = lo[-1], hi[0]
            w = (i - a) / (b - a)
            x[i] = (1 - w) * x[a] + w * x[b]
        else:
            x[i] = x[good[0] if not lo.size else good[-1]]
    return x, int(bad.sum())


def load_day(row, fs_target=FS_COMMON):
    f = corrected_path(row.file)
    if not f.is_file():
        return None
    with h5py.File(f, "r") as h:
        g = h["Acquisition/Raw[0]"]
        fs = float(g.attrs.get("OutputDataRate", 500.0))
        dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 1.0))
        n = int(min(SECONDS * fs, g["RawData"].shape[0]))
        x = g["RawData"][:n, :N_CH].astype(np.float32).T.astype(np.float64)
    x, dropped = edit_traces(x)
    x = np.diff(x, axis=1)
    if abs(fs - fs_target) > 1e-6:
        factor = int(round(fs / fs_target))
        if factor >= 1 and abs(fs / factor - fs_target) < 1e-6:
            x = resample_poly(x, 1, factor, axis=1)
            fs = fs_target
    x = detrend(x, axis=1)
    x = sosfiltfilt(butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos"),
                    x, axis=1)
    return dict(x=x, fs=fs, dx=dx, dropped=dropped, file=f.name)


def leading_vectors(x, rank=5):
    """First `rank` left singular vectors (spatial patterns) and their energies."""
    xc = x - x.mean(axis=1, keepdims=True)
    u, s, _ = np.linalg.svd(xc, full_matrices=False)
    frac = (s ** 2) / float(np.sum(s ** 2))
    return u[:, :rank], s[:rank], frac[:rank]


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    say("Is the interrogator to blame for the static fixed-k pattern?")
    say("  2024-25 raw HDF5 via h5py only (readFile_HDF's median=True would")
    say("  delete the quantity under test)")
    say("  band %g-%g Hz | %.0f Hz | channels 0-%d (lead-in INCLUDED) | %.0f s/day"
        % (FMIN, FMAX, FS_COMMON, N_CH - 1, SECONDS))
    say("")

    db = pd.read_csv(CSV, sep=r"\s+")
    db = db[db.nSamples == 30000].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    db = db.dropna(subset=["t"]).sort_values("t").reset_index(drop=True)
    db["day"] = db.t.dt.strftime("%Y-%m-%d")
    days = sorted(db.day.unique())
    if len(days) < 2:
        raise SystemExit("need at least two days in the manifest")
    picks = [days[int(round(i * (len(days) - 1) / (N_DAYS - 1)))]
             for i in range(N_DAYS)]
    picks = sorted(set(picks))
    say("  days sampled across the archive: %s" % ", ".join(picks))
    say("")

    loaded = []
    for d in picks:
        sub = db[db.day == d]
        rec = load_day(sub.iloc[len(sub) // 2])
        if rec is None:
            say("  %s unreadable, skipped" % d)
            continue
        u, s, frac = leading_vectors(rec["x"])
        rec.update(day=d, u=u, frac=frac)
        loaded.append(rec)
        say("  %s  %-58s u1 holds %5.2f %% of variance, %d channels edited"
            % (d, rec["file"], 100 * frac[0], rec["dropped"]))
    if len(loaded) < 2:
        raise SystemExit("fewer than two days loaded; cannot test stationarity")
    say("")

    # ---------------- T1: stationarity of the spatial pattern ----------------
    say("=== T1  stationarity of the dominant spatial pattern (the strong test) ===")
    n = len(loaded)
    same = []
    for i in range(n):
        for j in range(i + 1, n):
            c = abs(float(np.corrcoef(loaded[i]["u"][:, 0], loaded[j]["u"][:, 0])[0, 1]))
            same.append(c)
            say("  |corr(u1)| %s vs %s = %.4f" % (loaded[i]["day"], loaded[j]["day"], c))
    # control: u1 of one day against u2..u5 of another -- unrelated patterns
    ctrl = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for r in range(1, loaded[j]["u"].shape[1]):
                ctrl.append(abs(float(np.corrcoef(
                    loaded[i]["u"][:, 0], loaded[j]["u"][:, r])[0, 1])))
    same = np.array(same); ctrl = np.array(ctrl)
    say("")
    say("  cross-day  |corr(u1,u1)| : median %.4f, min %.4f, max %.4f  (n=%d)"
        % (np.median(same), same.min(), same.max(), same.size))
    say("  control    |corr(u1,u2+)|: median %.4f, 95th %.4f              (n=%d)"
        % (np.median(ctrl), np.percentile(ctrl, 95), ctrl.size))
    t1_stable = float(np.median(same)) > max(0.7, float(np.percentile(ctrl, 95)))
    t1_random = float(np.median(same)) <= float(np.percentile(ctrl, 95))
    say("  -> %s" % ("STABLE across months: points to an INSTRUMENTAL origin"
                     if t1_stable else
                     ("indistinguishable from unrelated patterns: points AWAY from a"
                      " fixed instrumental fingerprint" if t1_random else
                      "partially stable: does not decide")))
    say("")

    # ---------------- T2: presence in the surface lead-in ----------------
    say("=== T2  is the pattern present in the uncoupled surface lead-in? ===")
    say("  channel %d is the wellhead (G0); 0-%d is surface lead-in"
        % (WELLHEAD_CH, WELLHEAD_CH - 1))
    t2 = []
    for rec in loaded:
        u1 = rec["u"][:, 0]
        lead = float(np.mean(u1[:WELLHEAD_CH] ** 2))
        deep = float(np.mean(u1[WELLHEAD_CH:] ** 2))
        ratio = lead / deep if deep > 0 else np.inf
        t2.append(ratio)
        say("  %s  mean u1^2 lead-in %.3e vs deep %.3e -> ratio %.2f"
            % (rec["day"], lead, deep, ratio))
    t2 = np.array(t2)
    say("  median lead-in/deep power ratio in u1: %.2f" % np.median(t2))
    say("  NOTE: the lead-in is at the surface and records real surface ground")
    say("  motion, so a high ratio alone is NOT proof of an instrumental origin.")
    say("  This test is reported for completeness and does not decide on its own.")
    say("")

    # ---------------- T3: spectral character of the common component ----------
    say("=== T3  spectral character of the projection onto u1 ===")
    t3 = []
    for rec in loaded:
        u1 = rec["u"][:, 0]
        ts = u1 @ (rec["x"] - rec["x"].mean(axis=1, keepdims=True))
        fq, pxx = welch(ts, fs=rec["fs"], nperseg=min(2048, ts.size))
        band = (fq >= FMIN) & (fq <= FMAX)
        p = pxx[band]
        if p.size < 8:
            continue
        line = float(p.max() / np.median(p))
        flat = float(np.exp(np.mean(np.log(p + 1e-300))) / np.mean(p))
        t3.append((line, flat))
        say("  %s  strongest line / median = %6.2f | spectral flatness = %.3f"
            % (rec["day"], line, flat))
    if t3:
        lines = np.array([a for a, _ in t3])
        say("  median line-to-median ratio: %.2f" % np.median(lines))
        say("  -> %s" % ("a strong narrow line is present: positive evidence for an"
                         " instrumental artefact" if np.median(lines) > 20 else
                         "no strong narrow line; broadband. NOT evidence against an"
                         " instrumental origin (a drifting laser term is broadband"
                         " too), so this test does not decide."))
    say("")

    # ---------------- combined reading ----------------
    say("=== combined reading ===")
    say("  T1 stationarity : %s" % ("instrument" if t1_stable else
                                    ("not instrument" if t1_random else "undecided")))
    say("  T2 lead-in      : reported only, cannot decide (surface motion is real)")
    say("  T3 spectral     : %s" % ("instrument" if t3 and np.median(
        np.array([a for a, _ in t3])) > 20 else "undecided"))
    say("")
    if t1_stable:
        say("  VERDICT: the dominant spatial pattern is essentially the SAME across")
        say("  months. An earth or site noise pattern does not hold its shape over")
        say("  that span while noise conditions change; a fixed optical/electronic")
        say("  response does. The interrogator (or the fixed per-channel response")
        say("  of this installation) is the leading explanation for the static")
        say("  fixed-k contaminant, and the remedy follows directly: estimate that")
        say("  per-channel response ONCE and divide it out, which is the spatial")
        say("  calibration C2_PERMEABILITY_FOLLOWUP.md already proposes.")
    elif t1_random:
        say("  VERDICT: the dominant spatial pattern is NOT stable across months,")
        say("  so a fixed instrumental fingerprint is not supported. The static")
        say("  fixed-k energy is then more likely a site/coupling effect that")
        say("  varies with conditions, and a single static calibration would not")
        say("  remove it -- it would have to be re-estimated per interval.")
    else:
        say("  VERDICT: undecided by T1, which is the only test here able to")
        say("  decide. Do not attribute the contaminant to the interrogator on")
        say("  this evidence.")
    say("")
    say("  LIMITS: %d days x %.0f s each, one file per day, so a pattern that is")
    say("  stable within a day but drifts over weeks could read either way;")
    say("  channels 0-%d are included, which mixes the surface lead-in into u1."
        % (WELLHEAD_CH - 1))

    fig, ax = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    for rec in loaded:
        ax[0].plot(rec["u"][:, 0], lw=0.9, label=rec["day"])
    ax[0].axvline(WELLHEAD_CH, color="k", ls="--", lw=1.2, label="wellhead")
    ax[0].set(xlabel="channel", ylabel="u1 (dominant spatial pattern)",
              title="T1: is the spatial pattern the same across months?")
    ax[0].legend(fontsize=6); ax[0].grid(alpha=.3)
    ax[1].hist(ctrl, bins=30, color="0.7", alpha=.85, density=True,
               label="control: u1 vs u2+")
    ax[1].hist(same, bins=15, color="tab:red", alpha=.75, density=True,
               label="cross-day u1 vs u1")
    ax[1].set(xlabel="|correlation|", ylabel="density",
              title="T1: cross-day similarity vs unrelated patterns")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    for rec in loaded:
        u1 = rec["u"][:, 0]
        ts = u1 @ (rec["x"] - rec["x"].mean(axis=1, keepdims=True))
        fq, pxx = welch(ts, fs=rec["fs"], nperseg=min(2048, ts.size))
        m = (fq >= FMIN) & (fq <= FMAX)
        ax[2].semilogy(fq[m], pxx[m], lw=0.9, label=rec["day"])
    ax[2].set(xlabel="frequency (Hz)", ylabel="PSD of the u1 projection",
              title="T3: broadband, or narrow lines?")
    ax[2].legend(fontsize=6); ax[2].grid(alpha=.3)
    fig.savefig(str(STEM) + ".png", dpi=190)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
