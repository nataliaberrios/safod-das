#!/usr/bin/env python3
"""Scan the archive for ILLUMINATED windows, instead of stacking everything.

THE RECOVERY IDEA, and why it is the one route never tried.

interrogator_and_illumination_v2.py established, with a working positive control,
that the ingredient Lellouch et al. (2019) relied on is missing from our data:
the downgoing/upgoing asymmetry in the body-wave fan that they cite as evidence
for surface sources is significant in their 2017 records (|A| = 0.348, p = 0.005)
and absent from a 2024-25 window (|A| = 0.04, p = 0.73, at every rank tested).

But that was ONE 30 s window from ONE day, and illumination by surface activity
is anthropogenic: it varies by hour of day and day of week. Meanwhile Behm (2016,
Geophys. Prospect. 10.1111/1365-2478.12424) found that under good illumination
"ambient noise from time periods as short as 30 seconds is sufficient to obtain
robust interferograms".

Put those together and our whole stacking strategy looks wrong. Every stack in
this tree pools indiscriminately -- 6 days, then a coherent 96 h stack -- which
would DILUTE an illuminated minority of windows into a majority of unilluminated
ones. That is consistent with the otherwise strange observation that our results
got WORSE with more data (per-day min p = 0.1345, 96 h coherent stack p = 0.9184).
If even a few per cent of windows are illuminated, selecting them and correlating
only those is a completely different experiment from stacking all of them.

WHAT THIS DOES
  1. Samples many 30 s windows spread across days AND hours of day.
  2. For each, measures |A| = |E(+k) - E(-k)| / (E(+k) + E(-k)) in the frozen
     2500-4000 m/s fan at positive frequencies, after projecting out the leading
     RANK_REMOVE spatial patterns (the static separable pattern would otherwise
     force |A| to zero by algebra -- see interrogator_and_illumination_v2.py).
  3. Nulls each window by shuffling fan cell powers between the +k and -k halves.
  4. Reports how many windows clear their own null, and whether those windows
     cluster by hour of day or day of week -- which is the signature of
     anthropogenic surface sources rather than of estimation noise.

DECISION RULE, FIXED BEFORE THE RUN.  With N windows tested at alpha = 0.05,
about 5 % clear by chance. Illumination is claimed only if the observed count
exceeds the 95th percentile of Binomial(N, 0.05), reported explicitly. Clustering
in hour of day is reported as supporting evidence but is NOT part of the gate.

If illuminated windows exist, the follow-up is defined in advance: rerun
ambient_lellouch2019_exact_stack.py over ONLY those windows and compare against
an equal number of randomly chosen unilluminated windows, so the comparison is
matched on window count and the selection is testable rather than assumed.

If none exist, the ambient route is closed on this archive at 5-20 Hz, and the
honest write-up is that negative plus the per-velocity result (0 of 181).

Reads raw HDF5 with h5py only. Never DASutils.readFile_HDF, whose median=True
default would delete the k = 0 structure this measurement depends on.

Output: illumination_window_scan.{npz,png,txt}
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
from scipy.stats import binom

HERE = Path(__file__).resolve().parent
STEM = HERE / "illumination_window_scan"
CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")

FMIN, FMAX = 5.0, 20.0
FS_COMMON = 250.0
WELLHEAD_CH = 23
APERTURE_CH = 686
SECONDS = 30.0
WIN_S = 2.0
V_LO, V_HI = 2500.0, 4000.0
RANK_REMOVE = 2
EDIT_LO, EDIT_HI = 0.2, 5.0
NULL_COUNT = 300
ALPHA = 0.05
N_WINDOWS = int(os.environ.get("SCAN_WINDOWS", "240"))
SEED = 20260819


def corrected_path(p):
    return Path(str(p).replace("/data/SAFODAS1-harddrive-transfer",
                               "/data/SAFOD/SAFODAS1-harddrive-transfer"))


def edit_traces(x):
    rms = np.sqrt(np.mean(x ** 2, axis=1))
    med = np.median(rms[rms > 0]) if np.any(rms > 0) else 0.0
    bad = (rms < EDIT_LO * med) | (rms > EDIT_HI * med) | ~np.isfinite(rms)
    good = np.flatnonzero(~bad)
    if good.size < 8:
        return None, -1
    for i in np.flatnonzero(bad):
        lo, hi = good[good < i], good[good > i]
        if lo.size and hi.size:
            a, b = lo[-1], hi[0]
            w = (i - a) / (b - a)
            x[i] = (1 - w) * x[a] + w * x[b]
        else:
            x[i] = x[good[0] if not lo.size else good[-1]]
    return x, int(bad.sum())


def prepare(x, fs):
    x, dropped = edit_traces(np.asarray(x, dtype=np.float64))
    if x is None:
        return None, fs, dropped
    x = np.diff(x, axis=1)
    if abs(fs - FS_COMMON) > 1e-6:
        factor = int(round(fs / FS_COMMON))
        if factor >= 1 and abs(fs / factor - FS_COMMON) < 1e-6:
            x = resample_poly(x, 1, factor, axis=1)
            fs = FS_COMMON
    x = detrend(x, axis=1)
    x = sosfiltfilt(butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos"),
                    x, axis=1)
    return x, fs, dropped


def remove_rank(x, rank):
    if rank <= 0:
        return x
    xc = x - x.mean(axis=1, keepdims=True)
    u, _, _ = np.linalg.svd(xc, full_matrices=False)
    u = u[:, :rank]
    return xc - u @ (u.T @ xc)


def asymmetry(x, fs, dx, rng):
    nw = int(WIN_S * fs)
    nch = x.shape[0]
    if x.shape[1] < nw:
        return None
    wt, wx = np.hanning(nw)[None, :], np.hanning(nch)[:, None]
    P = None
    for s in range(0, x.shape[1] - nw + 1, nw):
        F = np.fft.fftshift(np.fft.fft2(x[:, s:s + nw] * wt * wx))
        P = np.abs(F) ** 2 if P is None else P + np.abs(F) ** 2
    k = np.fft.fftshift(np.fft.fftfreq(nch, dx))
    f = np.fft.fftshift(np.fft.fftfreq(nw, 1.0 / fs))
    K, Fm = np.meshgrid(k, f, indexing="ij")
    pos_f = (Fm >= FMIN) & (Fm <= FMAX)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(np.abs(K) > 0, Fm / np.abs(K), np.inf)
    fan = pos_f & (v >= V_LO) & (v <= V_HI)
    up, dn = fan & (K > 0), fan & (K < 0)
    if up.sum() < 4 or dn.sum() < 4:
        return None
    a, b = float(P[up].sum()), float(P[dn].sum())
    if a + b <= 0:
        return None
    obs = abs(a - b) / (a + b)
    signed = (a - b) / (a + b)
    pool = np.concatenate([P[up].ravel(), P[dn].ravel()])
    n_up = int(up.sum())
    nulls = np.empty(NULL_COUNT)
    for i in range(NULL_COUNT):
        s = rng.permutation(pool)
        x_, y_ = float(s[:n_up].sum()), float(s[n_up:].sum())
        nulls[i] = abs(x_ - y_) / (x_ + y_) if (x_ + y_) > 0 else 0.0
    return dict(asym=obs, signed=signed,
                p=float((np.sum(nulls >= obs) + 1) / (NULL_COUNT + 1)),
                null95=float(np.percentile(nulls, 95)),
                fan_share=100.0 * (a + b) / float(P[pos_f].sum()))


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    rng = np.random.default_rng(SEED)

    say("Scanning the archive for ILLUMINATED windows")
    say("  |A| in the %.0f-%.0f m/s fan, positive f, after rank-%d spatial removal"
        % (V_LO, V_HI, RANK_REMOVE))
    say("  channels %d-%d | %.0f s windows | %d windows requested"
        % (WELLHEAD_CH, WELLHEAD_CH + APERTURE_CH - 1, SECONDS, N_WINDOWS))
    say("  DECISION RULE fixed before the run: illumination is claimed only if the")
    say("  count of windows with p < %.2f exceeds the 95th percentile of" % ALPHA)
    say("  Binomial(N, %.2f). Hour-of-day clustering is supporting, not gating." % ALPHA)
    say("")

    db = pd.read_csv(CSV, sep=r"\s+")
    db = db[db.nSamples == 30000].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    db = db.dropna(subset=["t"]).sort_values("t").reset_index(drop=True)
    if db.empty:
        raise SystemExit("no continuous files in the manifest")

    # spread the sample over days AND hours so a diurnal pattern is visible
    idx = np.unique(np.linspace(0, len(db) - 1, N_WINDOWS).astype(int))
    say("  manifest spans %s to %s (%d continuous files)"
        % (db.t.iloc[0], db.t.iloc[-1], len(db)))
    say("")

    rows = []
    for n, i in enumerate(idx):
        row = db.iloc[int(i)]
        f = corrected_path(row.file)
        if not f.is_file():
            continue
        try:
            with h5py.File(f, "r") as h:
                g = h["Acquisition/Raw[0]"]
                fs = float(g.attrs.get("OutputDataRate", 500.0))
                dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 1.0))
                ns = int(min(SECONDS * fs, g["RawData"].shape[0]))
                hi = WELLHEAD_CH + APERTURE_CH
                x = g["RawData"][:ns, WELLHEAD_CH:hi].astype(np.float32).T
        except Exception as exc:
            say("  [%3d] %s unreadable: %s" % (n, f.name, exc))
            continue
        if x.shape[0] < APERTURE_CH:
            continue
        x, fs2, dropped = prepare(x, fs)
        if x is None:
            continue
        res = asymmetry(remove_rank(x, RANK_REMOVE), fs2, dx, rng)
        if res is None:
            continue
        rows.append(dict(t=row.t, hour=int(row.t.hour), dow=int(row.t.dayofweek),
                         dropped=dropped, **res))
        if n % 20 == 0:
            say("  [%3d/%d] %s |A| %.4f p %.4f" % (n, len(idx), row.t, res["asym"], res["p"]))

    if not rows:
        raise SystemExit("no windows measured")
    df = pd.DataFrame(rows)
    N = len(df)
    hits = df[df.p < ALPHA]
    thresh = int(binom.ppf(0.95, N, ALPHA))
    say("")
    say("=== result ===")
    say("  windows measured            : %d" % N)
    say("  |A| median %.4f, 90th %.4f, max %.4f"
        % (df.asym.median(), df.asym.quantile(0.90), df.asym.max()))
    say("  windows with p < %.2f        : %d" % (ALPHA, len(hits)))
    say("  expected by chance           : %.1f (95th percentile of Binomial = %d)"
        % (N * ALPHA, thresh))
    say("")

    detected = len(hits) > thresh
    if detected:
        say("  ILLUMINATED WINDOWS EXIST: %d observed against a chance ceiling of %d."
            % (len(hits), thresh))
        say("  Strongest windows:")
        for _, r in hits.nlargest(min(10, len(hits)), "asym").iterrows():
            say("    %s  |A| %.4f  p %.4f  signed %+.4f  fan %.2f%%"
                % (r.t, r.asym, r.p, r.signed, r.fan_share))
        same_sign = float(np.mean(np.sign(hits.signed) == np.sign(hits.signed.iloc[0])))
        say("")
        say("  directional consistency: %.0f %% of hits share one sign" % (100 * same_sign))
        say("  (a real one-sided illumination should be CONSISTENT in sign; a mix")
        say("   of both signs is more consistent with estimation noise)")
        say("")
        say("  hour-of-day distribution of hits (UTC):")
        for h in range(0, 24, 3):
            c = int(((hits.hour >= h) & (hits.hour < h + 3)).sum())
            a = int(((df.hour >= h) & (df.hour < h + 3)).sum())
            say("    %02d-%02d h  %3d hits / %3d windows%s"
                % (h, h + 3, c, a, "  " + "#" * c if c else ""))
        say("")
        say("  NEXT STEP, predefined: rerun ambient_lellouch2019_exact_stack.py over")
        say("  ONLY these windows and compare against an equal number of randomly")
        say("  chosen non-hit windows, matched on count so the selection itself is")
        say("  testable. Behm (2016) obtained robust interferograms from 30 s under")
        say("  good illumination, so %d windows is not obviously too few." % len(hits))
    else:
        say("  NO ILLUMINATED WINDOWS. %d windows clear p < %.2f against a chance"
            % (len(hits), ALPHA))
        say("  ceiling of %d, so the count is consistent with pure chance." % thresh)
        say("")
        say("  This closes the ambient route on this archive at %g-%g Hz. The field"
            % (FMIN, FMAX))
        say("  carries no net directional preference in the body-wave fan in ANY")
        say("  window sampled, and Lellouch's downgoing P is an inference FROM that")
        say("  asymmetry, so there is no arrival of that kind to recover. This is a")
        say("  property of the recording, not of the processing: the same")
        say("  measurement finds the asymmetry in his 2017 data (|A| = 0.348,")
        say("  p = 0.005) from only ~5 s.")
        say("")
        say("  The honest write-up is this negative plus the per-velocity result")
        say("  (0 of 181 velocities clear the per-velocity null). The velocity model")
        say("  itself remains reachable via earthquakes -- Lellouch's Figure 9 was")
        say("  mainly earthquake-derived, 206 events are cached, and G0 already")
        say("  matched the 2005 check shot to 0.2 %.")
    say("")
    say("  LIMITS: %.0f s per window and one window per sampled file, so a burst"
        % SECONDS)
    say("  shorter than %.0f s could be missed; rank-%d removal is conservative"
        % (SECONDS, RANK_REMOVE))
    say("  against a low-rank plane-wave arrival; the fan is the frozen")
    say("  %.0f-%.0f m/s selection; and |A| is along-fibre, near-vertical but not"
        % (V_LO, V_HI))
    say("  exactly vertical in this borehole.")

    fig, ax = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    ax[0].plot(df.t, df.asym, ".", ms=3, color="0.5", label="all windows")
    if len(hits):
        ax[0].plot(hits.t, hits.asym, "o", ms=5, color="crimson", label="p < %.2f" % ALPHA)
    ax[0].set(xlabel="time (UTC)", ylabel="|A| in the fan",
              title="Illumination across the archive")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    ax[0].tick_params(axis="x", rotation=30, labelsize=7)
    ax[1].hist(df.asym, bins=40, color="steelblue", alpha=.85, label="observed |A|")
    ax[1].axvline(df.null95.median(), color="crimson", ls="--",
                  label="median null 95th")
    ax[1].set(xlabel="|A|", ylabel="windows", title="Distribution vs null")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    hr = df.groupby("hour").asym.median()
    ax[2].bar(hr.index, hr.values, color="tab:green", alpha=.85)
    ax[2].set(xlabel="hour of day (UTC)", ylabel="median |A|",
              title="Diurnal pattern\n(surface activity would show one)")
    ax[2].grid(alpha=.3, axis="y")
    fig.savefig(str(STEM) + ".png", dpi=190)
    np.savez(str(STEM) + ".npz",
             t=df.t.astype("int64").to_numpy(), asym=df.asym.to_numpy(),
             signed=df.signed.to_numpy(), p=df.p.to_numpy(),
             null95=df.null95.to_numpy(), hour=df.hour.to_numpy(),
             dow=df.dow.to_numpy(), fan_share=df.fan_share.to_numpy(),
             n_hits=len(hits), chance_ceiling=thresh, alpha=ALPHA,
             rank_removed=RANK_REMOVE)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
