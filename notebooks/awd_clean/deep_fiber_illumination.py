#!/usr/bin/env python3
"""Is the DEEP fibre illuminated?  With a genuine active-source positive control.

WHY THIS DATASET IS THE BEST CANDIDATE IN THE ARCHIVE.

Everything in the Figure 7c work so far used the 2024-2025 continuous archive on
the main-hole fibre, recorded unattended. `interrogator_and_illumination_v2.py`
and `illumination_window_scan.py` showed that field carries no net downgoing
component in the body-wave fan (11 significant windows out of 240 against 12.0
expected by chance), and `AMBIENT_CC_LITERATURE_REVIEW.md` section 4a argues the
reason is that a vertical array needs stationary-phase sources in a narrow cone
directly above the wellhead -- illumination that is local and usually cultural.

The June 2026 AWD acquisition on the **Deep** fibre is the one dataset here that
should have had exactly that:

  - ~47 h continuous, 2026-06-15T23:06 to 2026-06-17T22:20, 2,835 x 60 s files;
  - **gauge length 10.209 m, essentially the 10 m Lellouch et al. (2019) used**,
    against 16.335 m on the 2024-25 fibre;
  - 1000 Hz, 2.0419 m channel spacing, 3,200 channels = 6.5 km of fibre;
  - and, decisively, it was recorded **during a manned field campaign** in which
    989 weight drops were performed at the surface. A crew, vehicles and a
    repeating near-vertical-incidence surface source is the illumination
    condition the analysis says is required.

WHAT MAKES THIS A REAL TEST RATHER THAN A HOPE: the weight drops give a
**positive control on this very fibre**. A weight drop IS a near-vertical surface
source, so windows containing one must show a strong directional asymmetry. If
they do, the measurement is validated on this fibre, this interrogator and this
geometry -- and then the ambient windows between drops become an interpretable
test rather than an unanchored null.

    drop windows      -> positive control; expect large |A|
    ambient windows   -> the actual question; no drop within GUARD_S

If ambient windows also show asymmetry, the Deep fibre is illuminated and Figure
7c-style processing should be run on it. If only the drop windows do, then even a
manned campaign did not leave a usable ambient downgoing field, which is a much
stronger version of the negative than the 2024-25 archive alone can support.

NO DEPTH REGISTRATION EXISTS FOR THIS FIBRE, and that is handled rather than
ignored. `CLAUDE.md` records that the AWD Deep branch intersections are *fibre
coordinates, not depth*. So instead of assuming where the borehole section is,
this slides a 700 m aperture along the whole 6.5 km of fibre and reports |A| for
each position. An illuminated borehole section would appear as a localised peak.
Consequence to keep in mind when reading the numbers: |A| measures a preferred
direction ALONG THE FIBRE, and where the fibre is not near-vertical that
direction is not "downgoing".

Processing is matched to the earlier work so numbers are comparable: differentiate,
decimate to 250 Hz, detrend, 5-20 Hz, 700 m aperture, Hann in space and time,
rank-2 spatial removal (required -- a separable static pattern forces |A| = 0 by
algebra), |A| over the frozen 2500-4000 m/s fan at positive frequencies, and the
same +k/-k shuffle null.

Raw HDF5 via h5py only, never DASutils.readFile_HDF, whose median=True default
would delete the quantity under test.

Output: deep_fiber_illumination.{npz,png,txt}
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
STEM = HERE / "deep_fiber_illumination"
DEEP = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/"
            "01_--_recording_2026-06-15T230629Z_--_active_source")
MANIFEST = HERE / "awd_manifest.csv"

FMIN, FMAX = 5.0, 20.0
FS_COMMON = 250.0
APERTURE_M = 700.0
WIN_S = 2.0
SECONDS = 20.0                 # per analysis window
GUARD_S = 45.0                 # ambient windows must be this far from any drop
V_LO, V_HI = 2500.0, 4000.0
RANK_REMOVE = 2
EDIT_LO, EDIT_HI = 0.2, 5.0
NULL_COUNT = 300
ALPHA = 0.05
N_DROP = int(os.environ.get("DEEP_N_DROP", "12"))
N_AMBIENT = int(os.environ.get("DEEP_N_AMBIENT", "24"))
APERTURE_STEP = int(os.environ.get("DEEP_APERTURE_STEP", "3"))
SEED = 20260819


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
        else:
            return None, fs, dropped
    x = detrend(x, axis=1)
    x = sosfiltfilt(butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos"),
                    x, axis=1)
    return x, fs, dropped


def remove_rank(x, rank):
    """eigh on the spatial covariance: same left singular vectors, far more robust."""
    if rank <= 0:
        return x
    xc = x - x.mean(axis=1, keepdims=True)
    c = xc @ xc.T
    c = 0.5 * (c + c.T)
    if not np.all(np.isfinite(c)):
        raise FloatingPointError("non-finite covariance")
    _, v = np.linalg.eigh(c)
    u = v[:, -rank:]
    return xc - u @ (u.T @ xc)


def asymmetry(x, fs, dx, rng):
    nw = int(WIN_S * fs)
    nch = x.shape[0]
    if x.shape[1] < nw or nch < 16:
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
    pool = np.concatenate([P[up].ravel(), P[dn].ravel()])
    n_up = int(up.sum())
    nulls = np.empty(NULL_COUNT)
    for i in range(NULL_COUNT):
        s = rng.permutation(pool)
        x_, y_ = float(s[:n_up].sum()), float(s[n_up:].sum())
        nulls[i] = abs(x_ - y_) / (x_ + y_) if (x_ + y_) > 0 else 0.0
    return dict(asym=obs, signed=(a - b) / (a + b),
                p=float((np.sum(nulls >= obs) + 1) / (NULL_COUNT + 1)),
                null95=float(np.percentile(nulls, 95)),
                fan_share=100.0 * (a + b) / float(P[pos_f].sum()))


def scan_apertures(x, fs, dx, rng, say, label):
    """Slide a 700 m aperture along the fibre; return the per-position results."""
    nch = int(round(APERTURE_M / dx))
    out = []
    starts = list(range(0, x.shape[0] - nch + 1, max(1, nch // APERTURE_STEP)))
    for c0 in starts:
        seg = x[c0:c0 + nch]
        try:
            seg = remove_rank(seg, RANK_REMOVE)
        except (np.linalg.LinAlgError, FloatingPointError):
            continue
        res = asymmetry(seg, fs, dx, rng)
        if res is None:
            continue
        out.append(dict(ch0=c0, fibre_m=c0 * dx, **res))
    return out


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    rng = np.random.default_rng(SEED)

    say("DEEP fibre illumination test, with an active-source positive control")
    say("  band %g-%g Hz | %.0f Hz | %.0f m aperture | %.0f s windows"
        % (FMIN, FMAX, FS_COMMON, APERTURE_M, SECONDS))
    say("  rank-%d spatial removal | fan %.0f-%.0f m/s | positive f only"
        % (RANK_REMOVE, V_LO, V_HI))
    say("  ambient windows must be >= %.0f s from any weight drop" % GUARD_S)
    say("")
    say("  NOTE: no depth registration exists for this fibre (AWD Deep branch")
    say("  intersections are fibre coordinates, not depth), so a %.0f m aperture" % APERTURE_M)
    say("  is slid along the whole fibre and |A| reported per position. Where the")
    say("  fibre is not near-vertical, 'along fibre' is not 'downgoing'.")
    say("")

    if not MANIFEST.is_file():
        raise SystemExit("missing %s (build_manifest.py)" % MANIFEST)
    man = pd.read_csv(MANIFEST)
    man = man[man.deep_available == 1].copy()
    man["t"] = pd.to_datetime(man.utc_time, utc=True, errors="coerce")
    man = man.dropna(subset=["t"]).sort_values("t").reset_index(drop=True)
    say("  %d drops with Deep coverage, %s to %s"
        % (len(man), man.t.iloc[0], man.t.iloc[-1]))

    files = sorted(DEEP.glob("SAFOD-Deep-*.h5"))
    if not files:
        raise SystemExit("no Deep h5 files under %s" % DEEP)
    import re
    def file_time(p):
        m = re.search(r"_(\d{4}-\d{2}-\d{2}T\d{6})Z\.h5$", p.name)
        return pd.to_datetime(m.group(1), format="%Y-%m-%dT%H%M%S", utc=True) if m else None
    ftimes = [(p, file_time(p)) for p in files]
    ftimes = [(p, t) for p, t in ftimes if t is not None]
    say("  %d Deep files, %s to %s" % (len(ftimes), ftimes[0][1], ftimes[-1][1]))
    say("")

    drops = man.t.to_numpy()

    def read_window(path, offset_s, fs_hint=1000.0):
        with h5py.File(path, "r") as h:
            g = h["Acquisition/Raw[0]"]
            fs = float(g.attrs.get("OutputDataRate", fs_hint))
            dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 2.0419))
            n_tot = g["RawData"].shape[0]
            lo = int(max(0, min(n_tot - int(SECONDS * fs), offset_s * fs)))
            hi = int(min(n_tot, lo + int(SECONDS * fs)))
            arr = g["RawData"][lo:hi, :].astype(np.float32).T
        return arr, fs, dx

    results = {"drop": [], "ambient": []}

    # ---- positive control: windows containing a weight drop ----
    say("=== positive control: windows containing a weight drop ===")
    picks = man.iloc[np.linspace(0, len(man) - 1, N_DROP).astype(int)]
    for _, row in picks.iterrows():
        path = DEEP / str(row.deep_file)
        if not path.is_file():
            continue
        off = max(0.0, float(row.deep_offset_s) - 2.0)
        try:
            arr, fs, dx = read_window(path, off)
        except Exception as exc:
            say("  %s unreadable: %s" % (path.name, exc))
            continue
        x, fs2, _ = prepare(arr, fs)
        if x is None:
            continue
        scan = scan_apertures(x, fs2, dx, rng, say, "drop")
        if not scan:
            continue
        best = max(scan, key=lambda r: r["asym"])
        results["drop"].append(dict(t=row.t, best=best, scan=scan))
        say("  %s  best |A| %.4f at fibre %.0f m (p %.4f)"
            % (row.t, best["asym"], best["fibre_m"], best["p"]))

    # ---- the actual test: ambient windows far from any drop ----
    say("")
    say("=== ambient windows, >= %.0f s from any drop ===" % GUARD_S)
    cand = []
    for p, t in ftimes:
        mid = t + pd.Timedelta(seconds=30)
        if drops.size:
            gap = np.min(np.abs((drops - np.datetime64(mid)).astype("timedelta64[s]").astype(float)))
        else:
            gap = 1e9
        if gap >= GUARD_S:
            cand.append((p, t, gap))
    say("  %d of %d files have no drop within %.0f s of their midpoint"
        % (len(cand), len(ftimes), GUARD_S))
    if not cand:
        raise SystemExit("no ambient windows available")
    sel = [cand[i] for i in np.linspace(0, len(cand) - 1, min(N_AMBIENT, len(cand))).astype(int)]
    for p, t, gap in sel:
        try:
            arr, fs, dx = read_window(p, 20.0)
        except Exception as exc:
            say("  %s unreadable: %s" % (p.name, exc))
            continue
        x, fs2, _ = prepare(arr, fs)
        if x is None:
            continue
        scan = scan_apertures(x, fs2, dx, rng, say, "ambient")
        if not scan:
            continue
        best = max(scan, key=lambda r: r["asym"])
        results["ambient"].append(dict(t=t, best=best, scan=scan))
        say("  %s  best |A| %.4f at fibre %.0f m (p %.4f, nearest drop %.0f s)"
            % (t, best["asym"], best["fibre_m"], best["p"], gap))

    # ---------------- reading ----------------
    say("")
    say("=== reading ===")
    for key in ("drop", "ambient"):
        rs = results[key]
        if not rs:
            say("  %-8s no windows measured" % key)
            continue
        a = np.array([r["best"]["asym"] for r in rs])
        pv = np.array([r["best"]["p"] for r in rs])
        n_sig = int(np.sum(pv < ALPHA))
        say("  %-8s n=%2d | best |A| median %.4f max %.4f | %d of %d windows p<%.2f"
            % (key, len(rs), np.median(a), a.max(), n_sig, len(rs), ALPHA))

    d, amb = results["drop"], results["ambient"]
    ok_control = bool(d) and np.mean([r["best"]["p"] < ALPHA for r in d]) >= 0.5
    say("")
    if not d:
        say("  NO POSITIVE CONTROL: no drop window measured, so an ambient null")
        say("  here would be uninterpretable. No conclusion.")
    elif not ok_control:
        say("  THE POSITIVE CONTROL FAILED. Fewer than half of the windows that")
        say("  CONTAIN a weight drop show significant fan asymmetry. A weight drop")
        say("  is a near-vertical surface source, so the measurement should detect")
        say("  it; that it does not means this configuration lacks sensitivity on")
        say("  this fibre -- most likely the 2500-4000 m/s fan or the 700 m aperture")
        say("  is wrong for this geometry, which has no depth registration. NO")
        say("  conclusion about ambient illumination may be drawn.")
    else:
        n_amb_sig = int(np.sum([r["best"]["p"] < ALPHA for r in amb])) if amb else 0
        say("  Positive control PASSES: the drop windows show the asymmetry a")
        say("  near-vertical surface source must produce, so the measurement works")
        say("  on this fibre, this interrogator and this geometry.")
        say("")
        if amb and n_amb_sig > max(1, int(0.2 * len(amb))):
            say("  AMBIENT WINDOWS ARE ALSO ILLUMINATED: %d of %d significant."
                % (n_amb_sig, len(amb)))
            say("  This is the first positive illumination result in this project and")
            say("  it makes the Deep fibre the right target for Figure 7c-style")
            say("  processing. NEXT: run the exact-stack geometry on these windows,")
            say("  restricted to the fibre positions where |A| peaks.")
        else:
            say("  AMBIENT WINDOWS ARE NOT ILLUMINATED: %d of %d significant."
                % (n_amb_sig, len(amb) if amb else 0))
            say("  With a working positive control on the same fibre, this is a")
            say("  stronger negative than the 2024-25 archive alone supports: even")
            say("  during a manned field campaign, the wavefield BETWEEN drops")
            say("  carries no net directional energy in the body-wave fan. The")
            say("  illumination that makes borehole ambient body-wave interferometry")
            say("  work is supplied by the active source itself, not by the ambient")
            say("  field at this site.")
    say("")
    say("  LIMITS: %.0f s windows; %d drop and %d ambient windows; the fan and"
        % (SECONDS, len(d), len(amb)))
    say("  aperture are the frozen selections from the main-hole work and may not")
    say("  suit this fibre; no depth registration, so fibre position is not depth;")
    say("  rank-%d removal is conservative against a low-rank arrival." % RANK_REMOVE)

    # ---------------- figure ----------------
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6), constrained_layout=True)
    for key, col in (("drop", "#D55E00"), ("ambient", "#0072B2")):
        for r in results[key]:
            s = r["scan"]
            ax[0].plot([q["fibre_m"] for q in s], [q["asym"] for q in s],
                       "-", color=col, alpha=.35, lw=0.9)
        if results[key]:
            ax[0].plot([], [], "-", color=col, label="%s windows" % key)
    ax[0].set(xlabel="position along fibre (m)", ylabel="|A| in the fan",
              title="|A| vs fibre position\n(700 m aperture slid along 6.5 km)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)

    data = [[r["best"]["asym"] for r in results[k]] for k in ("drop", "ambient")]
    if any(data):
        ax[1].boxplot([d_ for d_ in data if d_], labels=[k for k, d_ in
                      zip(("drop", "ambient"), data) if d_])
        ax[1].set(ylabel="best |A| per window",
                  title="Positive control vs ambient")
        ax[1].grid(alpha=.3, axis="y")
    for key, col in (("drop", "#D55E00"), ("ambient", "#0072B2")):
        rs = results[key]
        if not rs:
            continue
        ax[2].semilogy([r["t"] for r in rs], [r["best"]["p"] for r in rs],
                       "o", color=col, ms=5, label=key)
    ax[2].axhline(ALPHA, color="k", ls="--", lw=1.2, label="alpha = %.2f" % ALPHA)
    ax[2].set(xlabel="time (UTC)", ylabel="p of best aperture",
              title="Significance over the survey")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
    ax[2].tick_params(axis="x", rotation=30, labelsize=7)
    fig.savefig(str(STEM) + ".png", dpi=190)

    np.savez(str(STEM) + ".npz",
             drop_asym=np.array([r["best"]["asym"] for r in d]),
             drop_p=np.array([r["best"]["p"] for r in d]),
             drop_fibre_m=np.array([r["best"]["fibre_m"] for r in d]),
             amb_asym=np.array([r["best"]["asym"] for r in amb]),
             amb_p=np.array([r["best"]["p"] for r in amb]),
             amb_fibre_m=np.array([r["best"]["fibre_m"] for r in amb]),
             aperture_m=APERTURE_M, rank_removed=RANK_REMOVE, guard_s=GUARD_S,
             alpha=ALPHA, fan=(V_LO, V_HI))
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
