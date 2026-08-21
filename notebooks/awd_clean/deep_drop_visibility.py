#!/usr/bin/env python3
"""Is a weight drop VISIBLE on the Deep fibre, and at what apparent velocity?

WHY THIS EXISTS. `deep_fiber_illumination.py` (v1, v2) asks whether the Deep
fibre carries a directional body-wave field, using the 989 weight drops as a
positive control. Both versions returned a weak or physically impossible
control, and both then argued about the statistic -- the rank removed, the null,
the correction for selecting over aperture positions. That is arguing downstream
of an unchecked premise. Neither version ever established that

    (a) a weight drop produces any measurable amplitude on this fibre at all, or
    (b) the frozen 2500-4000 m/s fan is where this fibre's energy actually is.

The manifest was built by intersecting weight-drop GPS times with Deep file
coverage. Coverage is not detection. This script tests (a) and (b) directly, with
no velocity fan, no aperture selection, no rank removal and no permutation null
in the detection step -- so it cannot be defeated by any of the choices that the
illumination test was being blamed on.

PART A -- DETECTABILITY, MODEL FREE.
For N drops spread across the survey: read +-PRE/POST s around the GPS drop time
with h5py directly (never DASutils.readFile_HDF: its median=True default removes
the per-sample cross-channel median and would delete common-mode energy), convert
to strain rate, decimate to 250 Hz, band-pass 5-20 Hz, take the analytic
envelope, and divide each channel by ITS OWN pre-drop RMS. That per-channel
normalisation is what makes the result independent of the static amplitude
response that `C2_PHASE0_RESULTS.md` records as uncalibrated on this array, and it
neutralises dead or screaming channels without an editing step.

Stacking those SNR envelopes over drops gives an image in (fibre position, time
since drop). A visible drop is a bright feature at t > 0. The control for that
image is a SHAM arm: the identical procedure at random times with no drop within
GUARD_S, with a fake drop time at the same offset in the window. Sham gives the
level the stack reaches from noise alone, so "visible" means visible ABOVE SHAM,
not above 1.

PART B -- WHERE THE ENERGY IS IN VELOCITY, SIGNED.
The same windows, but binned by signed apparent velocity v = f/k along the fibre
instead of tested inside a pre-chosen fan. Two aperture lengths are used on
purpose. 700 m is the frozen main-hole aperture; at 5-20 Hz it gives
dk = 1/(343*2.0419) = 0.00143 /m, so the whole 2500-4000 m/s fan is only a couple
of wavenumber bins wide -- the same under-resolution that
"Why no filter can work: contaminant and target are unresolved at this aperture"
(commit 8fa2266) recorded for the main hole. 2800 m is included so the fan is
resolved and we can see whether 2500-4000 m/s was ever the right place to look.

WHAT THIS SCRIPT CAN AND CANNOT SETTLE. It can settle whether the drops are in
the data and what apparent velocities they excite, which is a prerequisite for
any illumination claim. It cannot settle illumination itself: |A| in a fan is a
directional statistic and needs its own null.

DEPTH REGISTRATION, added 2026-08-20. When this script was first written none
existed, and CLAUDE.md's warning that AWD Deep branch intersections are FIBRE
COORDINATES, NOT DEPTH applied to every number it produced.
`safod_geometry.py` now supplies the mapping from
`SAFOD_Phase2_GeoReferenced_Channels.xlsx`. Positions are still reported as
metres along fibre so runs stay comparable, but each responding section also
carries its TVD and limb, and sections lying in the SURFACE LEAD-IN -- fibre
that never enters the hole -- are labelled and kept out of the in-hole summary.
This matters here: the strongest section of the first run began inside that
lead-in. Note also that "along fibre" is only "downgoing" where the fibre is
near-vertical, which the mapping now makes checkable rather than unknown.

Output: deep_drop_visibility.{npz,png,txt}
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, detrend, hilbert, resample_poly, sosfiltfilt

HERE = Path(__file__).resolve().parent
STEM = HERE / "deep_drop_visibility"
DEEP = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/"
            "01_--_recording_2026-06-15T230629Z_--_active_source")
MANIFEST = HERE / "awd_manifest.csv"

FMIN, FMAX = 5.0, 20.0
FS_COMMON = 250.0
PRE_S = 1.5                     # read this much before the drop
POST_S = 4.0                    # and this much after
PAD_S = 0.5                     # extra either side, trimmed after filtering
NOISE_LO, NOISE_HI = -1.4, -0.3  # per-channel reference RMS window, s from drop
N_DROP = int(os.environ.get("VIS_N_DROP", "30"))
N_SHAM = int(os.environ.get("VIS_N_SHAM", "30"))
GUARD_S = 45.0
# Part B
FK_WIN_S = 2.0
FK_LEAD_S = 0.2                 # f-k window starts this long before the drop
APERTURES_M = (700.0, 2800.0)
N_POS = int(os.environ.get("VIS_N_POS", "8"))
V_EDGES = np.geomspace(300.0, 20000.0, 25)
SEED = 20260819


import safod_geometry as GEO


def band(x, fs):
    """Strain rate -> 250 Hz -> detrend -> 5-20 Hz.  x is (channels, time)."""
    x = np.diff(np.asarray(x, dtype=np.float64), axis=1)
    if abs(fs - FS_COMMON) > 1e-6:
        factor = int(round(fs / FS_COMMON))
        if factor < 1 or abs(fs / factor - FS_COMMON) > 1e-6:
            return None, fs
        x = resample_poly(x, 1, factor, axis=1)
        fs = FS_COMMON
    x = detrend(x, axis=1)
    x = sosfiltfilt(butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos"),
                    x, axis=1)
    return x, fs


def signed_velocity_spectrum(x, fs, dx, aperture_m, n_pos):
    """Sum band power into signed apparent-velocity bins, per aperture position.

    Returns (positions_m, P_plus, P_minus) with P_* of shape
    (n_positions, len(V_EDGES) - 1).  No fan, no selection: every cell in
    5-20 Hz with k != 0 lands in the bin its own f/|k| puts it in.
    """
    nch = int(round(aperture_m / dx))
    nw = int(FK_WIN_S * fs)
    if x.shape[0] < nch or x.shape[1] < nw:
        return None
    starts = np.unique(np.linspace(0, x.shape[0] - nch, n_pos).astype(int))
    k = np.fft.fftshift(np.fft.fftfreq(nch, dx))
    f = np.fft.fftshift(np.fft.fftfreq(nw, 1.0 / fs))
    K, Fm = np.meshgrid(k, f, indexing="ij")
    inband = (Fm >= FMIN) & (Fm <= FMAX) & (np.abs(K) > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        V = np.where(inband, Fm / np.abs(K), np.nan)
    nb = len(V_EDGES) - 1
    ib = np.digitize(V, V_EDGES) - 1
    up = inband & (K > 0)
    dn = inband & (K < 0)
    wt, wx = np.hanning(nw)[None, :], np.hanning(nch)[:, None]
    Pp = np.zeros((len(starts), nb))
    Pm = np.zeros((len(starts), nb))
    for i, c0 in enumerate(starts):
        seg = x[c0:c0 + nch, :nw]
        seg = seg - seg.mean(axis=1, keepdims=True)
        P = np.abs(np.fft.fftshift(np.fft.fft2(seg * wt * wx))) ** 2
        for b in range(nb):
            m = ib == b
            Pp[i, b] = P[m & up].sum()
            Pm[i, b] = P[m & dn].sum()
    return starts * dx, Pp, Pm


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    rng = np.random.default_rng(SEED)

    say("Deep fibre: is a weight drop VISIBLE, and at what apparent velocity?")
    say("  h5py direct (never readFile_HDF -- median=True would delete common mode)")
    say("  strain rate | %.0f Hz | %g-%g Hz | envelope / per-channel pre-drop RMS"
        % (FS_COMMON, FMIN, FMAX))
    say("  window %.1f s before to %.1f s after the GPS drop time" % (PRE_S, POST_S))
    say("  sham arm = identical processing at times with no drop within %.0f s" % GUARD_S)
    say("")
    say("  Depth registration NOW EXISTS (safod_geometry.py, from")
    say("  SAFOD_Phase2_GeoReferenced_Channels.xlsx, 2026-08-20). Positions are")
    say("  still quoted as metres along fibre for continuity with earlier runs,")
    say("  but responding sections are also given in TVD, and any section in the")
    say("  SURFACE LEAD-IN -- fibre that never enters the hole -- is labelled and")
    say("  excluded from the in-hole summary.")
    say("")

    man = pd.read_csv(MANIFEST)
    man = man[man.deep_available == 1].copy()
    man["t"] = pd.to_datetime(man.utc_time, utc=True, errors="coerce")
    man = man.dropna(subset=["t"]).sort_values("t").reset_index(drop=True)
    # the whole window must fit inside the 60 s file, so no clamping is needed;
    # clamping would silently misalign the stack, which is the failure this avoids
    fits = (man.deep_offset_s >= PRE_S + PAD_S + 1.0) & \
           (man.deep_offset_s <= 60.0 - POST_S - PAD_S - 1.0)
    man = man[fits].reset_index(drop=True)
    say("  %d drops with Deep coverage and a fully-contained window, %s to %s"
        % (len(man), man.t.iloc[0], man.t.iloc[-1]))
    drops_ns = man.t.astype("int64").to_numpy()

    files = sorted(DEEP.glob("SAFOD-Deep-*.h5"))
    def ftime(p):
        m = re.search(r"_(\d{4}-\d{2}-\d{2}T\d{6})Z\.h5$", p.name)
        return pd.to_datetime(m.group(1), format="%Y-%m-%dT%H%M%S",
                              utc=True) if m else None
    ftimes = [(p, ftime(p)) for p in files]
    ftimes = [(p, t) for p, t in ftimes if t is not None]
    say("  %d Deep files, %s to %s" % (len(ftimes), ftimes[0][1], ftimes[-1][1]))

    with h5py.File(ftimes[len(ftimes) // 2][0], "r") as h:
        g = h["Acquisition/Raw[0]"]
        fs_raw = float(g.attrs.get("OutputDataRate", 1000.0))
        n_loci = int(g["RawData"].shape[1])
        dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 2.0419))
        gl = h["Acquisition"].attrs.get("GaugeLength", np.nan)
    say("  %d channels x %.4f m = %.2f km | %.0f Hz | gauge length %s m"
        % (n_loci, dx, n_loci * dx / 1e3, fs_raw, gl))
    say("")

    def read_around(path, offset_s):
        """Read [offset - PRE - PAD, offset + POST + PAD] as (channels, time)."""
        with h5py.File(path, "r") as h:
            g = h["Acquisition/Raw[0]"]
            fs = float(g.attrs.get("OutputDataRate", fs_raw))
            n_tot = g["RawData"].shape[0]
            lo = int(round((offset_s - PRE_S - PAD_S) * fs))
            hi = int(round((offset_s + POST_S + PAD_S) * fs))
            if lo < 0 or hi > n_tot:
                return None, fs
            return g["RawData"][lo:hi, :].astype(np.float32).T, fs

    def snr_envelope(arr, fs):
        """Envelope / per-channel pre-drop RMS, on the trimmed time axis."""
        x, fs2 = band(arr, fs)
        if x is None:
            return None, None, None
        npad = int(round(PAD_S * fs2))
        x = x[:, npad:x.shape[1] - npad] if npad > 0 else x
        t = -PRE_S + np.arange(x.shape[1]) / fs2
        env = np.abs(hilbert(x, axis=1))
        ref = (t >= NOISE_LO) & (t <= NOISE_HI)
        rms = np.sqrt(np.mean(x[:, ref] ** 2, axis=1))
        rms = np.where(rms > 0, rms, np.nan)
        return env / rms[:, None], t, (x, fs2)

    # ---------------- pass over drop windows and sham windows ----------------
    cand = []
    for p, t in ftimes:
        mid = t + pd.Timedelta(seconds=30)
        gap = float(np.min(np.abs(drops_ns - mid.value)) / 1e9) if drops_ns.size else 1e9
        if gap >= GUARD_S:
            cand.append(p)
    say("  %d of %d files are drop-free within %.0f s of their midpoint"
        % (len(cand), len(ftimes), GUARD_S))

    picks = man.iloc[np.linspace(0, len(man) - 1, min(N_DROP, len(man))).astype(int)]
    sham_files = [cand[i] for i in
                  np.linspace(0, len(cand) - 1, min(N_SHAM, len(cand))).astype(int)]

    acc = {}
    for arm in ("drop", "sham"):
        acc[arm] = dict(n=0, env=None, tvec=None,
                        fk={a: None for a in APERTURES_M},
                        pos={a: None for a in APERTURES_M})

    def accumulate(arm, path, offset_s):
        arr, fs = read_around(path, offset_s)
        if arr is None:
            return "window does not fit"
        s, t, xf = snr_envelope(arr, fs)
        if s is None:
            return "decimation factor not integral"
        a = acc[arm]
        a["env"] = s if a["env"] is None else a["env"] + s
        a["tvec"] = t
        # Part B on the same window: the f-k window starts FK_LEAD_S before t=0
        x, fs2 = xf
        i0 = int(round((PRE_S - FK_LEAD_S) * fs2))
        seg = x[:, i0:i0 + int(FK_WIN_S * fs2)]
        for ap in APERTURES_M:
            out = signed_velocity_spectrum(seg, fs2, dx, ap, N_POS)
            if out is None:
                continue
            pos, Pp, Pm = out
            a["pos"][ap] = pos
            if a["fk"][ap] is None:
                a["fk"][ap] = [Pp, Pm]
            else:
                a["fk"][ap][0] += Pp
                a["fk"][ap][1] += Pm
        a["n"] += 1
        return None

    say("")
    say("=== reading %d drop windows ===" % len(picks))
    for _, row in picks.iterrows():
        path = DEEP / str(row.deep_file)
        if not path.is_file():
            say("  missing %s" % path.name); continue
        why = accumulate("drop", path, float(row.deep_offset_s))
        if why:
            say("  %s skipped: %s" % (row.t, why))
    say("  stacked %d drop windows" % acc["drop"]["n"])

    say("")
    say("=== reading %d sham windows (no drop within %.0f s) ==="
        % (len(sham_files), GUARD_S))
    for p in sham_files:
        why = accumulate("sham", p, 30.0)
        if why:
            say("  %s skipped: %s" % (p.name, why))
    say("  stacked %d sham windows" % acc["sham"]["n"])

    if not acc["drop"]["n"] or not acc["sham"]["n"]:
        raise SystemExit("one arm is empty; nothing to compare")

    for arm in ("drop", "sham"):
        acc[arm]["env"] = acc[arm]["env"] / acc[arm]["n"]
        for ap in APERTURES_M:
            if acc[arm]["fk"][ap] is not None:
                acc[arm]["fk"][ap] = [q / acc[arm]["n"] for q in acc[arm]["fk"][ap]]

    # ---------------- Part A reading ----------------
    t = acc["drop"]["tvec"]
    post = t > 0.05
    say("")
    say("=== PART A: is the drop visible? ===")
    say("  Stacked SNR envelope, peak over t > 0.05 s, per fibre position.")
    peak = {}
    for arm in ("drop", "sham"):
        peak[arm] = np.nanmax(acc[arm]["env"][:, post], axis=1)
    fibre = np.arange(peak["drop"].size) * dx
    ratio = peak["drop"] / np.where(peak["sham"] > 0, peak["sham"], np.nan)
    # Depth registration, available since 2026-08-20. Channels in the surface
    # lead-in are NOT in the hole, so a "response" there is not a formation
    # response and must be reported separately rather than mixed in.
    g = GEO.load()
    nch = peak["drop"].size
    chan = np.arange(nch)
    j = np.clip(np.searchsorted(g["channel"], chan), 0, g["channel"].size - 1)
    tvd_m = g["tvd_m"][j]
    in_hole = g["in_hole"][j]
    outbound = g["outbound"][j]
    say("  geometry: %d of %d channels in the hole; surface lead-in %d channels"
        % (int(in_hole.sum()), nch, int((~in_hole).sum())))
    say("  arm    median peak SNR   90th pct   max      argmax fibre position")
    for arm in ("drop", "sham"):
        pk = peak[arm]
        k = int(np.nanargmax(pk))
        say("  %-6s %14.3f %10.3f %8.3f %14.0f m  [%s]"
            % (arm, np.nanmedian(pk), np.nanpercentile(pk, 90), np.nanmax(pk),
               fibre[k], GEO.describe(int(chan[k]))))
    say("")
    say("  drop / sham peak-SNR ratio: median %.3f, max %.3f at fibre %.0f m"
        % (np.nanmedian(ratio), np.nanmax(ratio), fibre[int(np.nanargmax(ratio))]))
    thr = 1.5
    hot = np.flatnonzero(ratio > thr)
    if hot.size:
        say("  %d of %d channels have drop/sham > %.1f, spanning fibre %.0f-%.0f m"
            % (hot.size, ratio.size, thr, fibre[hot[0]], fibre[hot[-1]]))
        # contiguous runs, so we can say whether the response LOCALISES
        runs, s0 = [], hot[0]
        for a_, b_ in zip(hot[:-1], hot[1:]):
            if b_ - a_ > 5:
                runs.append((s0, a_)); s0 = b_
        runs.append((s0, hot[-1]))
        runs = [r for r in runs if (r[1] - r[0]) * dx >= 20.0]
        say("  contiguous responding sections (>= 20 m), with depth:")
        for a_, b_ in sorted(runs, key=lambda r: -(r[1] - r[0]))[:8]:
            frac_air = float((~in_hole[a_:b_ + 1]).mean())
            if frac_air > 0.99:
                where = "SURFACE LEAD-IN -- not in the hole"
            elif frac_air > 0.01:
                where = ("STRADDLES lead-in and hole (%.0f%% lead-in), "
                         "TVD up to %.0f m" % (100 * frac_air,
                                               np.nanmax(tvd_m[a_:b_ + 1])))
            else:
                limb = "outbound" if outbound[a_:b_ + 1].mean() > 0.5 else "return limb"
                where = "TVD %.0f-%.0f m, %s" % (np.nanmin(tvd_m[a_:b_ + 1]),
                                                 np.nanmax(tvd_m[a_:b_ + 1]), limb)
            say("    fibre %6.0f - %6.0f m (%4.0f m, %3d ch), max ratio %5.2f  | %s"
                % (fibre[a_], fibre[b_], fibre[b_] - fibre[a_], b_ - a_ + 1,
                   np.nanmax(ratio[a_:b_ + 1]), where))
        say("")
        say("  in-hole channels only (surface lead-in excluded):")
        rh = np.where(in_hole, ratio, np.nan)
        nh = int(np.nansum(rh > thr))
        say("    %d of %d in-hole channels have drop/sham > %.1f; median ratio %.3f"
            % (nh, int(in_hole.sum()), thr, np.nanmedian(rh)))
        if nh:
            kk = int(np.nanargmax(rh))
            say("    strongest in-hole: fibre %.0f m, TVD %.0f m, ratio %.2f"
                % (fibre[kk], tvd_m[kk], rh[kk]))
    else:
        say("  NO channel reaches drop/sham > %.1f. The weight drops are NOT" % thr)
        say("  detectable on this fibre in 5-20 Hz strain rate at this stack depth.")

    # time of the peak, which is what carries the moveout
    tp = t[np.nanargmax(np.where(post, acc["drop"]["env"], -np.inf), axis=1)]
    if hot.size:
        say("")
        say("  arrival time of the stacked peak, responding channels only:")
        say("    median %.3f s, 10-90 pct %.3f-%.3f s"
            % (np.nanmedian(tp[hot]), np.nanpercentile(tp[hot], 10),
               np.nanpercentile(tp[hot], 90)))

    # ---------------- Part B reading ----------------
    say("")
    say("=== PART B: where is the energy in signed apparent velocity? ===")
    vmid = np.sqrt(V_EDGES[:-1] * V_EDGES[1:])
    fanbin = (vmid >= 2500.0) & (vmid <= 4000.0)
    for ap in APERTURES_M:
        nch = int(round(ap / dx))
        say("")
        say("  aperture %.0f m (%d ch): dk = %.5f /m; at 10 Hz the 2500-4000 m/s"
            % (ap, nch, 1.0 / (nch * dx)))
        say("  fan spans k = %.5f-%.5f /m, i.e. %.1f wavenumber bins wide"
            % (10.0 / 4000, 10.0 / 2500,
               (10.0 / 2500 - 10.0 / 4000) * nch * dx))
        for arm in ("drop", "sham"):
            if acc[arm]["fk"][ap] is None:
                continue
            Pp, Pm = acc[arm]["fk"][ap]
            tot = Pp + Pm
            share = 100.0 * tot[:, fanbin].sum() / max(tot.sum(), 1e-300)
            # velocity of the band-power maximum, averaged over positions
            prof = tot.sum(axis=0)
            say("    %-5s fan share of 5-20 Hz power %5.2f %% | band-power peak at"
                " %.0f m/s" % (arm, share, vmid[int(np.argmax(prof))]))
            # signed asymmetry inside the frozen fan, per aperture position
            a_ = Pp[:, fanbin].sum(axis=1); b_ = Pm[:, fanbin].sum(axis=1)
            asym = np.abs(a_ - b_) / np.maximum(a_ + b_, 1e-300)
            say("          |A| in the frozen fan, per position: median %.3f, max"
                " %.3f at fibre %.0f m"
                % (np.median(asym), asym.max(),
                   acc[arm]["pos"][ap][int(np.argmax(asym))]))

    # ---------------- figure ----------------
    fig = plt.figure(figsize=(16.0, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    ext = [t[0], t[-1], fibre[-1] / 1e3, 0.0]
    for i, arm in enumerate(("drop", "sham")):
        axi = fig.add_subplot(gs[i, 0])
        im = axi.imshow(acc[arm]["env"], aspect="auto", extent=ext,
                        vmin=0.6, vmax=2.2, cmap="magma", interpolation="nearest")
        axi.axvline(0, color="c", lw=0.9, ls="--")
        axi.set(xlabel="time from drop (s)", ylabel="fibre position (km)",
                title="%s: stacked SNR envelope (n=%d)" % (arm, acc[arm]["n"]))
        fig.colorbar(im, ax=axi, label="envelope / pre-drop RMS")

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(peak["drop"], fibre / 1e3, lw=0.7, color="#D55E00", label="drop")
    ax2.plot(peak["sham"], fibre / 1e3, lw=0.7, color="#0072B2", label="sham")
    ax2.invert_yaxis()
    ax2.set(xlabel="peak stacked SNR (t > 0)", ylabel="fibre position (km)",
            title="Peak per channel")
    ax2.legend(fontsize=8); ax2.grid(alpha=.3)

    ax3 = fig.add_subplot(gs[1, 1])
    ax3.plot(ratio, fibre / 1e3, lw=0.7, color="k")
    ax3.axvline(1.0, color="0.5", lw=1.0)
    ax3.axvline(thr, color="r", lw=1.0, ls="--", label="ratio %.1f" % thr)
    ax3.invert_yaxis()
    ax3.set(xlabel="drop / sham peak SNR", ylabel="fibre position (km)",
            title="Where the fibre responds to the drops")
    ax3.legend(fontsize=8); ax3.grid(alpha=.3)

    for j, ap in enumerate(APERTURES_M):
        axv = fig.add_subplot(gs[j, 2])
        for arm, col in (("drop", "#D55E00"), ("sham", "#0072B2")):
            if acc[arm]["fk"][ap] is None:
                continue
            Pp, Pm = acc[arm]["fk"][ap]
            prof = (Pp + Pm).sum(axis=0)
            axv.loglog(vmid, prof / max(prof.max(), 1e-300), color=col, label=arm)
        axv.axvspan(2500, 4000, color="0.85", zorder=0, label="frozen fan")
        axv.set(xlabel="apparent velocity along fibre (m/s)",
                ylabel="normalised band power",
                title="%.0f m aperture: 5-20 Hz power vs |v|" % ap)
        axv.legend(fontsize=8); axv.grid(alpha=.3, which="both")
    fig.savefig(str(STEM) + ".png", dpi=170)

    out = dict(fibre_m=fibre, t=t, dx=dx, n_drop=acc["drop"]["n"],
               n_sham=acc["sham"]["n"], v_mid=vmid,
               peak_drop=peak["drop"], peak_sham=peak["sham"], ratio=ratio,
               env_drop=acc["drop"]["env"].astype(np.float32),
               env_sham=acc["sham"]["env"].astype(np.float32))
    for arm in ("drop", "sham"):
        for ap in APERTURES_M:
            if acc[arm]["fk"][ap] is None:
                continue
            out["%s_%d_Pp" % (arm, int(ap))] = acc[arm]["fk"][ap][0]
            out["%s_%d_Pm" % (arm, int(ap))] = acc[arm]["fk"][ap][1]
            out["%s_%d_pos" % (arm, int(ap))] = acc[arm]["pos"][ap]
    np.savez_compressed(str(STEM) + ".npz", **out)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
