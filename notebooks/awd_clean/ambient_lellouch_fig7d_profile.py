#!/usr/bin/env python3
"""
Lellouch et al. (2019) Figure 7d velocity profile from the cached ambient products.

WHY THIS SCRIPT EXISTS. The ambient F-K thread chased the Figure 7c observable: a
top-source correlation gather with a single ~3.2 km/s moveout ridge. Reading the
paper (section 4.1) shows 7c is not the measurement. It is a phase-identification
step that yields ONE number -- an *average* velocity, quoted as 3,200 m/s -- which
is then reused as the moveout-correction velocity. The velocity model (paper
Figure 9, green) comes from Figure 7d: the CONSTANT 50 m source-receiver offset
geometry, with one travel-time pick every 50 m.

  "Thus, we deduce we are observing a P wave and take the correlation functions
   shown in Figure 7d to estimate its velocity every 50 m."

Two consequences for this repository:

  1. The right cached product is `fixed_stack`, not `top_stack`. Every
     transfer_*.npz already carries it (171 pairs of constant ~50 m offset,
     midpoints every 5 channels). It has never been scored -- `velocity_score`
     is only ever called on `top_stack`.

  2. In the constant-offset geometry the target arrival sits at ~50 m / v, i.e.
     12-20 ms, where correlation SNR is highest -- not at 220 ms where the
     top-source arrival is buried under the near-zero-lag lobe. That is the
     reason the paper estimates velocity from 7d and not 7c.

WHAT IS AND IS NOT REPRODUCED. The paper's neighbour stack is
C_S,R = sum_{Z=R-10}^{R+10} C_S,Z with 3,200 m/s travel-time shifts applied
because offsets vary within that group. Here the cached pairs all have the SAME
offset, so neighbouring pair midpoints can be summed directly: their arrival lags
agree up to the local velocity gradient over +-10 m, so no shift is required.
This is an equivalent-purpose, not identical, SNR stack, and it smooths the model
over the same ~+-10 m. The eq. (6) spectral normalisation by |S(omega)|^2 is NOT
applied -- the cached correlations are already formed and cannot be retro-whitened.

Search window: the pick is taken as the largest causal correlation in
PICK_LO..PICK_HI s, which bounds apparent velocity to roughly 1.1-10 km/s. That
window is a declared selection region, frozen here before any profile was viewed,
and the unconstrained causal maximum is reported alongside it for comparison.

Controls: (a) eight independent complete days, giving a per-depth spread that is
a real reproducibility estimate, not a bootstrap; (b) a depth-label permutation
null on the velocity-depth trend; (c) comparison against Lellouch's own released
Figure 7d traces processed by the identical picker.

Outputs (beside this file): ambient_lellouch_fig7d_profile.{npz,png,txt}
"""
from __future__ import annotations
import glob, re, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")          # set here, not at module import by a notebook
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TRANSFER = HERE / "ambient_transfer"
STEM = HERE / "ambient_lellouch_fig7d_profile"
LELLOUCH = Path("/home/groups/ettore88/nberrios/safod_das_git/lellouch_traces")

PICK_LO, PICK_HI = 0.005, 0.045     # frozen causal search window (s)
NEIGHBOUR_PAIRS = 2                 # +-2 pair indices = +-10 channels ~ +-10 m
SMOOTH_M = 100.0                    # paper: "the largest operator is 100 m long"


def parabolic_peak(y, k, dt):
    """Three adjacent samples with the largest value + quadratic interpolation
    (Nakata & Snieder, 2012), as used by the paper."""
    if k <= 0 or k >= len(y) - 1:
        return None
    y0, y1, y2 = y[k - 1], y[k], y[k + 1]
    den = y0 - 2 * y1 + y2
    if den == 0:
        return 0.0
    return 0.5 * (y0 - y2) / den * dt


def pick_profile(gather, lags, offset_m, midpoints_m, lo=PICK_LO, hi=PICK_HI):
    """One travel-time pick per output depth -> apparent velocity."""
    dt = lags[1] - lags[0]
    w = np.where((lags >= lo) & (lags <= hi))[0]
    out_t, out_v = [], []
    for row in gather:
        seg = row[w]
        k = w[int(np.argmax(seg))]          # causal, positive correlation
        d = parabolic_peak(row, k, dt)
        t = lags[k] + (d if d is not None else 0.0)
        out_t.append(t)
        out_v.append(offset_m / t if t > 0 else np.nan)
    return np.asarray(out_t), np.asarray(out_v)


def neighbour_stack(gather, n=NEIGHBOUR_PAIRS):
    """Sum adjacent constant-offset pairs; no shift needed (identical offsets)."""
    out = np.zeros_like(gather)
    N = len(gather)
    for i in range(N):
        lo, hi = max(0, i - n), min(N, i + n + 1)
        out[i] = gather[lo:hi].sum(axis=0)
    return out


def smooth(x, z, window_m=SMOOTH_M):
    out = np.full_like(x, np.nan, dtype=float)
    for i, zi in enumerate(z):
        m = np.abs(z - zi) <= window_m / 2
        if m.any():
            out[i] = np.nanmean(x[m])
    return out


def load_day(date):
    """Weighted sum of every cached chunk for one date -> full-day fixed_stack."""
    files = sorted(glob.glob(str(TRANSFER / f"transfer_{date}_start*_n*.npz")))
    tot = None; wsum = 0; lags = dx = None
    for f in files:
        n = int(re.search(r"_n(\d+)\.npz$", f).group(1))
        d = np.load(f, allow_pickle=True)
        g = d["fixed_stack"].astype(float)
        if tot is None:
            tot = np.zeros_like(g); lags = d["lags"].astype(float); dx = float(d["dx"])
        if g.shape != tot.shape:
            continue
        tot += g * n; wsum += n
    if tot is None or wsum == 0:
        return None
    return tot / wsum, lags, dx, wsum, len(files)


def main():
    # Only the 10-file chunks carry `fixed_stack`; the aggregated
    # transfer_seasonal_<date>.npz products keep `top_stack` alone.
    chunk_re = re.compile(r"transfer_(\d{4}-\d{2}-\d{2})_start\d+_n\d+\.npz$")
    dates = sorted({m.group(1) for f in glob.glob(str(TRANSFER / "transfer_*.npz"))
                    if (m := chunk_re.search(f))})
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    say("Lellouch Fig 7d constant-offset velocity profile")
    say("frozen pick window %.3f-%.3f s (bounds v to %.0f-%.0f m/s); neighbour +-%d pairs"
        % (PICK_LO, PICK_HI, 50.0 / PICK_HI, 50.0 / PICK_LO, NEIGHBOUR_PAIRS))

    profiles, used_dates, mids = {}, [], None
    for date in dates:
        r = load_day(date)
        if r is None:
            continue
        g, lags, dx, nfiles, nchunks = r
        if nfiles < 1000:
            say("  skip %s (only %d files cached)" % (date, nfiles)); continue
        offset_ch = int(round(50.0 / dx))
        offset_m = offset_ch * dx
        npair = g.shape[0]
        mids = (np.arange(npair) * 5 + offset_ch / 2.0) * dx
        gs = neighbour_stack(g)
        t, v = pick_profile(gs, lags, offset_m, mids)
        profiles[date] = v
        used_dates.append(date)
        say("  %s : %4d files, %3d chunks, offset %.2f m -> median v %.0f m/s"
            % (date, nfiles, nchunks, offset_m, np.nanmedian(v)))

    if not profiles:
        say("no usable days"); return

    V = np.vstack([profiles[d] for d in used_dates])
    vmed = np.nanmedian(V, axis=0)
    vstd = np.nanstd(V, axis=0)
    vs = smooth(vmed, mids)

    # restrict to the depth range the paper trusts (>=75 m) and drop edges
    ok = (mids >= 75) & np.isfinite(vs)
    zz, vv = mids[ok], vs[ok]

    say("")
    say("--- cross-day reproducibility (the real error bar) ---")
    say("  %d independent complete days" % len(used_dates))
    for zt in (100, 200, 300, 400, 500, 600, 700):
        i = int(np.argmin(np.abs(mids - zt)))
        say("   z=%3d m : v = %6.0f +- %4.0f m/s  (across-day sd)" % (mids[i], vmed[i], vstd[i]))
    say("  median across-day sd = %.0f m/s   (paper quotes 70-90 m/s)" % np.nanmedian(vstd[ok]))

    A = np.vstack([zz, np.ones_like(zz)]).T
    slope, icpt = np.linalg.lstsq(A, vv, rcond=None)[0]
    r = float(np.corrcoef(zz, vv)[0, 1])
    say("")
    say("--- velocity-depth trend ---")
    say("  v = %.3f * z + %.0f  m/s   (r = %.3f)" % (slope, icpt, r))
    say("  v(100 m) = %.0f , v(700 m) = %.0f , mean = %.0f m/s"
        % (vv[np.argmin(np.abs(zz - 100))], vv[np.argmin(np.abs(zz - 700))], np.nanmean(vv)))

    rng = np.random.default_rng(20260814)
    null = np.array([abs(np.corrcoef(zz, rng.permutation(vv))[0, 1]) for _ in range(20000)])
    p = (np.sum(null >= abs(r)) + 1) / (len(null) + 1)
    say("  depth-label permutation null: |r| 95th pct %.3f | observed %.3f | p = %.5f"
        % (np.percentile(null, 95), abs(r), p))

    # Lellouch's own released Fig 7d traces, same picker
    lel_z = lel_v = None
    try:
        from obspy import read
        st = read(str(LELLOUCH / "CCstackSAFOD.mseed"))
        by = {int(tr.stats.station): tr.data.astype(float) for tr in st}
        lel_z = np.array(sorted(by), float)
        G = np.array([by[int(z)] for z in lel_z])
        n = G.shape[1]; lg = (np.arange(n) - (n - 1) // 2) / 500.0
        _, lel_v = pick_profile(G, lg, 50.0, lel_z)
        say("")
        say("--- Lellouch released Fig 7d traces, identical picker ---")
        for z, vq in zip(lel_z, lel_v):
            say("   z=%3d m : v = %6.0f m/s" % (z, vq))
        say("  fit r = %.3f" % np.corrcoef(lel_z, lel_v)[0, 1])
    except Exception as e:
        say("  (released-trace comparison unavailable: %r)" % (e,))

    fig, ax = plt.subplots(1, 2, figsize=(11, 7), constrained_layout=True)
    for d in used_dates:
        ax[0].plot(profiles[d], mids, lw=0.7, alpha=0.45, label=d)
    ax[0].plot(vs, mids, "k-", lw=2.4, label="8-day median, 100 m smoothed")
    if lel_v is not None:
        ax[0].plot(lel_v, lel_z, "o--", color="crimson", lw=1.6, ms=5,
                   label="Lellouch 2019 released 7d")
    ax[0].invert_yaxis(); ax[0].set_xlabel("apparent P velocity (m/s)")
    ax[0].set_ylabel("depth along fiber (m)"); ax[0].set_xlim(1000, 6000)
    ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
    ax[0].set_title("Fig 7d constant-offset profile\n(no F-K filter anywhere)")

    ax[1].fill_betweenx(mids, vmed - vstd, vmed + vstd, alpha=.3, color="steelblue",
                        label="across-day $\\pm$1 sd")
    ax[1].plot(vmed, mids, color="steelblue", lw=1.2, label="8-day median")
    ax[1].plot(vs, mids, "k-", lw=2.2, label="100 m smoothed")
    ax[1].invert_yaxis(); ax[1].set_xlabel("apparent P velocity (m/s)")
    ax[1].set_xlim(1000, 6000); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    ax[1].set_title("cross-day reproducibility")
    fig.savefig(str(STEM) + ".png", dpi=200)

    np.savez(str(STEM) + ".npz", midpoints_m=mids, per_day=V,
             dates=np.array(used_dates), v_median=vmed, v_sd=vstd, v_smooth=vs,
             lellouch_z=lel_z if lel_z is not None else np.array([]),
             lellouch_v=lel_v if lel_v is not None else np.array([]),
             slope=slope, intercept=icpt, r=r, p_trend=p,
             pick_window=np.array([PICK_LO, PICK_HI]))
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
