"""
N1 -- the measurement only the fiber can make: repeaters or neighbours, geometrically.

THE IDEA. A repeater pair from the SAME patch must arrive at the same ray parameter.
Two patches one source-radius apart arrive at slightly different angles, and on a
900-channel vertical array that appears as a differential moveout. For each pair,
per channel, measure the delay dt(z) between the two events and fit

    dt(z) = a + b * z

The intercept a absorbs the origin-time difference (a nuisance, 0.1-0.5 s of catalog
error). THE SLOPE b IS THE DIFFERENTIAL RAY PARAMETER, which maps to the relative
source offset perpendicular to the ray.

WHY NO SEISMOMETER NETWORK CAN DO THIS. A single 3-component sensor gets incidence
angle from particle motion to a few degrees; this needs 0.3 degrees. A network
constrains relative location through differential travel times, but that is
velocity-model-limited -- it is exactly the dedicated relative relocation Nadeau
1995 had to perform to turn 200 m of routine scatter into 10-20 m. The array
measures the angle LOCALLY, along a line, with 800 samples of it. And for a pair the
measurement is model-free in the differential sense: both events traverse identical
structure, so velocity errors cancel and only the source offset survives.

IT ANSWERS THE QUESTION ELLSWORTH ASKED: "close in magnitude and location, which
will then need to be verified as either repeaters or neighbors."

THE ONE REAL CONFOUND. A velocity change between the two occurrences also produces a
dt gradient linear in depth. Two things separate it, both reported:
  (a) magnitude -- the independent HRSN measurement bounds |dv/v| < 0.1% at
      251-284 m, which over 864 m at 3000 m/s is 0.29 ms against 1.3 ms for a
      one-radius source offset, a factor of 4 below signal;
  (b) common mode -- a medium change is shared by all pairs spanning the same
      interval, a source offset is pair-specific. Regress slope on elapsed time.

THE CONTROL THAT DECIDES EVERYTHING (N1.3). The clusters split by prior expectation:
long-interval candidates (440 d, 272 d) should give slope ~ 0 if one patch;
same-day bursts are almost certainly different patches and should give resolvably
non-zero slope. If those groups do not separate, the method is wrong and nothing
downstream should be written.

Frozen before evaluation: channels 23-896 (G0), band, P window.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dvv_core import sub_sample_delay, bulk_align                  # noqa: E402
from moveout_test import slant_scan, PRE_S, DX                     # noqa: E402

CH_LO, CH_HI = 23, 896        # G0: wellhead at 23, coupling flat to 896
WELLHEAD_CH = 23
BAND = (5.0, 20.0)
P_WIN = (-0.20, 1.00)         # about the picked P, per channel, follows the moveout
MIN_COH = 0.5                 # per-channel phase-stability floor
MIN_CH = 100                  # channels needed for a slope
CLIP_MS = 20.0                # coarse guard before robust rejection
NSIG = 3.5                    # robust residual rejection, iterative
MU_PA, DSIGMA_PA = 30.0e9, 3.0e6


def prep_abs(tag):
    """Load an event over the G0 channel range, keeping ABSOLUTE channel indices.

    moveout_test.prep cannot be reused directly: it is hardwired to channels
    100-800 (pre-G0) and returns compressed arrays after dropping bad channels, so
    two events of a pair come back with different lengths. Indexing those by
    position would compare different depths and the error would be invisible.
    Here the surviving absolute channel numbers are returned so the pair can be
    intersected on channel identity.
    """
    f = os.path.join(HERE, 'cache_all', f'{tag}.npz')
    if not os.path.exists(f):
        return None
    d = np.load(f)
    X, fs = d['X'], float(d['fs'])
    if X.shape[0] < CH_HI:
        return None
    sos = butter(4, list(BAND), btype='band', fs=fs, output='sos')
    A = sosfiltfilt(sos, X[CH_LO:CH_HI].astype(np.float64), axis=-1)
    A -= A.mean(axis=1, keepdims=True)
    rms = np.sqrt((A ** 2).mean(axis=1))
    good = np.isfinite(rms) & (rms > 0)
    if good.sum() > 20:
        med = np.median(rms[good])
        good &= (rms > 0.05 * med) & (rms < 20.0 * med)
    ch = np.arange(CH_LO, CH_HI)[good]
    return A[good] / rms[good, None], ch, fs


def source_radius_m(mag):
    m0 = 10.0 ** (1.5 * float(mag) + 9.1)
    return (7.0 * m0 / (16.0 * DSIGMA_PA)) ** (1.0 / 3.0)


def robust_line(z, dt, nsig=NSIG, iters=4):
    """Least-squares line with iterative rejection on the residual about the fit.

    Rejection must be on the RESIDUAL, not on |dt|: the delays legitimately trend
    with depth (that trend is the signal), so a fixed |dt| cut either keeps cycle
    skips or discards real signal at the ends of the interval. g5_shallow_dvv made
    exactly that mistake and produced -21% velocity changes that were skips.
    """
    m = np.isfinite(z) & np.isfinite(dt)
    z, dt = z[m], dt[m]
    if z.size < MIN_CH:
        return np.nan, np.nan, np.nan, 0, np.nan
    keep = np.ones(z.size, bool)
    for _ in range(iters):
        if keep.sum() < MIN_CH:
            break
        b, a = np.polyfit(z[keep], dt[keep], 1)
        r = dt - (a + b * z)
        s = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        if not np.isfinite(s) or s <= 0:
            break
        new = np.abs(r) < nsig * s
        if new.sum() == keep.sum():
            keep = new
            break
        keep = new
    if keep.sum() < MIN_CH:
        return np.nan, np.nan, np.nan, int(keep.sum()), np.nan
    b, a = np.polyfit(z[keep], dt[keep], 1)
    r = dt[keep] - (a + b * z[keep])
    sig = float(np.std(r))
    # slope error for a least-squares fit over an aperture
    zs = z[keep]
    se = sig / (np.sqrt(keep.sum()) * np.std(zs)) if np.std(zs) > 0 else np.nan
    return float(b), float(a), se, int(keep.sum()), sig


def pair_delays(tag_a, tag_b):
    """Per-channel delay of B relative to A on the direct-P window.

    The window follows the moveout: for channel at fibre distance z the P arrival
    sits at tP + p*z, using the slowness fitted to event A. Both events are windowed
    on A's moveout so the comparison is like-for-like; a difference in their own
    moveouts is precisely the signal being measured.
    """
    ga, gb = prep_abs(tag_a), prep_abs(tag_b)
    if ga is None or gb is None:
        return None
    Aa, cha, fs = ga
    Ab, chb, _ = gb
    # intersect on ABSOLUTE channel index so both events give the same depths
    common = np.intersect1d(cha, chb)
    if common.size < MIN_CH:
        return None
    Aa = Aa[np.searchsorted(cha, common)]
    Ab = Ab[np.searchsorted(chb, common)]
    # slant_scan wants fibre distance from the first retained channel
    zfib = (common - common[0]) * DX
    p, tp, sem, _, _ = slant_scan(Aa, zfib, fs)
    if not np.isfinite(sem) or sem <= 0:
        return None
    # BULK-ALIGN THE PAIR FIRST. Catalog origin times carry 0.1-0.5 s of error.
    # Placing the P window from event A's pick and applying it unshifted to event B
    # puts B's arrival 100-500 ms outside the window for a months-apart pair; every
    # per-channel delay then exceeds the 20 ms clip and the pair returns ZERO
    # channels. That is exactly what happened on the first run: 11 of 13 pairs came
    # back empty and only same-day pairs -- whose origin times happen to be
    # consistent -- survived. The bulk shift is bookkeeping; the signal is the
    # RESIDUAL variation of delay with depth after it is removed.
    sa = Aa.mean(axis=0)
    sb = Ab.mean(axis=0)
    _, lag_s = bulk_align(sa, sb, fs, max_lag_s=2.0)
    kshift = int(round(lag_s * fs))

    n = common.size
    dt = np.full(n, np.nan)
    coh = np.zeros(n)
    for k in range(n):
        t_arr = tp + p * zfib[k]
        i0 = int((t_arr + P_WIN[0]) * fs)
        i1 = int((t_arr + P_WIN[1]) * fs)
        j0, j1 = i0 + kshift, i1 + kshift          # B's window, shifted
        if i0 < 0 or i1 > Aa.shape[1] or j0 < 0 or j1 > Ab.shape[1]:
            continue
        d, c = sub_sample_delay(Aa[k, i0:i1], Ab[k, j0:j1], fs, BAND,
                                max_lag_s=0.25)
        if np.isfinite(d) and abs(d) < CLIP_MS / 1000.0:
            dt[k], coh[k] = d, c
    ok = coh > MIN_COH
    dt[~ok] = np.nan
    depth = (common - WELLHEAD_CH) * DX          # G0 registration
    return depth, dt, p, tp, sem, int(np.isfinite(dt).sum())


def main():
    ev = pd.read_csv(os.path.join(HERE, 'phaseA_events.csv'))
    ev['time'] = pd.to_datetime(ev.time, utc=True, format='mixed')
    ev = ev[ev.cov_full].reset_index(drop=True)
    B = pd.read_csv(os.path.join(HERE, 'beta_similarity.csv'))
    for c in ('t_i', 't_j'):
        B[c] = pd.to_datetime(B[c], utc=True, format='mixed')

    cand = B[B.beta >= 0.90].copy()
    print(f'{len(cand)} pairs at beta >= 0.90\n', flush=True)

    rows, curves = [], {}
    for _, r in cand.iterrows():
        i, j = int(r.i), int(r.j)
        got = pair_delays(ev.tag[i], ev.tag[j])
        if got is None:
            print(f'  {ev.tag[i]}/{ev.tag[j]}: not cached'); continue
        depth, dt, p, tp, sem, nch = got
        b, a, se, nfit, sig = robust_line(depth, dt)
        hyp = np.sqrt((r.sep_km * 0 + 1) * 0 +
                      (ev.depth[i]) ** 2 + (ev.km_safod[i]) ** 2)
        rad = source_radius_m(0.5 * (ev.mag[i] + ev.mag[j]))
        # slope (s/m) -> source offset perpendicular to the ray
        # dp = cos(i)/V * dtheta, dtheta = offset / hypocentral distance
        V = 3000.0
        offset_m = (abs(b) * V * hyp * 1000.0 / max(np.cos(np.radians(10.0)), 1e-9)
                    if np.isfinite(b) else np.nan)
        kind = ('burst' if r.days < 5 else
                'long' if r.days > 200 else 'mid')
        rows.append(dict(i=i, j=j, t_i=r.t_i, t_j=r.t_j, days=r.days,
                         beta=r.beta, dmag=r.dmag, kind=kind,
                         nch=nch, nfit=nfit, sigma_ms=1000 * sig if sig == sig else np.nan,
                         slope_s_per_m=b, slope_err=se,
                         slope_ms_across=1000 * b * 864.0 if b == b else np.nan,
                         slope_err_ms=1000 * se * 864.0 if se == se else np.nan,
                         nsigma=abs(b) / se if (se and se > 0) else np.nan,
                         hyp_km=hyp, radius_m=rad, offset_m=offset_m,
                         semblance=sem))
        curves[(i, j)] = (depth, dt)
        print(f'  {r.t_i:%Y-%m-%d}/{r.t_j:%Y-%m-%d} {kind:5s} '
              f'{nfit:4d} ch  sigma {1000*sig:5.2f} ms  '
              f'slope {1000*b*864:+7.3f} +/- {1000*se*864:5.3f} ms  '
              f'{abs(b)/se if se and se>0 else np.nan:5.1f}s', flush=True)

    R = pd.DataFrame(rows)
    if R.empty:
        print('nothing measured'); return
    R.to_csv(os.path.join(HERE, 'ray_parameter_test.csv'), index=False)

    print(f'\nACHIEVED PRECISION')
    print(f'  per-channel timing sigma : {R.sigma_ms.median():.2f} ms '
          f'(plan assumed 2.0 from g5_shallow_dvv)')
    print(f'  slope error across 864 m : {R.slope_err_ms.median():.3f} ms')

    print('\nN1.3  THE CONTROL -- do bursts separate from long-interval pairs?')
    for k in ('long', 'mid', 'burst'):
        g = R[R.kind == k].dropna(subset=['slope_ms_across'])
        if len(g):
            print(f'  {k:5s} (n={len(g)}): |slope| median '
                  f'{g.slope_ms_across.abs().median():6.3f} ms   '
                  f'median significance {g.nsigma.median():4.1f} sigma')
    lo = R[(R.kind == 'long')].dropna(subset=['slope_ms_across'])
    bu = R[(R.kind == 'burst')].dropna(subset=['slope_ms_across'])
    print('\n  VERDICT')
    if len(lo) and len(bu):
        if bu.slope_ms_across.abs().median() > 2 * lo.slope_ms_across.abs().median():
            print('  -> bursts show larger differential moveout than long-interval')
            print('     pairs, as different patches should. The discriminant works;')
            print('     proceed to N2.')
        else:
            print('  -> bursts and long-interval pairs are NOT separated. The')
            print('     discriminant does not work as posed. Report the achieved')
            print('     slope precision as an upper bound on resolvable source')
            print('     separation and stop -- do not reclassify anything on it.')
    else:
        print(f'  -> insufficient sample (long n={len(lo)}, burst n={len(bu)});')
        print('     cannot run the control, so nothing downstream is justified.')

    print('\nN1.2  CONFOUND -- is the slope common-mode (medium) or pair-specific?')
    g = R.dropna(subset=['slope_s_per_m', 'days'])
    if len(g) > 3:
        rr = np.corrcoef(g.days, g.slope_s_per_m)[0, 1]
        print(f'  r(slope, elapsed days) = {rr:+.2f}')
        print('  a shared velocity change would trend with elapsed time; a source')
        print('  offset would not. |dv/v| < 0.1% (HRSN) is 0.29 ms across 864 m,')
        print('  versus 1.3 ms for a one-source-radius offset.')

    # ---- figure: the poster panel -------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for a_, k, col in ((ax[0], 'long', 'C0'), (ax[1], 'burst', 'C3')):
        sub = R[R.kind == k]
        for _, r in sub.iterrows():
            c = curves.get((int(r.i), int(r.j)))
            if c is None:
                continue
            d, dt = c
            m = np.isfinite(dt)
            a_.plot(1000 * dt[m], d[m], '.', ms=2, alpha=.35, color=col)
            if np.isfinite(r.slope_s_per_m):
                zz = np.array([np.nanmin(d[m]), np.nanmax(d[m])])
                a_.plot(1000 * (r.slope_s_per_m * zz +
                        (r.slope_ms_across / 864.0 * 0)), zz, '-', lw=1.4, color='k')
        a_.invert_yaxis()
        a_.set(xlabel='delay B-A (ms)', ylabel='depth (m)', xlim=(-8, 8),
               title=f'{"A  long-interval pairs" if k=="long" else "B  same-day bursts"}\n'
                     f'(one patch expected)' if k == 'long' else
                     'B  same-day bursts\n(different patches expected)')
        a_.grid(alpha=.3)
    a_ = ax[2]
    for k, col in (('long', 'C0'), ('mid', '0.5'), ('burst', 'C3')):
        g = R[R.kind == k].dropna(subset=['slope_ms_across'])
        a_.errorbar(g.days, g.slope_ms_across.abs(), yerr=g.slope_err_ms,
                    fmt='o', ms=6, capsize=3, color=col, label=k)
    a_.set(xscale='symlog', xlabel='days between events',
           ylabel='|differential moveout| across 864 m (ms)',
           title='C  The discriminant')
    a_.legend(fontsize=8); a_.grid(alpha=.3)
    fig.suptitle('Differential ray parameter: repeaters or neighbours, '
                 'measured geometrically by the fibre', fontsize=12)
    fig.tight_layout()
    p_ = os.path.join(HERE, 'ray_parameter_test.png')
    fig.savefig(p_, dpi=140)
    print(f'\nwrote ray_parameter_test.csv and {os.path.basename(p_)}')


if __name__ == '__main__':
    main()
