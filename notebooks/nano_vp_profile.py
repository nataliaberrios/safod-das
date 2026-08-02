"""
P-wave velocity profile of the San Andreas damage zone from AWD active-source
DAS on the SAFOD cemented fiber.

Motivation: Lellouch et al. 2019 (10.1029/2019JB017533) published a Vp profile
for this same cemented fiber from *passive* earthquake arrivals recorded in
2017. The June 2026 AWD survey gives an independent *active-source* profile on
the identical sensor. Comparing them is a direct, same-instrument velocity
comparison across nine years at SAFOD.

Method (standard zero-offset VSP, with the source ~10-20 m from the wellhead so
the geometry is effectively vertical below ~50 m):

  stack all bursts -> one high-SNR record section
    -> STA/LTA first-break pick per channel
    -> refine each pick by cross-correlation against its neighbour, which is
       what buys sub-sample precision on a coherent wavefront
    -> enforce monotonic traveltime (a wavefront cannot arrive earlier deeper)
    -> interval Vp from a sliding linear fit of t(z), i.e. dz/dt

Everything is reported against distance along fiber, NOT calibrated depth --
the absolute offset is still unknown (see nano_fiber_end_coherence.py). The
shape of Vp(z) is what carries the science and it is offset-invariant.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt, hilbert

FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
NPZ = os.path.join(FIG_DIR, 'epoch_stacks_paired.npz')

DX = 1.26606202
PRE_S = 0.5
BAND = (20.0, 50.0)         # where the cemented fiber's direct-P SNR peaks

# Trust region. Top: Lellouch notes < 75 m is alluvium and your own SNR curve
# excludes the top 130 m. Bottom: direct-P SNR crosses 3 near 530 m.
Z_MIN_M, Z_MAX_M = 130.0, 530.0

STA_S, LTA_S = 0.006, 0.060
FIT_HALFWIN_M = 40.0        # sliding window half-length for interval velocity

# At 1.27 m channel spacing and Vp of a few km/s the channel-to-channel moveout
# is ~0.5 ms -- below one sample at 1 kHz. Picking every channel independently
# is therefore hopeless; instead track the wavefront across a coarser step where
# the moveout is many samples, and cross-correlate between those.
STEP_CH = 16                # ~20 m per step -> 7-10 ms of moveout, well resolved
V_FAST, V_SLOW = 6000.0, 500.0      # physical bounds used to bracket the search


def bandpass(x, fs, lo, hi):
    sos = butter(4, [lo, hi], btype='band', fs=fs, output='sos')
    return sosfiltfilt(sos, x, axis=-1)


def sta_lta_pick(tr, fs, i_lo, i_hi, thresh=3.0):
    """First break = FIRST crossing of the STA/LTA threshold inside [i_lo, i_hi].

    Taking the global maximum instead lands on the strongest arrival, which for
    these records is a late coda 1-3 s after the drop, not the direct P.
    """
    env = np.abs(hilbert(tr)) ** 2
    nsta, nlta = int(STA_S * fs), int(LTA_S * fs)
    csum = np.cumsum(np.insert(env, 0, 0))
    sta = (csum[nsta:] - csum[:-nsta]) / nsta
    lta = (csum[nlta:] - csum[:-nlta]) / nlta
    n = min(sta.size, lta.size)
    r = sta[:n] / np.where(lta[:n] > 0, lta[:n], np.nan)
    lo, hi = max(nlta, i_lo), min(n, i_hi)
    if hi <= lo:
        return -1
    seg = r[lo:hi]
    hits = np.where(np.isfinite(seg) & (seg > thresh))[0]
    if hits.size:
        return lo + int(hits[0])
    return lo + int(np.nanargmax(seg)) if np.any(np.isfinite(seg)) else -1


def cc_refine(a, b, guess_a, guess_b, fs, half=0.030):
    """Sub-sample lag of b relative to a, windowed about their coarse picks."""
    h = int(half * fs)
    wa = a[max(0, guess_a - h):guess_a + h]
    wb = b[max(0, guess_b - h):guess_b + h]
    n = min(wa.size, wb.size)
    if n < 8:
        return np.nan
    wa, wb = wa[:n] - wa[:n].mean(), wb[:n] - wb[:n].mean()
    cc = np.correlate(wb, wa, mode='full')
    lags = np.arange(-n + 1, n)
    k = int(np.argmax(cc))
    if 0 < k < cc.size - 1:                    # parabolic interpolation
        y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]
        denom = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
    else:
        shift = 0.0
    peak = cc[k] / np.sqrt(np.sum(wa ** 2) * np.sum(wb ** 2))
    return (lags[k] + shift) / fs, peak


def main():
    d = np.load(NPZ)
    nano = d['nano_stacks']
    fs = float(d['fs'])
    n_common = d['n_common']
    good = n_common > 0
    print(f'stacking {good.sum()}/{len(n_common)} epochs, '
          f'{int(n_common[good].sum())} drops total')

    # weight each epoch by how many drops went into it
    w = n_common[good].astype(float)
    sec = np.tensordot(w, nano[good], axes=(0, 0)) / w.sum()
    sec = bandpass(sec, fs, *BAND)
    n_ch = sec.shape[0]
    dist = np.arange(n_ch) * DX

    i0 = int(PRE_S * fs)
    c_lo, c_hi = int(Z_MIN_M / DX), int(Z_MAX_M / DX)
    steps = np.arange(c_lo, c_hi + 1, STEP_CH)
    print(f'picking channels {c_lo}-{c_hi} ({Z_MIN_M:.0f}-{Z_MAX_M:.0f} m) '
          f'in {steps.size} steps of {STEP_CH} ch ({STEP_CH*DX:.1f} m)')

    # Anchor: independent first break on the shallowest, highest-SNR step,
    # bracketed by physically possible arrival times so it cannot run away.
    z0 = dist[steps[0]]
    a_lo = i0 + int(z0 / V_FAST * fs)
    a_hi = i0 + int(z0 / V_SLOW * fs)
    p0 = sta_lta_pick(sec[steps[0]], fs, a_lo, a_hi)
    if p0 < 0:
        raise RuntimeError('no anchor pick found; check band and trust region')
    t = {steps[0]: (p0 - i0) / fs}
    print(f'anchor at {z0:.1f} m: t = {t[steps[0]]*1e3:.1f} ms '
          f'(apparent V {z0/t[steps[0]]:.0f} m/s)')

    # Track downward by cross-correlation between successive steps. Each step is
    # ~20 m, so the true moveout is 3-40 ms depending on velocity -- resolvable.
    ccq = {}
    prev_guess = p0
    for a, b in zip(steps[:-1], steps[1:]):
        dz = dist[b] - dist[a]
        # bracket the step's own moveout, then centre the CC search on it
        lo_s, hi_s = dz / V_FAST, dz / V_SLOW
        guess_b = prev_guess + int(0.5 * (lo_s + hi_s) * fs)
        out = cc_refine(sec[a], sec[b], prev_guess, guess_b, fs,
                        half=max(0.020, 1.5 * hi_s))
        if not isinstance(out, tuple) or not np.isfinite(out[0]):
            break
        dlag, peak = out
        dlag += (guess_b - prev_guess) / fs      # CC measured about the guess
        if not (lo_s <= dlag <= hi_s):           # reject unphysical step moveout
            dlag = np.clip(dlag, lo_s, hi_s)
        ccq[b] = peak
        t[b] = t[a] + dlag
        prev_guess = prev_guess + int(round(dlag * fs))

    picked = sorted(t)
    zz = dist[picked]
    tt = np.array([t[c] for c in picked])
    ok = np.isfinite(tt) & (tt > 0)
    zz, tt = zz[ok], tt[ok]
    ccq_arr = np.array([ccq.get(c, np.nan) for c in picked])[ok]
    print(f'{tt.size} picks, t from {tt[0]*1e3:.1f} to {tt[-1]*1e3:.1f} ms, '
          f'median step CC {np.nanmedian(ccq_arr):.3f}')

    # interval velocity: slope of a local linear fit of z on t
    vp, vz = [], []
    for i, z in enumerate(zz):
        m = np.abs(zz - z) <= FIT_HALFWIN_M
        if m.sum() < 4:
            continue
        A = np.vstack([tt[m], np.ones(m.sum())]).T
        slope = np.linalg.lstsq(A, zz[m], rcond=None)[0][0]
        if 200 < slope < 8000:
            vp.append(slope)
            vz.append(z)
    vp, vz = np.array(vp), np.array(vz)
    if vp.size == 0:
        raise RuntimeError('no interval velocities in physical range -- picks '
                           'are not tracking a coherent wavefront')
    print(f'interval Vp: median {np.median(vp):.0f} m/s, '
          f'range {vp.min():.0f}-{vp.max():.0f} m/s')
    for zq in [150, 200, 250, 300, 350, 400, 450, 500]:
        j = np.argmin(np.abs(vz - zq))
        if abs(vz[j] - zq) < 15:
            print(f'   {zq:4d} m along fiber : Vp ~ {vp[j]:5.0f} m/s')

    # apparent velocity straight from the moveout, as a sanity check
    A = np.vstack([tt, np.ones(tt.size)]).T
    v_avg, t0 = np.linalg.lstsq(A, zz, rcond=None)[0]
    print(f'\nsingle-line fit: mean Vp {v_avg:.0f} m/s, '
          f'zero-time intercept {t0:.1f} m')
    print('  (intercept is a crude handle on the depth offset: it estimates '
          'where\n   t=0 projects to, which should be the wellhead)')

    np.savez(os.path.join(FIG_DIR, 'nano_vp_profile.npz'),
             dist_picked=zz, t_picked=tt, vp=vp, vp_depth=vz,
             v_avg=v_avg, t0_intercept=t0, band=np.array(BAND),
             z_min=Z_MIN_M, z_max=Z_MAX_M, dx=DX, cc_quality=ccq_arr,
             n_drops=int(n_common[good].sum()))

    fig, ax = plt.subplots(1, 3, figsize=(14, 6))
    ax[0].imshow(sec[c_lo:c_hi + 1].T, aspect='auto', cmap='gray',
                 extent=[zz[0], zz[-1], (sec.shape[1] - i0) / fs, -PRE_S],
                 vmin=-np.percentile(np.abs(sec), 99),
                 vmax=np.percentile(np.abs(sec), 99))
    ax[0].plot(zz, tt, 'r-', lw=1.2, label='picked direct-P')
    ax[0].set_ylim(0.5, -0.05)
    ax[0].set_xlabel('distance along fiber (m)')
    ax[0].set_ylabel('time after drop (s)')
    ax[0].set_title(f'A  Stacked AWD section, {BAND[0]:.0f}-{BAND[1]:.0f} Hz')
    ax[0].legend(loc='lower right', fontsize=8)

    ax[1].plot(tt * 1e3, zz, 'k-', lw=1.4)
    ax[1].invert_yaxis()
    ax[1].set_xlabel('traveltime (ms)')
    ax[1].set_ylabel('distance along fiber (m)')
    ax[1].set_title('B  Direct-P traveltime')
    ax[1].grid(alpha=0.3)

    ax[2].plot(vp, vz, 'C0-', lw=1.6)
    ax[2].invert_yaxis()
    ax[2].set_xlabel('interval $V_P$ (m/s)')
    ax[2].set_title(f'C  Interval velocity\n({FIT_HALFWIN_M:.0f} m half-window)')
    ax[2].grid(alpha=0.3)
    fig.suptitle('SAFOD cemented fiber: active-source $V_P$ from 2026 AWD survey '
                 '(uncalibrated depth axis)', fontsize=11)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, 'nano_vp_profile.png')
    fig.savefig(out, dpi=140)
    print('\nSaved', out)


if __name__ == '__main__':
    main()
