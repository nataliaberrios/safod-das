"""
D0: can this array measure a depth-localized dv/v at all? Gate before any science.

THE MEASUREMENT BEING GATED. For a repeater pair, the per-channel delay dt(z)
accumulates along the ray, so inside an interval of uniform fractional velocity
change dv/v the delay is linear in one-way travel time tau:

    dt(z) = dt_0 - (dv/v) * tau(z)        =>   dv/v = -d(dt)/d(tau)

The intercept dt_0 absorbs the origin-time difference, so this is immune to the
unknown source instant. tau(z) = p * z is taken from the event's OWN measured
apparent slowness (moveout_test.slant_scan), so no velocity model is assumed.

WHY THIS IS WORTH MEASURING. Every published dv/v at Parkfield -- including the
2002-2022 record in doi:10.1029/2023jb028084 -- is volume-averaged over an unknown
scattering volume, and depth is inferred only from frequency dependence. A vertical
aperture localises the change to a known interval. Takagi & Okada 2012
(doi:10.1029/2012gl051342) did this on KiK-net vertical arrays, but with one or two
sensors per borehole, so they had to compare five stations of differing borehole
depth to get any depth information at all. This fiber has 900 receiver depths in one
hole.

TWO THINGS THAT MUST NOT BE GOT WRONG, both learned the hard way here.

1. BOTH EVENTS GET THE SAME ALIGNMENT SHIFT. Aligning each event with its own
   fitted slowness would remove exactly the differential moveout being measured.
   Applying one common shift to both preserves dt(z) identically while still
   flattening the wavefront for any per-bin stacking. The mean of the two fitted
   slownesses is used.

2. THE BOOTSTRAP IS BLOCKED AT ONE GAUGE LENGTH. 16.335 m / 1.021 m = 16 channels
   share the same averaged strain, so resampling channels singly understates the
   error by ~4x. dvv_core.slope_dvv(..., block=16). ray_parameter_test.py made
   precisely this mistake and its error bars were 4x too small.

--------------------------------------------------------------------------------
PASS CRITERIA, FIXED BEFORE RUNNING.

  G1  sigma(dv/v) from the random-pair null must be <= 0.1 %.
      Requirement comes from tau ~ 864 m / 2500 m/s = 0.35 s: resolving 0.1 %
      needs 0.35 ms of timing, and the predicted budget is 0.23 ms (Cramer-Rao
      calibrated against the measured 5.31 ms at the old 100 Hz configuration).
      The margin is only ~1.5x, which is why this is measured rather than assumed.

  G2  SYNTHETIC RECOVERY UNBIASED. Inject a known uniform dv/v of 0.1, 0.3 and
      1.0 % into one event of a pair and recover it within error. This is the
      single most important test in the plan: it validates the estimator against
      ground truth rather than against my expectations. Injection reuses
      moveout_test.align with slowness -eps*p, which applies exactly the
      shift -eps*p*z that a uniform fractional velocity change produces.

  G3  ACAUSAL FLOOR. Time-reverse one event of each repeater pair. A reversed
      record cannot contain a real interval arrival, so any apparent dv/v is pure
      method output and must be consistent with the null.

FAIL IS AN ACCEPTABLE OUTCOME. If sigma exceeds 0.1 % the deliverable becomes a
depth-localized UPPER BOUND on dv/v in 0-864 m, which no one has published for
Parkfield either. Report the achieved precision and stop; do not tune thresholds
until something passes.
--------------------------------------------------------------------------------
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.signal import decimate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import moveout_test as MT                                       # noqa: E402
from dvv_core import sub_sample_delay, slope_dvv                # noqa: E402

CACHE_HF = os.path.join(HERE, 'cache_hf')
FS_TARGET = 500.0
BAND = (5.0, 80.0)          # frozen: 5-80 Hz, the usable band per hf_snr_test
GAUGE_CH = 16               # 16.335 m gauge / 1.021 m spacing
WIN_P = (-0.2, 1.5)         # about the picked P, seconds
DELAY_MAX_LAG = 0.05        # per-channel residual search; the pair is pre-aligned
MIN_COH = 0.6               # phase-stability floor from sub_sample_delay
INJECT = [0.001, 0.003, 0.010]


def load(tag):
    """Load one event at the common rate. Some files are 5000 Hz, not 500."""
    f = os.path.join(CACHE_HF, f'{tag}.npz')
    if not os.path.exists(f):
        return None
    d = np.load(f)
    X, fs = d['X'].astype(np.float64), float(d['fs'])
    if fs > FS_TARGET * 1.5:
        q = int(round(fs / FS_TARGET))
        if abs(fs / q - FS_TARGET) > 1:
            return None
        X = decimate(X, q, axis=-1, ftype='fir', zero_phase=True)
        fs = fs / q
    if abs(fs - FS_TARGET) > 1:
        return None
    return X, fs


def prep_pair(tag_a, tag_b):
    """Both events, bandpassed, on common channels, with ONE shared alignment."""
    ra, rb = load(tag_a), load(tag_b)
    if ra is None or rb is None:
        return None
    MT.BAND = BAND
    MT.CACHE = CACHE_HF
    ga, _ = MT.prep(tag_a)
    gb, _ = MT.prep(tag_b)
    if ga is None or gb is None:
        return None
    Aa, za, fsa = ga
    Ab, zb, fsb = gb
    # prep() drops dead channels independently per event, so intersect on depth
    zc, ia, ib = np.intersect1d(za, zb, return_indices=True)
    if zc.size < 100:
        return None
    Aa, Ab = Aa[ia], Ab[ib]
    if abs(fsa - fsb) > 1:
        return None

    pa, tpa, _sa, _S, _ps = MT.slant_scan(Aa, zc, fsa)
    pb, tpb, _sb, _S2, _p2 = MT.slant_scan(Ab, zc, fsb)
    pbar = 0.5 * (pa + pb)
    # SAME shift to both -- see docstring point 1
    Ba = MT.align(Aa, zc, fsa, pbar)
    Bb = MT.align(Ab, zc, fsb, pbar)
    return dict(A=Ba, B=Bb, z=zc, fs=fsa, p=pbar, tp=0.5 * (tpa + tpb),
                pa=pa, pb=pb, raw_A=Aa, raw_B=Ab)


def per_channel_delays(P, invert_eps=None, reverse_b=False):
    """dt(z) between the two events, channel by channel. Returns (dt, tau, coh)."""
    A, B, z, fs, p = P['A'], P['B'], P['z'], P['fs'], P['p']
    if invert_eps is not None:
        # inject a uniform fractional velocity change: shift -eps*p*z, which is
        # exactly align() with slowness -eps*p
        B = MT.align(B, z, fs, -invert_eps * p)
    if reverse_b:
        B = B[:, ::-1]
    i0 = int((P['tp'] + WIN_P[0]) * fs)
    i1 = int((P['tp'] + WIN_P[1]) * fs)
    if i0 < 0 or i1 > A.shape[1]:
        return None, None, None
    dt = np.full(A.shape[0], np.nan)
    coh = np.zeros(A.shape[0])
    for k in range(A.shape[0]):
        d_, c_ = sub_sample_delay(A[k, i0:i1], B[k, i0:i1], fs, BAND,
                                  max_lag_s=DELAY_MAX_LAG)
        dt[k], coh[k] = d_, c_
    tau = p * z                      # one-way travel time from measured slowness
    return dt, tau, coh


def measure(P, **kw):
    dt, tau, coh = per_channel_delays(P, **kw)
    if dt is None:
        return dict(dvv=np.nan, err=np.nan, n=0, rms=np.nan)
    m = np.isfinite(dt) & (coh > MIN_COH)
    if m.sum() < 20:
        return dict(dvv=np.nan, err=np.nan, n=int(m.sum()), rms=np.nan)
    dvv, err, n, rms = slope_dvv(dt[m], tau[m], weights=coh[m],
                                 block=GAUGE_CH, n_boot=400)
    return dict(dvv=dvv, err=err, n=n, rms=rms, tau_span=float(np.ptp(tau[m])))


def main():
    mt = pd.read_csv(os.path.join(HERE, 'moveout_test.csv'))
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    have = {f[:-4] for f in os.listdir(CACHE_HF) if f.endswith('.npz')}
    mt = mt[mt.i.map(lambda k: ev.tag[k] in have) &
            mt.j.map(lambda k: ev.tag[k] in have)].copy()
    rep = mt[mt.is_cand]
    ctl = mt[~mt.is_cand]
    print(f'{len(rep)} repeater pairs, {len(ctl)} control pairs usable at 500 Hz')
    print(f'band {BAND[0]:.0f}-{BAND[1]:.0f} Hz, bootstrap block {GAUGE_CH} ch\n',
          flush=True)

    rows = []
    cache = {}

    def get(i, j):
        key = (i, j)
        if key not in cache:
            cache[key] = prep_pair(ev.tag[i], ev.tag[j])
        return cache[key]

    print('--- repeater pairs, measured ---', flush=True)
    print(f'{"events":>24}{"dv/v %":>10}{"err %":>9}{"nch":>6}'
          f'{"rms ms":>9}{"tau s":>8}')
    for _, r in rep.iterrows():
        P = get(int(r.i), int(r.j))
        if P is None:
            continue
        m = measure(P)
        lab = f'{ev.tag[int(r.i)][3:11]}/{ev.tag[int(r.j)][3:11]}'
        print(f'{lab:>24}{100*m["dvv"]:10.4f}{100*m["err"]:9.4f}{m["n"]:6d}'
              f'{1e3*m.get("rms",np.nan):9.3f}{m.get("tau_span",np.nan):8.3f}',
              flush=True)
        rows.append(dict(kind='repeater', i=int(r.i), j=int(r.j), lab=lab, **m))

    print('\n--- G3 acausal floor (event B time-reversed) ---', flush=True)
    for _, r in rep.iterrows():
        P = get(int(r.i), int(r.j))
        if P is None:
            continue
        m = measure(P, reverse_b=True)
        rows.append(dict(kind='acausal', i=int(r.i), j=int(r.j),
                         lab=f'{ev.tag[int(r.i)][3:11]}R', **m))
    ac = [x['dvv'] for x in rows if x['kind'] == 'acausal' and np.isfinite(x['dvv'])]
    if ac:
        print(f'  n={len(ac)}  median {100*np.median(ac):+.4f} %  '
              f'scatter {100*np.std(ac):.4f} %')

    print('\n--- G1 random-pair null ---', flush=True)
    for _, r in ctl.iterrows():
        P = get(int(r.i), int(r.j))
        if P is None:
            continue
        m = measure(P)
        rows.append(dict(kind='null', i=int(r.i), j=int(r.j), lab='ctl', **m))
    nu = [x['dvv'] for x in rows if x['kind'] == 'null' and np.isfinite(x['dvv'])]
    sigma_null = float(np.std(nu)) if len(nu) > 3 else np.nan
    if nu:
        mad = 1.4826 * np.median(np.abs(np.array(nu) - np.median(nu)))
        print(f'  n={len(nu)}  median {100*np.median(nu):+.4f} %  '
              f'std {100*sigma_null:.4f} %  MAD {100*mad:.4f} %')

    print('\n--- G2 synthetic recovery (inject into event B) ---', flush=True)
    print(f'{"injected %":>12}{"recovered %":>14}{"err %":>9}{"bias %":>9}'
          f'{"npairs":>8}')
    recov = {}
    for eps in INJECT:
        got, ers = [], []
        for _, r in rep.iterrows():
            P = get(int(r.i), int(r.j))
            if P is None:
                continue
            m = measure(P, invert_eps=eps)
            base = next((x['dvv'] for x in rows if x['kind'] == 'repeater'
                         and x['i'] == int(r.i) and x['j'] == int(r.j)), np.nan)
            if np.isfinite(m['dvv']) and np.isfinite(base):
                got.append(m['dvv'] - base)      # differential removes any real dv/v
                ers.append(m['err'])
        if got:
            recov[eps] = (float(np.median(got)), float(np.median(ers)), len(got))
            print(f'{100*eps:12.3f}{100*recov[eps][0]:14.4f}'
                  f'{100*recov[eps][1]:9.4f}{100*(recov[eps][0]-eps):9.4f}'
                  f'{recov[eps][2]:8d}')

    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(HERE, 'interval_dvv_gate.csv'), index=False)

    print('\n' + '=' * 74)
    print('VERDICT against criteria fixed before running')
    g1 = np.isfinite(sigma_null) and sigma_null <= 0.001
    print(f'  G1 null sigma <= 0.100 % : {100*sigma_null:.4f} %   '
          f'{"PASS" if g1 else "FAIL"}')
    g2 = True
    for eps in INJECT:
        if eps not in recov:
            g2 = False; continue
        bias, er, _ = recov[eps]
        ok = abs(bias - eps) <= max(2 * er, 0.2 * eps)
        g2 &= ok
        print(f'  G2 recover {100*eps:.1f} %      : bias '
              f'{100*(bias-eps):+.4f} %  {"PASS" if ok else "FAIL"}')
    g3 = bool(ac) and abs(np.median(ac)) <= 2 * (sigma_null if np.isfinite(sigma_null) else np.inf)
    print(f'  G3 acausal ~ null        : {"PASS" if g3 else "FAIL"}')
    print()
    if g1 and g2 and g3:
        print('  -> GATE PASSED. Proceed to D1/D2 with band and channel range frozen.')
    else:
        print('  -> GATE NOT PASSED. The deliverable is a depth-localized UPPER')
        print('     BOUND on dv/v over 0-864 m. Report the achieved precision.')
        print('     Do not retune thresholds to force a pass.')

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for kind, c in [('repeater', 'C3'), ('null', '0.6'), ('acausal', 'C0')]:
        v = D[(D.kind == kind) & np.isfinite(D.dvv)].dvv * 100
        if len(v):
            ax[0].hist(v, bins=18, alpha=.6, color=c, label=f'{kind} (n={len(v)})')
    ax[0].set(xlabel='dv/v (%)', ylabel='pairs', title='A  signal vs null vs acausal')
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    if recov:
        xs = [100 * e for e in recov]
        ys = [100 * recov[e][0] for e in recov]
        es = [100 * recov[e][1] for e in recov]
        ax[1].errorbar(xs, ys, yerr=es, fmt='o-', color='C2')
        lim = [0, max(xs) * 1.1]
        ax[1].plot(lim, lim, 'k--', lw=1, label='y = x')
        ax[1].set(xlabel='injected dv/v (%)', ylabel='recovered (%)',
                  title='B  G2 synthetic recovery')
        ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    rp = D[(D.kind == 'repeater') & np.isfinite(D.dvv)]
    if len(rp):
        ax[2].errorbar(range(len(rp)), 100 * rp.dvv, yerr=100 * rp.err,
                       fmt='o', color='C3')
        ax[2].axhline(0, color='k', lw=1)
        if np.isfinite(sigma_null):
            ax[2].axhspan(-100 * sigma_null, 100 * sigma_null, color='0.85',
                          label='null 1 sigma')
            ax[2].legend(fontsize=8)
        ax[2].set(xlabel='pair', ylabel='dv/v (%)', title='C  per-pair')
        ax[2].grid(alpha=.3)
    fig.tight_layout()
    p = os.path.join(HERE, 'interval_dvv_gate.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
