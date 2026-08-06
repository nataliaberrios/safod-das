"""
Re-check dv/v at 500 Hz with moveout-aligned stacks. The -3.6% verdict is stale.

WHY THIS EXISTS. dv/v from repeaters was closed as a negative result at +/-3.6%,
"70x worse than coda stretching". That number was measured on 100 Hz caches with
FLAT channel stacking, i.e. through both defects since fixed:

    flat stack across a 0.31 s moveout   ~27 dB of array gain destroyed
    decimation to 100 Hz + 40 Hz lowpass  80% of the band discarded

Repeater CC under that processing was 0.680; corrected it is 0.956, at parity with
HRSN's 0.954. Stretching precision scales with sqrt((1-C)/C) and with bandwidth, so
+/-3.6% is not the current answer to anything. It is simply unmeasured at the
working configuration.

The estimator is dvv_core.stretch_dvv_bootstrap, unchanged, so this is a like-for-
like re-measurement and not a new method.

PREDICTION, REGISTERED BEFORE RUNNING. The (1-C) term alone gives
sqrt(0.32/0.044) = 2.7x. Bandwidth 5-20 -> 5-80 Hz gives up to sqrt(5) = 2.2x if
the coda carries signal there, which hf_moveout_test shows it partly does
(repeater CC 0.66-0.87 in 20-80 Hz alone). Combined ceiling ~6x, so
+/-3.6% -> ~+/-0.6%. That is still well short of the ~0.1% needed for a seasonal
signal, so the honest expectation is IMPROVED BUT STILL INSUFFICIENT. If it lands
near 0.6% the negative result stands with a better bound; if it lands near 0.1%
the direction reopens.

Stating that in advance so a mediocre number is not talked up afterwards.

CONTROLS, same as G3 required them:
  * random non-repeating pairs -> the method's own noise floor
  * coda split into halves -> a real medium change affects both
  * seasonal phase vs elapsed time -> reported per pair
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
import moveout_test as MT                                        # noqa: E402
from dvv_core import stretch_dvv_bootstrap, bulk_align           # noqa: E402

CACHE_HF = os.path.join(HERE, 'cache_hf')
FS_TARGET = 500.0
BANDS = [(5.0, 20.0), (5.0, 40.0), (5.0, 80.0)]
PRE_S = 5.0
LAPSE = [(1.0, 6.0), (2.0, 12.0), (4.0, 14.0)]   # same windows G3 used
EPS_MAX = 0.04
# BOOTSTRAP SAMPLE SIZE. The first run used 4 subarrays, so stretch_dvv_bootstrap
# resampled 4 traces and the resulting +/-0.033 % error bar came from 4 samples --
# not trustworthy, and flagged as provisional in METHODS_STATUS 18.2. 24 subarrays
# of ~29 channels each still exceeds one 16-channel gauge length per subarray, so
# the samples stay quasi-independent while the bootstrap becomes meaningful.
SUBARRAYS = 24


def load(tag):
    f = os.path.join(CACHE_HF, f'{tag}.npz')
    if not os.path.exists(f):
        return None
    d = np.load(f)
    X, fs = d['X'].astype(np.float64), float(d['fs'])
    if fs > FS_TARGET * 1.5:
        q = int(round(fs / FS_TARGET))
        X = decimate(X, q, axis=-1, ftype='fir', zero_phase=True)
        fs = fs / q
    return (X, fs) if abs(fs - FS_TARGET) < 1 else None


def aligned_subarray_stacks(tag, band, nsub=SUBARRAYS):
    """Moveout-aligned stacks, one per depth subarray. This is the whole point:
    the old measurement stacked flat and lost ~27 dB."""
    if load(tag) is None:
        return None, None
    MT.BAND = band
    MT.CACHE = CACHE_HF
    got, _ = MT.prep(tag)
    if got is None:
        return None, None
    A, z, fs = got
    p, tp, _s, _S, _ps = MT.slant_scan(A, z, fs)
    B = MT.align(A, z, fs, p)
    edges = np.linspace(0, B.shape[0], nsub + 1).astype(int)
    stacks = []
    for k in range(nsub):
        seg = B[edges[k]:edges[k + 1]]
        if seg.shape[0] < 20:
            stacks.append(None); continue
        s = seg.mean(axis=0)
        s -= s.mean()
        n = np.sqrt(np.sum(s ** 2))
        stacks.append(s / n if np.isfinite(n) and n > 0 else None)
    return stacks, fs


def pair_dvv(sa, sb, fs, win, half=None):
    ta, tb = [], []
    for a, b in zip(sa, sb):
        if a is None or b is None:
            continue
        b, _lag = bulk_align(a, b, fs, max_lag_s=2.0)   # NECESSARY, NOT OPTIONAL
        i0, i1 = int((PRE_S + win[0]) * fs), int((PRE_S + win[1]) * fs)
        if i1 > min(a.size, b.size):
            continue
        aw, bw = a[i0:i1], b[i0:i1]
        n = min(aw.size, bw.size)
        aw, bw = aw[:n], bw[:n]
        if half == 'first':
            aw, bw = aw[:n // 2], bw[:n // 2]
        elif half == 'second':
            aw, bw = aw[n // 2:], bw[n // 2:]
        ta.append(aw); tb.append(bw)
    if not ta:
        return np.nan, np.nan, 0, np.nan
    return stretch_dvv_bootstrap(ta, tb, fs, eps_max=EPS_MAX)


def main():
    mt = pd.read_csv(os.path.join(HERE, 'moveout_test.csv'))
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    have = {f[:-4] for f in os.listdir(CACHE_HF) if f.endswith('.npz')}
    mt = mt[mt.i.map(lambda k: ev.tag[k] in have) &
            mt.j.map(lambda k: ev.tag[k] in have)]
    rep, ctl = mt[mt.is_cand], mt[~mt.is_cand]
    print(f'{len(rep)} repeater pairs, {len(ctl)} control pairs at 500 Hz')
    print(f'{SUBARRAYS} depth subarrays, moveout-aligned before stacking')
    print('reference: 100 Hz flat-stack result was +/-3.6 %\n', flush=True)

    rows = []
    for band in BANDS:
        cache = {}

        def stacks(k):
            if k not in cache:
                cache[k] = aligned_subarray_stacks(ev.tag[k], band)
            return cache[k]

        print(f'=== band {band[0]:.0f}-{band[1]:.0f} Hz ===', flush=True)
        for kind, D in (('repeater', rep), ('control', ctl)):
            for _, r in D.iterrows():
                sa, fsa = stacks(int(r.i))
                sb, fsb = stacks(int(r.j))
                if sa is None or sb is None:
                    continue
                for win in LAPSE:
                    d_, e_, n_, cc_ = pair_dvv(sa, sb, fsa, win)
                    h1 = pair_dvv(sa, sb, fsa, win, 'first')[0]
                    h2 = pair_dvv(sa, sb, fsa, win, 'second')[0]
                    rows.append(dict(kind=kind, band=f'{band[0]:.0f}-{band[1]:.0f}',
                                     lapse=f'{win[0]:.0f}-{win[1]:.0f}',
                                     i=int(r.i), j=int(r.j), dvv=d_, err=e_,
                                     nsub=n_, coda_cc=cc_, h1=h1, h2=h2,
                                     dt_days=r.dt_days, hrsn=r.hrsn))
        cur = pd.DataFrame([x for x in rows
                            if x['band'] == f'{band[0]:.0f}-{band[1]:.0f}'])
        if cur.empty:
            print('  nothing\n'); continue
        print(f'  {"lapse":>8}{"repN":>6}{"|dvv| med %":>13}{"err med %":>11}'
              f'{"codaCC":>8}{"ctrlN":>7}{"ctrl |dvv| %":>14}')
        for lp in sorted(cur.lapse.unique()):
            R = cur[(cur.lapse == lp) & (cur.kind == 'repeater')].dropna(subset=['dvv'])
            C = cur[(cur.lapse == lp) & (cur.kind == 'control')].dropna(subset=['dvv'])
            print(f'  {lp:>8}{len(R):6d}{100*np.median(np.abs(R.dvv)):13.4f}'
                  f'{100*np.median(R.err):11.4f}{np.median(R.coda_cc):8.3f}'
                  f'{len(C):7d}'
                  f'{100*np.median(np.abs(C.dvv)) if len(C) else np.nan:14.4f}')
        print(flush=True)

    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(HERE, 'recheck_dvv_500.csv'), index=False)

    print('=' * 72)
    print('VERDICT against the prediction registered in the docstring')
    R = D[(D.kind == 'repeater')].dropna(subset=['err'])
    C = D[(D.kind == 'control')].dropna(subset=['dvv'])
    if R.empty:
        print('  no measurements'); return
    best = R.loc[R.err.idxmin()]
    prec = float(R.err.median())
    floor = float(np.median(np.abs(C.dvv))) if len(C) else np.nan
    print(f'  median bootstrap error       : {100*prec:.4f} %')
    print(f'  best configuration           : {best.band} Hz, lapse {best.lapse} s,'
          f' err {100*best.err:.4f} %')
    print(f'  control (random-pair) floor  : {100*floor:.4f} %')
    print(f'  100 Hz flat-stack reference  : 3.6000 %')
    print(f'  improvement factor           : {3.6/(100*prec):.1f}x'
          if prec > 0 else '')
    print()
    if floor <= 0.15e-2 * 100:
        pass
    if 100 * floor <= 0.15:
        print('  -> control floor at or below 0.15 %: the direction REOPENS.')
    elif 100 * floor <= 0.8:
        print('  -> improved but still insufficient for a ~0.1 % seasonal signal,')
        print('     which is what the docstring predicted. The negative result')
        print('     stands, now with a much tighter bound.')
    else:
        print('  -> no material improvement. Negative result stands as before.')

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    for band in D.band.unique():
        s = D[(D.band == band) & (D.kind == 'repeater')].dropna(subset=['err'])
        if len(s):
            ax[0].scatter([band] * len(s), 100 * s.err, alpha=.6, label=band)
    ax[0].axhline(3.6, color='C3', ls='--', label='100 Hz flat stack')
    ax[0].set(ylabel='bootstrap error (%)', yscale='log',
              title='A  precision vs band')
    ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
    for kind, c in (('repeater', 'C3'), ('control', '0.6')):
        v = D[(D.kind == kind)].dropna(subset=['dvv']).dvv * 100
        if len(v):
            ax[1].hist(v, bins=25, alpha=.6, color=c, label=f'{kind} n={len(v)}')
    ax[1].set(xlabel='dv/v (%)', title='B  signal vs random-pair floor')
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    fig.tight_layout()
    p = os.path.join(HERE, 'recheck_dvv_500.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
