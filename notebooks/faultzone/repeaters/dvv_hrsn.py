"""
G3: does a velocity change exist at all? Ask the better instrument first.

THE LOGIC OF RUNNING THIS BEFORE ANY DAS WORK. HRSN borehole seismometers have far
better SNR than DAS, three components, 250 Hz, and are the array the Parkfield
repeater catalogs were built on. If the best available instrument sees no velocity
change between two occurrences of a repeating earthquake, DAS will not either, and
the depth-profile work has nothing to resolve. One day spent here beats three weeks
spent downstream. Data is already cached in hrsn_cache/ -- no downloads.

WHAT A NULL WOULD MEAN, which is why this is informative either way. HRSN sensors
sit at 251-284 m depth (station_geometry.py). Li & Ben-Zion 2023 put the seasonal
signal in the top few tens of metres with peak sensitivity ~17 m. So HRSN is BELOW
the layer where the signal is expected. A null at HRSN is therefore consistent with
-- even supporting -- shallow confinement, and would make the shallow DAS channels
the whole point rather than a bonus. A positive result at 250 m depth would instead
mean the change penetrates further than Li & Ben-Zion's sensitivity reaches.

Either outcome is a result. What would kill the project is a null everywhere.

CONTROLS RUN HERE, NOT LATER.
  (a) random non-repeating pairs must give incoherent stretch -- they share no
      source, so any apparent dv/v is the method's noise floor
  (b) the coda split into halves must agree -- a real medium change affects both
  (c) seasonal phase vs elapsed time -- reported per pair so the discriminator is
      available even at this stage

SELECTION CAUTION. Pairs were selected on direct-arrival similarity, not coda
similarity, deliberately: a velocity change decorrelates the coda, so selecting on
coda CC would preferentially discard the pairs carrying the largest signal.
"""
import os
import sys
import numpy as np
import pandas as pd
from obspy import read
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dvv_core import (bandpass, stretch_dvv, stretch_dvv_bootstrap,   # noqa: E402
                      seasonal_phase_days, bulk_align)

WF = os.path.join(HERE, 'hrsn_cache')
BAND = (5.0, 20.0)
PRE_S = 5.0                      # cache starts 5 s before catalog origin
# lapse windows relative to catalog origin. The direct arrival is at ~1-3 s for
# these distances, so all three are coda. Overlapping on purpose: a real medium
# change is stable across lapse time, a processing artifact usually is not.
# chosen from coda_window_survey.py: coda CC 0.85-0.96 and SNR 7-254
# across the whole 0.5-14 s range once pairs are aligned
LAPSE = [(1.0, 6.0), (2.0, 12.0), (4.0, 14.0)]
EPS_MAX = 0.04                   # +/-4%, covers Li & Ben-Zion's amplitude
MIN_CC = 0.5                     # per-station coda CC floor to be counted


def load(tag):
    f = os.path.join(WF, f'{tag}.mseed')
    if not os.path.exists(f):
        return None
    try:
        return read(f)
    except Exception:
        return None


def full(st):
    """Per-station bandpassed full traces. Windowing happens after alignment."""
    out = {}
    if st is None:
        return out
    for tr in st:
        try:
            fs = tr.stats.sampling_rate
            d = bandpass(tr.data.astype(float) - float(np.mean(tr.data)), fs, BAND)
            if np.all(d == 0) or not np.all(np.isfinite(d)):
                continue
            out[tr.stats.station] = (d, fs)
        except Exception:
            continue
    return out


def _win(d, fs, w):
    """Extract one lapse window, normalised."""
    try:
        i0, i1 = int((PRE_S + w[0]) * fs), int((PRE_S + w[1]) * fs)
        if i1 > d.size or i0 < 0:
            return None
        seg = d[i0:i1]
        return seg if (np.any(seg) and np.all(np.isfinite(seg))) else None
    except Exception:
        return None


def pair_dvv(A, B, w, half=None):
    """dv/v for one pair over lapse window w, bootstrapped across stations."""
    common = sorted(set(A) & set(B))
    ta, tb, fs = [], [], None
    for s in common:
        a, fa = A[s]; b, fb = B[s]
        if fa != fb:
            continue
        # ALIGN ON THE FULL TRACE FIRST. Catalog origin times carry 0.1-0.5 s of
        # error; at 5-20 Hz that drives a windowed correlation to zero. The first
        # pass of this script omitted it and returned coda CC ~0.2 with ~1% dv/v
        # on pairs that correlate at 0.99 once aligned.
        b, _lag = bulk_align(a, b, fa, max_lag_s=2.0)
        aw, bw = _win(a, fa, w), _win(b, fa, w)
        if aw is None or bw is None:
            continue
        n = min(aw.size, bw.size)
        aw, bw = aw[:n], bw[:n]
        if half == 'first':
            aw, bw = aw[:n // 2], bw[:n // 2]
        elif half == 'second':
            aw, bw = aw[n // 2:], bw[n // 2:]
        ta.append(aw); tb.append(bw); fs = fa
    if not ta:
        return np.nan, np.nan, 0, np.nan
    return stretch_dvv_bootstrap(ta, tb, fs, eps_max=EPS_MAX)


def main():
    R = pd.read_csv(os.path.join(HERE, 'hrsn_control.csv'))
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    ev['time'] = pd.to_datetime(ev['time'], utc=True, format='mixed')
    R['t_i'] = pd.to_datetime(R['t_i'], utc=True, format='mixed')
    R['t_j'] = pd.to_datetime(R['t_j'], utc=True, format='mixed')

    cand = R[R.is_cand & (R.hrsn > 0.78)].copy()
    ctrl = R[~R.is_cand].copy()
    print(f'{len(cand)} confirmed repeater pairs (HRSN CC > 0.78)')
    print(f'{len(ctrl)} random control pairs\n', flush=True)

    rows = []
    for label, D, is_rep in [('repeater', cand, True), ('control', ctrl, False)]:
        for _, r in D.iterrows():
            ti, tj = ev.tag[int(r.i)], ev.tag[int(r.j)]
            sa, sb = load(ti), load(tj)
            if sa is None or sb is None:
                continue
            A, B = full(sa), full(sb)
            ph = seasonal_phase_days(r.t_i, r.t_j)
            for w in LAPSE:
                d, e, n, cc = pair_dvv(A, B, w)
                h1 = pair_dvv(A, B, w, 'first')[0]
                h2 = pair_dvv(A, B, w, 'second')[0]
                rows.append(dict(kind=label, is_rep=is_rep, i=int(r.i), j=int(r.j),
                                 lapse=f'{w[0]:.0f}-{w[1]:.0f}',
                                 dvv=d, err=e, nsta=n, coda_cc=cc,
                                 dvv_h1=h1, dvv_h2=h2,
                                 dt_days=r.dt_days, phase_days=ph,
                                 hrsn_cc=r.hrsn, das_cc=r.das,
                                 t_i=r.t_i, t_j=r.t_j))
            if is_rep:
                print(f'  {r.t_i:%Y-%m-%d}/{r.t_j:%Y-%m-%d} '
                      f'dt={r.dt_days:5.0f}d phase={ph:5.1f}d  done', flush=True)

    D = pd.DataFrame(rows)
    if D.empty:
        print('no measurements'); return
    D.to_csv(os.path.join(HERE, 'dvv_hrsn.csv'), index=False)

    rep = D[D.is_rep]; con = D[~D.is_rep]
    print('\nREPEATER PAIRS -- dv/v by lapse window')
    print(f'{"events":>24}{"lapse":>9}{"dv/v %":>9}{"err %":>8}{"nsta":>6}'
          f'{"codaCC":>8}{"half1%":>8}{"half2%":>8}')
    for _, r in rep.sort_values(['t_i', 'lapse']).iterrows():
        print(f'{r.t_i:%Y-%m-%d}/{r.t_j:%Y-%m-%d}{r.lapse:>9}'
              f'{100*r.dvv:9.3f}{100*r.err:8.3f}{int(r.nsta):6d}{r.coda_cc:8.3f}'
              f'{100*r.dvv_h1:8.3f}{100*r.dvv_h2:8.3f}')

    print('\n(a) CONTROL -- random non-repeating pairs (the method noise floor)')
    for lp in sorted(con.lapse.unique()):
        c = con[con.lapse == lp].dvv.dropna()
        r_ = rep[rep.lapse == lp].dvv.dropna()
        print(f'  lapse {lp:>7}s: control |dv/v| median '
              f'{100*np.median(np.abs(c)):.3f}%  spread {100*c.std():.3f}%  '
              f'(n={c.size})   repeaters |dv/v| median '
              f'{100*np.median(np.abs(r_)):.3f}%')

    print('\n(b) CODA-HALVES CONSISTENCY (a real change affects both halves)')
    ok = rep.dropna(subset=['dvv_h1', 'dvv_h2'])
    if len(ok) > 2:
        dh = np.abs(ok.dvv_h1 - ok.dvv_h2)
        print(f'  median |half1 - half2| = {100*np.median(dh):.3f}% '
              f'vs median |dv/v| = {100*np.median(np.abs(ok.dvv)):.3f}%')
        rr = np.corrcoef(ok.dvv_h1, ok.dvv_h2)[0, 1]
        print(f'  r(half1, half2) = {rr:+.2f}  '
              f'({"consistent" if rr > 0.5 else "NOT consistent"})')

    print('\n(c) SEASONAL vs SECULAR')
    m = rep[rep.lapse == f'{LAPSE[1][0]:.0f}-{LAPSE[1][1]:.0f}'].dropna(subset=['dvv'])
    if len(m) >= 4:
        for xk, nm in [('phase_days', 'seasonal phase'), ('dt_days', 'elapsed time')]:
            rr = np.corrcoef(m[xk], np.abs(m.dvv))[0, 1]
            print(f'  r(|dv/v|, {nm:>14}) = {rr:+.2f}')
        print('  seasonal predicts phase; drift predicts elapsed time')

    print('\n(d) SELECTION-BIAS CHECK (control 8)')
    if len(m) >= 4:
        rr = np.corrcoef(m.hrsn_cc, np.abs(m.dvv))[0, 1]
        print(f'  r(selection CC, |dv/v|) = {rr:+.2f}  '
              f'{"<- negative would mean the pair set is censored against signal" if rr < -0.3 else "(no evidence of censoring)"}')

    print('\nG3 VERDICT')
    mid = f'{LAPSE[1][0]:.0f}-{LAPSE[1][1]:.0f}'
    rm = rep[rep.lapse == mid].dropna(subset=['dvv'])
    cm = con[con.lapse == mid].dropna(subset=['dvv'])
    if rm.empty:
        print('  no repeater measurements converged -- estimator or windows wrong')
    else:
        floor = np.percentile(np.abs(cm.dvv), 90) if len(cm) > 3 else np.nan
        n_sig = int(np.sum(np.abs(rm.dvv) > np.maximum(2 * rm.err, floor)))
        print(f'  control-derived floor (90th pct of |dv/v| on random pairs): '
              f'{100*floor:.3f}%')
        print(f'  repeater pairs exceeding max(2*err, floor): {n_sig}/{len(rm)}')
        if n_sig >= 2:
            print('  -> PASS. A velocity change is resolvable at HRSN (251-284 m). '
                  'It\n     penetrates below the top-tens-of-metres layer; proceed '
                  'to G4/G5.')
        else:
            print('  -> NULL at HRSN depth. Consistent with the change being '
                  'confined\n     ABOVE 250 m, which makes the shallow DAS channels '
                  'the whole point.\n     Not a failure: proceed to G5, but the '
                  'expectation is shallow-only.')

    # ------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    a = ax[0]
    for lp, mk in zip(sorted(rep.lapse.unique()), 'os^'):
        g = rep[rep.lapse == lp].dropna(subset=['dvv'])
        a.errorbar(g.dt_days, 100 * g.dvv, yerr=100 * g.err, fmt=mk, ms=6,
                   capsize=3, label=f'lapse {lp} s')
    if len(cm):
        a.axhspan(-100 * np.percentile(np.abs(cm.dvv), 90),
                  100 * np.percentile(np.abs(cm.dvv), 90),
                  color='0.85', label='control floor')
    a.axhline(0, color='k', lw=0.8)
    a.set(xlabel='elapsed time (days)', ylabel='dv/v (%)',
          title='A  dv/v vs elapsed time\n(drift would trend here)')
    a.legend(fontsize=7); a.grid(alpha=0.3)

    a = ax[1]
    for lp, mk in zip(sorted(rep.lapse.unique()), 'os^'):
        g = rep[rep.lapse == lp].dropna(subset=['dvv'])
        a.errorbar(g.phase_days, 100 * g.dvv, yerr=100 * g.err, fmt=mk, ms=6,
                   capsize=3, label=f'lapse {lp} s')
    a.axhline(0, color='k', lw=0.8)
    a.set(xlabel='seasonal phase separation (days)', ylabel='dv/v (%)',
          title='B  dv/v vs seasonal phase\n(seasonal would trend here)')
    a.legend(fontsize=7); a.grid(alpha=0.3)

    a = ax[2]
    a.hist(100 * con.dvv.dropna(), bins=20, alpha=0.6, color='0.6',
           label='random pairs (null)')
    for v in rep[rep.lapse == mid].dvv.dropna():
        a.axvline(100 * v, color='C3', lw=1.4)
    a.set(xlabel='dv/v (%)', ylabel='count',
          title='C  Repeaters (lines) vs null (bars)')
    a.legend(fontsize=8); a.grid(alpha=0.3)

    fig.suptitle('G3: is a velocity change resolvable on HRSN borehole '
                 'seismometers (251-284 m)?', fontsize=12)
    fig.tight_layout()
    p = os.path.join(HERE, 'dvv_hrsn.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
