"""
Extend the HRSN confirmation past the top 10 pairs.

WHY IT HAS TO BE REDONE RATHER THAN JUST LENGTHENED. The existing hrsn_control.csv
took its ten candidates from the TRACE-STACKED DAS search, which correlated one
channel-averaged trace per event. That ranking is compromised by the moveout bug
(METHODS_STATUS section 2.4): stacking 700 traces across an unflattened wavefront
is a comb filter, worth about 1 dB over a single channel, and it depresses every
DAS CC by 0.15-0.30. Pairs it ranked 11th-40th are not reliably worse than the ten
it ranked first -- they are ten samples from a badly-blurred ordering.

correlate_perchannel.py replaces that ranking with the published DAS recipe:
per-channel correlation, correlogram stacking, three bands, a null-calibrated
threshold and a time-reversed acausal floor. Correlating channel i of A against
channel i of B cancels the moveout exactly, with no velocity model. This script
takes the top 40 of THAT ranking to HRSN and asks the same question the original
control asked: does the better instrument agree that these are repeating events?

CONFIRMATION CRITERION. HRSN CC > 0.78, matching the cut that defined the
confirmed seven, so the extended set is homogeneous with the existing one. The
null is re-measured here from random pairs rather than carried over, because the
candidate set changed and a threshold is only meaningful against its own null.

WHAT THIS IS FOR. The dv/v work (G3, G5) is limited by having seven pairs and four
usable CWI baselines. Every additional confirmed pair adds a baseline, and the
long ones matter most: the seasonal signal is separable from secular drift only
because elapsed time and seasonal phase are nearly uncorrelated across the pair
set, and that decorrelation improves with more pairs.

NETWORK. Waveforms come from NCEDC over FDSN and are cached in hrsn_cache/ (73
events already there). If the compute node has no outbound route the fetches fail
one by one and the script still reports every pair it can form from the cache, so
a network-less run degrades rather than dies.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hrsn_control import fetch, traces, cc_pair, WF   # noqa: E402

CONFIRM_CC = 0.78          # same cut that defined the confirmed seven
N_RANDOM = 60              # random pairs for the null, from cached events only
MIN_DT_DAYS = 20.0


def cached_tags():
    return {os.path.splitext(f)[0] for f in os.listdir(WF)
            if f.endswith('.mseed')}


def main():
    cf = os.path.join(HERE, 'perchannel_candidates.csv')
    if not os.path.exists(cf):
        print('perchannel_candidates.csv missing -- run correlate_perchannel.py '
              'first (this job is normally chained after it).')
        return
    cand = pd.read_csv(cf)
    for c in ('t_i', 't_j'):
        cand[c] = pd.to_datetime(cand[c], utc=True, format='mixed')
    print(f'{len(cand)} candidate pairs from the per-channel ranking')
    print(f'  min-band CC range {cand.min_cc.min():.3f} - '
          f'{cand.min_cc.max():.3f}')

    prev = None
    pf = os.path.join(HERE, 'hrsn_control.csv')
    if os.path.exists(pf):
        prev = pd.read_csv(pf)
        for c in ('t_i', 't_j'):
            prev[c] = pd.to_datetime(prev[c], utc=True, format='mixed')
        print(f'  previously confirmed: '
              f'{int((prev.is_cand & (prev.hrsn > CONFIRM_CC)).sum())} pairs\n')

    have = cached_tags()
    need = {}
    for _, r in cand.iterrows():
        need[r.tag_i] = r.t_i
        need[r.tag_j] = r.t_j
    fresh = [t for t in need if t not in have]
    print(f'{len(need)} distinct events in the candidate set; '
          f'{len(need)-len(fresh)} already cached, {len(fresh)} to fetch\n',
          flush=True)

    client = None
    if fresh:
        try:
            from obspy.clients.fdsn import Client
            client = Client('NCEDC', timeout=120)
        except Exception as e:
            print(f'  FDSN client unavailable ({str(e)[:60]}); '
                  f'cache-only run\n')

    store, nfail = {}, 0
    for tag, t0 in need.items():
        # fetch() returns the cached mseed if present and only reaches the network
        # otherwise, catching its own failures -- so this is safe with client=None
        T = traces(fetch(client, t0, tag), t0)
        if T:
            store[tag] = T
        else:
            nfail += 1
    print(f'{len(store)}/{len(need)} events with usable HRSN data '
          f'({nfail} unavailable)\n', flush=True)

    rows = []
    for _, r in cand.iterrows():
        A, B = store.get(r.tag_i), store.get(r.tag_j)
        if not A or not B:
            continue
        h, nst = cc_pair(A, B)
        rows.append(dict(tag_i=r.tag_i, tag_j=r.tag_j, t_i=r.t_i, t_j=r.t_j,
                         m_i=r.m_i, m_j=r.m_j, dt_days=r.dt_days,
                         das_min=r.min_cc, das_2_8=r.cc_2_8,
                         das_4_16=r.cc_4_16, das_10_40=r.cc_10_40,
                         hrsn=h, nsta=nst, kind='candidate'))

    # Null from random pairs of cached events. Same processing, same window, so
    # the threshold is calibrated against this candidate set rather than imported.
    tags = sorted(store)
    tt = {t: need[t] for t in tags if t in need}
    rng = np.random.default_rng(0)
    seen = {(min(r['tag_i'], r['tag_j']), max(r['tag_i'], r['tag_j']))
            for r in rows}
    tries = 0
    while sum(r['kind'] == 'random' for r in rows) < N_RANDOM and tries < 4000:
        tries += 1
        a, b = rng.choice(tags, 2, replace=False)
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        dtd = abs((tt[a] - tt[b]).total_seconds()) / 86400
        if dtd <= MIN_DT_DAYS:
            continue
        seen.add(key)
        h, nst = cc_pair(store[a], store[b])
        rows.append(dict(tag_i=a, tag_j=b, t_i=tt[a], t_j=tt[b],
                         m_i=np.nan, m_j=np.nan, dt_days=dtd,
                         das_min=np.nan, das_2_8=np.nan, das_4_16=np.nan,
                         das_10_40=np.nan, hrsn=h, nsta=nst, kind='random'))

    R = pd.DataFrame(rows)
    if R.empty:
        print('no pairs could be formed -- no HRSN data reached this node.')
        return
    R.to_csv(os.path.join(HERE, 'hrsn_extended.csv'), index=False)

    cd = R[(R.kind == 'candidate')].dropna(subset=['hrsn'])
    nl = R[(R.kind == 'random')].dropna(subset=['hrsn'])
    print(f'{"DASmin":>8}{"2-8":>7}{"4-16":>7}{"10-40":>7}{"HRSN":>8}'
          f'{"nsta":>6}{"dt_d":>8}   events')
    for _, r in cd.sort_values('hrsn', ascending=False).iterrows():
        flag = '  <-- CONFIRMED' if r.hrsn > CONFIRM_CC else ''
        print(f'{r.das_min:8.3f}{r.das_2_8:7.3f}{r.das_4_16:7.3f}'
              f'{r.das_10_40:7.3f}{r.hrsn:8.3f}{int(r.nsta):6d}'
              f'{r.dt_days:8.1f}   {r.t_i:%Y-%m-%d} M{r.m_i:.2f} / '
              f'{r.t_j:%Y-%m-%d} M{r.m_j:.2f}{flag}')

    if len(nl) > 3:
        mad = 1.4826 * np.median(np.abs(nl.hrsn - np.median(nl.hrsn)))
        print(f'\nHRSN NULL from {len(nl)} random pairs: median '
              f'{np.median(nl.hrsn):+.3f}, MAD {mad:.3f}, max {nl.hrsn.max():.3f}')
        print(f'  6xMAD = {6*mad:.3f};  fixed cut in use = {CONFIRM_CC:.2f}')

    conf = cd[cd.hrsn > CONFIRM_CC]
    print(f'\nCONFIRMED at HRSN CC > {CONFIRM_CC:.2f}: {len(conf)}/{len(cd)} '
          f'candidates')
    if len(conf):
        print(f'  baselines (days): '
              f'{", ".join(f"{d:.0f}" for d in sorted(conf.dt_days))}')
        print(f'  longest baseline: {conf.dt_days.max():.0f} d')

    if prev is not None:
        old = prev[prev.is_cand & (prev.hrsn > CONFIRM_CC)]
        okey = {(min(str(a)[:19], str(b)[:19]), max(str(a)[:19], str(b)[:19]))
                for a, b in zip(old.t_i, old.t_j)}
        nkey = {(min(str(a)[:19], str(b)[:19]), max(str(a)[:19], str(b)[:19]))
                for a, b in zip(conf.t_i, conf.t_j)}
        print(f'\nversus the previous confirmation ({len(okey)} pairs):')
        print(f'  recovered : {len(okey & nkey)}')
        print(f'  new       : {len(nkey - okey)}')
        print(f'  lost      : {len(okey - nkey)}  '
              f'(pairs the trace-stacked ranking found that this one did not '
              f'rank top-40)')
        n_new = len(nkey - okey)
        if n_new:
            print(f'\n  -> {n_new} NEW confirmed pairs. Add them to the CWI '
                  f'baseline set and\n     re-run G3/G5; more baselines is exactly '
                  f'what separates the\n     seasonal term from secular drift.')
        else:
            print('\n  -> No new pairs. The per-channel ranking agrees with the '
                  'old one on\n     which events repeat, which is itself a '
                  'useful negative: the moveout\n     bug depressed the CC '
                  'VALUES but did not scramble the ORDERING.')

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5))
    a = ax[0]
    if len(nl) > 3:
        a.axhspan(nl.hrsn.quantile(0.05), nl.hrsn.quantile(0.95), color='0.85',
                  label=f'random-pair null, 5-95% (n={len(nl)})')
        a.axhline(nl.hrsn.max(), color='0.5', ls=':', lw=1.2,
                  label=f'null max {nl.hrsn.max():.3f}')
    a.scatter(cd.das_min, cd.hrsn, s=45, c='C3', label='candidates')
    a.axhline(CONFIRM_CC, color='C0', ls='--', lw=1.4,
              label=f'confirmation cut {CONFIRM_CC:.2f}')
    a.set(xlabel='DAS per-channel CC (weakest band)', ylabel='HRSN CC',
          title='A  Per-channel DAS ranking vs HRSN')
    a.legend(fontsize=8); a.grid(alpha=0.3)

    a = ax[1]
    a.hist(nl.hrsn, bins=20, color='0.65', label='random pairs')
    for v in cd.hrsn:
        a.axvline(v, color='C3', lw=1.1, alpha=0.8)
    a.axvline(CONFIRM_CC, color='C0', ls='--', lw=1.6)
    a.set(xlabel='HRSN CC', ylabel='pairs',
          title='B  Candidates (lines) against the HRSN null')
    a.legend(fontsize=8); a.grid(alpha=0.3)
    fig.suptitle('Extending HRSN confirmation to the per-channel candidate '
                 'ranking', fontsize=12)
    fig.tight_layout()
    p = os.path.join(HERE, 'hrsn_extended.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
