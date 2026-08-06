"""
Stage HF-2: does the discarded 80% of the band change any of the moveout verdict?

WHAT IS ALREADY SETTLED (METHODS_STATUS 17, at 100 Hz / 5-20 Hz):
  * moveout correction takes repeater CC 0.680 -> 0.956, HRSN parity. Confirmed.
  * it drives the random-pair null the WRONG way, +0.124 -> +0.376. Refuted.
  * detection therefore gets WORSE: d' 2.18 -> 1.74.

All of that was measured through a cache that was decimated to 100 Hz and
low-passed at 40 Hz. This script asks whether bandwidth changes it.

--------------------------------------------------------------------------------
PREDICTIONS, REGISTERED BEFORE RUNNING. The last four directions in this project
died because a physics argument was committed to before the derivation was
checked, so the falsifiable content goes here, in advance.

P1. The null must NARROW with bandwidth. Two independent signals correlated over
    a time-bandwidth product BT give CC scattered as ~1/sqrt(BT). The origin
    window is 11 s, so BT goes 165 (5-20 Hz) -> 385 (5-40) -> 825 (5-80), and the
    null MAD should fall roughly 0.255 -> 0.17 -> 0.11.
    NOTE this is a weak model: at 100 Hz the P window (2 s, BT=30) had a SMALLER
    null MAD (0.165) than the origin window (BT=165, MAD 0.283), which 1/sqrt(BT)
    does not predict. So the scaling is checked, not assumed.

P2. Repeater CC should hold up with bandwidth ONLY where there is real signal.
    hf_snr_test.py measured that directly: M3.17 at 1.8 km has SNR>3 on 855
    channels at 80-100 Hz, while M0.65 at 2.4 km is dead by ~70 Hz. So the
    high-band CC must be MAGNITUDE-DEPENDENT.

P3. Therefore d' should IMPROVE with bandwidth -- the opposite of what alignment
    did. If the null narrows (P1) while repeater CC holds (P2), discrimination
    improves for free. This is the one route back to detection below completeness
    after 17.3 closed off moveout correction as that route.

THE DECISIVE CONTROL is the 20-80 Hz band, which uses ONLY frequencies the old
cache threw away. If repeaters correlate there at all, the discarded band
demonstrably carried signal. If they do not, the 100 Hz cache lost nothing and
this whole line is closed -- which is a clean, cheap negative.

FAILURE MODE TO WATCH. More bandwidth also means more noise samples. Both real
signal and pure noise narrow the null, so a narrower null on its own proves
nothing. Only P2's magnitude dependence separates them. Do not read a d'
improvement as real unless the large events outperform the small ones.
--------------------------------------------------------------------------------

Reuses moveout_test's prep/slant_scan/align/stack/cc by overriding its module
constants, so the alignment physics has exactly one implementation.
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
import moveout_test as MT                                   # noqa: E402

CACHE_HF = os.path.join(HERE, 'cache_hf')

# (5,20) reproduces the published configuration at the new rate, isolating
# sample-rate from bandwidth. (20,80) is the previously-discarded band alone.
BANDS = [(5.0, 20.0), (5.0, 40.0), (5.0, 80.0), (20.0, 80.0)]


def robust(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan
    m = np.median(x)
    return m, 1.4826 * np.median(np.abs(x - m))


def run_band(mt, ev, band, fs_expect):
    """All four variants for one band. Returns a tidy frame."""
    MT.BAND = band
    MT.CACHE = CACHE_HF
    need = sorted(set(mt.i) | set(mt.j))

    prepped, slow = {}, {}
    for k in need:
        tag = ev.tag[k]
        got, n = MT.prep(tag)
        if got is None:
            continue
        A, z, fs = got
        if fs_expect and abs(fs - fs_expect) > 1:
            print(f'    WARNING {tag}: fs={fs:g}, expected {fs_expect:g}')
        p, tp, sem, _S, _ps = MT.slant_scan(A, z, fs)
        prepped[k] = (A, z, fs)
        slow[k] = (p, tp, sem)

    rows = []
    for _, r in mt.iterrows():
        i, j = int(r.i), int(r.j)
        if i not in prepped or j not in prepped:
            continue
        Ai, zi, fs = prepped[i]
        Aj, zj, _ = prepped[j]
        pi, tpi, _ = slow[i]
        pj, tpj, _ = slow[j]

        oi = int(MT.PRE_S * fs)
        wo = (oi + int(MT.WIN_ORIGIN[0] * fs), oi + int(MT.WIN_ORIGIN[1] * fs))
        wpi = (int((tpi + MT.WIN_P[0]) * fs), int((tpi + MT.WIN_P[1]) * fs))
        wpj = (int((tpj + MT.WIN_P[0]) * fs), int((tpj + MT.WIN_P[1]) * fs))

        Bi, Bj = MT.align(Ai, zi, fs, pi), MT.align(Aj, zj, fs, pj)
        out = dict(i=i, j=j, is_cand=bool(r.is_cand), hrsn=r.hrsn,
                   band=f'{band[0]:.0f}-{band[1]:.0f}', fs=fs,
                   mag_min=float(min(ev.mag[i], ev.mag[j])),
                   mag_max=float(max(ev.mag[i], ev.mag[j])))
        out['flat_o'] = MT.cc(MT.stack(Ai, fs, *wo), MT.stack(Aj, fs, *wo), fs)
        out['algn_o'] = MT.cc(MT.stack(Bi, fs, *wo), MT.stack(Bj, fs, *wo), fs)
        out['flat_p'] = MT.cc(MT.stack(Ai, fs, *wpi), MT.stack(Aj, fs, *wpj), fs)
        out['algn_p'] = MT.cc(MT.stack(Bi, fs, *wpi), MT.stack(Bj, fs, *wpj), fs)
        rows.append(out)
    return pd.DataFrame(rows)


def main():
    mt = pd.read_csv(os.path.join(HERE, 'moveout_test.csv'))
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    have = {f[:-4] for f in os.listdir(CACHE_HF) if f.endswith('.npz')}
    need = sorted(set(mt.i) | set(mt.j))
    miss = [k for k in need if ev.tag[k] not in have]
    print(f'{len(need)} events needed, {len(need)-len(miss)} present in cache_hf')
    if miss:
        print(f'  missing {len(miss)}: {[ev.tag[k] for k in miss][:6]}')
    mt = mt[mt.i.isin([k for k in need if ev.tag[k] in have]) &
            mt.j.isin([k for k in need if ev.tag[k] in have])].copy()
    print(f'{len(mt)} pairs usable '
          f'({int(mt.is_cand.sum())} repeater, {int((~mt.is_cand).sum())} control)\n',
          flush=True)

    allr = []
    for band in BANDS:
        print(f'=== band {band[0]:.0f}-{band[1]:.0f} Hz ===', flush=True)
        D = run_band(mt, ev, band, 500.0)
        if D.empty:
            print('  nothing\n'); continue
        allr.append(D)
        r, c = D[D.is_cand], D[~D.is_cand]
        bt = (band[1] - band[0]) * (MT.WIN_ORIGIN[1] - MT.WIN_ORIGIN[0])
        print(f'  BT (origin window) = {bt:.0f}')
        print(f'  {"variant":>14}{"rep med":>9}{"null med":>10}{"null MAD":>10}'
              f'{"d-prime":>9}{"gap":>8}')
        for k in ['flat_o', 'algn_o', 'flat_p', 'algn_p']:
            rm, _ = robust(r[k]); nm, nmad = robust(c[k])
            dp = (rm - nm) / nmad if nmad and np.isfinite(nmad) else np.nan
            gap = np.nanmin(r[k]) - np.nanmax(c[k])
            print(f'  {k:>14}{rm:9.3f}{nm:10.3f}{nmad:10.3f}{dp:9.2f}{gap:+8.3f}')
        print(flush=True)

    if not allr:
        print('no results'); return
    A = pd.concat(allr, ignore_index=True)
    A.to_csv(os.path.join(HERE, 'hf_moveout_test.csv'), index=False)

    print('\n' + '=' * 78)
    print('P1  does the null narrow with bandwidth?  (flat / origin window)')
    print(f'  {"band":>10}{"BT":>7}{"null MAD":>10}{"pred 1/sqrtBT":>15}')
    base = None
    for band in BANDS:
        s = A[(A.band == f'{band[0]:.0f}-{band[1]:.0f}') & (~A.is_cand)]
        if s.empty:
            continue
        bt = (band[1] - band[0]) * (MT.WIN_ORIGIN[1] - MT.WIN_ORIGIN[0])
        _, nmad = robust(s.flat_o)
        if base is None:
            base = (bt, nmad)
        pred = base[1] * np.sqrt(base[0] / bt)
        print(f'  {band[0]:.0f}-{band[1]:.0f}'.rjust(10)
              + f'{bt:7.0f}{nmad:10.3f}{pred:15.3f}')

    print('\nP2  is high-band repeater CC magnitude-dependent?  (aligned / P)')
    hi = A[(A.band == '20-80') & (A.is_cand)]
    lo = A[(A.band == '5-20') & (A.is_cand)]
    if not hi.empty:
        m = hi.merge(lo[['i', 'j', 'algn_p']], on=['i', 'j'],
                     suffixes=('_hi', '_lo'))
        m = m.sort_values('mag_max')
        print(f'  {"Mmax":>6}{"Mmin":>6}{"CC 5-20":>10}{"CC 20-80":>10}'
              f'{"retained":>10}')
        for _, x in m.iterrows():
            ret = x.algn_p_hi / x.algn_p_lo if x.algn_p_lo else np.nan
            print(f'  {x.mag_max:6.2f}{x.mag_min:6.2f}{x.algn_p_lo:10.3f}'
                  f'{x.algn_p_hi:10.3f}{ret:10.2f}')
        if len(m) > 3:
            rr = np.corrcoef(m.mag_max, m.algn_p_hi)[0, 1]
            print(f'\n  r(Mmax, CC in 20-80 Hz) = {rr:+.3f}')
            print('  P2 predicts this to be clearly POSITIVE. If it is ~0, the'
                  '\n  high band is noise and any d-prime gain is an artifact'
                  ' of\n  counting more noise samples.')

    print('\nP3  does d-prime improve with bandwidth?')
    print(f'  {"band":>10}' + ''.join(f'{k:>10}' for k in
                                      ['flat_o', 'algn_o', 'flat_p', 'algn_p']))
    for band in BANDS:
        lab = f'{band[0]:.0f}-{band[1]:.0f}'
        s = A[A.band == lab]
        if s.empty:
            continue
        row = f'  {lab:>10}'
        for k in ['flat_o', 'algn_o', 'flat_p', 'algn_p']:
            rm, _ = robust(s[s.is_cand][k]); nm, nmad = robust(s[~s.is_cand][k])
            dp = (rm - nm) / nmad if nmad and np.isfinite(nmad) else np.nan
            row += f'{dp:10.2f}'
        print(row)
    print('\n  reference, 100 Hz / 5-20 Hz (METHODS_STATUS 17):'
          '\n     flat_o 2.45   algn_o 1.67   flat_p 2.61   algn_p 1.74')

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    for k, mk in [('flat_o', 'o-'), ('algn_o', 's-'),
                  ('flat_p', 'o--'), ('algn_p', 's--')]:
        xs, dps, nm_ = [], [], []
        for band in BANDS:
            s = A[A.band == f'{band[0]:.0f}-{band[1]:.0f}']
            if s.empty:
                continue
            rm, _ = robust(s[s.is_cand][k]); n0, nmad = robust(s[~s.is_cand][k])
            xs.append(band[1] - band[0])
            dps.append((rm - n0) / nmad if nmad else np.nan)
            nm_.append(nmad)
        ax[0].plot(xs, dps, mk, label=k)
        ax[1].plot(xs, nm_, mk, label=k)
    ax[0].set(xlabel='bandwidth (Hz)', ylabel="d'", title="A  discrimination")
    ax[0].axhline(2.45, color='0.6', ls=':', label='100 Hz flat_o')
    ax[1].set(xlabel='bandwidth (Hz)', ylabel='null MAD', title='B  null width')
    for a in ax[:2]:
        a.legend(fontsize=7); a.grid(alpha=.3)
    if not hi.empty:
        ax[2].scatter(hi.mag_max, hi.algn_p, c='C3')
        ax[2].set(xlabel='max magnitude in pair', ylabel='CC, 20-80 Hz',
                  title='C  P2: is the discarded band real signal?')
        ax[2].grid(alpha=.3)
    fig.tight_layout()
    p = os.path.join(HERE, 'hf_moveout_test.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
