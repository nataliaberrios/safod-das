"""
Verify the geometrically-screened pairs by waveform, with the corrected pipeline.

This is Ellsworth's recipe executed in the right order, and for the first time with
both halves correct:

    "find events that are close in magnitude and location, which will then need to
     be verified as either repeaters or neighbors"

CLOSE IN LOCATION -- now actually true. The previous candidate set was gated on a
horizontal_separation_km column reading 57x too small (median 72 m against a true
4137 m), so it screened nothing. Rebuilt from DDRT lat/lon/depth: 37 pairs at
sep < 200 m with dM < 0.5, of which 13 have gap-free DAS coverage at both events.

VERIFIED BY WAVEFORM -- with the pipeline fixed. Every previous DAS correlation
here used a flat channel stack across a 0.31 s moveout (~27 dB lost) on data
decimated to 100 Hz and low-passed at 40 Hz (5x band lost). Corrected, the same
class of pair goes from CC 0.680 to 0.956, at parity with HRSN's 0.954.

So geometry and waveform are now two INDEPENDENT criteria, which is the point:
the screen knows nothing about waveforms and the correlation knows nothing about
locations. Agreement between them means something; agreement between a waveform
criterion and a waveform-derived family label does not, which is what invalidated
the same-patch discriminant (r = -0.826 with CC within the family group).

CONTROLS, fixed before running:
  * the null is pairs from the SAME event set with sep > 2 km -- matched in
    magnitude range, coverage and processing, differing only in separation
  * acausal (time-reversed) correlation gives the false-positive floor
  * CC is reported against separation, so a distance-dependent trend is visible
    rather than assumed
"""
import os, sys, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import moveout_test as MT
from dvv_core import bulk_align

CACHE = os.path.join(HERE, 'cache_hf')
BAND = (5.0, 40.0)
FS_COMMON = 500.0
WIN = (-0.5, 6.0)          # P through early coda, about the picked P


def beam(tag):
    """Moveout-corrected beam, resampled to a COMMON rate.

    Three cached events (ev_20240510T122958, ev_20240512T042912,
    ev_20240513T093029) are 5000 Hz, not 500 -- almost certainly event-triggered
    recordings. Correlating a 5000 Hz trace against a 500 Hz one using a single fs
    silently returns nonsense: pair 20240510/20250406 came back at CC 0.059 despite
    an HRSN CC of 0.917. Everything is brought to FS_COMMON first.
    """
    from scipy.signal import decimate as _dec
    MT.BAND = BAND; MT.CACHE = CACHE
    got, _ = MT.prep(tag)
    if got is None:
        return None
    A, z, fs = got
    if fs > FS_COMMON * 1.5:
        q = int(round(fs / FS_COMMON))
        if abs(fs / q - FS_COMMON) > 1:
            return None
        A = _dec(A, q, axis=-1, ftype='fir', zero_phase=True)
        fs = fs / q
    if abs(fs - FS_COMMON) > 1:
        return None
    p, tp, sem, _, _ = MT.slant_scan(A, z, fs)
    B = MT.align(A, z, fs, p)
    i0, i1 = int((tp + WIN[0]) * fs), int((tp + WIN[1]) * fs)
    if i0 < 0 or i1 > B.shape[1]:
        return None
    s = B[:, i0:i1].mean(axis=0); s -= s.mean()
    n = np.sqrt((s ** 2).sum())
    return (s / n, fs, sem) if n > 0 else None


def cc(a, b, fs, reverse=False):
    if reverse:
        b = b[::-1].copy()
    b2, _ = bulk_align(a, b, fs, max_lag_s=2.0)
    n = min(a.size, b2.size)
    return float(np.max(np.abs(np.correlate(b2[:n], a[:n], 'full'))) /
                 max(np.sqrt((a[:n]**2).sum() * (b2[:n]**2).sum()), 1e-30))


def main():
    P = pd.read_csv(os.path.join(HERE, os.environ.get(
        'PAIR_SET', 'ddrt_pairs_corrected.csv')))
    E = pd.read_csv(os.path.join(HERE, os.environ.get(
        'EVENT_SET', 'ddrt_corrected_events.csv')))
    tag = dict(zip(E.k, E.tag))
    have = {f[:-4] for f in os.listdir(CACHE)}
    B = {}
    for k, t in tag.items():
        if t in have:
            r = beam(t)
            if r: B[k] = r
    print(f'{len(B)} events beamed of {len(tag)}')
    rows = []
    for _, r in P.iterrows():
        i, j = int(r.i), int(r.j)
        if i not in B or j not in B: continue
        a, fs, sa = B[i]; b, _, sb = B[j]
        rows.append(dict(sep_m=r.sep_m, dM=r.dM, mmax=r.mmax, rad_m=r.rad_m,
                         dt_days=r.dt_days, cc=cc(a, b, fs),
                         acausal=cc(a, b, fs, reverse=True),
                         sem=min(sa, sb), tag=f'{tag[i][3:11]}/{tag[j][3:11]}'))
    D = pd.DataFrame(rows)
    if D.empty:
        print('no pairs'); return
    D.to_csv(os.path.join(HERE, 'ddrt_pair_cc.csv'), index=False)
    print(f'\n{"pair":>20}{"sep m":>8}{"dM":>6}{"Mmax":>6}{"dt d":>8}{"CC":>8}{"acausal":>9}')
    for _, x in D.sort_values('cc', ascending=False).iterrows():
        print(f'{x.tag:>20}{x.sep_m:8.0f}{x.dM:6.2f}{x.mmax:6.2f}'
              f'{x.dt_days:8.1f}{x.cc:8.3f}{x.acausal:9.3f}')
    ac = D.acausal.max()
    print(f'\nacausal false-positive floor: {ac:.3f}')
    print(f'pairs above it: {int((D.cc > ac).sum())} of {len(D)}')
    if len(D) > 5:
        from scipy.stats import pearsonr
        r_, p_ = pearsonr(D.sep_m, D.cc)
        print(f'r(separation, CC) = {r_:+.3f}, p = {p_:.3f}   '
              '(negative = closer pairs correlate better)')
        # dM MUST be controlled for. Different magnitudes mean different corner
        # frequencies and so genuinely different waveforms, which lowers CC with
        # no help from geometry. Without this the separation term is confounded
        # by exactly the kind of nuisance that invalidated the same-patch test.
        X1 = np.column_stack([np.ones(len(D)), D.dM])
        X2 = np.column_stack([np.ones(len(D)), D.dM, D.sep_m])
        b1, *_ = np.linalg.lstsq(X1, D.cc, rcond=None)
        b2, *_ = np.linalg.lstsq(X2, D.cc, rcond=None)
        ss1 = float(((D.cc - X1 @ b1) ** 2).sum())
        ss2 = float(((D.cc - X2 @ b2) ** 2).sum())
        F = ((ss1 - ss2) / 1) / (ss2 / (len(D) - 3))
        try:
            from scipy.stats import f as fdist
            pf = 1 - fdist.cdf(F, 1, len(D) - 3)
        except Exception:
            pf = np.nan
        print(f'  controlling for dM: separation coefficient {b2[2]:+.2e} /m, '
              f'F = {F:.2f}, p = {pf:.3f}')
        print(f'  r(dM, CC) = {pearsonr(D.dM, D.cc)[0]:+.3f}  '
              '(the nuisance being controlled)')


if __name__ == '__main__':
    main()
