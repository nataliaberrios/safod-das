"""
Stage 3: correlate every covered event against every other, and let similarity
define the sequences.

This is the test the previous three attempts were a poor substitute for. Those
used catalog LOCATION as the prior for which events might repeat; three of four
co-located, magnitude-matched pairs then came back at CC ~ 0, because location at
+/-150-250 m does not predict whether two events rupture the same patch.

The decisive advantage of doing all pairs: **the CC distribution is its own
null**. Of 23,005 pairs, essentially all are unrelated events, so the bulk of the
histogram measures what "not a repeater" looks like on THIS instrument, in THIS
band, with THIS processing. Repeaters are the upper-tail outliers. No threshold
needs to be imported from studies that used single seismometers at different
frequencies -- which was the weakness of the 0.9 criterion used earlier.

Method:
  * channel-stack each event over channels 100-800 (skips the noisy top and the
    degraded zone below 800 m that G1 identified), giving one high-SNR trace per
    event
  * correlate every pair with a lag search, since catalog origin times carry
    0.1-0.5 s error and zero-lag correlation is meaningless at that offset
  * characterise the null, set the threshold from it, then cluster
  * only THEN look at where the clustered events sit in space -- location as a
    check on the answer, not a filter on the input
"""
import os
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt, correlate
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache_all')

CH_LO, CH_HI = 100, 800
BAND = (5.0, 20.0)
PRE_S = 5.0
WIN = (-1.0, 10.0)       # relative to origin: P, S and early coda
MAX_LAG_S = 2.0


def main():
    ev = pd.read_csv(os.path.join(HERE, 'all_events.csv'))
    ev['time'] = pd.to_datetime(ev['time'], utc=True, format='mixed')

    traces, keep = [], []
    fs = None
    for _, e in ev.iterrows():
        f = os.path.join(CACHE, f'{e["tag"]}.npz')
        if not os.path.exists(f):
            continue
        d = np.load(f)
        X, fs = d['X'], float(d['fs'])
        sos = butter(4, list(BAND), btype='band', fs=fs, output='sos')
        Xb = sosfiltfilt(sos, X[CH_LO:CH_HI].astype(np.float64), axis=-1)
        # DO NOT subtract the across-channel median here. Median subtraction is the
        # standard DAS trick for common-mode NOISE, but following it with
        # mean(axis=0) computes (mean - median) ~ 0, so an earlier run correlated
        # residual noise and produced a spurious null.
        #
        # THE FLAT STACK BELOW IS DELIBERATE, BUT NOT FOR THE REASON ORIGINALLY
        # WRITTEN HERE. This comment used to assert the arrival is flat across the
        # 900 channels because incidence is near-vertical. That is backwards: for a
        # VERTICAL fiber, near-vertical incidence gives the steepest moveout the
        # geometry allows (~0.31 s end to end). See METHODS_STATUS 2.4.
        #
        # moveout_test.py measured both ways on the 10 HRSN-confirmed pairs
        # (METHODS_STATUS 17). Aligning before stacking raises repeater CC from
        # 0.680 to 0.956 -- HRSN parity -- so the moveout is unambiguously real.
        # But it raises the random-pair null further still, 0.124 -> 0.376, because
        # an aligned stack is a clean P wavelet and all P wavelets at similar
        # distance resemble each other. Detection margin therefore DROPS:
        # d' 2.18 -> 1.74, and repMIN - nullMAX 0.299 -> 0.177.
        #
        # This script's job is DISCRIMINATION, so the flat stack is kept. Use
        # aligned stacking for anything needing waveform FIDELITY -- dv/v, spectra,
        # differential timing -- where the null is irrelevant. Do not "fix" this
        # line without re-reading section 17.
        i0 = int((PRE_S + WIN[0]) * fs)
        i1 = int((PRE_S + WIN[1]) * fs)
        s = Xb[:, i0:i1].mean(axis=0)          # channel stack -> one trace
        s -= s.mean()
        n = np.sqrt(np.sum(s ** 2))
        if not np.isfinite(n) or n <= 0:
            continue
        traces.append(s / n)
        keep.append(e)
        del X, Xb
    ev = pd.DataFrame(keep).reset_index(drop=True)
    T = np.array(traces)
    N = len(ev)
    print(f'{N} events loaded, fs={fs}, trace {T.shape[1]} samples '
          f'({T.shape[1]/fs:.1f} s)')
    print(f'pairs: {N*(N-1)//2}\n', flush=True)

    pad = int(MAX_LAG_S * fs)
    C = np.eye(N)
    L = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            cc = correlate(T[j], T[i], mode='full', method='fft')
            mid = len(T[i]) - 1
            lo, hi = mid - pad, mid + pad + 1
            seg = cc[max(lo, 0):min(hi, len(cc))]
            k = int(np.argmax(np.abs(seg)))
            C[i, j] = C[j, i] = seg[k]
            L[i, j] = (k + max(lo, 0) - mid) / fs
            L[j, i] = -L[i, j]
        if i % 40 == 0:
            print(f'  row {i}/{N}', flush=True)

    iu = np.triu_indices(N, 1)
    cc = C[iu]
    print(f'\nCC distribution over {cc.size} pairs (the empirical null):')
    for q in [50, 90, 99, 99.9]:
        print(f'  {q:5.1f}th percentile: {np.percentile(cc, q):.3f}')
    print(f'  max: {cc.max():.3f}')
    mu, sd = np.median(cc), 1.4826 * np.median(np.abs(cc - np.median(cc)))
    print(f'  robust centre {mu:.3f}, robust sigma {sd:.3f}')

    # The null is BROAD (robust sigma ~0.25), so a median+8sigma rule lands above
    # 1.0 and admits nothing -- correlation is bounded, so sigma rules break here.
    # Count the tail explicitly instead and look for a separated population.
    print('\nupper tail:')
    for t_ in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]:
        n_ = int(np.sum(cc > t_))
        print(f'  CC > {t_:.2f}: {n_:5d} pairs  ({100*n_/cc.size:7.4f}%)   '
              f'expected under null if tail were smooth: '
              f'{cc.size*(1-np.searchsorted(np.sort(cc), t_)/cc.size):.0f}')

    # Duplicate origin times would produce CC = 1 spuriously; check before trusting
    dup = ev.time.duplicated(keep=False)
    if dup.any():
        print(f'\nWARNING: {int(dup.sum())} events share an origin time -- '
              f'these self-pairs inflate the tail')
        print(ev.loc[dup, ['tag', 'time', 'mag']].to_string(index=False))

    print('\ntop 15 pairs by CC:')
    order = np.argsort(cc)[::-1][:15]
    for k in order:
        i, j = iu[0][k], iu[1][k]
        dt_ = abs((ev.time[j] - ev.time[i]).total_seconds()) / 86400
        P = np.array([[ev.t[i], ev.r[i], ev.depth[i]],
                      [ev.t[j], ev.r[j], ev.depth[j]]]) * 1000
        sep = float(np.linalg.norm(P[0] - P[1]))
        print(f'  CC {cc[k]:6.3f}  {ev.time[i]:%Y-%m-%d %H:%M} M{ev.mag[i]:.2f} / '
              f'{ev.time[j]:%Y-%m-%d %H:%M} M{ev.mag[j]:.2f}  '
              f'dt={dt_:7.1f}d  sep={sep:6.0f}m  lag={L[i,j]:+.2f}s')

    thr = 0.7      # set after inspecting the tail; see printout above
    n_out = int(np.sum(cc > thr))
    print(f'\nusing threshold {thr:.2f}: {n_out} pairs '
          f'({100*n_out/cc.size:.4f}% of all)')

    # cluster on similarity
    D = 1.0 - np.clip(C, -1, 1)
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2
    Z = linkage(squareform(D, checks=False), method='average')
    lab = fcluster(Z, t=1.0 - thr, criterion='distance')
    ev['cl'] = lab
    fam = [(c, g) for c, g in ev.groupby('cl') if len(g) >= 2]
    fam.sort(key=lambda x: -len(x[1]))
    print(f'\nsimilarity families (>=2 events): {len(fam)}')
    if fam:
        print(f'{"cl":>5}{"n":>4}{"maxCC":>8}{"span_d":>9}{"sep_m":>8}'
              f'{"dmag":>7}   dates')
    for c, g in fam[:15]:
        ii = g.index.values
        sub = C[np.ix_(ii, ii)]
        mx = sub[np.triu_indices(len(ii), 1)].max()
        span = (g.time.max() - g.time.min()).total_seconds() / 86400
        klat = 111.19e3; klon = 111.19e3 * np.cos(np.radians(35.982))
        P = np.column_stack([g.t * 1000, g.r * 1000, g.depth * 1000])
        sep = float(np.max([np.linalg.norm(P[a] - P[b])
                            for a in range(len(P)) for b in range(a + 1, len(P))]))
        ds = ', '.join(f'{t:%Y-%m-%d}' for t in sorted(g.time))
        print(f'{c:>5}{len(g):>4}{mx:8.3f}{span:9.0f}{sep:8.0f}'
              f'{g.mag.max()-g.mag.min():7.2f}   {ds[:44]}')

    np.savez(os.path.join(HERE, 'correlate_all.npz'),
             C=C, L=L, tags=ev.tag.values, times=ev.time.astype(str).values,
             mag=ev.mag.values, depth=ev.depth.values, t=ev.t.values,
             r=ev.r.values, cl=ev.cl.values, thr=thr, band=np.array(BAND))
    ev.to_csv(os.path.join(HERE, 'correlate_all_events.csv'), index=False)

    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    ax[0].hist(cc, bins=200, color='0.6')
    ax[0].axvline(thr, color='C3', lw=2, label=f'threshold {thr:.2f}')
    ax[0].set(xlabel='cross-correlation', ylabel='pairs', yscale='log',
              title=f'A  Null from {cc.size} pairs')
    ax[0].legend(fontsize=8)

    im = ax[1].imshow(C, vmin=-0.2, vmax=1, cmap='magma')
    ax[1].set(title='B  CC matrix', xlabel='event', ylabel='event')
    ax[1].grid(False)
    plt.colorbar(im, ax=ax[1], fraction=0.046)

    ax[2].scatter(ev.t, -ev.depth, s=6, c='0.8')
    for c, g in fam:
        ax[2].scatter(g.t, -g.depth, s=45, facecolors='none', edgecolors='C3',
                      lw=1.3)
    ax[2].set(xlabel='distance along strike (km)', ylabel='depth (km)',
              title=f'C  {len(fam)} similarity families')
    ax[2].grid(alpha=0.3)
    fig.suptitle('Similarity-first repeater search: all covered events, all pairs',
                 fontsize=12)
    fig.tight_layout()
    p = os.path.join(HERE, 'correlate_all.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
