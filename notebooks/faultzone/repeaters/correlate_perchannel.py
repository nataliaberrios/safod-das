"""
Repeater identification following the published DAS recipe, not my improvisation.

WHAT WAS WRONG BEFORE. I channel-stacked 700 DAS traces into one and then
correlated. Nobody does that. Every reference -- seismometer networks and DAS
alike -- correlates PER SENSOR and combines afterwards:

  Waldhauser & Schaff 2008   mean CC >= 0.9 at >= 5 common stations
  Schaff & Waldhauser 2005   CC >= 0.7 at >= 4 stations
  Li & Zhan 2018 (DAS)       per-channel CC, "stacks all individual cross
                             correlations into a network mean"
  Lellouch et al. 2021 (DAS) "N independent cross-correlations for the N
                             channels... correlograms are stacked, using weights
                             based on average channel-SNR"

The reason is geometric and decisive: two co-located events have IDENTICAL
moveout, so correlating channel i of A against channel i of B cancels the moveout
exactly -- no velocity model needed. Stacking correlograms is moveout-free.
Stacking traces is not: across a 919 m aperture the moveout is 0.15-0.46 s, while
coherent stacking needs it well under half a period (0.025 s at 20 Hz). Trace
stacking acts as a comb filter.

THRESHOLD. Do not import 0.9 -- that number is calibrated on NCSN vertical
short-period seismometers, 1 s P windows, 1.5-15 Hz, averaged over >=5 stations.
Lellouch et al. 2021, on downhole DAS, calibrate against their own null instead:
6 x MAD of the CC distribution (~0.20 for them), validated with a TIME-REVERSED
acausal template that can only produce random CC. Both are computed here.

MULTI-BAND. Igarashi 2020 selects passband by source size. For M0-M2 the Brune
corner is 23-228 Hz, so 5-20 Hz sits below corner for everything -- forgiving of
magnitude differences but blind to fine differences. Run 2-8, 4-16 and 10-40 Hz:
a true repeater survives all three; a merely-nearby pair decorrelates upward.
"""
import os
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt, correlate
from scipy.fft import next_fast_len
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache_all')

CH_LO, CH_HI = 100, 800      # Lellouch et al. exclude the degraded deep channels
BANDS = [(2.0, 8.0), (4.0, 16.0), (10.0, 40.0)]
PRE_S = 5.0
WIN = (-0.5, 6.0)            # P through S plus a little coda
MAX_LAG_S = 1.0              # Schaff & Waldhauser: +/-1 s covers pick error
MIN_CH = 200                 # DAS analogue of ">= 5 stations" (Li et al. 2021)


def load(tag, band):
    f = os.path.join(CACHE, f'{tag}.npz')
    if not os.path.exists(f):
        return None, None
    d = np.load(f)
    X, fs = d['X'][CH_LO:CH_HI].astype(np.float64), float(d['fs'])
    sos = butter(4, list(band), btype='band', fs=fs, output='sos')
    Xb = sosfiltfilt(sos, X, axis=-1)
    i0 = int((PRE_S + WIN[0]) * fs)
    i1 = int((PRE_S + WIN[1]) * fs)
    W = Xb[:, i0:i1]
    W = W - W.mean(axis=1, keepdims=True)
    # per-channel L2 normalisation (Lellouch: "trace-by-trace L2 normalization")
    nrm = np.sqrt(np.sum(W ** 2, axis=1, keepdims=True))
    good = (nrm[:, 0] > 0) & np.isfinite(nrm[:, 0])
    W[good] /= nrm[good]
    return W, (fs, good)


def stacked_correlogram(A, B, fs, good, reverse=False):
    """Per-channel CC, then stack the CORRELOGRAMS (not the traces).

    Returns peak of the stack, its lag, and the contributing channel count.
    reverse=True time-flips A, giving an acausal template that can only produce
    random correlation -- the empirical false-positive floor.
    """
    pad = int(MAX_LAG_S * fs)
    n = A.shape[1]
    acc = None
    cnt = 0
    for k in np.where(good)[0]:
        a = A[k][::-1] if reverse else A[k]
        c = correlate(B[k], a, mode='full', method='fft')
        mid = n - 1
        seg = c[max(mid - pad, 0):mid + pad + 1]
        if acc is None:
            acc = np.zeros_like(seg)
        if seg.size == acc.size:
            acc += seg
            cnt += 1
    if acc is None or cnt == 0:
        return np.nan, np.nan, 0
    acc /= cnt
    j = int(np.argmax(np.abs(acc)))
    return float(acc[j]), (j - pad) / fs, cnt


def precompute_fft(Ws, meta, nfft):
    """rfft of every event's channel block, with unusable channels zeroed.

    Stacking correlograms is linear, so the whole 700-channel stack for a pair is
    one frequency-domain product summed over channels -- sum_k conj(F_i[k]) F_j[k]
    -- followed by a single inverse transform. The explicit loop in
    stacked_correlogram recomputes an FFT of both events for every one of the
    21,115 pairs, i.e. 206 transforms of work done 21,115 times. Precomputing turns
    a ~1 hour band into about 2 minutes, which is the difference between this
    running and this not running.

    Zeroing bad channels rather than masking them makes the pairwise intersection
    automatic: a channel unusable in either event contributes exactly zero to the
    product, which is what the loop's `gi & gj` does.
    """
    F = []
    for W, m in zip(Ws, meta):
        if W is None:
            F.append(None)
            continue
        Wz = np.where(m[1][:, None], W, 0.0)
        F.append(np.fft.rfft(Wz, nfft, axis=1).astype(np.complex64))
    return F


def stacked_correlogram_fft(Fi, Fj, gi, gj, fs, nfft, pad):
    """Frequency-domain equivalent of stacked_correlogram. Same value, same lag.

    Accumulation is forced to complex128: summing 700 complex64 terms in float32
    would leave a relative error near 1e-5, and while that is far below the ~0.02
    MAD that sets the threshold, the acausal floor is compared against the same
    numbers and there is no reason to spend accuracy here.
    """
    cnt = int(np.count_nonzero(gi & gj))
    if cnt == 0:
        return np.nan, np.nan, 0
    S = np.einsum('kf,kf->f', Fi.conj(), Fj, dtype=np.complex128)
    r = np.fft.irfft(S, nfft)
    # r[l] is lag +l, negative lags wrap to the end; reorder to -pad..+pad
    seg = np.concatenate([r[nfft - pad:], r[:pad + 1]]) / cnt
    j = int(np.argmax(np.abs(seg)))
    return float(seg[j]), (j - pad) / fs, cnt


def main():
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    ev['time'] = pd.to_datetime(ev['time'], utc=True, format='mixed')

    # DUPLICATE-TAG FIX. Cache tags are second-resolution, so two catalog entries
    # inside one second share a cache file and correlate with themselves at
    # CC = 1.000 -- error #7 in METHODS_STATUS. There is exactly one such
    # collision here, 2024-12-30 04:39:33, and checking the catalog settles what it
    # is: nc75109596 (M0.98, 18 stations, gap 60, rms 0.07) and nc75109601
    # (M1.02, 16 stations, gap 84, rms 0.09) carry the SAME origin time to the
    # millisecond, 0.3 km apart. That is the NCSN double-listing one earthquake
    # with two solutions, not two earthquakes. So the collapse was physically
    # correct and no re-extraction is needed; what is needed is to drop the second
    # listing so it cannot pair with itself. Keep the better-constrained solution.
    dup = ev.tag.duplicated(keep=False)
    if dup.any():
        print(f'duplicate cache tags: {int(dup.sum())} rows')
        print(ev.loc[dup, ['idx', 'tag', 'time', 'mag', 'lat', 'lon']]
              .to_string(index=False))
        # Both listings point at the same cache file, so which row survives has no
        # effect on the waveforms -- only on the metadata carried alongside. Keep
        # the earlier row; for the record that is nc75109596, the better-constrained
        # of the two solutions (18 stations vs 16, gap 60 vs 84, rms 0.07 vs 0.09).
        ev = ev.drop_duplicates(subset='tag', keep='first').reset_index(drop=True)
        print(f'  -> same earthquake double-listed; kept one, '
              f'{len(ev)} events remain\n')
    N = len(ev)
    print(f'{N} unique events\n', flush=True)

    out = {}
    for band in BANDS:
        print(f'=== band {band[0]:.0f}-{band[1]:.0f} Hz ===', flush=True)
        Ws = meta = F = None      # release the previous band before allocating
        Ws, meta = [], []
        for _, e in ev.iterrows():
            W, m = load(e['tag'], band)
            Ws.append(W); meta.append(m)
        ok = [i for i in range(N) if Ws[i] is not None]
        fs = meta[ok[0]][0]
        pad = int(MAX_LAG_S * fs)
        npts = Ws[ok[0]].shape[1]
        nfft = next_fast_len(2 * npts - 1)
        F = precompute_fft(Ws, meta, nfft)
        print(f'  {len(ok)} events, {npts} samples, nfft {nfft}', flush=True)

        # Verify the fast path against the explicit loop before trusting it on
        # 21,115 pairs. Cheap, and it is the only thing standing between a
        # transcription error and a plausible-looking wrong threshold.
        worst = 0.0
        for i, j in [(ok[0], ok[1]), (ok[0], ok[len(ok) // 2]),
                     (ok[1], ok[-1]), (ok[len(ok) // 3], ok[-2])]:
            g = meta[i][1] & meta[j][1]
            if g.sum() < MIN_CH:
                continue
            c0, l0, n0 = stacked_correlogram(Ws[i], Ws[j], fs, g)
            c1, l1, n1 = stacked_correlogram_fft(F[i], F[j], meta[i][1],
                                                 meta[j][1], fs, nfft, pad)
            worst = max(worst, abs(c0 - c1))
            assert n0 == n1 and abs(l0 - l1) < 1e-9, (n0, n1, l0, l1)
        print(f'  fast-path check vs explicit loop: max |dCC| = {worst:.2e}',
              flush=True)
        assert worst < 1e-5, 'FFT correlogram disagrees with the loop'

        C = np.full((N, N), np.nan)
        L = np.full((N, N), np.nan)
        NC = np.zeros((N, N), int)
        AC = []
        for a_i, i in enumerate(ok):
            gi = meta[i][1]
            for j in ok[a_i + 1:]:
                gj = meta[j][1]
                if np.count_nonzero(gi & gj) < MIN_CH:
                    continue
                c, lg, n_ = stacked_correlogram_fft(F[i], F[j], gi, gj, fs,
                                                    nfft, pad)
                C[i, j] = C[j, i] = c
                L[i, j] = lg; L[j, i] = -lg
                NC[i, j] = NC[j, i] = n_
            if a_i % 40 == 0:
                print(f'  row {a_i}/{len(ok)}', flush=True)

        # acausal null: time-reversed templates, a random subset is enough
        rng = np.random.default_rng(0)
        for _ in range(300):
            i, j = rng.choice(ok, 2, replace=False)
            g = meta[i][1] & meta[j][1]
            if g.sum() < MIN_CH:
                continue
            c, _, _ = stacked_correlogram(Ws[i], Ws[j], fs, g, reverse=True)
            if np.isfinite(c):
                AC.append(abs(c))

        iu = np.triu_indices(N, 1)
        cc = C[iu][np.isfinite(C[iu])]
        mad = 1.4826 * np.median(np.abs(cc - np.median(cc)))
        thr = 6 * mad
        acmax = max(AC) if AC else np.nan
        print(f'  pairs {cc.size}, median {np.median(cc):.4f}, MAD {mad:.4f}')
        print(f'  6xMAD threshold      : {thr:.4f}')
        print(f'  acausal max (n={len(AC)}): {acmax:.4f}   <- false-positive floor')
        print(f'  pairs above threshold: {int(np.sum(cc > thr))}')
        print(f'  max CC               : {cc.max():.4f}\n', flush=True)
        out[band] = dict(C=C, L=L, NC=NC, thr=thr, acmax=acmax, mad=mad)

    # a repeater must survive every band
    surv = np.ones((N, N), bool)
    for band in BANDS:
        surv &= (out[band]['C'] > max(out[band]['thr'], out[band]['acmax']))
    iu = np.triu_indices(N, 1)
    hits = [(i, j) for i, j in zip(*iu) if surv[i, j]]
    dt = lambda i, j: abs((ev.time[j] - ev.time[i]).total_seconds()) / 86400
    hits.sort(key=lambda p: -min(out[b]['C'][p] for b in BANDS))

    print(f'pairs above threshold in ALL {len(BANDS)} bands: {len(hits)}')
    print(f'\n{"2-8":>7}{"4-16":>8}{"10-40":>8}{"nch":>6}{"dt_d":>8}{"dM":>6}   events')
    for i, j in hits[:30]:
        cs = [out[b]['C'][i, j] for b in BANDS]
        print(f'{cs[0]:7.3f}{cs[1]:8.3f}{cs[2]:8.3f}'
              f'{out[BANDS[0]]["NC"][i,j]:6d}{dt(i,j):8.1f}'
              f'{abs(ev.mag[i]-ev.mag[j]):6.2f}   '
              f'{ev.time[i]:%Y-%m-%d} M{ev.mag[i]:.2f} / '
              f'{ev.time[j]:%Y-%m-%d} M{ev.mag[j]:.2f}')

    # Hand-off to the HRSN confirmation, which currently covers only the top 10
    # pairs from the old trace-stacked search. Ranked by the WEAKEST band, so a
    # pair that survives every band outranks one that is spectacular in 2-8 Hz and
    # gone by 10-40 -- the multi-band criterion, applied to the ordering as well as
    # to the cut. dt > 20 d keeps out same-aftershock-sequence pairs, which are
    # similar for reasons that have nothing to do with a repeating source.
    rank = []
    for i, j in zip(*iu):
        cs = [out[b]['C'][i, j] for b in BANDS]
        if not all(np.isfinite(c) for c in cs) or dt(i, j) <= 20:
            continue
        rank.append((min(cs), i, j, cs))
    rank.sort(key=lambda x: -x[0])
    cdf = pd.DataFrame([dict(i=int(i), j=int(j), min_cc=float(mn),
                             cc_2_8=cs[0], cc_4_16=cs[1], cc_10_40=cs[2],
                             dt_days=dt(i, j), tag_i=ev.tag[i], tag_j=ev.tag[j],
                             t_i=ev.time[i], t_j=ev.time[j],
                             m_i=ev.mag[i], m_j=ev.mag[j])
                        for mn, i, j, cs in rank[:40]])
    cdf.to_csv(os.path.join(HERE, 'perchannel_candidates.csv'), index=False)
    print(f'\nwrote perchannel_candidates.csv: top {len(cdf)} pairs for HRSN '
          f'confirmation (currently confirmed: 10)')

    np.savez(os.path.join(HERE, 'correlate_perchannel.npz'),
             tags=ev.tag.values, times=ev.time.astype(str).values,
             mag=ev.mag.values, depth=ev.depth.values,
             **{f'C_{b[0]:.0f}_{b[1]:.0f}': out[b]['C'] for b in BANDS},
             **{f'L_{b[0]:.0f}_{b[1]:.0f}': out[b]['L'] for b in BANDS},
             thr=np.array([out[b]['thr'] for b in BANDS]),
             acmax=np.array([out[b]['acmax'] for b in BANDS]))

    fig, ax = plt.subplots(1, len(BANDS), figsize=(5.2 * len(BANDS), 4.6))
    for a, b in zip(np.atleast_1d(ax), BANDS):
        c = out[b]['C'][iu]; c = c[np.isfinite(c)]
        a.hist(c, bins=150, color='0.6')
        a.axvline(out[b]['thr'], color='C3', lw=2,
                  label=f'6xMAD {out[b]["thr"]:.3f}')
        a.axvline(out[b]['acmax'], color='C0', lw=1.6, ls='--',
                  label=f'acausal max {out[b]["acmax"]:.3f}')
        a.set(yscale='log', xlabel='stacked-correlogram CC', ylabel='pairs',
              title=f'{b[0]:.0f}-{b[1]:.0f} Hz')
        a.legend(fontsize=7)
    fig.suptitle('Per-channel CC, stacked correlograms, null-calibrated threshold',
                 fontsize=12)
    fig.tight_layout()
    p = os.path.join(HERE, 'correlate_perchannel.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
