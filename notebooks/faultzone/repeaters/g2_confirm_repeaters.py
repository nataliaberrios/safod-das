"""
GATE G2: are the catalog clusters actually REPEATING earthquakes?

G1 passed -- the fiber sees coda (785/900 channels above SNR 3 for M1.6-1.9 at
~3 km). That establishes the measurement is possible. It says nothing about
whether the events repeat.

The candidate list in doublet_candidates.csv came from catalog COLOCATION, which
is an upper bound and has two known defects I found by inspecting it:
  * single-link clustering chains, so seq 0 spans 1 km of depth despite a 0.5 km
    threshold -- these are over-merged;
  * some "repeats" are seconds apart (2025-03-29 06:43:47 / 06:43:58), which give
    no temporal baseline for CWI at all.

Only waveform similarity settles it. True repeaters rupture the same patch the
same way, so their waveforms are near-identical; merely nearby events are not.

PASS: >= 3 sequences with median cross-channel CC > 0.9 for at least one pair
      separated by > 30 days.
FAIL: the catalog clusters are not repeaters. The plan says stop.

Caches every extracted snippet to disk -- week 2 must not re-extract, since
extraction is the expensive step (~5-9 min/event against a 392k-row manifest).
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
from g1_coda_snr import load_manifest, normalize, DX          # noqa: E402

for p in ['/home/groups/edunham/nberrios/safod_das/DAS-utilities/python',
          '/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python']:
    if os.path.exists(p) and p not in sys.path:
        sys.path.insert(0, p)
import DASutils                                                # noqa: E402

CACHE = os.path.join(HERE, 'cache')
os.makedirs(CACHE, exist_ok=True)

# Shorter window than G1: repeater confirmation needs the direct arrivals and a
# little coda, not 60 s. Cuts extraction cost roughly threefold.
PRE_S, POST_S = 5.0, 15.0
FMIN, FMAX = 5.0, 40.0
CC_BAND = (5.0, 20.0)
CC_WIN = (0.0, 5.0)          # lapse window for similarity, relative to origin
CH_LO, CH_HI = 100, 800      # skip the noisy top and the sub-800 m failure zone

# Candidate source is now selectable. The DDRT (double-difference relocated)
# list supersedes the standard-catalog one: DD relative locations are accurate to
# tens of metres, which is the scale repeaters live at, whereas the standard
# catalog's 0.2-0.5 km error guaranteed the first attempt clustered different
# patches together.
CAND_FILE = os.environ.get('CAND_FILE', 'ddrt_candidates.csv')
SEQS = [int(x) for x in os.environ.get('SEQS', '0,1').split(',')]
MAX_EV = 8                   # cap per sequence to bound runtime
PASS_CC, PASS_DAYS, PASS_NSEQ = 0.9, 30.0, 3


def extract(db, origin, tag):
    """Snippet with on-disk caching. Returns (X, fs) or (None, None)."""
    f = os.path.join(CACHE, f'{tag}.npz')
    if os.path.exists(f):
        d = np.load(f)
        return d['X'], float(d['fs'])
    a = origin - pd.Timedelta(seconds=PRE_S)
    b = origin + pd.Timedelta(seconds=POST_S)
    sel = db[(db['t0'] < b) & (db['t1'] > a)].copy()
    # Check emptiness BEFORE masking: a boolean index built from .map() on an
    # empty frame loses the columns, so the later sort_values('t0') raises
    # KeyError rather than returning an empty result.
    if sel.empty:
        return None, None
    sel = sel[sel['fn'].map(os.path.exists)]
    if sel.empty:
        return None, None
    sel = sel.sort_values('t0')
    acc = hits = None
    fs = dt = npts = None
    for _, r in sel.iterrows():
        try:
            D, info = DASutils.readFile_HDF(
                [r['fn']], FMIN, FMAX, verbose=0, preproc=True, diff=True,
                taper=False, desampling=True, nChbuffer=900, system='OptaSense')
        except Exception:
            return None, None
        if fs is None:
            fs = info['fs']; dt = 1.0 / fs
            npts = int(round((PRE_S + POST_S) * fs))
            acc = np.zeros((D.shape[0], npts)); hits = np.zeros(npts, int)
        o0, o1 = max(r['t0'], a), min(r['t1'], b)
        if o1 <= o0:
            continue
        s0 = max(int(round((o0 - r['t0']).total_seconds() / dt)), 0)
        d0 = max(int(round((o0 - a).total_seconds() / dt)), 0)
        n = min(int(round((o1 - o0).total_seconds() / dt)),
                D.shape[1] - s0, npts - d0)
        if n > 0:
            acc[:, d0:d0 + n] += D[:, s0:s0 + n]
            hits[d0:d0 + n] += 1
        del D
    if hits is None or np.any(hits == 0):
        return None, None
    X = acc / hits[None, :]
    np.savez_compressed(f, X=X.astype(np.float32), fs=fs)
    return X, fs


MAX_LAG_S = 2.0     # covers catalog origin-time error for M~1 events


def cc_pair(A, B, fs):
    """Median CC across channels, MAXIMISED over lag.

    Zero-lag correlation is wrong here and was the first version's bug. Catalog
    origin times for M~1 events carry ~0.1-0.5 s uncertainty, which at 100 Hz is
    10-50 samples; two *identical* waveforms offset by that much give near-zero
    zero-lag correlation. Repeater identification is always done by taking the
    maximum over lags.

    The moveout is flat across the array (near-vertical incidence, confirmed in
    G1), so one lag applies to every channel: estimate it once from the
    channel-stacked traces, then evaluate per-channel CC at that lag. Cheaper
    than a per-channel search and less prone to each channel picking its own
    spurious peak.
    """
    i0 = int((PRE_S + CC_WIN[0]) * fs)
    i1 = int((PRE_S + CC_WIN[1]) * fs)
    pad = int(MAX_LAG_S * fs)
    a = A[CH_LO:CH_HI, i0:i1].astype(np.float64)
    # widen B's window so the shifted comparison stays inside the record
    b_wide = B[CH_LO:CH_HI, max(i0 - pad, 0):i1 + pad].astype(np.float64)
    a -= a.mean(axis=1, keepdims=True)
    b_wide -= b_wide.mean(axis=1, keepdims=True)

    sa, sb = a.mean(axis=0), b_wide.mean(axis=0)
    cc_full = np.correlate(sb, sa, mode='valid')
    den = np.sqrt(np.sum(sa ** 2) * np.sum(sb ** 2))
    if den <= 0 or cc_full.size == 0:
        return np.nan, np.full(a.shape[0], np.nan), 0.0
    k = int(np.argmax(np.abs(cc_full)))
    lag = (k - min(pad, i0)) / fs

    b = b_wide[:, k:k + a.shape[1]]
    if b.shape[1] != a.shape[1]:
        return np.nan, np.full(a.shape[0], np.nan), lag
    b = b - b.mean(axis=1, keepdims=True)
    num = np.sum(a * b, axis=1)
    d = np.sqrt(np.sum(a * a, axis=1) * np.sum(b * b, axis=1))
    cc = num / np.where(d > 0, d, np.nan)
    return float(np.nanmedian(cc)), cc, lag


def main():
    cand = pd.read_csv(os.path.join(HERE, CAND_FILE))
    print(f'candidates: {CAND_FILE}, sequences {SEQS}', flush=True)
    cand['time'] = pd.to_datetime(cand['time'], utc=True, format='mixed')
    db = load_manifest()
    print(f'manifest: {len(db)} rows\n', flush=True)

    summary, passing = [], []
    for seq in SEQS:
        ev = cand[cand['seq'] == seq].sort_values('time').reset_index(drop=True)
        # prefer events spread in time over near-duplicates seconds apart
        keep, last = [], None
        for i, r in ev.iterrows():
            if last is None or (r['time'] - last).total_seconds() > 3600:
                keep.append(i); last = r['time']
        ev = ev.loc[keep].reset_index(drop=True).head(MAX_EV)
        print(f'--- sequence {seq}: {len(ev)} events after removing '
              f'near-duplicates ---', flush=True)

        recs = []
        for _, e in ev.iterrows():
            tag = f'seq{seq}_{e["time"]:%Y%m%dT%H%M%S}'
            X, fs = extract(db, e['time'], tag)
            if X is None:
                print(f'   {e["time"]:%Y-%m-%d %H:%M}  M{e["mag"]:.2f}  skip',
                      flush=True)
                continue
            Xb = DASutils.bandpass2D_c(X - np.median(X, axis=0, keepdims=True),
                                       CC_BAND[0], CC_BAND[1], 1.0 / fs,
                                       zerophase=True)
            recs.append(dict(t=e['time'], m=e['mag'], X=Xb, fs=fs))
            print(f'   {e["time"]:%Y-%m-%d %H:%M}  M{e["mag"]:.2f}  ok',
                  flush=True)
        if len(recs) < 2:
            print('   too few events\n', flush=True)
            continue

        n = len(recs)
        M = np.full((n, n), np.nan)
        best = None
        for i in range(n):
            M[i, i] = 1.0
            for j in range(i + 1, n):
                c, _, lg = cc_pair(recs[i]['X'], recs[j]['X'], recs[i]['fs'])
                M[i, j] = M[j, i] = c
                dd = abs((recs[j]['t'] - recs[i]['t']).total_seconds()) / 86400
                if c > PASS_CC and dd > PASS_DAYS and (best is None or c > best[0]):
                    best = (c, dd, recs[i]['t'], recs[j]['t'], lg)
        print(f'   CC matrix: max off-diagonal {np.nanmax(M[~np.eye(n,dtype=bool)]):.3f}',
              flush=True)
        if best:
            print(f'   PASS: CC {best[0]:.3f} over {best[1]:.0f} days '
                  f'({best[2]:%Y-%m-%d} / {best[3]:%Y-%m-%d}), lag {best[4]:+.3f} s\n',
                  flush=True)
            passing.append(seq)
        else:
            print('   no pair with CC > %.2f separated by > %.0f days\n'
                  % (PASS_CC, PASS_DAYS), flush=True)
        summary.append(dict(seq=seq, n=n, M=M,
                            times=[r['t'] for r in recs],
                            mags=[r['m'] for r in recs], best=best))

    verdict = 'PASS' if len(passing) >= PASS_NSEQ else 'FAIL'
    print(f'G2 RESULT: {verdict} -- {len(passing)}/{len(SEQS)} sequences '
          f'confirmed (need >= {PASS_NSEQ}): {passing}', flush=True)
    if verdict == 'FAIL':
        print('  -> The plan says stop: the catalog clusters are not repeaters.',
              flush=True)

    np.savez(os.path.join(HERE, 'g2_confirm_repeaters.npz'),
             seqs=np.array([s['seq'] for s in summary]),
             passing=np.array(passing), verdict=verdict,
             **{f'M_{s["seq"]}': s['M'] for s in summary})

    k = len(summary)
    if k:
        fig, ax = plt.subplots(1, k, figsize=(4 * k, 4.4))
        ax = np.atleast_1d(ax)
        for a, s in zip(ax, summary):
            im = a.imshow(s['M'], vmin=0, vmax=1, cmap='viridis')
            a.set_title(f'seq {s["seq"]}  (n={s["n"]})'
                        + ('\nCONFIRMED' if s['seq'] in passing else ''),
                        fontsize=9)
            a.set_xticks(range(s['n']))
            a.set_yticks(range(s['n']))
            a.set_xticklabels([f'{t:%m-%d}' for t in s['times']], rotation=90,
                              fontsize=6)
            a.set_yticklabels([f'{t:%m-%d}' for t in s['times']], fontsize=6)
            plt.colorbar(im, ax=a, fraction=0.046)
        fig.suptitle(f'G2: waveform similarity between candidate repeats '
                     f'({CC_BAND[0]:.0f}-{CC_BAND[1]:.0f} Hz, '
                     f'ch {CH_LO}-{CH_HI}) — {verdict}', fontsize=11)
        fig.tight_layout()
        out = os.path.join(HERE, 'g2_confirm_repeaters.png')
        fig.savefig(out, dpi=140)
        print(f'\nSaved {out}', flush=True)


if __name__ == '__main__':
    main()
