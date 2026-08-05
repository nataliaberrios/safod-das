"""
GATE G1: is there usable CODA SNR for M0-M1 earthquakes on this fiber?

This is the assumption the whole repeater project rests on and nobody has tested
it. Lellouch et al. 2019 showed earthquake ARRIVALS on this cable. Coda is a much
harder bar -- it is the scattered tail, tens of dB below the direct arrival, and
coda-wave interferometry needs it well above noise for many seconds of lapse time.
At 16.34 m gauge length, for M~1 events 3-5 km away, nobody has demonstrated it.

PASS: coda SNR > 3 over at least 200 channels, for a 2-10 s lapse window.
FAIL: stop the repeater project. No amount of processing recovers coda that is
      not there, and the plan says so explicitly.

Targets: sequence 4 from doublet_candidates.csv -- the largest magnitudes in the
candidate list (M1.58-1.86) at ~3 km from the borehole. If the coda is not visible
for these, it will not be visible for anything.

Extraction reuses the approach in event_cc_20s_exact_snippet.py: manifest load,
path rewriting, and a time-aligned snippet assembled across file boundaries with a
coverage check. fmax is raised from that script's 24 Hz to 40 Hz so the desampled
rate is ~100 Hz rather than 60 -- CWI stretching wants sample rate to spare.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

for p in ['/home/groups/edunham/nberrios/safod_das/DAS-utilities/python',
          '/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python']:
    if os.path.exists(p) and p not in sys.path:
        sys.path.insert(0, p)
import DASutils

HERE = os.path.dirname(os.path.abspath(__file__))
# BOTH manifests for this cable. Using only SAFODAS1 meant events could be
# SELECTED from the merged coverage list but not EXTRACTED -- the second member
# of both candidate pairs was silently skipped. Selection and extraction must see
# the same file universe.
CSVS = [
    ('/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/'
     'SAFOD_2024_2025.csv'),
    ('/oak/stanford/groups/ettore88/data/SAFOD/SAFOD-Harvest-2026-01-28/'
     'SAFOD_vertical_2026_01_28.csv'),
]
OLD = '/oak/stanford/groups/ettore88/data/SAFODAS1-harddrive-transfer'
NEW = '/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer'

DX = 1.0209523439407349
FMIN, FMAX = 5.0, 40.0        # read band; coda measured in 5-20 Hz
CODA_BAND = (5.0, 20.0)
PRE_S, POST_S = 20.0, 40.0    # snippet extent around origin time
NOISE_WIN = (-18.0, -3.0)     # pre-event noise, relative to origin
CODA_WIN = (2.0, 10.0)        # lapse window, relative to origin

PASS_SNR, PASS_NCH = 3.0, 200


def normalize(path):
    path = str(path)
    if os.path.exists(path):
        return path
    if path.startswith(OLD):
        alt = path.replace(OLD, NEW, 1)
        if os.path.exists(alt):
            return alt
    return path


def load_manifest():
    """Parse the file manifest, with an on-disk cache.

    Two costs were killing job walltime. `sep=r'\\s+'` is a regex separator,
    which forces pandas onto its slow pure-Python parser -- on a 108 MB, 392k-row
    file that alone ran past a 25 minute limit. And every job re-parsed from
    scratch. Use the C parser via delim_whitespace, then pickle the result so
    subsequent jobs load in seconds.
    """
    cache = os.path.join(HERE, 'manifest_cache.pkl')
    if os.path.exists(cache):
        return pd.read_pickle(cache)
    parts = []
    for c in CSVS:
        if os.path.exists(c):
            d = pd.read_csv(c, delim_whitespace=True)
            print(f'  {os.path.basename(c)}: {len(d)} rows', flush=True)
            parts.append(d)
    db = pd.concat(parts, ignore_index=True).drop_duplicates()
    db = db[db['nSamples'] > 0].reset_index(drop=True)
    db['t0'] = pd.to_datetime(db['startTime'], errors='coerce', utc=True)
    db['t1'] = pd.to_datetime(db['endTime'], errors='coerce', utc=True)
    db = db.dropna(subset=['t0', 't1']).reset_index(drop=True)
    db = db[(db['t1'] - db['t0']).dt.total_seconds() > 1.0].reset_index(drop=True)
    db['fn'] = db['file'].map(normalize)
    db.to_pickle(cache)
    return db


def snippet(db, origin):
    """Time-aligned window around `origin`, assembled across files.

    Same construction as build_event_snippet in event_cc_20s_exact_snippet.py:
    average overlapping samples, and refuse to return a partially covered window.
    """
    a = origin - pd.Timedelta(seconds=PRE_S)
    b = origin + pd.Timedelta(seconds=POST_S)
    sel = db[(db['t0'] < b) & (db['t1'] > a)].copy()
    # Emptiness must be checked before masking; a boolean index from .map() on
    # an empty frame drops the columns and sort_values('t0') then raises.
    if sel.empty:
        return None, None, 'no files'
    sel = sel[sel['fn'].map(os.path.exists)]
    if sel.empty:
        return None, None, 'no existing files'
    sel = sel.sort_values('t0')

    acc = hits = None
    fs = dt = npts = None
    for _, r in sel.iterrows():
        try:
            D, info = DASutils.readFile_HDF(
                [r['fn']], FMIN, FMAX, verbose=0, preproc=True, diff=True,
                taper=False, desampling=True, nChbuffer=900, system='OptaSense')
        except Exception as e:
            return None, None, f'read failed: {e}'
        if fs is None:
            fs = info['fs']; dt = 1.0 / fs
            npts = int(round((PRE_S + POST_S) * fs))
            acc = np.zeros((D.shape[0], npts))
            hits = np.zeros(npts, int)
        o0, o1 = max(r['t0'], a), min(r['t1'], b)
        if o1 <= o0:
            continue
        s0 = int(round((o0 - r['t0']).total_seconds() / dt))
        d0 = int(round((o0 - a).total_seconds() / dt))
        n = int(round((o1 - o0).total_seconds() / dt))
        s0, d0 = max(s0, 0), max(d0, 0)
        n = min(n, D.shape[1] - s0, npts - d0)
        if n <= 0:
            continue
        acc[:, d0:d0 + n] += D[:, s0:s0 + n]
        hits[d0:d0 + n] += 1
        del D
    if hits is None or np.any(hits == 0):
        miss = int(np.sum(hits == 0)) if hits is not None else -1
        return None, None, f'incomplete coverage ({miss} samples)'
    return acc / hits[None, :], fs, 'ok'


def main():
    cand = pd.read_csv(os.path.join(HERE, 'doublet_candidates.csv'))
    seq = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    ev = cand[cand['seq'] == seq].sort_values('time').reset_index(drop=True)
    # format='mixed': the saved timestamps vary in fractional-second precision,
    # which the default single-format parser rejects.
    ev['time'] = pd.to_datetime(ev['time'], utc=True, format='mixed')
    print(f'GATE G1 -- sequence {seq}: {len(ev)} events, '
          f'M{ev.mag.min():.2f}-{ev.mag.max():.2f}, '
          f'{ev.km_from_safod.mean():.1f} km from SAFOD\n')

    db = load_manifest()
    print(f'manifest: {len(db)} usable rows\n')

    results = []
    for _, e in ev.iterrows():
        X, fs, status = snippet(db, e['time'])
        if X is None:
            print(f'  {e["time"]:%Y-%m-%d %H:%M:%S}  M{e["mag"]:.2f}  SKIP: {status}')
            continue
        Xb = DASutils.bandpass2D_c(X - np.median(X, axis=0, keepdims=True),
                                   CODA_BAND[0], CODA_BAND[1], 1.0 / fs,
                                   zerophase=True)
        i_org = int(PRE_S * fs)

        def rms(w):
            i0 = i_org + int(w[0] * fs); i1 = i_org + int(w[1] * fs)
            i0, i1 = max(i0, 0), min(i1, Xb.shape[1])
            return np.sqrt(np.mean(Xb[:, i0:i1] ** 2, axis=1))

        noise, coda = rms(NOISE_WIN), rms(CODA_WIN)
        snr = coda / np.where(noise > 0, noise, np.nan)
        nch = int(np.sum(snr > PASS_SNR))
        print(f'  {e["time"]:%Y-%m-%d %H:%M:%S}  M{e["mag"]:.2f}  '
              f'median coda SNR {np.nanmedian(snr):6.2f}  '
              f'channels>{PASS_SNR:.0f}: {nch:4d}/{snr.size}'
              f'{"   <-- PASS" if nch >= PASS_NCH else ""}')
        results.append(dict(time=e['time'], mag=e['mag'], snr=snr, sec=Xb, fs=fs))

    if not results:
        print('\nG1 RESULT: FAIL -- no events could be extracted.')
        return

    S = np.vstack([r['snr'] for r in results])
    npass = int(np.sum(np.nanmedian(S, axis=0) > PASS_SNR))
    print(f'\nacross {len(results)} events: median-SNR channels above '
          f'{PASS_SNR:.0f} = {npass}/{S.shape[1]}')
    verdict = 'PASS' if npass >= PASS_NCH else 'FAIL'
    print(f'G1 RESULT: {verdict}  (need >= {PASS_NCH} channels)')
    if verdict == 'FAIL':
        print('  -> The plan says stop here. Coda that is not there cannot be '
              'processed into existence.')

    z = np.arange(S.shape[1]) * DX
    np.savez(os.path.join(HERE, f'g1_coda_snr_seq{seq}.npz'),
             snr=S, depth=z, times=np.array([str(r['time']) for r in results]),
             mags=np.array([r['mag'] for r in results]),
             coda_win=np.array(CODA_WIN), noise_win=np.array(NOISE_WIN),
             band=np.array(CODA_BAND), verdict=verdict)

    fig, ax = plt.subplots(1, 2, figsize=(14, 7))
    for r, s in zip(results, S):
        ax[0].semilogx(s, z, lw=0.8, alpha=0.6,
                       label=f'{r["time"]:%Y-%m-%d} M{r["mag"]:.2f}')
    ax[0].semilogx(np.nanmedian(S, axis=0), z, 'k-', lw=2.2, label='median')
    ax[0].axvline(PASS_SNR, ls='--', c='C3', lw=1.4)
    ax[0].invert_yaxis()
    ax[0].set(xlabel=f'coda SNR ({CODA_WIN[0]:.0f}-{CODA_WIN[1]:.0f} s lapse)',
              ylabel='distance along fiber (m)',
              title=f'A  Coda SNR vs depth — G1 {verdict}')
    ax[0].legend(fontsize=7)
    ax[0].grid(alpha=0.3)

    r = results[int(np.argmax([x['mag'] for x in results]))]
    Xb, fs = r['sec'], r['fs']
    t = (np.arange(Xb.shape[1]) - PRE_S * fs) / fs
    clip = np.percentile(np.abs(Xb), 99)
    ax[1].imshow(Xb.T, aspect='auto', cmap='gray_r',
                 extent=[0, z[-1], t[-1], t[0]], vmin=-clip, vmax=clip)
    for w, c, lab in [(NOISE_WIN, 'C0', 'noise'), (CODA_WIN, 'C3', 'coda')]:
        ax[1].axhspan(w[0], w[1], color=c, alpha=0.15)
        ax[1].text(z[-1] * 0.02, np.mean(w), lab, color=c, fontsize=9)
    ax[1].set(ylim=(15, -5), xlabel='distance along fiber (m)',
              ylabel='time after origin (s)',
              title=f'B  {r["time"]:%Y-%m-%d} M{r["mag"]:.2f}, '
                    f'{CODA_BAND[0]:.0f}-{CODA_BAND[1]:.0f} Hz')
    ax[1].grid(False)
    fig.suptitle('G1: can this fiber see earthquake coda?', fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, f'g1_coda_snr_seq{seq}.png')
    fig.savefig(out, dpi=140)
    print(f'\nSaved {out}')


if __name__ == '__main__':
    main()
