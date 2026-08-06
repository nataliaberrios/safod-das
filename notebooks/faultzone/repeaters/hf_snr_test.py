"""
Is there usable signal at 80-100 Hz? The one thing that kills the stress-drop idea.

WHY THIS MATTERS. Corner frequency is where the source spectrum bends, and it is
what stress drop is computed from. For these events (3 MPa Brune, Vs 3.5 km/s):

    M0.6 -> 114 Hz     M1.0 -> 72 Hz     M1.5 -> 40 Hz     M2.0 -> 23 Hz

Everything cached so far was extracted at 100 Hz, Nyquist 50 Hz, so the corner of
almost every one of these events was invisible. The files are actually 500 Hz
(OutputDataRate in the HDF5 header; the manifest's "fs 10000" is the interrogator
pulse rate, not the output rate). Nyquist is therefore 250 Hz.

The binding limit is then the 16.335 m gauge length, which averages strain over
that distance and rolls off as |sinc(pi f G / V)|:

    V = 3000 m/s : -3 dB at ~81 Hz, first null 184 Hz
    V = 2000 m/s : -3 dB at ~54 Hz, first null 122 Hz

So the question is empirical: at 60-120 Hz, on these magnitudes, at these
distances, is the earthquake above the noise floor on enough channels to fit a
spectrum? If not, the stress-drop route is closed and nothing else needs building.

This reads raw HDF5 directly at full rate rather than the 100 Hz cache. It touches
only a handful of files, so it does not hammer Lustre.
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
from g1_coda_snr import load_manifest                              # noqa: E402

for p in ['/home/groups/edunham/nberrios/safod_das/DAS-utilities/python',
          '/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python']:
    if os.path.exists(p) and p not in sys.path:
        sys.path.insert(0, p)
import DASutils                                                    # noqa: E402

CH_LO, CH_HI = 23, 896          # G0 usable range
WELLHEAD = 23
DX = 1.0210
GAUGE = 16.335
PRE_S, POST_S = 5.0, 15.0
FMIN, FMAX = 1.0, 240.0         # read wide; Nyquist is 250
SIG_WIN = (0.0, 3.0)            # about the origin: P, S and early coda
NOISE_WIN = (-4.5, -0.5)        # pre-event
BANDS = [(5, 10), (10, 20), (20, 40), (40, 60), (60, 80),
         (80, 100), (100, 140), (140, 200)]


def gauge_response(f, v):
    """|sinc(pi f G / V)| -- the strain-averaging rolloff."""
    return np.abs(np.sinc(f * GAUGE / v))


def read_event(db, origin):
    a = origin - pd.Timedelta(seconds=PRE_S)
    b = origin + pd.Timedelta(seconds=POST_S)
    sel = db[(db['t0'] < b) & (db['t1'] > a)].sort_values('t0')
    if sel.empty:
        return None, None
    acc = hits = None
    fs = dt = npts = None
    for _, r in sel.iterrows():
        if not os.path.exists(r['fn']):
            continue
        try:
            D, info = DASutils.readFile_HDF(
                [r['fn']], FMIN, FMAX, verbose=0, preproc=True, diff=True,
                taper=False, desampling=False, nChbuffer=900,
                system='OptaSense')
        except Exception as e:
            print(f'    read fail: {str(e)[:50]}'); continue
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
        return None, fs
    return acc / hits[None, :], fs


def band_snr(X, fs):
    """Per-channel signal/noise amplitude ratio in each band."""
    i0 = int((PRE_S + SIG_WIN[0]) * fs); i1 = int((PRE_S + SIG_WIN[1]) * fs)
    j0 = int((PRE_S + NOISE_WIN[0]) * fs); j1 = int((PRE_S + NOISE_WIN[1]) * fs)
    S = X[:, i0:i1]; N = X[:, j0:j1]
    fS = np.fft.rfftfreq(S.shape[1], 1 / fs)
    fN = np.fft.rfftfreq(N.shape[1], 1 / fs)
    PS = np.abs(np.fft.rfft(S * np.hanning(S.shape[1]), axis=-1)) ** 2
    PN = np.abs(np.fft.rfft(N * np.hanning(N.shape[1]), axis=-1)) ** 2
    out = {}
    for lo, hi in BANDS:
        ms = (fS >= lo) & (fS < hi); mn = (fN >= lo) & (fN < hi)
        if ms.sum() < 2 or mn.sum() < 2:
            continue
        # power spectral density, so divide by the number of bins and duration
        ps = PS[:, ms].mean(axis=1) / S.shape[1]
        pn = PN[:, mn].mean(axis=1) / N.shape[1]
        out[(lo, hi)] = np.sqrt(ps / np.maximum(pn, 1e-300))
    return out


def main():
    ev = pd.read_csv(os.path.join(HERE, 'phaseA_events.csv'))
    ev['time'] = pd.to_datetime(ev.time, utc=True, format='mixed')
    ev = ev[ev.cov_full].reset_index(drop=True)
    # test the largest covered events plus two of the confirmed-sequence members
    want = list(ev.nlargest(3, 'mag').index)
    for d in ('2024-07-08', '2025-04-06'):
        m = ev[ev.time.dt.strftime('%Y-%m-%d') == d]
        if len(m):
            want.append(m.index[0])
    want = list(dict.fromkeys(want))[:5]

    db = load_manifest()
    print(f'manifest {len(db)} intervals\n', flush=True)
    print('gauge-length rolloff |sinc(pi f G / V)|, G = 16.335 m')
    ff = np.array([40, 60, 80, 100, 120, 140])
    for v in (2000, 3000, 4000):
        print(f'  V={v:5d} m/s: ' +
              '  '.join(f'{f:3.0f}Hz {gauge_response(f, v):.2f}' for f in ff))
    print()

    results = {}
    for k in want:
        e = ev.loc[k]
        print(f'{e.tag}  M{e.mag:.2f}  depth {e.depth:.2f} km  '
              f'{e.km_safod:.2f} km from SAFOD', flush=True)
        X, fs = read_event(db, e.time)
        if X is None:
            print('   could not assemble a gap-free window\n'); continue
        print(f'   fs = {fs:g} Hz  (Nyquist {fs/2:g})  shape {X.shape}')
        X = X[CH_LO:CH_HI]
        snr = band_snr(X, fs)
        results[e.tag] = (e.mag, snr)
        print(f'   {"band":>10}{"med SNR":>10}{"ch>3":>7}{"ch>5":>7}')
        for (lo, hi), s in snr.items():
            print(f'   {f"{lo}-{hi}":>10}{np.median(s):10.1f}'
                  f'{int((s > 3).sum()):7d}{int((s > 5).sum()):7d}')
        print(flush=True)

    if not results:
        print('no events read'); return

    print('VERDICT -- can corner frequencies be measured?')
    ok = True
    for tag, (mag, snr) in results.items():
        m0 = 10 ** (1.5 * mag + 9.1)
        r = (7 * m0 / (16 * 3e6)) ** (1 / 3)
        fc = 0.37 * 3500 / r
        band = next((b for b in BANDS if b[0] <= fc < b[1]), None)
        if band and band in snr:
            n3 = int((snr[band] > 3).sum())
            print(f'  {tag} M{mag:.2f}: corner ~{fc:.0f} Hz, band {band[0]}-{band[1]} Hz, '
                  f'{n3} channels with SNR>3')
            if n3 < 100:
                ok = False
        else:
            print(f'  {tag} M{mag:.2f}: corner ~{fc:.0f} Hz is outside the '
                  f'tested bands')
            ok = False
    print()
    if ok:
        print('  -> usable signal at the corner frequency on enough channels.')
        print('     Re-extract at 500 Hz and proceed to spectral fitting.')
    else:
        print('  -> NOT enough high-frequency signal at the corner for at least')
        print('     one event. Either restrict to larger events (lower corner),')
        print('     or the stress-drop route is closed. Report which.')

    fig, ax = plt.subplots(figsize=(8, 5))
    for tag, (mag, snr) in results.items():
        f = [0.5 * (lo + hi) for lo, hi in snr]
        s = [np.median(v) for v in snr.values()]
        ax.loglog(f, s, 'o-', label=f'{tag[3:11]} M{mag:.2f}')
    ax.axhline(3, color='C3', ls='--', label='SNR = 3')
    for v, c in ((2000, '0.7'), (3000, '0.4')):
        fg = np.logspace(0.7, 2.3, 100)
        ax.loglog(fg, 30 * gauge_response(fg, v), ':', color=c,
                  label=f'gauge rolloff, V={v}')
    ax.set(xlabel='frequency (Hz)', ylabel='median channel SNR (amplitude)',
           title='High-frequency SNR: is the corner frequency observable?')
    ax.legend(fontsize=8); ax.grid(alpha=.3, which='both')
    fig.tight_layout()
    p = os.path.join(HERE, 'hf_snr_test.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
