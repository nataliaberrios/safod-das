"""
Ambient-noise correlation between SF.MH029 and every DAS channel.

WHAT THIS IS FOR, AND WHY IT IS NOT MORE AMBIENT INTERFEROMETRY. awd_clean already
holds an extensive, versioned ambient F-K pipeline (v48-v50, injection-recovery,
directional audits, ground-truth validation, a Lellouch 2019 reproduction). Every
correlation in it -- and in stack_daily.py, sanity/sanity_cc.py and
awd_virtual_source.py -- has a DAS CHANNEL AT BOTH ENDS.

That is the gap. Both endpoints of a DAS-DAS correlation share one interrogator,
one laser and one clock, so a common-mode instrumental artefact is present at both
and is indistinguishable from a coherent arrival BY CONSTRUCTION. Injection-recovery
testing does not close this: the injected signal traverses the same shared path as
the data, so it validates the processing chain rather than the instrument.

SF.MH029 has a different digitiser and a different clock. A coherent MH029-DAS
arrival cannot be an interrogator artefact. That is the one thing this buys, and no
amount of work on the DAS-DAS side can substitute for it.

GEOMETRY, which is the best available at SAFOD:

    MH029   3-component, 1000 Hz, surveyed depth 2555.1 m, INSIDE the SAFOD hole
    fibre   873 usable channels, 0-864 m depth, same hole
    offset  0.00 km lateral -- the path runs straight up the borehole

For comparison the nearest HRSN station, CCRB, is 1.89 km away, which is 6.3
wavelengths at 10 Hz; station_geometry.py already concluded that at that separation
a station-channel correlation measures path, not instrument.

COVERAGE, checked against NCEDC. MH029 continuous data runs 2022-05-21 to about
2024-09-18 and then stops (metadata says operational to 3000, so this is an
archiving gap, not a dead station). The DAS spans 2024-05-07 to 2026-01-27, so the
overlap is 2024-05-07 to ~2024-09-18, about 4.3 months. Enough for a stable noise
Green's function over a ~2 km path; not enough for a time-lapse, and it does not
cover the confirmed repeater pairs, nearly all of which straddle 2024 and 2025.

--------------------------------------------------------------------------------
THE PREDICTION, REGISTERED BEFORE RUNNING, AND IT IS SHARP.

The path is vertical and its length is known: 2555.1 - z for a channel at depth z.
So the correlation arrival must appear at

    lag(z) = (2555.1 - z) / V

i.e. it must move LINEARLY with channel depth, with slope 1/V. Across the 864 m
aperture that is a lag change of

    V = 3000 m/s -> 0.288 s      V = 4000 -> 0.216 s      V = 5000 -> 0.173 s

and the absolute lag at the wellhead is 0.51-0.85 s.

PASS  a coherent arrival whose lag decreases linearly with depth, with a fitted
      slope giving V in a physically sensible 2000-6000 m/s, and an intercept
      consistent with 2555.1 m at that V.
FAIL  a flat, incoherent, or zero-lag-peaked gather. Any of those means there is no
      interpretable path and the measurement stops here.

A flat arrival is the specific failure to watch: it would indicate correlated noise
common to both instruments (mains, GPS timing, shared telemetry) rather than a
propagating wavefield. That is why the linear moveout, not the mere presence of a
peak, is the criterion.

CONTROL. The same processing with MH029 taken from a DIFFERENT DAY than the DAS,
which preserves every statistical property and destroys the physical link. Any
feature surviving that is instrumental or numerical.
--------------------------------------------------------------------------------

Preprocessing follows the project's existing ambient convention exactly -- 5-20 Hz
bandpass, 5 s running-absolute-mean normalisation -- so this is comparable to the
DAS-DAS results rather than a different recipe. awd_clean is imported, never
modified.
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'awd_clean'))
from g1_coda_snr import load_manifest                            # noqa: E402
import h5py                                                      # noqa: E402
from scipy.signal import butter, sosfiltfilt, detrend, decimate  # noqa: E402
from scipy.ndimage import uniform_filter1d                       # noqa: E402

OUT = os.environ.get('MH_OUT', '/scratch/groups/ettore88/nberrios/safod_mh029')
os.makedirs(OUT, exist_ok=True)

CH_LO, CH_HI = 23, 896
DX = 1.0209523
FS = 500.0                      # DAS output rate; MH029 is 1000 Hz -> decimate
BAND = (5.0, 20.0)              # the project's ambient band
NORM_S = 5.0                    # running-abs-mean window, project convention
MAX_LAG = 1.5                   # must span the predicted 0.34-0.85 s
MH029_DEPTH = 2555.1
STA = ('SF', 'MH029', '01', 'GP1')      # GP1 first; all three saved
DAY = os.environ.get('MH_DAY', '2024-06-01')
NULL_OFFSET_DAYS = int(os.environ.get('MH_NULL_DAYS', 0))   # 0 = causal, >0 = control


def preprocess(x, fs, band=BAND, norm_seconds=NORM_S):
    """Identical to awd_clean/ambient_transfer_test.preprocess, reimplemented here
    so this script carries no write dependency on that tree."""
    x = detrend(np.atleast_2d(x), axis=1, type='linear')
    sos = butter(4, list(band), btype='bandpass', fs=fs, output='sos')
    x = sosfiltfilt(sos, x, axis=1)
    nwin = max(3, int(norm_seconds * fs))
    m = uniform_filter1d(np.abs(x), size=nwin, axis=1, mode='nearest')
    floor = np.percentile(m, 5, axis=1, keepdims=True) * 0.1 + 1e-12
    return x / np.maximum(m, floor)


def fetch_mh029(t0, t1):
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    c = Client('NCEDC', timeout=180)
    st = c.get_waveforms(STA[0], STA[1], STA[2], 'GP?',
                         UTCDateTime(str(t0)), UTCDateTime(str(t1)))
    st.merge(method=1, fill_value=0)
    st.sort()
    return st


def das_window(fn):
    with h5py.File(fn, 'r') as h:
        a = h['Acquisition'].attrs
        if abs(float(a.get('MaximumFrequency', 250.0)) * 2 - FS) > 1:
            return None                      # skip 5000 Hz event-triggered files
        g = [k for k in h['Acquisition'].keys() if 'Raw' in k][0]
        X = np.asarray(h[f'Acquisition/{g}/RawData'][:, :], dtype=np.float64)
    if X.shape[0] > X.shape[1]:
        X = X.T
    if X.shape[0] < CH_HI:
        return None
    A = X[CH_LO:CH_HI]
    return np.diff(A, axis=1, prepend=A[:, :1])       # phase -> strain rate


def main():
    day = pd.Timestamp(DAY, tz='UTC')
    lab = 'null' if NULL_OFFSET_DAYS else 'causal'
    print(f'MH029 x DAS, {day:%Y-%m-%d}, mode={lab}'
          f'{f" (MH029 shifted {NULL_OFFSET_DAYS} d)" if NULL_OFFSET_DAYS else ""}',
          flush=True)

    mh_day = day + pd.Timedelta(days=NULL_OFFSET_DAYS)
    try:
        st = fetch_mh029(mh_day, mh_day + pd.Timedelta(days=1))
    except Exception as e:
        print(f'  MH029 fetch failed: {str(e)[:120]}'); return
    print(f'  MH029: {len(st)} traces, {[t.stats.channel for t in st]}, '
          f'{st[0].stats.sampling_rate:g} Hz, '
          f'{st[0].stats.endtime - st[0].stats.starttime:.0f} s', flush=True)

    db = load_manifest()
    sel = db[(db['t0'] >= day) & (db['t0'] < day + pd.Timedelta(days=1))]
    sel = sel[sel['fn'].map(os.path.exists)].sort_values('t0')
    print(f'  DAS: {len(sel)} files', flush=True)
    if len(sel) < 20 or not len(st):
        print('  insufficient data'); return

    ml = int(MAX_LAG * FS)
    lags = np.arange(-ml, ml + 1) / FS
    nch = CH_HI - CH_LO
    acc = {tr.stats.channel: np.zeros((nch, lags.size)) for tr in st}
    nstack = 0

    for k, (_, r) in enumerate(sel.iterrows()):
        A = das_window(r['fn'])
        if A is None:
            continue
        n = A.shape[1]
        Ap = preprocess(A, FS)
        Fa = np.fft.rfft(Ap, n=2 * n, axis=-1)
        # matching MH029 window; when NULL_OFFSET_DAYS != 0 this is a different day
        w0 = pd.Timestamp(r['t0']) + pd.Timedelta(days=NULL_OFFSET_DAYS)
        for tr in st:
            seg = tr.slice(starttime=__import__('obspy').UTCDateTime(str(w0)),
                           endtime=__import__('obspy').UTCDateTime(
                               str(w0 + pd.Timedelta(seconds=n / FS))))
            d = np.asarray(seg.data, dtype=np.float64)
            if d.size < 10:
                continue
            if abs(seg.stats.sampling_rate - FS) > 1:
                q = int(round(seg.stats.sampling_rate / FS))
                if q > 1:
                    d = decimate(d, q, ftype='fir', zero_phase=True)
            if d.size < n:
                d = np.pad(d, (0, n - d.size))
            d = preprocess(d[:n], FS)[0]
            Fb = np.fft.rfft(d, n=2 * n)
            cc = np.fft.irfft(Fa * np.conj(Fb), n=2 * n, axis=-1)
            cc = np.concatenate([cc[:, -ml:], cc[:, :ml + 1]], axis=1)
            nrm = np.sqrt((Ap ** 2).sum(axis=1) * (d ** 2).sum())
            acc[tr.stats.channel] += cc / np.maximum(nrm[:, None], 1e-30)
        nstack += 1
        if k % 120 == 0:
            print(f'    {k}/{len(sel)}  stacked {nstack}', flush=True)

    if nstack == 0:
        print('  nothing stacked'); return
    z = (np.arange(CH_LO, CH_HI) - CH_LO) * DX
    out = os.path.join(OUT, f'mh029_das_{day:%Y%m%d}_{lab}.npz')
    np.savez_compressed(out, lags=lags, z=z, nstack=nstack,
                        **{f'cc_{c}': v / nstack for c, v in acc.items()})
    print(f'\n  wrote {out}  ({nstack} windows stacked)')

    # first look: where is the peak, and does it move with depth?
    for c, v in acc.items():
        g = v / nstack
        pk = lags[np.argmax(np.abs(g), axis=1)]
        good = np.isfinite(pk)
        if good.sum() > 100:
            A_ = np.column_stack([np.ones(good.sum()), z[good]])
            b, *_ = np.linalg.lstsq(A_, pk[good], rcond=None)
            V = -1.0 / b[1] if b[1] != 0 else np.inf
            print(f'  {c}: peak lag {pk[0]:+.3f} s at z=0 -> {pk[-1]:+.3f} s at '
                  f'z={z[-1]:.0f} m; fitted slope gives V = {V:.0f} m/s')
            print(f'      predicted intercept for 2555.1 m at that V: '
                  f'{MH029_DEPTH/abs(V):.3f} s' if np.isfinite(V) and V != 0 else '')


if __name__ == '__main__':
    main()
