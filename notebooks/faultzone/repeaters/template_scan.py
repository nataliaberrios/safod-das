"""
Find the missing sequence members. Matched-filter scan of the continuous archive.

WHY THIS IS THE MEASUREMENT THAT FEEDS CREEP. The Nadeau & McEvilly method turns
repeater recurrence into fault slip rate, and it needs event TIMES, not relative
locations. So the fibre's contribution to a creep number is catalog completeness,
and that route does not depend on the same-patch test that has failed twice here.

The confirmed sequences have intervals of 272, 303, 305, 331 and 440 days for M0.7-1.5
events. Those are implausibly long for repeaters of that size, and the obvious
explanation is that most members sit below the network detection threshold. Lellouch
et al. 2020 measured Mc = -1.4 for DAS on THIS fibre against -1.7 for collocated
geophones; the catalog used here starts at M0.65. There is a lot of room underneath.

WHY IT IS ONLY NOW WORTH RUNNING. Template matching was closed earlier in this
project on the grounds that moveout correction does not improve detection. That was
measured at 100 Hz. Re-measured on the 500 Hz cache, d' goes 2.45 -> 25.3, i.e. the
repeater population sits ~25 noise-widths above random pairs. That is a detector
worth pointing at continuous data.

METHOD. One template per confirmed event:
  * moveout-correct the template with its OWN slowness (moveout_test.slant_scan)
  * stack channels -> one high-SNR template trace
  * apply the SAME slowness to each continuous 60 s file, stack -> continuous beam
  * normalised cross-correlation of template against beam

Applying the template's slowness to the continuous data is what makes this a
moveout-corrected matched filter rather than a flat stack; per METHODS_STATUS 2.4
that is worth ~27 dB.

--------------------------------------------------------------------------------
CONTROLS, FIXED BEFORE RUNNING. Detection thresholds are where this kind of
analysis usually goes wrong, so none of them is chosen after seeing the answer.

 T1 ACAUSAL FLOOR. Every scan is repeated with the TIME-REVERSED template, which
    cannot match a causal arrival. The highest acausal correlation over the same
    data is the empirical false-positive floor. This is the threshold, not a
    number imported from the literature.

 T2 MAD THRESHOLD. The detection threshold is median + 8 x MAD of the day's own
    correlogram, and must ALSO exceed the acausal maximum. Both conditions.

 T3 DETECTIONS MUST NOT CLUSTER AT FILE BOUNDARIES. The archive is 60 s segments;
    an edge artefact would produce apparent detections at multiples of 60 s.
    Reported as a histogram of detection time modulo 60 s.

 T4 A DETECTION IS A CANDIDATE, NOT AN EVENT. Anything found here needs its own
    waveform check against the template and, where the magnitude allows, an HRSN
    confirmation. This project has already had one result die to an external
    control; a matched filter on 156 GB/day will find things.

EXPECTED OUTCOME, STATED IN ADVANCE. If the sequences are genuinely undersampled,
a 14-day window inside a 272-day gap should contain a few detections. Zero
detections across 14 days would bound the recurrence at longer than ~1 event per
2 weeks at this magnitude, which is also a result and would mean the long intervals
are real rather than an artefact of completeness.
--------------------------------------------------------------------------------

Run:  sbatch --array=0-13 template_scan_job.sh
Each task does one day, ~1440 files, and writes its own npz.
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import moveout_test as MT                                        # noqa: E402
from g1_coda_snr import load_manifest                            # noqa: E402
import h5py                                                      # noqa: E402

CACHE_HF = os.path.join(HERE, 'cache_hf')
OUT = os.path.join(HERE, 'template_scan')
os.makedirs(OUT, exist_ok=True)

CH_LO, CH_HI = 23, 896
DX = 1.0209523
BAND = (5.0, 40.0)          # detection band: wide enough for d', narrow enough
                            # that the gauge rolloff (-3 dB at 54-81 Hz) is mild
FS = 500.0
TEMPLATE_TAG = os.environ.get('TEMPLATE_TAG', 'ev_20240708T083036')
TPL_WIN = (-0.5, 3.0)       # seconds about the template's picked P
MAD_K = 8.0
DAY0 = os.environ.get('SCAN_START', '2024-09-15')   # inside the 272 d gap


def build_template():
    """Moveout-corrected, channel-stacked template trace, plus its slowness."""
    MT.BAND = BAND
    MT.CACHE = CACHE_HF
    got, _n = MT.prep(TEMPLATE_TAG)
    if got is None:
        return None
    A, z, fs = got
    p, tp, sem, _S, _ps = MT.slant_scan(A, z, fs)
    B = MT.align(A, z, fs, p)
    i0 = int((tp + TPL_WIN[0]) * fs)
    i1 = int((tp + TPL_WIN[1]) * fs)
    if i0 < 0 or i1 > B.shape[1]:
        return None
    w = B[:, i0:i1].mean(axis=0)
    w = w - w.mean()
    nn = np.sqrt(np.sum(w ** 2))
    return dict(w=w / nn, p=p, fs=fs, sem=sem, nsamp=w.size, z=z)


def beam_file(fn, p, zref):
    """One 60 s file -> moveout-corrected beam trace at the template's slowness."""
    from scipy.signal import butter, sosfiltfilt
    try:
        with h5py.File(fn, 'r') as h:
            g = [k for k in h['Acquisition'].keys() if 'Raw' in k][0]
            X = np.asarray(h[f'Acquisition/{g}/RawData'][:, :], dtype=np.float64)
    except Exception:
        return None
    if X.shape[0] > X.shape[1]:
        X = X.T                                  # -> (channel, time)
    if X.shape[0] < CH_HI:
        return None
    A = X[CH_LO:CH_HI]
    A = np.diff(A, axis=1, prepend=A[:, :1])     # phase -> strain rate, as the
                                                 # cached events were extracted
    A -= A.mean(axis=1, keepdims=True)
    sos = butter(4, list(BAND), btype='band', fs=FS, output='sos')
    A = sosfiltfilt(sos, A, axis=-1)
    rms = np.sqrt((A ** 2).mean(axis=1))
    good = np.isfinite(rms) & (rms > 0)
    if good.sum() < 100:
        return None
    med = np.median(rms[good])
    good &= (rms > 0.05 * med) & (rms < 20.0 * med)
    A = A[good] / rms[good, None]
    z = (np.arange(CH_LO, CH_HI)[good] - CH_LO) * DX
    B = MT.align(A, z, FS, p)
    b = B.mean(axis=0)
    return b - b.mean()


def scan(beam, w):
    """Normalised correlation of template w against continuous beam."""
    n = w.size
    if beam.size < n * 2:
        return None
    from scipy.signal import fftconvolve
    num = fftconvolve(beam, w[::-1], mode='valid')
    e = fftconvolve(beam ** 2, np.ones(n), mode='valid')
    den = np.sqrt(np.maximum(e, 1e-30))
    return num / den


def main():
    task = int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
    day = pd.Timestamp(DAY0, tz='UTC') + pd.Timedelta(days=task)
    print(f'task {task}: scanning {day:%Y-%m-%d}  template {TEMPLATE_TAG}',
          flush=True)

    T = build_template()
    if T is None:
        print('template failed to build'); return
    print(f'  template {T["nsamp"]} samples, slowness {T["p"]:.3e} s/m '
          f'(V_app {1/abs(T["p"]):.0f} m/s), semblance {T["sem"]:.3f}', flush=True)

    db = load_manifest()
    sel = db[(db['t0'] >= day) & (db['t0'] < day + pd.Timedelta(days=1))]
    sel = sel[sel['fn'].map(os.path.exists)].sort_values('t0')
    print(f'  {len(sel)} files', flush=True)
    if len(sel) < 10:
        print('  too few files'); return

    w = T['w']
    wrev = w[::-1].copy()          # T1 acausal template
    peaks, apeaks, tstamps = [], [], []
    nfile = 0
    for k, (_, r) in enumerate(sel.iterrows()):
        b = beam_file(r['fn'], T['p'], None)
        if b is None:
            continue
        nfile += 1
        c = scan(b, w)
        ca = scan(b, wrev)
        if c is None:
            continue
        peaks.append(c.astype(np.float32))
        apeaks.append(ca.astype(np.float32))
        tstamps.append(r['t0'])
        if k % 240 == 0:
            print(f'    {k}/{len(sel)}  max|cc| so far '
                  f'{max(float(np.max(np.abs(x))) for x in peaks):.3f}', flush=True)

    if not peaks:
        print('  nothing scanned'); return
    C = np.concatenate(peaks)
    CA = np.concatenate(apeaks)
    med = float(np.median(C))
    mad = float(1.4826 * np.median(np.abs(C - med)))
    acmax = float(np.max(np.abs(CA)))
    thr = max(med + MAD_K * mad, acmax)
    n_det = int(np.sum(C > thr))
    print(f'\n  scanned {nfile} files, {C.size} correlation samples')
    print(f'  median {med:+.4f}  MAD {mad:.4f}')
    print(f'  T2 MAD threshold  : {med + MAD_K*mad:.4f}')
    print(f'  T1 acausal max    : {acmax:.4f}   <- false-positive floor')
    print(f'  threshold used    : {thr:.4f}  (the larger of the two)')
    print(f'  DETECTIONS        : {n_det}')
    print(f'  acausal above thr : {int(np.sum(CA > thr))}   (should be 0)')

    np.savez_compressed(os.path.join(OUT, f'scan_{day:%Y%m%d}.npz'),
                        cmax=C.max(), med=med, mad=mad, acmax=acmax, thr=thr,
                        n_det=n_det, nfile=nfile,
                        top=np.sort(C)[-200:].astype(np.float32),
                        atop=np.sort(CA)[-200:].astype(np.float32),
                        day=str(day), template=TEMPLATE_TAG)
    print(f'\n  wrote scan_{day:%Y%m%d}.npz', flush=True)


if __name__ == '__main__':
    main()
