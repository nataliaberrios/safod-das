"""
Matched-filter scan of the ENTIRE continuous archive, all templates.

SCOPE. 574,301 file intervals, 2024-05-07 to 2026-01-27, 398.6 days of actual
recording, ~62 TB. Fourteen templates -- the unique events in the 9 pairs confirmed
above the acausal floor by geometry AND waveform.

WHY THIS IS THE MEASUREMENT THAT MATTERS. Recurrence-based creep needs event TIMES.
The confirmed sequences show intervals of 39-331 days for M1.1-1.9, and the derived
slip rates come out at 51-421 mm/yr against a geodetic 25-30 -- 2-14x too high. The
most likely reading is that most sequence members are missing: the catalog is
complete to M0.65 while this array, calibrated against a repeater of known size,
detects to M ~ -0.25. That is 0.9 magnitude units, ~8x more events at b = 1.

A 14-day pilot found nothing, which is expected at these recurrence intervals and
was never going to settle it. 398 days is the experiment.

--------------------------------------------------------------------------------
EFFICIENCY, which is what makes this feasible at all.

Each template has its OWN slowness, so a beam formed for one is wrong for another,
and 14 independent alignments per file would be 14x the cost. But alignment is a
per-channel phase shift, so the channel sum can be done in the FREQUENCY domain:

    beam_p(f) = (1/N) * sum_ch  F[ch, f] * exp(2i*pi*p*z_ch*f)

One forward FFT of the 700 x 30000 record serves all 14 slownesses; each then costs
one channel-sum and one inverse FFT. Roughly 1.5 s per file total, so ~239 core-hours
over the archive, or ~1.2 h per task at 200 tasks.

--------------------------------------------------------------------------------
THRESHOLDS. Same discipline as the pilot, which is the part most likely to go wrong
at this scale: with ~10^10 correlation samples, a threshold that is slightly too low
produces thousands of false detections and one that is slightly too high produces
none, and either can be made to look like a result.

  * The ACAUSAL (time-reversed) template is run on every file alongside the causal
    one. It cannot match a real arrival, so its maximum is the empirical
    false-positive floor. That is the threshold. No literature value is imported.
  * The floor is computed PER DAY PER TEMPLATE, because the pilot showed it varies
    (0.438-0.543 across days) and a fixed cut would be wrong in both directions.
  * Thresholding happens in POST, not here: this pass stores the per-day median,
    MAD, acausal maximum and the top 500 causal peaks with their times. That keeps
    the decision reversible and means the scan does not have to be repeated if the
    criterion changes.

WHAT A DETECTION IS NOT. A peak above the floor is a candidate. It is not an event
until it survives a waveform check against the template and, where the magnitude
allows, an HRSN confirmation. A matched filter over 62 TB will find things.

STATED IN ADVANCE:
  * If the sequences are undersampled, members should appear at intervals well
    below the catalog's 39-331 days, and the creep rates should FALL toward the
    geodetic 25-30 mm/yr as the recurrence series fills in. That is the prediction.
  * If 398 days of scanning at M ~ -0.25 yields no new members, the long intervals
    are real, the sequences are genuinely sparse, and the 2-14x creep discrepancy
    needs a different explanation -- most likely that these pairs are neighbours
    rather than repeaters, which their dM of 0.26-0.71 already suggests.
Both outcomes are informative and neither requires the same-patch test.
--------------------------------------------------------------------------------

Run:  sbatch --array=0-199 full_scan_job.sh
Each task takes every 200th file, so a dead task costs an even sample of the
archive rather than a contiguous block of it.
"""
import os
import sys
import numpy as np
import pandas as pd
import h5py
import time as _time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import moveout_test as MT                                        # noqa: E402
from g1_coda_snr import load_manifest                            # noqa: E402

CACHE_HF = os.path.join(HERE, 'cache_hf')
OUT = os.environ.get('SCAN_OUT', '/scratch/groups/ettore88/nberrios/safod_fullscan')
os.makedirs(OUT, exist_ok=True)

CH_LO, CH_HI = 23, 896
DX = 1.0209523
BAND = (5.0, 40.0)
FS = 500.0
TPL_WIN = (-0.5, 3.0)
KEEP_TOP = 500          # causal peaks retained per (day, template)


TPL_CACHE = os.path.join(HERE, 'full_scan_templates.npz')
# FROZEN CHANNEL SET. Per-file bad-channel rejection makes z vary file to file,
# which forces the 700 x 15001 phase matrix to be rebuilt for every file and every
# template -- the dominant cost. Freezing the channel set makes the phase matrices
# precomputable once per task. It is also the better scientific choice: a frozen
# selection region is what this project's conventions require, and a channel set
# that changes with the noise of each individual file is not one.
CH_KEEP = None          # filled by freeze_channels()


def freeze_channels(db, nsample=40, seed=0):
    """Channels usable in most files, decided ONCE from a sample of the archive."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(db), size=min(nsample, len(db)), replace=False)
    votes = np.zeros(CH_HI - CH_LO)
    n = 0
    for k in idx:
        fn = db.iloc[int(k)]['fn']
        try:
            with h5py.File(fn, 'r') as h:
                g = [q for q in h['Acquisition'].keys() if 'Raw' in q][0]
                X = np.asarray(h[f'Acquisition/{g}/RawData'][::10, :], dtype=np.float64)
        except Exception:
            continue
        if X.shape[0] > X.shape[1]:
            X = X.T
        if X.shape[0] < CH_HI:
            continue
        A = X[CH_LO:CH_HI]
        rms = np.sqrt((np.diff(A, axis=1) ** 2).mean(axis=1))
        ok = np.isfinite(rms) & (rms > 0)
        if ok.sum() < 100:
            continue
        med = np.median(rms[ok])
        votes += ((rms > 0.05 * med) & (rms < 20.0 * med) & ok).astype(float)
        n += 1
    if n == 0:
        return np.ones(CH_HI - CH_LO, bool)
    return (votes / n) > 0.8          # usable in >80% of sampled files


def build_templates():
    """One moveout-corrected, channel-stacked trace per confirmed-pair event.

    Cached to disk: this costs ~22 min (14 slant scans over 81 slownesses each) and
    is identical for every array task, so rebuilding it 200 times would waste ~73
    core-hours.

    Templates are resolved from the tag_i/tag_j columns, NOT by matching the 8-char
    date in the pair label. The date-matching version silently built 32 templates
    instead of 14, because several of those dates carry more than one catalogued
    event.
    """
    if os.path.exists(TPL_CACHE):
        d = np.load(TPL_CACHE, allow_pickle=True)
        return [dict(tag=str(t), w=w[np.isfinite(w)], p=float(p), sem=float(s))
                for t, w, p, s in zip(d['tag'], d['w'], d['p'], d['sem'])]
    D = pd.read_csv(os.path.join(HERE, 'ddrt_pair_cc.csv'))
    ac = D.acausal.max()
    C = D[D.cc > ac]
    if 'tag_i' in C:
        full = sorted(set(C.tag_i) | set(C.tag_j))
    else:
        raise SystemExit('ddrt_pair_cc.csv lacks tag_i/tag_j; rerun ddrt_pair_cc.py')
    MT.BAND = BAND
    MT.CACHE = CACHE_HF
    T = []
    for tag in full:
        got, _ = MT.prep(tag)
        if got is None:
            continue
        A, z, fs = got
        if abs(fs - FS) > 1:          # the three 5000 Hz event-triggered files
            continue
        p, tp, sem, _, _ = MT.slant_scan(A, z, fs)
        B = MT.align(A, z, fs, p)
        i0, i1 = int((tp + TPL_WIN[0]) * fs), int((tp + TPL_WIN[1]) * fs)
        if i0 < 0 or i1 > B.shape[1]:
            continue
        w = B[:, i0:i1].mean(axis=0)
        w -= w.mean()
        n = np.sqrt((w ** 2).sum())
        if n > 0:
            T.append(dict(tag=tag, w=w / n, p=p, sem=sem))
    if T:
        L = max(len(t['w']) for t in T)
        np.savez_compressed(
            TPL_CACHE,
            tag=np.array([t['tag'] for t in T]),
            w=np.array([np.pad(t['w'], (0, L - len(t['w'])),
                               constant_values=np.nan) for t in T]),
            p=np.array([t['p'] for t in T]),
            sem=np.array([t['sem'] for t in T]))
    return T


def beams(fn, slow, phase=None, keep=None, sos=None):
    """All beams for one file from ONE forward FFT. Returns (nslow, nt) or None.

    `phase` is the precomputed exp(2i.pi.p.z.f) stack, valid only because the
    channel set is frozen. Rebuilding it per file was the dominant cost.
    """
    from scipy.signal import butter, sosfiltfilt
    try:
        with h5py.File(fn, 'r') as h:
            g = [k for k in h['Acquisition'].keys() if 'Raw' in k][0]
            X = np.asarray(h[f'Acquisition/{g}/RawData'][:, :], dtype=np.float64)
    except Exception:
        return None
    if X.shape[0] > X.shape[1]:
        X = X.T
    if X.shape[0] < CH_HI:
        return None
    A = X[CH_LO:CH_HI][keep]
    A = np.diff(A, axis=1, prepend=A[:, :1])       # phase -> strain rate
    A -= A.mean(axis=1, keepdims=True)
    if sos is None:
        sos = butter(4, list(BAND), btype='band', fs=FS, output='sos')
    A = sosfiltfilt(sos, A, axis=-1)
    rms = np.sqrt((A ** 2).mean(axis=1))
    ok = np.isfinite(rms) & (rms > 0)
    if ok.sum() < 100:
        return None
    A = A[ok] / rms[ok, None]
    n = A.shape[1]
    F = np.fft.rfft(A, axis=-1).astype(np.complex64)   # one expensive transform
    out = np.empty((len(slow), n))
    for k in range(len(slow)):
        out[k] = np.fft.irfft((F * phase[k][ok]).mean(axis=0).astype(np.complex128), n=n)
    return out - out.mean(axis=1, keepdims=True)


def corr(beam, w):
    from scipy.signal import fftconvolve
    n = w.size
    if beam.size < 2 * n:
        return None
    num = fftconvolve(beam, w[::-1], mode='valid')
    den = np.sqrt(np.maximum(fftconvolve(beam ** 2, np.ones(n), mode='valid'), 1e-30))
    return num / den


def main():
    task = int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
    ntask = os.environ.get('SLURM_ARRAY_TASK_COUNT')
    if ntask:
        ntask = int(ntask)
    else:
        lo = int(os.environ.get('SLURM_ARRAY_TASK_MIN', 0))
        hi = int(os.environ.get('SLURM_ARRAY_TASK_MAX', 0))
        ntask = hi - lo + 1

    T = build_templates()
    if not T:
        print('no templates built'); return
    slow = [t['p'] for t in T]
    print(f'{len(T)} templates', flush=True)
    for t in T:
        print(f"  {t['tag']}  p={t['p']:+.3e} s/m  V_app={1/abs(t['p']):.0f} m/s"
              f"  sem={t['sem']:.3f}", flush=True)

    db = load_manifest().sort_values('t0').reset_index(drop=True)
    db = db[db['fn'].map(os.path.exists)]
    mine = db.iloc[task::ntask]
    print(f'\ntask {task}/{ntask}: {len(mine)} files of {len(db)}', flush=True)

    from scipy.signal import butter
    keep = freeze_channels(db)
    z = (np.arange(CH_LO, CH_HI)[keep] - CH_LO) * DX
    nfft = 30000                       # nominal 60 s at 500 Hz
    f = np.fft.rfftfreq(nfft, 1.0 / FS)
    # complex64: the stack is 14 x 700 x 15001, which is 2.4 GB at complex128 and
    # 1.2 GB here, and the per-file multiply is the inner loop. Max phase is
    # p*z*f ~ 2.5e-4 * 715 * 250 = 45 cycles, and float32 carries ~7 digits, so
    # the shift is accurate to ~1e-5 cycles -- far below the 2 ms timing target.
    phase = np.array([np.exp(2j * np.pi * np.outer(p * z, f)) for p in slow],
                     dtype=np.complex64)
    sos = butter(4, list(BAND), btype='band', fs=FS, output='sos')
    print(f'frozen channel set: {int(keep.sum())} of {CH_HI-CH_LO}; '
          f'phase stack {phase.nbytes/1e9:.2f} GB', flush=True)

    acc = {}          # (day, tag) -> dict of running stats
    nread = 0
    t_start = _time.time()
    for c, (_, r) in enumerate(mine.iterrows()):
        B = beams(r['fn'], slow, phase=phase, keep=keep, sos=sos)
        if B is None:
            continue
        nread += 1
        day = pd.Timestamp(r['t0']).strftime('%Y%m%d')
        for k, t in enumerate(T):
            cc = corr(B[k], t['w'])
            ca = corr(B[k], t['w'][::-1])
            if cc is None:
                continue
            key = (day, t['tag'])
            a = acc.setdefault(key, dict(n=0, s=0.0, s2=0.0, acmax=-9.0,
                                         peaks=[], times=[]))
            a['n'] += cc.size
            a['s'] += float(cc.sum()); a['s2'] += float((cc ** 2).sum())
            if ca is not None:
                a['acmax'] = max(a['acmax'], float(np.max(np.abs(ca))))
            m = int(np.argmax(cc))
            a['peaks'].append(float(cc[m]))
            a['times'].append(str(r['t0']) + f'+{m/FS:.3f}')
        if c % 250 == 0:
            el = _time.time() - t_start
            rate = el / max(nread, 1)
            eta = rate * (len(mine) - c) / 3600.0
            print(f'  {c}/{len(mine)}  read {nread}  {rate:.2f} s/file  '
                  f'ETA {eta:.1f} h  keys {len(acc)}', flush=True)

    rows = []
    for (day, tag), a in acc.items():
        pk = np.array(a['peaks'])
        o = np.argsort(pk)[::-1][:KEEP_TOP]
        mu = a['s'] / max(a['n'], 1)
        sd = np.sqrt(max(a['s2'] / max(a['n'], 1) - mu ** 2, 0.0))
        rows.append(dict(day=day, tag=tag, nsamp=a['n'], mean=mu, std=sd,
                         acmax=a['acmax'], nfile=len(pk),
                         top=pk[o].astype(np.float32),
                         toptime=np.array(a['times'])[o]))
    if not rows:
        print('nothing accumulated'); return
    out = os.path.join(OUT, f'full_task{task:03d}.npz')
    np.savez_compressed(
        out,
        day=np.array([r['day'] for r in rows]),
        tag=np.array([r['tag'] for r in rows]),
        nsamp=np.array([r['nsamp'] for r in rows]),
        mean=np.array([r['mean'] for r in rows]),
        std=np.array([r['std'] for r in rows]),
        acmax=np.array([r['acmax'] for r in rows]),
        nfile=np.array([r['nfile'] for r in rows]),
        top=np.array([np.pad(r['top'], (0, KEEP_TOP - r['top'].size),
                             constant_values=np.nan) for r in rows]),
        toptime=np.array([np.pad(r['toptime'], (0, KEEP_TOP - r['toptime'].size),
                                 constant_values='') for r in rows]))
    print(f'\nwrote {out}: {len(rows)} (day,template) cells, {nread} files read')


if __name__ == '__main__':
    main()
