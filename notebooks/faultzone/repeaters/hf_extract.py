"""
Stage HF-1: re-extract the moveout-test events at the NATIVE 500 Hz.

WHY. Every DAS measurement in this project reads cache_all/, which extract_all.py
wrote with `desampling=True` AND `FMAX = 40.0`. The cache is therefore band-limited
twice over: decimated to 100 Hz (Nyquist 50) and low-passed at 40 Hz. The files
actually carry 500 Hz (OutputDataRate in the HDF5 header; the manifest's "fs 10000"
is the interrogator PULSE rate, not the output rate -- that mislabel is what hid
this). Nyquist is 250 Hz, so ~80% of the available band was discarded before any
analysis ran.

SCOPE. Only the 73 events appearing in moveout_test.csv (10 HRSN-confirmed repeater
pairs + 40 random control pairs), not all 206. That is 3.3 GB rather than 9.5 GB and
keeps the null intact, which is the part that matters.

Writes to cache_hf/, leaving cache_all/ untouched so published numbers stay
reproducible.

Run:  sbatch hf_moveout_job.sh
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from g1_coda_snr import load_manifest                      # noqa: E402

for p in ['/home/groups/edunham/nberrios/safod_das/DAS-utilities/python',
          '/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python']:
    if os.path.exists(p) and p not in sys.path:
        sys.path.insert(0, p)
import DASutils                                             # noqa: E402

CACHE = os.path.join(HERE, 'cache_hf')
os.makedirs(CACHE, exist_ok=True)

PRE_S, POST_S = 5.0, 20.0        # identical window to cache_all, so the only
FMIN, FMAX = 1.0, 240.0          # difference is bandwidth. Nyquist at 500 Hz = 250.


def extract(db, origin, tag):
    f = os.path.join(CACHE, f'{tag}.npz')
    if os.path.exists(f):
        return 'cached'
    a = origin - pd.Timedelta(seconds=PRE_S)
    b = origin + pd.Timedelta(seconds=POST_S)
    sel = db[(db['t0'] < b) & (db['t1'] > a)].copy()
    if sel.empty:
        return 'no files'
    sel = sel[sel['fn'].map(os.path.exists)]
    if sel.empty:
        return 'files missing'
    sel = sel.sort_values('t0')

    acc = hits = None
    fs = dt = npts = None
    for _, r in sel.iterrows():
        try:
            D, info = DASutils.readFile_HDF(
                [r['fn']], FMIN, FMAX, verbose=0, preproc=True, diff=True,
                taper=False, desampling=False, nChbuffer=900, system='OptaSense')
        except Exception as ex:
            return f'read fail: {str(ex)[:40]}'
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
        miss = int(np.sum(hits == 0)) if hits is not None else -1
        return f'gap ({miss} samples)'
    np.savez_compressed(f, X=(acc / hits[None, :]).astype(np.float32), fs=fs,
                        pre_s=PRE_S, post_s=POST_S)
    return f'ok fs={fs:g}'


def main():
    task = int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
    # SLURM_ARRAY_TASK_COUNT is not set by every Slurm build; fall back to the
    # MIN/MAX pair. Getting this wrong silently skips events rather than failing.
    ntask = os.environ.get('SLURM_ARRAY_TASK_COUNT')
    if ntask:
        ntask = int(ntask)
    else:
        lo = int(os.environ.get('SLURM_ARRAY_TASK_MIN', 0))
        hi = int(os.environ.get('SLURM_ARRAY_TASK_MAX', 0))
        ntask = hi - lo + 1

    mt = pd.read_csv(os.path.join(HERE, 'moveout_test.csv'))
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    ev['time'] = pd.to_datetime(ev['time'], utc=True, format='mixed')
    need = sorted(set(mt.i) | set(mt.j))
    print(f'{len(mt)} pairs -> {len(need)} unique events; '
          f'task {task}/{ntask}', flush=True)

    db = load_manifest()
    mine = [k for n, k in enumerate(need) if n % ntask == task]
    print(f'this task: {len(mine)} events\n', flush=True)

    for k in mine:
        e = ev.loc[k]
        st = extract(db, e['time'], e['tag'])
        print(f'  [{k:3d}] {e["tag"]}  M{e["mag"]:.2f}  {st}', flush=True)


if __name__ == '__main__':
    main()
