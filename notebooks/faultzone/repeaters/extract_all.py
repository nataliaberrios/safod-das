"""
Stage 2: extract every covered event to the waveform cache. SLURM array task.

Extraction is the expensive step -- each event assembles a time-aligned window
across several HDF files. Splitting it over array tasks turns ~15 hours serial
into about an hour wall-clock.

Already-cached events are skipped, so the array is safe to re-run after a
partial failure, and a task that dies costs only its own slice.

Run as:  sbatch --array=0-19  ...  python -u extract_all.py
Each task handles events where (index % n_tasks) == task_id, so the work
interleaves rather than blocking on one dense time period.
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

CACHE = os.path.join(HERE, 'cache_all')
os.makedirs(CACHE, exist_ok=True)

PRE_S, POST_S = 5.0, 20.0        # 25 s: direct arrivals plus coda for CWI later
FMIN, FMAX = 5.0, 40.0


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
                taper=False, desampling=True, nChbuffer=900, system='OptaSense')
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
    return 'ok'


def main():
    task = int(os.environ.get('SLURM_ARRAY_TASK_ID', '0'))
    ntask = int(os.environ.get('N_TASKS', '20'))
    ev = pd.read_csv(os.path.join(HERE, 'all_events.csv'))
    ev['time'] = pd.to_datetime(ev['time'], utc=True, format='mixed')
    mine = ev[ev['idx'] % ntask == task].reset_index(drop=True)
    print(f'task {task}/{ntask}: {len(mine)} events', flush=True)

    db = load_manifest()
    print(f'manifest {len(db)} rows', flush=True)

    counts = {}
    for _, e in mine.iterrows():
        st = extract(db, e['time'], e['tag'])
        counts[st.split(' ')[0]] = counts.get(st.split(' ')[0], 0) + 1
        print(f'  {e["tag"]}  M{e["mag"]:.2f}  {st}', flush=True)
    print(f'task {task} done: {counts}', flush=True)


if __name__ == '__main__':
    main()
