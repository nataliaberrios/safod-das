# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Orientation

SAFOD borehole DAS processing, run on Stanford's Sherlock cluster. The repo holds
two distinct bodies of work:

| | Ambient-noise cross-correlation | AWD active-source coupling |
|---|---|---|
| Status | Tracked in git; stable | **Current active work**; untracked |
| Data | `SAFOD_2024_2025.csv` manifest → HDF5 | `/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/` |
| Reader | `DASutils.readFile_HDF` | `readFile_protobuf` for `Nano/*.pb`, `readFile_HDF` for `SAFOD-Deep-*.h5` |
| Entry points | `stack_daily.py`, `sanity/`, `event_cc_20s_exact_snippet.py` | `notebooks/*_fml_coupling.ipynb` + scripts lifted from them |
| SLURM | `-p gpu`, CUDA + PyTorch modules | CPU only, `conda activate das` |

`git ls-files` shows only the first column. The active AWD notebooks and their
job scripts are untracked, so do not assume the tracked file list reflects what
is being worked on.

### Which directory is the real one

`/home/groups/ettore88/nberrios/safod_das_git` is canonical.
`/home/groups/edunham/nberrios/safod_das_git_link` is a symlink to it — that path
exists because JupyterLab only exposes the Edunham group tree.

`/home/groups/ettore88/nberrios/safod-das` is a **second clone of the same
GitHub repo**. Its git history has been pushed to `origin`, but its working tree
has been idle since May 2026. Work in `safod_das_git`.

### Check `git fetch` before concluding a file is untracked

This checkout drifts behind `origin/main` — it has been 20+ commits behind, which
makes committed files (`cc_tools.py`, parts of `sanity/`) appear untracked and
leaves stale entries staged in the index for months. Run `git fetch` and compare
against `origin/main`, not local `HEAD`, before drawing conclusions about what is
or isn't in version control.

## Environment: two conventions, don't mix them

Nothing here is pip-installed. `DASutils` is imported off `PYTHONPATH` (or an
in-script `sys.path` insert) from a `DAS-utilities/python` checkout, which is
itself untracked and not a submodule.

### CC pipeline (GPU jobs)

Every `.sbatch` in the root and in `sanity/` repeats this block verbatim:

```bash
module --force purge
module load devel math
module load cuda/11.7.1
module load py-pytorch/2.0.0_py39
module unload py-scipy          # collides with the das env's scipy
export PYTHONNOUSERSITE=1
export PYTHONPATH=<das site-packages>:<py-pytorch site-packages>:<repo>/DAS-utilities/python
```

Invoke the interpreter by absolute path,
`/home/users/nberrios/miniconda3/envs/das/bin/python` — never bare `python`, which
resolves to the module's Python and will not find the env's packages.

### AWD / notebook work (CPU jobs)

Different and simpler: `ml gcc/12.4.0` then `conda activate das`. These scripts
bootstrap their own imports rather than relying on `PYTHONPATH`, using a fallback
list that tries both `DAS-utilities` checkouts (`notebooks/paired_stack_job.py`):

```python
SEARCH_DIRS = [
    '/home/groups/edunham/nberrios/safod_das/DAS-utilities/python',
    '/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python',
]
```

The matching Jupyter kernel is named `das`. For an interactive `sh_dev` shell,
`enable_das.sh` sets up the equivalent.

### Three `DAS-utilities` paths are live at once

`run_daily_array.sbatch:19` uses this repo's copy;
`run_event_cc_exact_snippet.sbatch:19` and `enable_das.sh` use the Edunham copy;
`sanity/run_sanity.sbatch:53` defaults `REPO_DIR` to the idle `safod-das` clone.
All three exist on disk, so pointing at the wrong one produces a completed job
with different code rather than an error. Set `REPO_DIR` explicitly when running
`sanity/`.

## CC pipeline architecture

One shared shape across `stack_daily.py`, `event_cc_20s_exact_snippet.py`, and
`sanity/sanity_cc.py`:

```
manifest CSV on $OAK  (whitespace-separated: file, startTime, endTime, nSamples)
  → filter rows, then rewrite paths through normalize_file_path()
  → DASutils.readFile_HDF(..., system="OptaSense")
  → slice channels [ch_start:ch_end]
  → 30 s windows at 50% overlap
  → subtract per-sample median across channels
  → DASutils.bandpass2D_c 5-20 Hz
  → optional running-absolute-mean temporal normalization
  → computeCC against one virtual-source channel  (torch FFT cross-correlation)
  → sum windows, divide by nstack
  → .npz on $OAK
```

`normalize_file_path` exists because some manifest rows point at a path that
predates a data move; it rewrites `.../data/SAFODAS1-harddrive-transfer` to
`.../data/SAFOD/SAFODAS1-harddrive-transfer` when the original does not resolve.
It is duplicated verbatim in all three scripts (`stack_daily.py:17`,
`event_cc_20s_exact_snippet.py:38`, `sanity/sanity_cc.py:103`).

### Array index → date, not a fixed file block

`SLURM_ARRAY_TASK_ID` indexes into the sorted list of *unique calendar dates*
present in the filtered manifest (`stack_daily.py:68`), not a fixed 1440-files-
per-day block. `submit_next.sh:24` rebuilds that same list independently. If the
manifest filter changes in one place and not the other, every day index shifts
and jobs silently process the wrong dates.

### Output layout and idempotency

Results always go to `$OAK`, never into the repo. `OUTPUT_VERSION` controls the
directory: `base` writes to `/oak/stanford/groups/ettore88/nberrios/`, anything
else writes to a `<version>/` subdirectory, so comparison runs cannot clobber the
baseline.

`stack_daily.py:89` exits early when both `<date>_day.npz` and `<date>_night.npz`
already exist, and `submit_next.sh:38` uses the same existence test to choose the
next N missing indices. Re-running a version therefore means deleting its `.npz`
files first — otherwise every job no-ops.

### Two deliberate divergences — do not "fix" them

**Temporal normalization window.** `stack_daily.py:32` defaults `TN_WINDOW` to
10.0 s. `sanity/sanity_cc.py:74` hardcodes 0.1 s, and `sanity/README.md` records
that 10 s is roughly 100× too long for body-wave cross-correlation. New scripts
should follow the sanity value.

**Manifest filter.** `stack_daily.py:59` keeps `nSamples > 0`, which includes the
~20 s event-triggered files. `sanity/sanity_cc.py:147` keeps `nSamples == 30000`,
continuous files only. The sanity pipeline is a strict Lellouch et al. 2019
reproduction and the event triggers do not belong in it.

### `cc_tools.py` defines everything twice

`nextpow2`, `torch_xcorr`, and `computeCC` each appear twice in the file. The
second set (from line 123) is what any importer actually gets: it adds the
`whitening_params` argument and shifts lag centering from `npts` to `(npts-1)`.
Editing the first set has no effect.

## AWD active-source pipeline

June 2026 active-source survey comparing two fibers in the same borehole. Survey
window 2026-06-16 23:45 → 2026-06-17 23:45 UTC.

- Weight-drop times come from `notebooks/p26.cc9.txt` (989 drops, `UTC_Date`
  column), picked by cross-correlation against a nearby seismometer.
- Two fibers: `ActiveJune2026/Nano/*.pb` via `readFile_protobuf`, and
  `ActiveJune2026/01_--_recording_.../SAFOD-Deep-*.h5` via `readFile_HDF`.
- Files are grouped 6 per epoch; each drop is windowed `PRE_S=0.5` to
  `POST_S=3.0` seconds.
- Outputs land in `notebooks/figures/awd_2026/`.

The two file types encode timestamps incompatibly, and both parsers are
reimplemented in several places:

```python
# .pb  — underscore-split, dotted time
os.path.basename(f).split('_')          # parts[1]=date, parts[2]=HH.MM.SS
# .h5  — regex
re.search(r'_(\d{4}-\d{2}-\d{2}T\d{6}Z)\.h5$', fname)
```

### Working pattern

Prototype in a `*_fml_coupling.ipynb` notebook, then lift the relevant cells into
a standalone script so it can run as a batch job without kernel state —
`notebooks/paired_stack_job.py` documents this in its header comment. Submit via
`notebooks/deep_stack_job.sh` or `notebooks/submit_paired_stacks.sh`.

## Configuration is entirely environment variables

Two disjoint namespaces, both delivered through `sbatch --export=ALL,...`:

| Pipeline | Variables |
|---|---|
| Main CC | `OUTPUT_VERSION`, `USE_TEMPORAL_NORMALIZATION`, `USE_SECOND_BANDPASS`, `TN_WINDOW` |
| `sanity/` | `SANITY_DATE`, `SANITY_DATES`, `SANITY_HOURS`, `SANITY_CH_START`, `SANITY_CH_END`, `SANITY_SOURCE_CH`, `SANITY_WHITEN`, `SANITY_BADCH_K`, `SANITY_OUT`, `STAGE`, `REPO_DIR` |

## Commands

Submit the next N days that have no output yet:

```bash
bash submit_next.sh 20
OUTPUT_VERSION=tn_bp2 USE_TEMPORAL_NORMALIZATION=true USE_SECOND_BANDPASS=true bash submit_next.sh 20
```

One explicit day index:

```bash
sbatch --array=64 run_daily_array.sbatch
```

Sanity (Lellouch reproduction) — `STAGE` selects the script:

```bash
sbatch --export=ALL,STAGE=preflight,REPO_DIR=/home/groups/ettore88/nberrios/safod_das_git run_sanity.sbatch
# STAGE = preflight | cc | cc_segy | cc_segy_cmn | plot | all
```

Single-event CC (positional args: UTC time, then label):

```bash
python event_cc_20s_exact_snippet.py 2024-05-25T04:56:28.150000Z EQ_20240525_045628
python plot_event_cc_clean.py --event-file /oak/.../EQ_20240525_045628_event20s.npz
```

AWD stacks:

```bash
sbatch notebooks/deep_stack_job.sh
sbatch notebooks/submit_paired_stacks.sh
```

Job inspection:

```bash
squeue --me
sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem
ls -lt logs | head          # logs/%x_%A_%a.out and .err
```

### There is no test suite

`import_test.py` and `synthetic_cc_test.py` are ad-hoc scripts, not a harness.
The real cheap check before committing hours of cluster time is
`STAGE=preflight`, which runs on one hour of data and produces an F-K spectrum
answering whether the wavefield contains body-wave energy at all. If it doesn't,
no amount of CC tuning will recover a Green's function — `sanity/README.md`
describes the decision point in detail.

## Git notes for this repo

- Never `git add -A`. Tens of GB of `.npz`, `.mseed`, and figure output sit
  untracked in the working tree. Stage by explicit path.
- Results belong on `$OAK`, not in git. `.npz` files reach multiple GB.
- Notebooks carry embedded outputs and get very large (one exceeds 50 MB). Clear
  outputs before committing.
- `git fetch` first — see the staleness note at the top.
- The repo root contains a self-referential symlink
  `safod_das_git -> /home/groups/ettore88/nberrios/safod_das_git`. Recursive
  `find` or `glob` from the root will loop.
- The system `git` on Sherlock login nodes is 1.8.3.1, which lacks `worktree`
  and `-C`. A modern build is at
  `/share/software/user/open/git/2.45.1/bin/git`.
