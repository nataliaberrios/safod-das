# SAFOD DAS Workflow Notes

This repository contains the working Sherlock copy of the SAFOD DAS cross-correlation workflow.

## Main workflow files

- `stack_daily.py`
- `submit_next.sh`
- `run_daily_array.sbatch`
- `find_dates.py`
- `notebooks/13_stacks_n_seis.ipynb`

## Current Sherlock repo location

Working git-controlled copy:

```bash
/home/groups/ettore88/nberrios/safod_das_git
```

Older non-git copy kept for safety:

```bash
/home/groups/edunham/nberrios/safod_das
```

## What was fixed

### 1. Date indexing no longer assumes 1440 files per day

`stack_daily.py` used to assume:

```python
day_idx -> file block day_idx * 1440
```

That failed once there were gaps in the recording sequence. It now:

- parses `startTime_dt`
- groups rows by actual calendar date
- processes all valid files for the selected date

`submit_next.sh` was updated to use the same date-based logic.

### 2. SAFOD file-path normalization

Some CSV rows referenced file paths under:

```bash
/oak/stanford/groups/ettore88/data/SAFODAS1-harddrive-transfer/...
```

but the actual files exist under:

```bash
/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/...
```

`stack_daily.py` now normalizes those paths before attempting to read the HDF5 files.

### 3. Batch script path update

`run_daily_array.sbatch` must point to the git-controlled repo location and its `DAS-utilities/python` directory.

## Running the daily CC workflow

Submit missing dates in chunks:

```bash
bash submit_next.sh 20
```

Submit a comparison run with temporal normalization and a second bandpass:

```bash
OUTPUT_VERSION=tn_bp2 USE_TEMPORAL_NORMALIZATION=true USE_SECOND_BANDPASS=true bash submit_next.sh 20
```

Test one day explicitly:

```bash
sbatch --array=64 run_daily_array.sbatch
```

Inspect logs:

```bash
ls -lt logs | head
cat logs/safod_daily_<JOBID>_<ARRAYIDX>.out
cat logs/safod_daily_<JOBID>_<ARRAYIDX>.err
```

## Baseline vs comparison processing

The main configuration options are at the top of `stack_daily.py`.

### Baseline run

```python
OUTPUT_VERSION = "base"
USE_TEMPORAL_NORMALIZATION = False
USE_SECOND_BANDPASS = False
```

This writes to:

```bash
/oak/stanford/groups/ettore88/nberrios/
```

### Comparison run: temporal normalization + second bandpass

```python
OUTPUT_VERSION = "tn_bp2"
USE_TEMPORAL_NORMALIZATION = True
USE_SECOND_BANDPASS = True
```

This writes to:

```bash
/oak/stanford/groups/ettore88/nberrios/tn_bp2/
```

So the original baseline outputs are not overwritten.

## Weekly CC + seismicity notebook

The main notebook is:

```bash
notebooks/13_stacks_n_seis.ipynb
```

It now supports:

- weekly CC plots with seismicity histograms
- one figure per week
- a chunking mode switch:

```python
CHUNK_MODE = 'calendar_strict'   # or 'available_files'
```

Recommended use:

- `calendar_strict` for cleaner figures
- `available_files` for more permissive analysis when gaps would discard too much data

## JupyterLab access

If JupyterLab only exposes the Edunham group directory tree, a symlink can make the Ettore repo visible:

```bash
ln -s /home/groups/ettore88/nberrios/safod_das_git /home/groups/edunham/nberrios/safod_das_git_link
```

Then open the notebook through the visible symlink path in JupyterLab.

## Git / GitHub

GitHub repo:

```text
https://github.com/nataliaberrios/safod-das
```

Typical workflow:

```bash
git status
git diff
git add <files>
git commit -m "Describe the change"
git push
```
