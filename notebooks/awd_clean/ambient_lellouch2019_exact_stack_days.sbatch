#!/bin/bash
#SBATCH --job-name=lel7c_days
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --array=0-95%28
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_days/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_days/logs/%x-%A_%a.err

# Extend the paper-faithful Figure 7c baseline (configuration 0 only) to the
# four EXACTLY CONTINUOUS days other than 2024-12-20.  Task = day x hour, 4 x 24 = 96.
#
# WHY.  The faithful operator has only ever been run on 2024-12-20.  The raw
# energy census (ambient_fk_energy_census.py) shows the 2.5-4 km/s fan share
# varies by day -- 0.3% on 2024-05-11 up to 4.9% on 2024-11-30 -- so the day
# that was tested is not the most favourable one available.  Days are ordered
# below by descending census fan share so the most promising day reports first.
#
# Writes to its OWN output directory so the products cited by
# AMBIENT_FIG7C_STATUS.md (jobs 38988141 / 38993456, replicated by 38995917)
# are not overwritten.
#
# Chunks store summed cross spectra, so per-day aggregates and an exact
# multi-day stack both follow from these files without touching raw data again.

set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

WORK=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUT="$WORK/ambient_transfer/lellouch2019_exact_stack_days"
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
mkdir -p "$OUT/logs"
cd "$WORK"

# Only days whose manifest is exactly continuous at 60.000 s are usable: the
# operator streams across file boundaries and its guard rightly refuses to splice
# a discontinuous day.  Measured over the manifest, 2024-11-30, 2024-10-28 and
# 2025-03-04 each carry two timing anomalies (longest continuous runs 1278, 1191
# and 826 files) and are therefore excluded here rather than silently truncated.
# Remaining days are ordered by census 2.5-4 km/s fan share, richest first.
DATES=(2024-06-17 2024-06-26 2025-02-24 2024-05-11)

DAY=$(( SLURM_ARRAY_TASK_ID / 24 ))
HOUR=$(( SLURM_ARRAY_TASK_ID % 24 ))
DATE=${DATES[$DAY]}
START=$(( HOUR * 60 ))

echo "day=$DAY date=$DATE hour=$HOUR start=$START (configuration 0, paper baseline)"
ARGS=(
  -u ambient_lellouch2019_exact_stack.py
  --action chunk
  --date "$DATE"
  --start "$START"
  --nfiles 60
  --block-files 10
  --source-channel 23
  --ram-seconds 0.1
  --spectral-mode cross_correlation
  --null-method ordered
  --realization 0
  --output-dir "$OUT"
)
srun "$PY" "${ARGS[@]}"
