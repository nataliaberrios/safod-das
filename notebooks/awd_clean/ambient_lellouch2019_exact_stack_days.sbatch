#!/bin/bash
#SBATCH --job-name=lel7c_days
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --array=0-119%28
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_days/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_days/logs/%x-%A_%a.err

# Extend the paper-faithful Figure 7c baseline (configuration 0 only) to the
# five days other than 2024-12-20.  Task = day x hour, 5 x 24 = 120.
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
# 2024-11-30 IS included here via --continuous-prefix, which analyses its leading
# 21.3 h contiguous block (1278 of 1441 rows) rather than dropping the day.  The
# ordering is arbitrary: the census that formerly ranked these days is withdrawn
# (see FIG7C_MULTIDAY_RESULT.md), so this is an unselected set.
DATES=(2024-06-17 2024-06-26 2025-02-24 2024-05-11 2024-11-30)

# FIG7C_CONFIG selects the branch.  Default 0 reproduces the previous behaviour
# exactly.  Config 3 subtracts the instantaneous median across all acquisition
# channels before RAM normalisation.
#
# WHY CONFIG 3 MATTERS.  The config-0 moveout statistic is ~97 % pedestal: its
# gates slide toward the dominant zero-lag lobe as trial velocity rises, giving
# corr(trial velocity, score) = +0.976 and an argmax pinned to the grid ceiling.
# Common-mode removal eliminates that -- corr = -0.381, curve flat at ~0.95 -- so
# config 3 is the only branch whose statistic actually measures moveout.  Config 3
# has never been stacked beyond a single day; every multi-day stack to date,
# including the coherent 96 h one, used config 0.
CONFIG=${FIG7C_CONFIG:-0}
case "$CONFIG" in
  0) LABEL=paper_baseline; COMMON_MODE=0 ;;
  3) LABEL=common_mode;    COMMON_MODE=1 ;;
  *) echo "unsupported FIG7C_CONFIG=$CONFIG (use 0 or 3)" >&2; exit 2 ;;
esac

# 2024-11-30 needs its leading contiguous block; the other days are exactly
# continuous and the flag is a no-op for them, so it is safe to pass always.
PREFIX_FLAG=(--continuous-prefix)

DAY=$(( SLURM_ARRAY_TASK_ID / 24 ))
HOUR=$(( SLURM_ARRAY_TASK_ID % 24 ))
DATE=${DATES[$DAY]}
START=$(( HOUR * 60 ))

# 2024-11-30's contiguous block is 1278 files, so only hours 0-20 (files 0-1259)
# lie wholly inside it.  Hours 21-23 would write ragged partial chunks that the
# aggregate's uniform-size expectation cannot combine, so skip them cleanly.
if [[ "$DATE" == "2024-11-30" && "$HOUR" -ge 21 ]]; then
  echo "skip: $DATE hour $HOUR lies beyond the 1278-file contiguous prefix"
  exit 0
fi

echo "day=$DAY date=$DATE hour=$HOUR start=$START config=$CONFIG ($LABEL)"
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
  "${PREFIX_FLAG[@]}"
)
if [[ "$COMMON_MODE" -eq 1 ]]; then
  ARGS+=(--common-mode)
fi
srun "$PY" "${ARGS[@]}"
