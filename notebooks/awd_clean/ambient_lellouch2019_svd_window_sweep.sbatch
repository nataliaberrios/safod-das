#!/bin/bash
#SBATCH --job-name=svdwin
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --array=0-95%28
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack/logs/%x-%A_%a.err

# WINDOWED coherent-subspace removal, 2024-12-20. 4 configs x 24 hours = 96.
#
# A whole-block rank-k projection was measured on 2026-08-14 NOT to suppress the
# pedestal: corr(trial velocity, score) stayed at +0.96..+0.98 for ranks 1/2/4/8
# against the untreated baseline's +0.976, while the per-sample median reached
# -0.381. A global truncation captures time-averaged modes; the pedestal comes
# from the INSTANTANEOUS common level. Windowing the projection approximates the
# instantaneous removal while still handling channel-dependent gain, which the
# median cannot -- and the median's own residual was measured to ACCUMULATE under
# stacking (corr -0.381 at 1 day, +0.951 at 4 days coherently stacked).
#
# Figure of merit: corr(trial velocity, score), the pedestal indicator. Target is
# near 0. Only then is a p-value on this statistic interpretable as moveout.

set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
WORK=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUT="$WORK/ambient_transfer/lellouch2019_exact_stack"
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
mkdir -p "$OUT/logs"; cd "$WORK"

CFG=$(( SLURM_ARRAY_TASK_ID / 24 ))
HOUR=$(( SLURM_ARRAY_TASK_ID % 24 ))
START=$(( HOUR * 60 ))

# rank, window seconds, and whether the median runs first
case "$CFG" in
  0) RANK=2; WIN=1.0;  MED=0 ;;   # projection alone, near-instantaneous
  1) RANK=2; WIN=5.0;  MED=1 ;;   # median then projection
  2) RANK=2; WIN=30.0; MED=1 ;;   # median then projection at the correlation window
  3) RANK=4; WIN=5.0;  MED=1 ;;
  *) exit 2 ;;
esac

ARGS=(-u ambient_lellouch2019_exact_stack.py --action chunk --date 2024-12-20
      --start "$START" --nfiles 60 --block-files 10 --source-channel 23
      --ram-seconds 0.1 --spectral-mode cross_correlation --null-method ordered
      --realization 0 --svd-rank "$RANK" --svd-window-s "$WIN" --output-dir "$OUT")
if [[ "$MED" -eq 1 ]]; then ARGS+=(--common-mode); fi

echo "cfg=$CFG rank=$RANK window=${WIN}s median=$MED hour=$HOUR"
srun "$PY" "${ARGS[@]}"
