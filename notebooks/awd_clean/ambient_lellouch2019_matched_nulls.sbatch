#!/bin/bash
#SBATCH --job-name=lellouch_null
#SBATCH --partition=serc
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-119%12
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_nulls/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_nulls/logs/%x-%A_%a.err

set -euo pipefail

WORK=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUT="$WORK/ambient_transfer/lellouch2019_exact_stack_nulls"
mkdir -p "$OUT/logs"
cd "$WORK"

MODE_BLOCK=$((SLURM_ARRAY_TASK_ID / 60))
INNER=$((SLURM_ARRAY_TASK_ID % 60))
METHOD_BLOCK=$((INNER / 20))
REALIZATION=$((INNER % 20))

case "$MODE_BLOCK" in
  0) FK_MODE=unfiltered ;;
  1) FK_MODE=downgoing_2p5_4p5 ;;
  *) exit 2 ;;
esac

case "$METHOD_BLOCK" in
  0) NULL_METHOD=white_noise ;;
  1) NULL_METHOD=channel_permutation ;;
  2) NULL_METHOD=circular_time_shift ;;
  *) exit 2 ;;
esac

/home/users/nberrios/miniconda3/envs/das/bin/python -u   ambient_lellouch2019_exact_stack.py   --date 2024-12-20   --start 0   --nfiles 300   --null-method "$NULL_METHOD"   --realization "$REALIZATION"   --fk-mode "$FK_MODE"   --spectral-mode cross_correlation   --alignment-velocity data   --output-dir "$OUT"
