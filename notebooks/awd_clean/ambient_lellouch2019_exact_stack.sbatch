#!/bin/bash
#SBATCH --job-name=lellouch_exact
#SBATCH --partition=serc
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-3
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack/logs/%x-%A_%a.err

set -euo pipefail

WORK=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUT="$WORK/ambient_transfer/lellouch2019_exact_stack"
mkdir -p "$OUT/logs"
cd "$WORK"

EXTRA=()
case "$SLURM_ARRAY_TASK_ID" in
  0)
    EXTRA=(--fk-mode unfiltered --spectral-mode cross_correlation)
    ;;
  1)
    EXTRA=(--fk-mode unfiltered --spectral-mode cross_correlation --common-mode)
    ;;
  2)
    EXTRA=(--fk-mode unfiltered --spectral-mode source_power_stabilized --waterlevel 0.001)
    ;;
  3)
    EXTRA=(--fk-mode downgoing_direction_only --spectral-mode cross_correlation)
    ;;
  *)
    echo "Unexpected array task $SLURM_ARRAY_TASK_ID" >&2
    exit 2
    ;;
esac

/home/users/nberrios/miniconda3/envs/das/bin/python -u   ambient_lellouch2019_exact_stack.py   --date 2024-12-20   --start 0   --nfiles 1440   --null-method ordered   --alignment-velocity data   --output-dir "$OUT"   "${EXTRA[@]}"
