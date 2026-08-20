#!/bin/bash
#SBATCH --job-name=lellouch_null
#SBATCH --partition=serc
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-79%12
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_nulls/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_nulls/logs/%x-%A_%a.err
#
# The NULL arm matching ambient_lellouch2019_matched_ordered.sbatch.  These are
# INPUT-LEVEL nulls -- the surrogate is built before the operator runs, so the
# operator is inside the null and cannot inflate the p-value.  That is the gate
# the F-K fan failed, and it is the gate that counts.
#
# UPDATED 2026-08-19.  Three things here were stale and made every task fail:
#   --fk-mode              removed from exact_stack (and the fan is QC-rejected)
#   --alignment-velocity   removed from exact_stack
#   circular_time_shift    never a valid --null-method for exact_stack; it exists
#                          only in ambient_fk_full_pipeline_null_v2.py.
#                          VALID_NULLS = (ordered, white_noise, channel_permutation)
# The array is now 2 branches x 2 methods x 20 realizations = 80 tasks.

set -euo pipefail

WORK=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUT="$WORK/ambient_transfer/lellouch2019_exact_stack_nulls"
mkdir -p "$OUT/logs"
cd "$WORK"

REALIZATIONS=20
BRANCH_BLOCK=$((SLURM_ARRAY_TASK_ID / 40))
INNER=$((SLURM_ARRAY_TASK_ID % 40))
METHOD_BLOCK=$((INNER / REALIZATIONS))
REALIZATION=$((INNER % REALIZATIONS))

case "$BRANCH_BLOCK" in
  0) BRANCH=() ;;                 # config 0: no common-mode removal
  1) BRANCH=(--common-mode) ;;    # config 3: per-sample median across channels
  *) exit 2 ;;
esac

case "$METHOD_BLOCK" in
  0) NULL_METHOD=white_noise ;;
  1) NULL_METHOD=channel_permutation ;;
  *) exit 2 ;;
esac

# bash 4.2 with `set -u` errors on "${BRANCH[@]}" when BRANCH is empty.
/home/users/nberrios/miniconda3/envs/das/bin/python -u \
  ambient_lellouch2019_exact_stack.py \
    --date 2024-12-20 \
    --start 0 \
    --nfiles 300 \
    --source-channel 23 \
    --ram-seconds 0.1 \
    --null-method "$NULL_METHOD" \
    --realization "$REALIZATION" \
    --spectral-mode cross_correlation \
    ${BRANCH[@]+"${BRANCH[@]}"} \
    --output-dir "$OUT"
