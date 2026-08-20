#!/bin/bash
#SBATCH --job-name=lellouch_ordered
#SBATCH --partition=serc
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-1
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_nulls/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_nulls/logs/%x-%A_%a.err
#
# The OBSERVED arm of the input-level null test.  Its matched null arm is
# ambient_lellouch2019_matched_nulls.sbatch -- every setting below must stay
# identical there, which is the whole point of the pair.
#
# UPDATED 2026-08-19.  This used to sweep --fk-mode {unfiltered,
# downgoing_2p5_4p5} and pass --alignment-velocity data.  Both options were
# removed from ambient_lellouch2019_exact_stack.py, so the file had been failing
# at argparse rather than running.  The F-K axis is not restored: the QC workflow
# rejects the 2.5-4.5 km/s fan because it fails the pre-filter channel-scramble
# gate, so sweeping it is not a comparison worth making.
#
# The array axis is now the common-mode branch, which is the comparison that
# matters -- config 0 (no removal) has a moveout statistic that is ~97 %
# pedestal, corr(v, score) = +0.976; config 3 (median removal) is the only one
# whose statistic actually measures moveout, corr = -0.381.

set -euo pipefail
WORK=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUT="$WORK/ambient_transfer/lellouch2019_exact_stack_nulls"
mkdir -p "$OUT/logs"
cd "$WORK"

if [[ "$SLURM_ARRAY_TASK_ID" -eq 0 ]]; then
  BRANCH=()                 # config 0: no common-mode removal
else
  BRANCH=(--common-mode)    # config 3: per-sample median across channels
fi

# bash 4.2 with `set -u` errors on "${BRANCH[@]}" when BRANCH is empty, so guard
# the expansion on the element count rather than expanding unconditionally.
/home/users/nberrios/miniconda3/envs/das/bin/python -u \
  ambient_lellouch2019_exact_stack.py \
    --date 2024-12-20 \
    --start 0 \
    --nfiles 300 \
    --source-channel 23 \
    --ram-seconds 0.1 \
    --null-method ordered \
    --realization 0 \
    --spectral-mode cross_correlation \
    ${BRANCH[@]+"${BRANCH[@]}"} \
    --output-dir "$OUT"
