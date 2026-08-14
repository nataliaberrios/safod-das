#!/bin/bash
#SBATCH --job-name=lel7c_agg
#SBATCH --partition=serc
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=0-6
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack/logs/%x-%A_%a.err

set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

WORK=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUT="$WORK/ambient_transfer/lellouch2019_exact_stack"
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
mkdir -p "$OUT/logs"
cd "$WORK"

SOURCE=23
RAM=0.1
SPECTRAL=cross_correlation
NULL_METHOD=ordered
REALIZATION=0
EXTRA=()

case "$SLURM_ARRAY_TASK_ID" in
  0) LABEL=paper_baseline ;;
  1) LABEL=source_ch0; SOURCE=0 ;;
  2) LABEL=ram_5s; RAM=5.0 ;;
  3) LABEL=common_mode; EXTRA=(--common-mode) ;;
  4) LABEL=equation6_stabilized; SPECTRAL=source_power_stabilized ;;
  5) LABEL=white_noise; NULL_METHOD=white_noise ;;
  6) LABEL=channel_permutation; NULL_METHOD=channel_permutation; REALIZATION=0 ;;
  *) echo "unexpected configuration $SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac

echo "aggregate label=$LABEL source=$SOURCE ram=$RAM spectral=$SPECTRAL null=$NULL_METHOD"
srun "$PY" -u ambient_lellouch2019_exact_stack.py \
  --action aggregate \
  --date 2024-12-20 \
  --nfiles 60 \
  --total-files 1440 \
  --source-channel "$SOURCE" \
  --ram-seconds "$RAM" \
  --spectral-mode "$SPECTRAL" \
  --null-method "$NULL_METHOD" \
  --realization "$REALIZATION" \
  --null-count 10000 \
  --output-dir "$OUT" \
  "${EXTRA[@]}"
