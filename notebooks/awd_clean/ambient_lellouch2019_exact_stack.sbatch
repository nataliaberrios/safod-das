#!/bin/bash
#SBATCH --job-name=lel7c_exact
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --array=0-167%28
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack/logs/%x-%A_%a.err

# Seven predeclared configurations x 24 one-hour chunks.  Configuration 0 is
# the paper-faithful baseline.  Configurations 1--6 are labelled sensitivities
# or nulls and must never be silently pooled with the baseline.

set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

WORK=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUT="$WORK/ambient_transfer/lellouch2019_exact_stack"
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
mkdir -p "$OUT/logs"
cd "$WORK"

CONFIG=$((SLURM_ARRAY_TASK_ID / 24))
HOUR=$((SLURM_ARRAY_TASK_ID % 24))
START=$((HOUR * 60))

SOURCE=23
RAM=0.1
SPECTRAL=cross_correlation
NULL_METHOD=ordered
REALIZATION=0
COMMON_MODE=0

case "$CONFIG" in
  0) LABEL=paper_baseline ;;
  1) LABEL=source_ch0; SOURCE=0 ;;
  2) LABEL=ram_5s; RAM=5.0 ;;
  3) LABEL=common_mode; COMMON_MODE=1 ;;
  4) LABEL=equation6_stabilized; SPECTRAL=source_power_stabilized ;;
  5) LABEL=white_noise; NULL_METHOD=white_noise ;;
  6) LABEL=channel_permutation; NULL_METHOD=channel_permutation; REALIZATION=0 ;;
  *) echo "unexpected configuration $CONFIG" >&2; exit 2 ;;
esac

echo "config=$CONFIG label=$LABEL hour=$HOUR start=$START source=$SOURCE ram=$RAM spectral=$SPECTRAL null=$NULL_METHOD"
ARGS=(
  -u ambient_lellouch2019_exact_stack.py
  --action chunk
  --date 2024-12-20
  --start "$START"
  --nfiles 60
  --block-files 10
  --source-channel "$SOURCE"
  --ram-seconds "$RAM"
  --spectral-mode "$SPECTRAL"
  --null-method "$NULL_METHOD"
  --realization "$REALIZATION"
  --output-dir "$OUT"
)
if [[ "$COMMON_MODE" -eq 1 ]]; then
  ARGS+=(--common-mode)
fi
srun "$PY" "${ARGS[@]}"
