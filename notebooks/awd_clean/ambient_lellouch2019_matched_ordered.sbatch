#!/bin/bash
#SBATCH --job-name=lellouch_ordered
#SBATCH --partition=serc
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-1
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_nulls/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack_nulls/logs/%x-%A_%a.err

set -euo pipefail
WORK=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUT="$WORK/ambient_transfer/lellouch2019_exact_stack_nulls"
mkdir -p "$OUT/logs"
cd "$WORK"

if [[ "$SLURM_ARRAY_TASK_ID" -eq 0 ]]; then
  FK_MODE=unfiltered
else
  FK_MODE=downgoing_2p5_4p5
fi

/home/users/nberrios/miniconda3/envs/das/bin/python -u   ambient_lellouch2019_exact_stack.py   --date 2024-12-20   --start 0   --nfiles 300   --null-method ordered   --realization 0   --fk-mode "$FK_MODE"   --spectral-mode cross_correlation   --alignment-velocity data   --output-dir "$OUT"
