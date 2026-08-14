#!/bin/bash
#SBATCH --job-name=fig7c
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24GB
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/logs/%x_%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/logs/%x_%A_%a.err

# Lellouch Figure 7c reproduction on the 2024-2025 archive.
# Array task i processes FILES_PER_TASK files starting at i*FILES_PER_TASK.
#   sbatch --array=0-23 ambient_lellouch_fig7c.sbatch                 # full day
#   sbatch --array=0 --export=ALL,FILES_PER_TASK=4 ...                # smoke test

set -euo pipefail
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean

DATE=${FIG7C_DATE:-2024-12-20}
FILES_PER_TASK=${FILES_PER_TASK:-60}
TASK=${SLURM_ARRAY_TASK_ID:-0}
START=$(( TASK * FILES_PER_TASK ))

# ablation switches, empty by default (all three corrections ON)
EXTRA=${FIG7C_EXTRA:-}

echo "date=$DATE task=$TASK start=$START nfiles=$FILES_PER_TASK extra='$EXTRA'"
srun $PY -u ambient_lellouch_fig7c.py \
    --date "$DATE" --start "$START" --nfiles "$FILES_PER_TASK" $EXTRA
