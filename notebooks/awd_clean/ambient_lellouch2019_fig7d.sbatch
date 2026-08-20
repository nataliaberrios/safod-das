#!/bin/bash
#SBATCH --job-name=fig7d
#SBATCH --partition=serc
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-23%24
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/logs/%x_%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/logs/%x_%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
START=$(( SLURM_ARRAY_TASK_ID * 60 ))
EXTRA=(); [[ "${FIG7D_CM:-0}" -eq 1 ]] && EXTRA=(--common-mode)
srun /home/users/nberrios/miniconda3/envs/das/bin/python -u ambient_lellouch2019_fig7d.py \
  --date 2024-12-20 --start "$START" --nfiles 60 ${EXTRA[@]+"${EXTRA[@]}"}
