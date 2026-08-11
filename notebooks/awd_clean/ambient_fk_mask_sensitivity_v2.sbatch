#!/bin/bash
#SBATCH --job-name=safod_fk_masks_v2
#SBATCH --partition=serc
#SBATCH --array=0-7
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_mask_sensitivity_v2/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_mask_sensitivity_v2/%x-%A_%a.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
/home/users/nberrios/miniconda3/envs/das/bin/python -u \
  awd_clean/run_ambient_fk_mask_sensitivity_v2.py \
  --day-index "$SLURM_ARRAY_TASK_ID"
