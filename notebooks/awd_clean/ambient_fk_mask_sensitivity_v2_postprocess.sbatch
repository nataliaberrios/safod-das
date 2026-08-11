#!/bin/bash
#SBATCH --job-name=safod_fk_masks_post
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_mask_sensitivity_v2/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_mask_sensitivity_v2/%x-%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
/home/users/nberrios/miniconda3/envs/das/bin/python -u \
  awd_clean/aggregate_ambient_fk_mask_sensitivity_v2.py \
  --permutations 2000
