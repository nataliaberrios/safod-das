#!/bin/bash
#SBATCH --job-name=safod_lellouch_fk_match
#SBATCH --partition=serc
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch_fk_matched_v2/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch_fk_matched_v2/%x-%j.err
set -euo pipefail

cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
mkdir -p awd_clean/ambient_transfer/lellouch_fk_matched_v2

/home/users/nberrios/miniconda3/envs/das/bin/python -u \
  awd_clean/ambient_lellouch_fk_matched_v2.py \
  --date 2024-12-20 \
  --start 0 \
  --nfiles 300 \
  --output-dir awd_clean/ambient_transfer/lellouch_fk_matched_v2
