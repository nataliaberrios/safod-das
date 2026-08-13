#!/bin/bash
#SBATCH --job-name=lellouch2019_day
#SBATCH --partition=serc
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_reproduction_v1/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_reproduction_v1/%x-%j.err

set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
/home/users/nberrios/miniconda3/envs/das/bin/python -u \
    ambient_lellouch2019_reproduction_v1.py \
    --date 2024-12-20
