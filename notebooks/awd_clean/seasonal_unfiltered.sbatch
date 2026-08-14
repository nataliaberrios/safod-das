#!/bin/bash
#SBATCH --partition=serc
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/%x-%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
#SBATCH --job-name=safod_unfiltered
/home/users/nberrios/miniconda3/envs/das/bin/python -u awd_clean/run_seasonal_unfiltered.py
