#!/bin/bash
#SBATCH --job-name=safod_signed_fk_v2_post
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/signed_lag_v2/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/signed_lag_v2/%x-%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
/home/users/nberrios/miniconda3/envs/das/bin/python -u awd_clean/aggregate_seasonal_signed_fk_v2.py
