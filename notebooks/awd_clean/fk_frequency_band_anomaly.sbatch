#!/bin/bash
#SBATCH --partition=serc
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=safod_fk_frequency_band
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/%x-%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
/home/users/nberrios/miniconda3/envs/das/bin/python -u awd_clean/fk_frequency_band_anomaly_test.py --blocks-per-day 4 --nfiles 10
