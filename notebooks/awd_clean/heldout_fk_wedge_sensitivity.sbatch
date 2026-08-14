#!/bin/bash
#SBATCH --job-name=safod_fk_heldout
#SBATCH --partition=serc
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/%x-%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
/home/users/nberrios/miniconda3/envs/das/bin/python -u awd_clean/heldout_fk_wedge_sensitivity.py "$@"
