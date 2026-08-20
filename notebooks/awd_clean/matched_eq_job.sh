#!/bin/bash
#SBATCH -J matched_eq
#SBATCH -p serc
#SBATCH -t 01:00:00
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
/home/users/nberrios/miniconda3/envs/das/bin/python -u matched_earthquake_census.py
