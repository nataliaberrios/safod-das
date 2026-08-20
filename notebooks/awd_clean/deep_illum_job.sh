#!/bin/bash
#SBATCH -J deep_illum
#SBATCH -p serc
#SBATCH -t 03:00:00
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
/home/users/nberrios/miniconda3/envs/das/bin/python -u deep_fiber_illumination.py
