#!/bin/bash
#SBATCH -J deepfast
#SBATCH -p serc
#SBATCH -t 00:45:00
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
DEEP_SECONDS=6 DEEP_N_DROP=8 DEEP_N_AMBIENT=10 DEEP_APERTURE_STEP=2 \
  /home/users/nberrios/miniconda3/envs/das/bin/python -u deep_fiber_illumination.py
