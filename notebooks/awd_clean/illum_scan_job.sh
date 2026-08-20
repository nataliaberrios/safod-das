#!/bin/bash
#SBATCH -J illum_scan
#SBATCH -p serc
#SBATCH -t 04:00:00
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
SCAN_WINDOWS=${SCAN_WINDOWS:-240} \
  /home/users/nberrios/miniconda3/envs/das/bin/python -u illumination_window_scan.py
