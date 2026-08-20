#!/bin/bash
#SBATCH -J fixedk2
#SBATCH -p serc
#SBATCH -t 00:20:00
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
/home/users/nberrios/miniconda3/envs/das/bin/python -u ambient_fixed_k_test.py
