#!/bin/bash
#SBATCH -J pervel
#SBATCH -p serc
#SBATCH -t 02:00:00
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
/home/users/nberrios/miniconda3/envs/das/bin/python -u deep_pervel.py
