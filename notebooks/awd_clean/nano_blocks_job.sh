#!/bin/bash
#SBATCH -J nanoblk
#SBATCH -p serc
#SBATCH -t 04:00:00
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
/home/users/nberrios/miniconda3/envs/das/bin/python -u nano_blocks.py --block-hours 2
