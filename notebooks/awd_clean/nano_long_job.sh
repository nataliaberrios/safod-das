#!/bin/bash
#SBATCH -J nano_long
#SBATCH -p serc
#SBATCH -t 08:00:00
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
# 300 x 5-minute records = 25 h, the full pre-drop Nano window.
# Cross-spectra accumulate batch by batch so peak memory is set by the 6-record
# batch, not by the 25 h total.
/home/users/nberrios/miniconda3/envs/das/bin/python -u nano_long_stack.py \
  --nfiles 300 --batch 6 --source 73
