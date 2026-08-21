#!/bin/bash
#SBATCH -J combC
#SBATCH -p serc
#SBATCH -t 02:00:00
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --array=0-1
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then B=5-20; else B=15-30; fi
$PY -u deep_timeseries.py --arm deepC --band "$B" --combine
