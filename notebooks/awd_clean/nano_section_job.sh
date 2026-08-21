#!/bin/bash
#SBATCH -J nano_sect
#SBATCH -p serc
#SBATCH -t 03:00:00
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH --array=0-1
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then EX=(--tag nano); else EX=(--common-mode --tag nanocm); fi
$PY -u deep_record_section.py --fibre nano --source 73 --nfiles 36 "${EX[@]}"
