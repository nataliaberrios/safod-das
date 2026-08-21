#!/bin/bash
#SBATCH -J nano_diag
#SBATCH -p serc
#SBATCH -t 01:00:00
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
$PY -u nano_diagnose.py
echo "=== Nano section with a VALID source (ch 73, the wellhead) ==="
$PY -u deep_record_section.py --fibre nano --source 73 --nfiles 24 --tag nano_valid
