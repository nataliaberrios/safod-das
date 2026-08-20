#!/bin/bash
#SBATCH -J blame
#SBATCH -p serc
#SBATCH -t 01:30:00
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
echo "############ IS THE INTERROGATOR TO BLAME ############"
$PY -u interrogator_blame_test.py
echo; echo "############ SURFACE ILLUMINATION (Lellouch's diagnostic) ############"
$PY -u ambient_directional_asymmetry.py
