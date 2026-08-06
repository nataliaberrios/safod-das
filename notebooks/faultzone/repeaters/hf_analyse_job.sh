#!/bin/bash
#SBATCH --job-name=hf_analyse
#SBATCH --partition=serc
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Stage HF-2: four bands x four stacking variants at 500 Hz.
# Memory is the binding constraint: prep() holds 700 ch x 12500 samples float64
# per event and slant_scan FFTs the whole array, so 73 events resident is ~5 GB
# before the 81-slowness scan working set.

ml gcc/12.4.0
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters

echo "start $(date) on $(hostname)"
$PY -u hf_moveout_test.py
echo "done $(date)"
