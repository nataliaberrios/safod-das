#!/bin/bash
#SBATCH --job-name=interval_dvv_gate
#SBATCH --partition=serc
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96GB
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# D0 precision gate for the depth-localized dv/v measurement.
#
# Memory: prep() holds ~700 ch x 12500 samples float64 per event (70 MB), and the
# pair cache keeps both events of every pair resident plus their aligned copies.
# 50 pairs x 4 arrays is the worst case, hence 96 GB.

ml gcc/12.4.0
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters

echo "start $(date) on $(hostname)"
$PY -u interval_dvv_gate.py
echo "done $(date)"
