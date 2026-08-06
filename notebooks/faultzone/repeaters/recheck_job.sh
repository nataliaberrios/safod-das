#!/bin/bash
#SBATCH --job-name=recheck
#SBATCH --partition=serc
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96GB
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
ml gcc/12.4.0
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
echo "### STRESS DROP 500 Hz ###"; $PY -u stress_drop_500.py
echo; echo "### DVV RECHECK 500 Hz ###"; $PY -u recheck_dvv_500.py
