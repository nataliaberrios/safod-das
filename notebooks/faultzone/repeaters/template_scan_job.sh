#!/bin/bash
#SBATCH --job-name=tscan
#SBATCH --partition=serc
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --array=0-13
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
# One day per task: ~1440 files x 108 MB read serially, only the beam retained.
ml gcc/12.4.0
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
TEMPLATE_TAG=${TEMPLATE_TAG:-ev_20240708T083036} SCAN_START=${SCAN_START:-2024-09-15} \
  /home/users/nberrios/miniconda3/envs/das/bin/python -u template_scan.py
