#!/bin/bash
#SBATCH --job-name=deep_creep
#SBATCH --partition=serc
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# 1249 LF files x ~15 MB; only the 1/min decimated product is retained.
ml gcc/12.4.0
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
FIBRE=${FIBRE:-deep} DEEP_DECIM=${DEEP_DECIM:-60} DEEP_NFILE=${DEEP_NFILE:-0} \
  /home/users/nberrios/miniconda3/envs/das/bin/python -u deep_creep_strain.py
