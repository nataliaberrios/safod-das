#!/bin/bash
#SBATCH --job-name=creep_gate
#SBATCH --partition=serc
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# 360 files x 30000 x 900 int32 read one at a time; only the 1 Hz decimated
# product is retained, so peak memory is one file (108 MB) plus 873 ch x 21600 s.
ml gcc/12.4.0
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
CREEP_NFILES=${CREEP_NFILES:-360} CREEP_START=${CREEP_START:-2024-06-01} \
  /home/users/nberrios/miniconda3/envs/das/bin/python -u creep_strain_gate.py
