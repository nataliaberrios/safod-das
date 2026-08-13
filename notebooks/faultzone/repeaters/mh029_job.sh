#!/bin/bash
#SBATCH --job-name=mh029_xc
#SBATCH --partition=serc
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
ml gcc/12.4.0
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
/home/users/nberrios/miniconda3/envs/das/bin/python -u mh029_das_xcorr.py
