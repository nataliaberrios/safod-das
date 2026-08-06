#!/bin/bash
#SBATCH --job-name=hrsn_reloc
#SBATCH --partition=serc
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
ml gcc/12.4.0
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
/home/users/nberrios/miniconda3/envs/das/bin/python -u hrsn_reloc.py
