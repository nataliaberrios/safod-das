#!/bin/bash
#SBATCH -J mh029cov
#SBATCH -p serc
#SBATCH -t 00:25:00
#SBATCH --mem=8G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/mh029cov_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/mh029cov_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
python -u mh029_coverage.py
