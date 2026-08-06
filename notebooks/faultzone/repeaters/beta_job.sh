#!/bin/bash
#SBATCH -J beta
#SBATCH -p serc
#SBATCH -t 08:00:00
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/beta_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/beta_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
python -u beta_similarity.py
