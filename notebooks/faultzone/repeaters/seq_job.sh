#!/bin/bash
#SBATCH -J seqcreep
#SBATCH -p serc
#SBATCH -t 00:30:00
#SBATCH --mem=8G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/seq_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/seq_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
python -u sequences_and_creep.py
