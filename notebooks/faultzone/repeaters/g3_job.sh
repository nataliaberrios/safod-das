#!/bin/bash
#SBATCH -J g3_hrsn
#SBATCH -p serc
#SBATCH -t 01:30:00
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/g3_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/g3_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
python -u dvv_hrsn.py
