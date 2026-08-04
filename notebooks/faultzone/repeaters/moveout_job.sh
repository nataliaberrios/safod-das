#!/bin/bash
#SBATCH -J moveout
#SBATCH -p serc
#SBATCH -t 03:00:00
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/moveout_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/moveout_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
python -u moveout_test.py
