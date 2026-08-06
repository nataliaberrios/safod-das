#!/bin/bash
#SBATCH -J chk2017
#SBATCH -p serc
#SBATCH -t 00:20:00
#SBATCH --mem=8G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/chk2017_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/chk2017_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
python -u check_2017_window.py
