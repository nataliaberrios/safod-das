#!/bin/bash
#SBATCH -J nbrun
#SBATCH -p serc
#SBATCH -t 00:20:00
#SBATCH --mem=8G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/nbrun_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/nbrun_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
python -u run_nb_cells.py REPEATERS_dashboard.ipynb 2>&1 | tail -30
