#!/bin/bash
#SBATCH -J g0_ref
#SBATCH -p serc
#SBATCH -t 00:20:00
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/g0ref_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/g0ref_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
python -u g0_refine.py
