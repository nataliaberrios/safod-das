#!/bin/bash
#SBATCH -J g0_reg
#SBATCH -p serc
#SBATCH -t 04:00:00
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/g0_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/g0_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
python -u channel_depth_registration.py
