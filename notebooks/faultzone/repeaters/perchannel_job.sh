#!/bin/bash
#SBATCH -J perchannel
#SBATCH -p serc
#SBATCH -t 06:00:00
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/perchannel_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/perchannel_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
python -u correlate_perchannel.py
