#!/bin/bash
#SBATCH -J fig7c_figs
#SBATCH -p serc
#SBATCH -t 00:40:00
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -o ../logs/%x_%j.out
#SBATCH -e ../logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/fig7c_negative
/home/users/nberrios/miniconda3/envs/das/bin/python -u make_figures.py
