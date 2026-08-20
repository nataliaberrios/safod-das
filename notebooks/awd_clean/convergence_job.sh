#!/bin/bash
#SBATCH -J stackconv
#SBATCH -p serc
#SBATCH -t 01:00:00
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
/home/users/nberrios/miniconda3/envs/das/bin/python -u ambient_stack_convergence.py
