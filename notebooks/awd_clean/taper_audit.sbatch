#!/bin/bash
#SBATCH -J taper_audit
#SBATCH -p serc
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/%x_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/%x_%j.err
set -euo pipefail
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
export MPLBACKEND=Agg
python -u awd_clean/taper_audit.py
