#!/bin/bash
#SBATCH -J scec_nb
#SBATCH -p serc
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24GB
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/%x_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/%x_%j.err
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
# execute to a scratch copy so the committed notebook stays output-free
export MPLBACKEND=Agg
python run_nb_cells.py SCEC_awd_resolving_power.ipynb
