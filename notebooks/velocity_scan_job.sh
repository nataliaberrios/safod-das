#!/bin/bash
#SBATCH -J nano_vscan
#SBATCH -p serc
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24GB
#SBATCH -o /home/groups/ettore88/nberrios/safod_das_git/logs/%x_%j.out
#SBATCH -e /home/groups/ettore88/nberrios/safod_das_git/logs/%x_%j.err

# Submitted with sbatch rather than srun so the run survives the client session.
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
python nano_velocity_scan.py
