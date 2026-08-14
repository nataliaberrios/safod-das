#!/bin/bash
#SBATCH --job-name=awd_canonical
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=awd_clean/rebuild_%j.out
#SBATCH --error=awd_clean/rebuild_%j.err

ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
python paired_stack_job_deep_all.py
