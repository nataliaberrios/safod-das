#!/bin/bash
#SBATCH --job-name=safod_daily
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-429          # 430 days ~ 14 months; adjust if needed
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

module --force purge
module load devel math
module load cuda/11.7.1
module load py-pytorch/2.0.0_py39
module unload py-scipy

export PYTHONNOUSERSITE=1
export PYTHONPATH=/home/users/nberrios/miniconda3/envs/das/lib/python3.9/site-packages:/share/software/user/open/py-pytorch/2.0.0_py39/lib/python3.9/site-packages:/home/groups/edunham/nberrios/safod_das/DAS-utilities/python

cd /home/groups/edunham/nberrios/safod_das || exit 1

mkdir -p logs

/home/users/nberrios/miniconda3/envs/das/bin/python stack_daily.py
# SLURM_ARRAY_TASK_ID is read inside stack_daily.py via os.environ
