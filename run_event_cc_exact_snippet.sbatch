#!/bin/bash
#SBATCH --job-name=event_cc_exact
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module --force purge
module load devel math
module load cuda/11.7.1
module load py-pytorch/2.0.0_py39
module unload py-scipy

P1=/home/users/nberrios/miniconda3/envs/das/lib/python3.9/site-packages
P2=/share/software/user/open/py-pytorch/2.0.0_py39/lib/python3.9/site-packages
P3=/home/groups/edunham/nberrios/safod_das/DAS-utilities/python

export PYTHONNOUSERSITE=1
export PYTHONPATH="$P1:$P2:$P3"

cd /home/groups/ettore88/nberrios/safod_das_git || exit 1
mkdir -p logs

EVENT_TIME=2024-05-25T04:56:28.150000Z
EVENT_LABEL=EQ_20240525_045628

OUTPUT_VERSION=event20_base \
USE_TEMPORAL_NORMALIZATION=false \
USE_SECOND_BANDPASS=false \
/home/users/nberrios/miniconda3/envs/das/bin/python \
event_cc_20s_exact_snippet.py "$EVENT_TIME" "$EVENT_LABEL"
