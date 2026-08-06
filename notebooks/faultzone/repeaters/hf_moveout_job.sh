#!/bin/bash
#SBATCH --job-name=hf_extract
#SBATCH --partition=serc
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --array=0-11
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Stage HF-1: re-extract the 73 moveout-test events at the native 500 Hz.
# Stage HF-2 (hf_analyse_job.sh) runs as a dependent job once this array clears.

ml gcc/12.4.0
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters

echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date)"
$PY -u hf_extract.py
echo "done at $(date)"
