#!/bin/bash
#SBATCH --job-name=fs_smoke
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# No pipe: piping python -u through grep re-buffers and left the last smoke test
# at a 0-byte log for 27 minutes with the job running fine.
ml gcc/12.4.0
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
export SLURM_ARRAY_TASK_ID=0 SLURM_ARRAY_TASK_MIN=0 SLURM_ARRAY_TASK_MAX=1999
export SCAN_OUT=/scratch/groups/ettore88/nberrios/safod_fullscan_smoke
/home/users/nberrios/miniconda3/envs/das/bin/python -u full_scan.py
