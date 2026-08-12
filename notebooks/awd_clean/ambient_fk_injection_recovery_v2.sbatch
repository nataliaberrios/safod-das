#!/bin/bash
#SBATCH --job-name=safod_fk_inject
#SBATCH --partition=serc
#SBATCH --array=0-7
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_injection_recovery_v2_logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_injection_recovery_v2_logs/%x-%A_%a.err
set -euo pipefail

cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean

VELOCITIES=(1800 1800 2750 2750 3200 3200 4000 4000)
DIRECTIONS=(1 -1 1 -1 1 -1 1 -1)
VELOCITY="${VELOCITIES[$SLURM_ARRAY_TASK_ID]}"
DIRECTION="${DIRECTIONS[$SLURM_ARRAY_TASK_ID]}"
OUTPUT="ambient_transfer/fk_injection_recovery_v2_n300"

mkdir -p "$OUTPUT" ambient_transfer/fk_injection_recovery_v2_logs
/home/users/nberrios/miniconda3/envs/das/bin/python -u \
  ambient_fk_injection_recovery_v2.py \
  --date 2024-12-20 \
  --start 0 \
  --nfiles 300 \
  --velocity "$VELOCITY" \
  --direction "$DIRECTION" \
  --output-dir "$OUTPUT"
