#!/bin/bash
#SBATCH --job-name=safod_fk_inject_pilot
#SBATCH --partition=serc
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_injection_recovery_v2_logs/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_injection_recovery_v2_logs/%x-%j.err
set -euo pipefail

cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUTPUT="ambient_transfer/fk_injection_recovery_v2_pilot"
mkdir -p "$OUTPUT" ambient_transfer/fk_injection_recovery_v2_logs
/home/users/nberrios/miniconda3/envs/das/bin/python -u \
  ambient_fk_injection_recovery_v2.py \
  --date 2024-12-20 \
  --start 0 \
  --nfiles 2 \
  --velocity 3200 \
  --direction 1 \
  --output-dir "$OUTPUT" \
  --stem fk_injection_v2_pilot_v3200_inc \
  --overwrite
