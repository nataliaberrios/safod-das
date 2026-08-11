#!/bin/bash
#SBATCH --job-name=safod_fk_prefilter
#SBATCH --partition=serc
#SBATCH --array=0-19
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_prefilter_energy_v1_logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_prefilter_energy_v1_logs/%x-%A_%a.err
set -euo pipefail

cd /home/groups/ettore88/nberrios/safod_das_git/notebooks

PREFILTER_DATE="${PREFILTER_DATE:-2024-12-20}"
PREFILTER_START="${PREFILTER_START:-0}"
PREFILTER_FILES="${PREFILTER_FILES:-300}"
PREFILTER_OUTPUT="${PREFILTER_OUTPUT:-awd_clean/ambient_transfer/fk_prefilter_energy_v1_n300_r20}"

mkdir -p "$PREFILTER_OUTPUT" awd_clean/ambient_transfer/fk_prefilter_energy_v1_logs
/home/users/nberrios/miniconda3/envs/das/bin/python -u \
  awd_clean/ambient_fk_prefilter_energy_v1.py \
  --date "$PREFILTER_DATE" \
  --start "$PREFILTER_START" \
  --nfiles "$PREFILTER_FILES" \
  --null-start "$SLURM_ARRAY_TASK_ID" \
  --nulls 1 \
  --output-dir "$PREFILTER_OUTPUT"
