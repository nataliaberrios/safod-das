#!/bin/bash
#SBATCH --job-name=safod_fk_null_v2
#SBATCH --partition=serc
#SBATCH --array=0-19
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_full_pipeline_null_v2/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_full_pipeline_null_v2/%x-%A_%a.err
set -euo pipefail

cd /home/groups/ettore88/nberrios/safod_das_git/notebooks

# Override these with sbatch --export when expanding beyond the default pilot.
NULL_DATE="${NULL_DATE:-2024-12-20}"
NULL_FILE_START="${NULL_FILE_START:-0}"
NULL_FILE_COUNT="${NULL_FILE_COUNT:-30}"
NULL_INDEX_OFFSET="${NULL_INDEX_OFFSET:-0}"
NULL_METHODS="${NULL_METHODS:-channel_permutation,circular_time_shift}"
NULL_OUTPUT_DIR="${NULL_OUTPUT_DIR:-awd_clean/ambient_transfer/fk_full_pipeline_null_v2}"
NULL_REALIZATION=$((NULL_INDEX_OFFSET + SLURM_ARRAY_TASK_ID))

mkdir -p "$NULL_OUTPUT_DIR"
/home/users/nberrios/miniconda3/envs/das/bin/python -u \
  awd_clean/ambient_fk_full_pipeline_null_v2.py \
  --date "$NULL_DATE" \
  --start "$NULL_FILE_START" \
  --nfiles "$NULL_FILE_COUNT" \
  --null-start "$NULL_REALIZATION" \
  --nulls 1 \
  --null-methods "$NULL_METHODS" \
  --output-dir "$NULL_OUTPUT_DIR"
