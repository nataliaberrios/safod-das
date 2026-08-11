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

NULL_DATE="${NULL_DATE:-2024-12-20}"
NULL_FILE_START="${NULL_FILE_START:-0}"
NULL_FILE_COUNT="${NULL_FILE_COUNT:-30}"
NULL_INDEX_OFFSET="${NULL_INDEX_OFFSET:-0}"
NULLS_PER_TASK="${NULLS_PER_TASK:-1}"
NULL_METHODS="${NULL_METHODS:-channel_permutation,circular_time_shift}"
NULL_OUTPUT_DIR="${NULL_OUTPUT_DIR:-awd_clean/ambient_transfer/fk_full_pipeline_null_v2}"

mkdir -p "$NULL_OUTPUT_DIR"
/home/users/nberrios/miniconda3/envs/das/bin/python -u \
  awd_clean/run_ambient_fk_full_pipeline_null_v2.py \
  --task-id "$SLURM_ARRAY_TASK_ID" \
  --date "$NULL_DATE" \
  --start "$NULL_FILE_START" \
  --nfiles "$NULL_FILE_COUNT" \
  --null-index-offset "$NULL_INDEX_OFFSET" \
  --nulls-per-task "$NULLS_PER_TASK" \
  --null-methods "$NULL_METHODS" \
  --output-dir "$NULL_OUTPUT_DIR"
