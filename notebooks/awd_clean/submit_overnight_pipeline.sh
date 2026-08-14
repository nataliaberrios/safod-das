#!/bin/bash
set -euo pipefail

parent_job="${1:?usage: submit_overnight_pipeline.sh PARENT_JOB_ID}"
repo_dir="/home/groups/ettore88/nberrios/safod_das_git/notebooks"

sbatch \
  --job-name=awd_overnight \
  --partition=serc \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem=8G \
  --time=00:30:00 \
  --dependency="afterany:${parent_job}" \
  --output="awd_clean/overnight_%j.out" \
  --error="awd_clean/overnight_%j.err" \
  --export="ALL,AWD_PARENT_JOB=${parent_job}" \
  --wrap="ml gcc/12.4.0; source /home/users/nberrios/miniconda3/etc/profile.d/conda.sh; conda activate das; cd ${repo_dir}; python awd_clean/overnight_audit.py"
