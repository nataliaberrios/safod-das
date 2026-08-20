#!/bin/bash
#SBATCH -J exec_qc_nb
#SBATCH -p serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# Execute the advisor-facing ambient F-K QC notebook in place.
# Build with SYSTEM python3 (it has nbformat); execute with the `das` KERNEL.
# The base miniconda env carries jupyter/nbconvert; `das` does not.
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
export PYDEVD_DISABLE_FILE_VALIDATION=1
export MPLBACKEND=Agg
/home/users/nberrios/miniconda3/bin/jupyter nbconvert --to notebook --execute --inplace \
  Ambient_FK_QC_workflow.ipynb \
  --ExecutePreprocessor.kernel_name=das \
  --ExecutePreprocessor.timeout=3600
echo "DONE $(date -u)"
