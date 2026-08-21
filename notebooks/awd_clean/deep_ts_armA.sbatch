#!/bin/bash
#SBATCH -J deep_tsA
#SBATCH -p serc
#SBATCH -t 06:00:00
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --array=0-3
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
# block_spectra now streams; refuse to spend node hours if chunking
# does not reproduce reading the block whole.
$PY -u test_block_spectra_chunking.py || exit 1
if [ "$SLURM_ARRAY_TASK_ID" -lt 2 ]; then BAND=5-20;  T=$SLURM_ARRAY_TASK_ID
else                                      BAND=15-30; T=$((SLURM_ARRAY_TASK_ID-2)); fi
# ARM A is only 24.6 h, so 2 h blocks every 2 h -> full coverage, no gaps
$PY -u deep_timeseries.py --arm deepA --band "$BAND" \
  --block-hours 2 --stride-hours 2 --task "$T" --ntasks 2
