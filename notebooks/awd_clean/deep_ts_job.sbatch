#!/bin/bash
#SBATCH -J deep_ts
#SBATCH -p serc
#SBATCH -t 12:00:00
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --array=0-15
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
# 16 tasks: 8 per band, both bands in parallel.
#   band 5-20 Hz  : Lellouch's, used everywhere else, kept for comparability
#   band 15-30 Hz : where the AWD active source says the Deep mode is STRONGEST
# ARM C = 917 h with no personnel on site. 4 h blocks every 12 h -> ~76 blocks,
# spanning 2026-05-01 to 06-08, giving both a time series and a grand stack.
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
# block_spectra now streams; refuse to spend node hours if chunking
# does not reproduce reading the block whole.
$PY -u test_block_spectra_chunking.py || exit 1
if [ "$SLURM_ARRAY_TASK_ID" -lt 8 ]; then BAND=5-20;  T=$SLURM_ARRAY_TASK_ID
else                                      BAND=15-30; T=$((SLURM_ARRAY_TASK_ID-8)); fi
$PY -u deep_timeseries.py --arm deepC --band "$BAND" \
  --block-hours 4 --stride-hours 12 --task "$T" --ntasks 8
