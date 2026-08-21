#!/bin/bash
#SBATCH -J regen_sec
#SBATCH -p serc
#SBATCH -t 03:00:00
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH --array=0-2
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
# Regenerate the Deep depth sections with:
#   - the lag window and plot xlim DERIVED from the aperture (0.35 s could not
#     reach 700 m below 1934 m/s, so the arrival was neither scored nor drawn)
#   - the velocity line taken from THIS source's own measurement
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
case "$SLURM_ARRAY_TASK_ID" in
  0) $PY -u deep_section_depth.py --source 211 --nfiles 300 --tag top211 ;;
  1) $PY -u deep_section_depth.py --source 400 --nfiles 300 --tag src400 ;;
  2) $PY -u deep_section_depth.py --source 308 --limb return --nfiles 300 --tag pair200ret ;;
esac
