#!/bin/bash
#SBATCH -J relabel2
#SBATCH -p serc
#SBATCH -t 06:00:00
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH --array=0-4
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
# Regenerate figures 03,05,06,07,08,09 so every one carries hours stacked and the
# LOCAL date range in the title, with UTC demoted to a footnote.
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
case "$SLURM_ARRAY_TASK_ID" in
  0) $PY -u deep_pervel.py ;;                                   # 09
  1) $PY -u nano_find_wellhead.py ;;                            # 05
  2) $PY -u ambient_stack_convergence.py ;;                     # 08
  3) $PY -u illumination_window_scan.py ;;                      # 07
  4) $PY -u interrogator_and_illumination_v2.py ;;              # 06
esac
