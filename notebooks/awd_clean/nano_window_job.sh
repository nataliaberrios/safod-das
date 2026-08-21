#!/bin/bash
#SBATCH -J nanowin
#SBATCH -p serc
#SBATCH -t 04:00:00
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH --array=0-2
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
# PAIRED test, matched to 63 records (5.25 h) so illumination is not confounded
# with stack depth:
#   0  crew active   Mon 15 Jun 15:48-20:58 local (recording starts 15:48, so
#                    the requested 09:00 start yields the evening only)
#   1  quiet night   Tue 16 Jun 00:00-05:15 local, nobody working
#   2  next daytime  Tue 16 Jun 09:00-14:15 local, before the survey begins
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
case "$SLURM_ARRAY_TASK_ID" in
  0) S=2026-06-15T09:00; E=2026-06-15T21:00; T=jun15_crew   ;;
  1) S=2026-06-16T00:00; E=2026-06-16T06:00; T=jun16_night  ;;
  2) S=2026-06-16T09:00; E=2026-06-16T15:00; T=jun16_day    ;;
esac
$PY -u nano_window.py --start "$S" --end "$E" --tag "$T" --max-records 63
