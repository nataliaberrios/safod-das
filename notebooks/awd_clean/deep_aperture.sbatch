#!/bin/bash
#SBATCH -J deepap
#SBATCH -p serc
#SBATCH -t 08:00:00
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH --array=0-1
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#
# DOUBLE THE APERTURE FROM THE WELLHEAD.
#
# The 50-700 m offsets are Lellouch's Figure 7c aperture, inherited by the
# reproduction and never revisited for this fibre. The Deep near-vertical
# section runs ch 211-949, TVD 2-1513 m, so from the wellhead 1400 m of offset
# is available -- twice what has been used -- and the far receiver still lands
# at ch 907, inside the near-vertical section (with room for the R+-10 sum).
#
# A longer lever arm means more moveout across the gather at a given velocity,
# so the velocity is better resolved. It is a change of selection region, so it
# is run ALONGSIDE the 700 m result rather than replacing it, and both are
# reported.
#
# The lag window scales with the aperture: 1400 m at the 300 m/s scan floor
# needs 4.67 s, so 5.5 s.
#
#   task 0: 1400 m aperture, wellhead
#   task 1: 700 m aperture, wellhead -- the control, same code path
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
OUT="ambient_transfer/deep_exact_stack_aperture"
mkdir -p "$OUT"
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then OFF=1400; LAG=5.5; else OFF=700; LAG=2.5; fi
COMMON=(--dataset deep --date deepA --nfiles 300 --block-files 10
        --source-channel 211 --ram-seconds 0.1
        --spectral-mode cross_correlation --null-method ordered --realization 0
        --max-offset "$OFF" --output-dir "$OUT")
$PY -u ambient_lellouch2019_exact_stack.py --action chunk --start 0 "${COMMON[@]}"
$PY -u ambient_lellouch2019_exact_stack.py --action aggregate --total-files 300 \
  --velocity-min 300 --velocity-max 6000 --velocity-step 25 --max-lag "$LAG" \
  "${COMMON[@]}"
