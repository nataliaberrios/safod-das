#!/bin/bash
#SBATCH -J deepwide
#SBATCH -p serc
#SBATCH -t 08:00:00
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH --array=0-3
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#
# RE-MEASURE THE DEEP ARRIVAL WITH A LAG WINDOW THAT HOLDS THE FULL APERTURE.
#
# The +-0.35 s window cannot hold a 700 m offset below 1934 m/s, so on the
# 300-6000 m/s grid 71 of 229 trial velocities were scored on whichever offsets
# still fitted -- the NEAR ones, which sit closest to the zero-lag lobe and
# therefore score high. That produces a score rising smoothly as velocity falls
# for a purely geometric reason, and the Deep arrival (1443-1675 m/s) sits
# entirely inside the affected range. Every low-velocity number from this scan
# was measured in that regime.
#
# 2.5 s holds 700 m down to 280 m/s, so the whole grid is scored on the same
# fixed offset set and the velocities become comparable to each other.
#
#   task 0,1: src 400 and 800 -- the sources the published statistics came from
#   task 2,3: src 211 -- the WELLHEAD, which is the geometry the figures show
#             but for which no per-velocity null has ever existed (needs chunking)
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
OUT="ambient_transfer/deep_exact_stack"
case "$SLURM_ARRAY_TASK_ID" in
  0) SRC=400; CM=() ;;
  1) SRC=800; CM=() ;;
  2) SRC=211; CM=() ;;
  3) SRC=211; CM=(--common-mode) ;;
esac
COMMON=(--dataset deep --date deepA --nfiles 300 --block-files 10
        --source-channel "$SRC" --ram-seconds 0.1
        --spectral-mode cross_correlation --null-method ordered --realization 0
        --output-dir "$OUT")
# src 211 has no chunks at all; 400/800 already do, so skip the expensive pass
if [ "$SRC" = "211" ]; then
  $PY -u ambient_lellouch2019_exact_stack.py --action chunk --start 0 \
    "${COMMON[@]}" ${CM[@]+"${CM[@]}"}
fi
$PY -u ambient_lellouch2019_exact_stack.py --action aggregate --total-files 300 \
  --velocity-min 300 --velocity-max 6000 --velocity-step 25 --max-lag 2.5 \
  "${COMMON[@]}" ${CM[@]+"${CM[@]}"}
