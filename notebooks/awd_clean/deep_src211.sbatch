#!/bin/bash
#SBATCH -J deep211
#SBATCH -p serc
#SBATCH -t 06:00:00
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --array=0-1
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#
# The per-velocity nulls that support the Deep arrival (p <= 0.0005,
# causal/acausal 2.34) were all computed at source channel 400, TVD 389 m.
# The figure the repo presents as "the Deep arrival" is sourced at channel 211,
# the WELLHEAD (TVD 2 m) -- which is the correct Lellouch Figure 7c geometry, a
# source at the top of the array. Those are different runs and the statistics
# from one do not license the figure from the other.
#
# There is no src211 aggregate at all, so this runs chunk AND aggregate, on the
# extended 300-6000 m/s grid (a fluid-guided mode near 1500-1700 m/s sits close
# to the old 1500 m/s grid floor, which is why the floor was dropped).
#   task 0: as processed everywhere else
#   task 1: common mode removed, since the Nano gather turned out to be 87 %
#           a single zero-moveout mode and Deep should be checked the same way
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
OUT="ambient_transfer/deep_exact_stack"
SRC=211
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then CM=(); TAG="no-cm"; else CM=(--common-mode); TAG="cm"; fi
echo "src $SRC  $TAG  grid 300-6000"
$PY -u ambient_lellouch2019_exact_stack.py --action chunk \
  --dataset deep --date deepA --start 0 --nfiles 300 --block-files 10 \
  --source-channel "$SRC" --ram-seconds 0.1 \
  --spectral-mode cross_correlation --null-method ordered --realization 0 \
  ${CM[@]+"${CM[@]}"} --output-dir "$OUT"
$PY -u ambient_lellouch2019_exact_stack.py --action aggregate \
  --dataset deep --date deepA --total-files 300 --nfiles 300 --block-files 10 \
  --source-channel "$SRC" --ram-seconds 0.1 \
  --velocity-min 300 --velocity-max 6000 --velocity-step 25 \
  --spectral-mode cross_correlation --null-method ordered --realization 0 \
  ${CM[@]+"${CM[@]}"} --output-dir "$OUT"
