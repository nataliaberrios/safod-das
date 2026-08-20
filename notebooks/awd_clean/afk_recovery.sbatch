#!/bin/bash
#SBATCH -J afk_rec
#SBATCH -p serc
#SBATCH -t 08:00:00
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --array=0-4
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#
# Adaptive f-k recovery attempt (Isken et al. 2022, 10.1093/gji/ggac229).
# AMBIENT_CC_LITERATURE_REVIEW.md section 1 predicts, BEFORE this run:
#   task 0 (AFK alone, no static removal) should be the WORST, because an
#          adaptive filter enhances the dominant coherent component and ours is
#          the static fixed-k pattern holding ~39 % of in-band energy;
#   tasks 2-4 (static removal FIRST, then AFK) are the ones with a chance.
# Task 1 is the config-3 baseline with no AFK, so the AFK effect is isolated.
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
OUT="ambient_transfer/lellouch2019_exact_stack_afk"
mkdir -p "$OUT"

case "$SLURM_ARRAY_TASK_ID" in
  0) ARGS=(--afk-alpha 1.0) ;;                                   # prediction: worst
  1) ARGS=(--common-mode) ;;                                     # baseline, no AFK
  2) ARGS=(--common-mode --afk-alpha 1.0) ;;
  3) ARGS=(--common-mode --afk-alpha 2.0) ;;
  4) ARGS=(--common-mode --svd-rank 2 --afk-alpha 1.0) ;;
  *) exit 2 ;;
esac
echo "config $SLURM_ARRAY_TASK_ID : ${ARGS[*]}"

$PY -u ambient_lellouch2019_exact_stack.py --action chunk \
  --date 2024-12-20 --start 0 --nfiles 300 --block-files 10 \
  --source-channel 23 --ram-seconds 0.1 --spectral-mode cross_correlation \
  --null-method ordered --realization 0 \
  "${ARGS[@]}" --output-dir "$OUT"

# --nfiles MUST match the chunk step: aggregate() builds expected chunk names as
# range(0, total_files, nfiles), so omitting it here made the aggregate look for
# five 60-file chunks when the chunk step had written one 300-file chunk.
$PY -u ambient_lellouch2019_exact_stack.py --action aggregate \
  --date 2024-12-20 --total-files 300 --nfiles 300 --block-files 10 \
  --source-channel 23 --ram-seconds 0.1 --spectral-mode cross_correlation \
  --null-method ordered --realization 0 \
  "${ARGS[@]}" --output-dir "$OUT"
