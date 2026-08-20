#!/bin/bash
#SBATCH -J afk_agg
#SBATCH -p serc
#SBATCH -t 02:00:00
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
OUT="ambient_transfer/lellouch2019_exact_stack_afk"
run() {  # $1 = label, rest = args
  local lab="$1"; shift
  echo "================ aggregate: $lab ================"
  $PY -u ambient_lellouch2019_exact_stack.py --action aggregate \
    --date 2024-12-20 --total-files 300 --nfiles 300 --block-files 10 \
    --source-channel 23 --ram-seconds 0.1 --spectral-mode cross_correlation \
    --null-method ordered --realization 0 "$@" --output-dir "$OUT" \
    || echo "  (chunks not ready or aggregate failed for $lab)"
}
run "afk1 only (predicted worst)"      --afk-alpha 1.0
run "cm baseline (no afk)"             --common-mode
run "cm + afk1"                        --common-mode --afk-alpha 1.0
run "cm + afk2"                        --common-mode --afk-alpha 2.0
run "cm + svd2 + afk1"                 --common-mode --svd-rank 2 --afk-alpha 1.0
