#!/bin/bash
#SBATCH -J fk_down
#SBATCH -p serc
#SBATCH -t 02:00:00
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --array=0-3
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
# Corrected 2026-08-20. Previously ran at V=1675 (retracted) on along-fibre
# sections, and the Nano arm used deep_record_section_nano.npz whose virtual
# source is channel 10 -- in the air. Now: the WELLHEAD depth gather at its own
# measured velocity, and a Nano gather from a downhole source at the in-band
# 3200 m/s rather than the 30-60 Hz value of 2950.
case "$SLURM_ARRAY_TASK_ID" in
  0) IN=deep_section_depth_top211.npz;      V=1350; T=deep_top211  ;;
  1) IN=deep_section_depth_src400.npz;      V=1525; T=deep_src400  ;;
  2) IN=deep_section_depth_pair200ret.npz;  V=1350; T=deep_return  ;;
  3) IN=deep_record_section_nano_valid.npz; V=3200; T=nano_ch73    ;;
esac
[ -f "$IN" ] || { echo "missing $IN"; exit 0; }
echo "=== $T : $IN at $V m/s ==="
$PY -u fk_downgoing.py --input "$IN" --v-ref "$V" --tag "$T"
