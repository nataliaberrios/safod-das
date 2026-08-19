#!/bin/bash
#SBATCH --job-name=svdsweep
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --array=0-95%28
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack/logs/%x-%A_%a.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/lellouch2019_exact_stack/logs/%x-%A_%a.err

# Rank-k coherent-subspace removal sweep on 2024-12-20. 4 ranks x 24 hours = 96.
#
# Median removal takes out only a uniform-gain common mode and leaves
# (a_i - median(a))*c(t), which is coherent at zero lag. Measured: the
# median-removed velocity curve is flat on one day, corr(v,score) = -0.381, but
# +0.951 once four days are coherently stacked -- the residual accumulates while
# the incoherent noise averages down. Rank-k handles any gain pattern.
#
# Diagnostic of merit is corr(trial velocity, score): the pedestal indicator.
# Baseline config 0 is +0.976, median is -0.381 at one day but +0.951 stacked.

set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
WORK=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
OUT="$WORK/ambient_transfer/lellouch2019_exact_stack"
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
mkdir -p "$OUT/logs"; cd "$WORK"

RANKS=(1 2 4 8)
RANK=${RANKS[$(( SLURM_ARRAY_TASK_ID / 24 ))]}
HOUR=$(( SLURM_ARRAY_TASK_ID % 24 ))
START=$(( HOUR * 60 ))

echo "rank=$RANK hour=$HOUR start=$START"
srun "$PY" -u ambient_lellouch2019_exact_stack.py \
  --action chunk --date 2024-12-20 --start "$START" --nfiles 60 --block-files 10 \
  --source-channel 23 --ram-seconds 0.1 --spectral-mode cross_correlation \
  --null-method ordered --realization 0 --svd-rank "$RANK" --output-dir "$OUT"
