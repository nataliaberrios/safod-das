#!/bin/bash
#SBATCH -J fixed_k
#SBATCH -p serc
#SBATCH -t 00:30:00
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
$PY -u _verify_prefix_tmp.py
echo "==================== fixed-k test ===================="
$PY -u ambient_fixed_k_test.py
