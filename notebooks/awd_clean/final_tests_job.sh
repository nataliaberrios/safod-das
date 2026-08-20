#!/bin/bash
#SBATCH -J fig7c_tests
#SBATCH -p serc
#SBATCH -t 01:00:00
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
echo "################ 2017 window audit (verdict-branch fix) ################"
$PY -u lellouch2017_window_audit.py
echo; echo "################ low-k centroid / fixed-k test ################"
$PY -u ambient_fixed_k_test.py
echo; echo "################ matched earthquake census ################"
$PY -u matched_earthquake_census.py
echo "(exit $? -- a non-zero exit here is the NOT TESTABLE guard, not a crash)"
