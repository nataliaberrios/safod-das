#!/bin/bash
#SBATCH --job-name=safod_fk_inject_post
#SBATCH --partition=serc
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_injection_recovery_v2_logs/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/fk_injection_recovery_v2_logs/%x-%j.err
set -euo pipefail

cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
/home/users/nberrios/miniconda3/envs/das/bin/python test_ambient_fk_injection_recovery_v2.py
/home/users/nberrios/miniconda3/envs/das/bin/python test_aggregate_ambient_fk_injection_recovery_v2.py
/home/users/nberrios/miniconda3/envs/das/bin/python -u aggregate_ambient_fk_injection_recovery_v2.py
/home/users/nberrios/miniconda3/envs/das/bin/python -u audit_ambient_fk_injection_recovery_v2.py
