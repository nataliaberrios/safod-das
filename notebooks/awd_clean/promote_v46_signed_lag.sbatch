#!/bin/bash
#SBATCH --job-name=safod_signed_fk_v46
#SBATCH --partition=serc
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/signed_lag_v2/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/signed_lag_v2/%x-%j.err
set -euo pipefail

cd /home/groups/ettore88/nberrios/safod_das_git/notebooks

/home/users/nberrios/miniconda3/envs/das/bin/python \
  awd_clean/test_audit_seasonal_signed_fk_v2.py
/home/users/nberrios/miniconda3/envs/das/bin/python \
  awd_clean/test_promote_v46_signed_lag.py
/home/users/nberrios/miniconda3/envs/das/bin/python \
  awd_clean/promote_v46_signed_lag.py
/home/users/nberrios/miniconda3/envs/das/bin/python \
  run_nb_cells.py awd_clean/AWD_results_dashboard.ipynb

git diff --check -- \
  awd_clean/AWD_results_dashboard.ipynb \
  awd_clean/AWD_advisor_figure_guide.tex \
  awd_clean/AWD_advisor_figure_guide_v46.tex \
  awd_clean/promote_v46_signed_lag.py \
  awd_clean/promote_v46_signed_lag.sbatch \
  awd_clean/test_audit_seasonal_signed_fk_v2.py \
  awd_clean/test_promote_v46_signed_lag.py

git add \
  awd_clean/AWD_results_dashboard.ipynb \
  awd_clean/AWD_advisor_figure_guide.tex \
  awd_clean/AWD_advisor_figure_guide_v46.tex \
  awd_clean/promote_v46_signed_lag.py \
  awd_clean/promote_v46_signed_lag.sbatch \
  awd_clean/test_audit_seasonal_signed_fk_v2.py \
  awd_clean/test_promote_v46_signed_lag.py
git add -f \
  awd_clean/ambient_transfer/signed_lag_v2/seasonal_signed_fk_v2_aggregate.json \
  awd_clean/ambient_transfer/signed_lag_v2/seasonal_signed_fk_v2_aggregate.png \
  awd_clean/ambient_transfer/signed_lag_v2/seasonal_signed_fk_v2_completion_audit.json

git diff --cached --check
git commit \
  -m "Promote corrected signed-lag comparison to v46" \
  -m "Record the independently audited eight-day physical-lag branch comparison, update stale status flags and detailed captions, and retain the conditional F-K interpretation. Synthetic audit and promotion tests plus the full dashboard smoke test pass. No remote push is performed."
