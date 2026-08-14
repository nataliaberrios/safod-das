#!/bin/bash
#SBATCH --job-name=safod_seasonal_aggregate
#SBATCH --partition=serc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/%x-%j.out
#SBATCH --error=/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/ambient_transfer/%x-%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
$PY -u awd_clean/aggregate_seasonal_fk.py
$PY -u awd_clean/aggregate_seasonal_unfiltered.py
$PY - <<'PY2'
from pathlib import Path
from datetime import datetime, timezone
Path("awd_clean/ambient_transfer/seasonal_postprocess_complete.txt").write_text(datetime.now(timezone.utc).isoformat()+"\n")
PY2
