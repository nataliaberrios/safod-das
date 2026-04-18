#!/bin/bash
# Usage: bash submit_next.sh [N]
# Checks which daily CC output files are missing and submits the next N day
# indices as a SLURM array job. Defaults to N=10.

N=${1:-10}

NEXT_INDICES=$(/home/users/nberrios/miniconda3/envs/das/bin/python - $N <<'EOF'
import os, sys, pandas as pd

CSV     = "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv"
OUT_DIR = "/oak/stanford/groups/ettore88/nberrios"
NFILES_PER_DAY = 1440
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

db = pd.read_csv(CSV, delim_whitespace=True).drop_duplicates()
db = db[db["nSamples"] > 0].reset_index(drop=True)
db["startTime_dt"] = pd.to_datetime(db["startTime"], errors="coerce", utc=True)
db = db.dropna(subset=["startTime_dt"]).reset_index(drop=True)

total_days = len(db) // NFILES_PER_DAY
missing = []
for day_idx in range(total_days):
    start = day_idx * NFILES_PER_DAY
    if start >= len(db):
        break
    date_str = db["startTime_dt"].iloc[start].strftime("%Y-%m-%d")
    day_file   = os.path.join(OUT_DIR, f"{date_str}_day.npz")
    night_file = os.path.join(OUT_DIR, f"{date_str}_night.npz")
    if not (os.path.exists(day_file) and os.path.exists(night_file)):
        missing.append(day_idx)
    if len(missing) == N:
        break

print(",".join(str(i) for i in missing))
EOF
)

if [ -z "$NEXT_INDICES" ]; then
    echo "All days are done!"
    exit 0
fi

echo "Submitting day indices: $NEXT_INDICES"
sbatch --array="$NEXT_INDICES" run_daily_array.sbatch
