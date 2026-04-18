import os
import sys
import numpy as np
import pandas as pd
import tqdm
import DASutils
from cc_tools import temporal_normalization, computeCC

# ================== USER PARAMETERS ==================
CSV     = "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv"
OUT_DIR = "/oak/stanford/groups/ettore88/nberrios"

DATA_ROOT_OLD = "/oak/stanford/groups/ettore88/data/SAFODAS1-harddrive-transfer"
DATA_ROOT_NEW = "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer"


def normalize_file_path(path):
    path = str(path)
    if os.path.exists(path):
        return path
    if path.startswith(DATA_ROOT_OLD):
        alt = path.replace(DATA_ROOT_OLD, DATA_ROOT_NEW, 1)
        if os.path.exists(alt):
            return alt
    return path

ch_start, ch_end         = 150, 800
virtual_source_original  = 200

fmin, fmax  = 5.0, 20.0
max_lag     = 1.0
TN_WINDOW   = 20.0
WINDOW_SEC  = 30.0
OVERLAP     = 0.5
DO_TEMP_NORM = False

DAY_START_HOUR = 8
DAY_END_HOUR   = 20
UTC_OFFSET     = -7
# =====================================================

# --- get day index from SLURM array task ID ---
day_idx = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1] if len(sys.argv) > 1 else 0))

# --- load CSV ---
print(f"Day index: {day_idx}")
db = pd.read_csv(CSV, delim_whitespace=True).drop_duplicates()
db = db[db["nSamples"] > 0].reset_index(drop=True)
db["startTime_dt"] = pd.to_datetime(db["startTime"], errors="coerce", utc=True)
db["endTime_dt"]   = pd.to_datetime(db["endTime"],   errors="coerce", utc=True)
db = db.dropna(subset=["startTime_dt", "endTime_dt"]).reset_index(drop=True)
dur_s = (db["endTime_dt"] - db["startTime_dt"]).dt.total_seconds()
db = db[dur_s > 1.0].reset_index(drop=True)
db["date"] = db["startTime_dt"].dt.strftime("%Y-%m-%d")
db["file_norm"] = db["file"].map(normalize_file_path)

unique_dates = db["date"].drop_duplicates().tolist()
if day_idx >= len(unique_dates):
    print(f"Day {day_idx}: only {len(unique_dates)} valid dates available, skipping.")
    sys.exit(0)

date_str = unique_dates[day_idx]
db_day = db[db["date"] == date_str].reset_index(drop=True)
selectedFiles = db_day["file_norm"].tolist()
selectedFiles = [f for f in selectedFiles if os.path.exists(f)]

if len(selectedFiles) == 0:
    print(f"Day {day_idx}: no files, skipping.")
    sys.exit(0)

# --- use the calendar date as output tag ---
print(f"Date: {date_str}  ({len(selectedFiles)} files)")

# --- skip if already done ---
os.makedirs(OUT_DIR, exist_ok=True)
outfile_day   = os.path.join(OUT_DIR, f"{date_str}_day.npz")
outfile_night = os.path.join(OUT_DIR, f"{date_str}_night.npz")
if os.path.exists(outfile_day) and os.path.exists(outfile_night):
    print(f"Output already exists for {date_str}, skipping.")
    sys.exit(0)

file_to_time = db.set_index("file_norm")["startTime_dt"].to_dict()

# --- get fs/dt from first file ---
DAS0, info0 = DASutils.readFile_HDF(
    [selectedFiles[0]], 0.05, 24.0, verbose=0,
    preproc=True, diff=True, taper=False,
    desampling=True, nChbuffer=900, system="OptaSense"
)
fs = info0["fs"]
dt = 1.0 / fs

nch          = DAS0[ch_start:ch_end, :].shape[0]
ch_buffer_in = nch
isource      = virtual_source_original - ch_start
assert 0 <= isource < nch, "Virtual source outside slice!"

win_npts  = int(WINDOW_SEC / dt)
step_npts = int(win_npts * (1 - OVERLAP))

print(f"fs={fs}  dt={dt}  window={WINDOW_SEC}s  overlap={OVERLAP}")

# --- stacks ---
variants = ["combo"]   # full run: only combo (median -> bandpass -> optional TN)
stacks_day    = {v: None for v in variants}
stacks_night  = {v: None for v in variants}
nstacks_day   = {v: 0 for v in variants}
nstacks_night = {v: 0 for v in variants}

for f in tqdm.tqdm(selectedFiles):
    try:
        file_time  = file_to_time.get(f)
        hour_local = (pd.Timestamp(file_time).hour + UTC_OFFSET) % 24
        is_day     = DAY_START_HOUR <= hour_local < DAY_END_HOUR
        stacks_use  = stacks_day   if is_day else stacks_night
        nstacks_use = nstacks_day  if is_day else nstacks_night

        DAS, info = DASutils.readFile_HDF(
            [f], 0.05, 24.0, verbose=0,
            preproc=True, diff=True, taper=False,
            desampling=True, nChbuffer=900, system="OptaSense"
        )
        X    = DAS[ch_start:ch_end, :].astype(np.float64, copy=False)
        npts = X.shape[1]

        for i0 in range(0, npts - win_npts + 1, step_npts):
            Xraw   = X[:, i0:i0 + win_npts].copy()
            Xcombo = Xraw - np.median(Xraw, axis=0, keepdims=True)
            Xcombo = DASutils.bandpass2D_c(Xcombo, fmin, fmax, dt, zerophase=True)
            if DO_TEMP_NORM:
                Xcombo = temporal_normalization(Xcombo.copy(), fs, window_time=TN_WINDOW)

            cc_i = computeCC(Xcombo, dt, max_lag,
                             isource=isource,
                             ch_buffer_in=ch_buffer_in,
                             whitening_params=None)
            if stacks_use["combo"] is None:
                stacks_use["combo"] = np.zeros_like(cc_i)
            stacks_use["combo"] += cc_i
            nstacks_use["combo"] += 1

    except Exception as e:
        print("Skipping file:", e)
        continue

# --- save ---
common = dict(
    date=date_str, day_idx=day_idx,
    fs=fs, dt=dt,
    ch_start=ch_start, ch_end=ch_end,
    src=virtual_source_original,
    fmin=fmin, fmax=fmax,
    max_lag=max_lag,
    window_sec=WINDOW_SEC, overlap=OVERLAP,
    do_temp_norm=DO_TEMP_NORM,
    day_start_hour=DAY_START_HOUR,
    day_end_hour=DAY_END_HOUR,
    utc_offset=UTC_OFFSET,
)

def make_save_dict(stacks, nstacks):
    d = dict(**common)
    for v in variants:
        d[f"cc_{v}"]     = stacks[v] / max(nstacks[v], 1) if stacks[v] is not None else np.array([])
        d[f"nstack_{v}"] = nstacks[v]
    return d

np.savez(outfile_day,   **make_save_dict(stacks_day,   nstacks_day))
np.savez(outfile_night, **make_save_dict(stacks_night, nstacks_night))

print(f"Saved: {outfile_day}   (nstack={nstacks_day['combo']})")
print(f"Saved: {outfile_night} (nstack={nstacks_night['combo']})")
print("Done.")
