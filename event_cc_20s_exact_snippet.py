import os
import sys
import glob
import numpy as np
import pandas as pd
import DASutils
from cc_tools import temporal_normalization, computeCC

# ================== USER PARAMETERS ==================
CSV = "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv"
BASE_OUTPUT_DIR = "/oak/stanford/groups/ettore88/nberrios/event_cc"
OUTPUT_VERSION = os.environ.get("OUTPUT_VERSION", "base")

DATA_ROOT_OLD = "/oak/stanford/groups/ettore88/data/SAFODAS1-harddrive-transfer"
DATA_ROOT_NEW = "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer"

ch_start, ch_end = 150, 800
virtual_source_original = 200

fmin, fmax = 5.0, 20.0
max_lag = 1.0

USE_TEMPORAL_NORMALIZATION = os.environ.get("USE_TEMPORAL_NORMALIZATION", "false").lower() == "true"
USE_SECOND_BANDPASS = os.environ.get("USE_SECOND_BANDPASS", "false").lower() == "true"
TN_WINDOW = float(os.environ.get("TN_WINDOW", "0.1"))

SNIPPET_SEC = 20.0
WINDOW_SEC = 20.0
OVERLAP = 0.0

# pass event origin time in UTC, e.g.
# python event_cc_20s_exact_snippet.py 2024-05-21T12:34:56Z EQ001
EVENT_TIME_STR = sys.argv[1]
EVENT_LABEL = sys.argv[2] if len(sys.argv) > 2 else EVENT_TIME_STR.replace(":", "").replace("-", "")
# =====================================================


def normalize_file_path(path):
    path = str(path)
    if os.path.exists(path):
        return path
    if path.startswith(DATA_ROOT_OLD):
        alt = path.replace(DATA_ROOT_OLD, DATA_ROOT_NEW, 1)
        if os.path.exists(alt):
            return alt
    return path


event_time = pd.Timestamp(EVENT_TIME_STR, tz="UTC")
t0 = event_time - pd.Timedelta(seconds=SNIPPET_SEC / 2)
t1 = event_time + pd.Timedelta(seconds=SNIPPET_SEC / 2)

print(f"Event time: {event_time}")
print(f"Snippet: {t0} to {t1}")

db = pd.read_csv(CSV, sep=r"\s+").drop_duplicates()
db = db[db["nSamples"] > 0].reset_index(drop=True)
db["startTime_dt"] = pd.to_datetime(db["startTime"], errors="coerce", utc=True)
db["endTime_dt"] = pd.to_datetime(db["endTime"], errors="coerce", utc=True)
db = db.dropna(subset=["startTime_dt", "endTime_dt"]).reset_index(drop=True)
dur_s = (db["endTime_dt"] - db["startTime_dt"]).dt.total_seconds()
db = db[dur_s > 1.0].reset_index(drop=True)
db["file_norm"] = db["file"].map(normalize_file_path)

# files overlapping the event snippet
sel = db[(db["startTime_dt"] < t1) & (db["endTime_dt"] > t0)].copy()
sel = sel[sel["file_norm"].map(os.path.exists)].copy()
sel = sel.sort_values("startTime_dt").reset_index(drop=True)
selectedFiles = sel["file_norm"].tolist()

if len(selectedFiles) == 0:
    raise ValueError("No files overlap the requested event snippet.")

print(f"Using {len(selectedFiles)} files")

DAS0, info0 = DASutils.readFile_HDF(
    [selectedFiles[0]], 0.05, 24.0, verbose=0,
    preproc=True, diff=True, taper=False,
    desampling=True, nChbuffer=900, system="OptaSense"
)

fs = info0["fs"]
dt = 1.0 / fs

nch = DAS0[ch_start:ch_end, :].shape[0]
ch_buffer_in = nch
isource = virtual_source_original - ch_start
assert 0 <= isource < nch

win_npts = int(WINDOW_SEC / dt)
step_npts = max(1, int(win_npts * (1 - OVERLAP)))
snippet_npts = int(round(SNIPPET_SEC / dt))

print(f"fs={fs} dt={dt} WINDOW_SEC={WINDOW_SEC} TN_WINDOW={TN_WINDOW}")


def build_event_snippet(rows, snippet_start, snippet_end, dt_seconds, n_channels):
    snippet_sum = np.zeros((n_channels, snippet_npts), dtype=np.float64)
    snippet_hits = np.zeros(snippet_npts, dtype=np.int32)

    for _, row in rows.iterrows():
        DAS, info = DASutils.readFile_HDF(
            [row["file_norm"]], 0.05, 24.0, verbose=0,
            preproc=True, diff=True, taper=False,
            desampling=True, nChbuffer=900, system="OptaSense"
        )

        X = DAS[ch_start:ch_end, :].astype(np.float64, copy=False)
        file_start = row["startTime_dt"]
        file_end = row["endTime_dt"]

        overlap_start = max(file_start, snippet_start)
        overlap_end = min(file_end, snippet_end)
        if overlap_end <= overlap_start:
            continue

        src_i0 = int(round((overlap_start - file_start).total_seconds() / dt_seconds))
        dst_i0 = int(round((overlap_start - snippet_start).total_seconds() / dt_seconds))
        nsamp = int(round((overlap_end - overlap_start).total_seconds() / dt_seconds))

        src_i0 = max(src_i0, 0)
        dst_i0 = max(dst_i0, 0)
        nsamp = min(nsamp, X.shape[1] - src_i0, snippet_npts - dst_i0)

        if nsamp <= 0:
            continue

        snippet_sum[:, dst_i0:dst_i0 + nsamp] += X[:, src_i0:src_i0 + nsamp]
        snippet_hits[dst_i0:dst_i0 + nsamp] += 1

    if np.any(snippet_hits == 0):
        missing = int(np.sum(snippet_hits == 0))
        raise ValueError(
            f"Event snippet is not fully covered by the selected files; "
            f"{missing} samples are missing."
        )

    snippet = snippet_sum / snippet_hits[np.newaxis, :]
    return snippet, snippet_hits


X_event, snippet_hits = build_event_snippet(sel, t0, t1, dt, nch)

cc_sum = None
nstack = 0

for i0 in range(0, X_event.shape[1] - win_npts + 1, step_npts):
    Xraw = X_event[:, i0:i0 + win_npts].copy()
    Xcombo = Xraw - np.median(Xraw, axis=0, keepdims=True)
    Xcombo = DASutils.bandpass2D_c(Xcombo, fmin, fmax, dt, zerophase=True)

    if USE_TEMPORAL_NORMALIZATION:
        Xcombo = temporal_normalization(Xcombo.copy(), fs, window_time=TN_WINDOW)

    if USE_SECOND_BANDPASS:
        Xcombo = DASutils.bandpass2D_c(Xcombo, fmin, fmax, dt, zerophase=True)

    cc_i = computeCC(
        Xcombo, dt, max_lag,
        isource=isource,
        ch_buffer_in=ch_buffer_in,
        whitening_params=None
    )

    if cc_sum is None:
        cc_sum = np.zeros_like(cc_i)

    cc_sum += cc_i
    nstack += 1

if cc_sum is None or nstack == 0:
    raise ValueError("No valid CC windows were produced for this event.")

cc_mean = cc_sum / nstack

out_dir = os.path.join(BASE_OUTPUT_DIR, OUTPUT_VERSION)
os.makedirs(out_dir, exist_ok=True)

outfile = os.path.join(out_dir, f"{EVENT_LABEL}_event20s.npz")

np.savez(
    outfile,
    event_time=str(event_time),
    snippet_start=str(t0),
    snippet_end=str(t1),
    output_version=OUTPUT_VERSION,
    fs=fs,
    dt=dt,
    ch_start=ch_start,
    ch_end=ch_end,
    src=virtual_source_original,
    fmin=fmin,
    fmax=fmax,
    max_lag=max_lag,
    snippet_sec=SNIPPET_SEC,
    window_sec=WINDOW_SEC,
    overlap=OVERLAP,
    nfiles_selected=len(selectedFiles),
    snippet_unique_coverage_samples=int(np.sum(snippet_hits > 0)),
    snippet_max_overlap_samples=int(np.max(snippet_hits)),
    use_temporal_normalization=USE_TEMPORAL_NORMALIZATION,
    use_second_bandpass=USE_SECOND_BANDPASS,
    tn_window=TN_WINDOW,
    cc_combo=cc_mean,
    nstack_combo=nstack,
)

print(f"Saved: {outfile} (nstack={nstack})")
