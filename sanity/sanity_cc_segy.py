"""
Lellouch et al. 2019 Figure 7c reproduction — 2017 SEG-Y data variant.

Same CC parameters as sanity_cc.py (5-20 Hz, 30s windows, running-AM, no whitening).
Reads PASSCAL SEG-Y files instead of HDF5.

Environment variables:
  SEGY_DIR        : directory containing .sgy/.gz files (required)
  SANITY_DATE     : YYYY-MM-DD  (default 2017-09-20 — change to your target day)
  SANITY_DATES    : YYYY-MM-DD:YYYY-MM-DD or comma list — overrides SANITY_DATE
  SANITY_HOURS    : all / day / night  (default all — Lellouch uses continuous)
  SANITY_CH_START : first channel index (default 0)
  SANITY_CH_END   : last channel exclusive (default 800)
  SANITY_SOURCE_CH: virtual source channel (default = CH_START)
  SANITY_DIFF     : true/false — differentiate data before CC (default true,
                    set false if data is already strain rate)
  SANITY_DECIMATE : integer decimation factor from raw fs to working fs
                    (default 40 → 2500 Hz → 62.5 Hz; set 1 to skip)
  SANITY_OUT      : output directory for npz
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import tqdm
from scipy.signal import decimate as scipy_decimate

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import DASutils
from cc_tools import temporal_normalization, computeCC


def _segy_nt_from_trace_header(infile):
    """Read nSamples from first trace header bytes 115-116 (SEG-Y standard location)."""
    try:
        with open(infile, 'rb') as f:
            f.seek(3600 + 114)
            return int.from_bytes(f.read(2), byteorder='big', signed=False)
    except Exception:
        return 0


def _segy_ns_ntr_from_filesize(infile, fs, max_seconds=300):
    """Recover (nSamples, nTraces) when header nSamples=0 (OptaSense: >65535 samples).

    Tries multiples of fs (round-second durations) since OptaSense files are
    always an integer number of seconds.  Returns (0, 0) on failure.
    """
    filesize = os.path.getsize(infile)
    data_bytes = filesize - 3600  # subtract text + binary file headers
    fs_int = int(round(fs))
    for n_sec in range(1, max_seconds + 1):
        ns = fs_int * n_sec
        record_size = 240 + ns * 4
        if data_bytes % record_size == 0:
            ntr = data_bytes // record_size
            if 1 <= ntr <= 5000:
                return ns, ntr
    return 0, 0


def _segy_end_time_from_samples(start_time, n_samples, fs):
    if fs <= 0 or n_samples <= 0:
        return start_time
    return start_time + pd.to_timedelta(float(n_samples) / float(fs), unit="s")

# ── parameters ────────────────────────────────────────────────────────────────
SEGY_DIR    = os.environ.get("SEGY_DIR", "")
TARGET_DATE = os.environ.get("SANITY_DATE",  "2017-09-20")
TARGET_DATES_ENV = os.environ.get("SANITY_DATES", "").strip()
HOUR_MODE   = os.environ.get("SANITY_HOURS", "all").strip().lower()
CH_START    = int(os.environ.get("SANITY_CH_START",  "0"))
CH_END      = int(os.environ.get("SANITY_CH_END",    "800"))
SOURCE_CH   = int(os.environ.get("SANITY_SOURCE_CH", str(CH_START)))
DO_DIFF     = os.environ.get("SANITY_DIFF",     "true").lower() == "true"
DECIMATE_Q  = int(os.environ.get("SANITY_DECIMATE", "40"))   # 2500 → 62.5 Hz
SEGY_SUFFIX = os.environ.get("SANITY_FILE_SUFFIX", "").strip()  # e.g. ".814+0000.sgy"
OUT_DIR     = os.environ.get("SANITY_OUT",
              "/oak/stanford/groups/ettore88/nberrios/sanity_segy")

# Lellouch-faithful CC parameters (identical to sanity_cc.py)
FMIN, FMAX   = 5.0, 20.0
WINDOW_SEC   = 30.0
OVERLAP      = 0.5
MAX_LAG      = 1.0
TN_WINDOW_S  = 0.1
RMS_MAD_K    = float(os.environ.get("SANITY_BADCH_K", "8.0"))

# Day/night hour sets (UTC). Oct 2024 PDT offset = −7 h.
# Adjust if 2017 data spans a different season.
DAY_UTC  = set(list(range(17, 24)) + list(range(0, 3)))
ALL_UTC  = set(range(24))
NIGHT_UTC = ALL_UTC - DAY_UTC
# ──────────────────────────────────────────────────────────────────────────────


def parse_target_dates():
    if not TARGET_DATES_ENV:
        return [TARGET_DATE]
    if ":" in TARGET_DATES_ENV:
        s, e = TARGET_DATES_ENV.split(":", 1)
        return [d.strftime("%Y-%m-%d") for d in pd.date_range(s.strip(), e.strip(), freq="D")]
    return [x.strip() for x in TARGET_DATES_ENV.split(",") if x.strip()]


def hour_set(mode):
    if mode == "all":   return ALL_UTC
    if mode == "day":   return DAY_UTC
    if mode == "night": return NIGHT_UTC
    raise ValueError(f"SANITY_HOURS must be all/day/night, got '{mode}'")


def build_manifest(segy_dir):
    """
    Scan directory for .sgy / .sgy.gz files, read headers, return DataFrame.
    Columns: file, startTime, endTime, fs, nTraces, nSamples
    """
    patterns = ["*.sgy", "*.sgy.gz", "*.SGY"]
    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(segy_dir, "**", pat), recursive=True))
        files.extend(glob.glob(os.path.join(segy_dir, pat)))
    files = sorted(set(files))
    if SEGY_SUFFIX:
        files = [f for f in files if os.path.basename(f).endswith(SEGY_SUFFIX)]

    if not files:
        print(f"No SEG-Y files found in {segy_dir}")
        sys.exit(1)

    rows = []
    print(f"Scanning {len(files)} SEG-Y files for headers...")
    for f in tqdm.tqdm(files):
        try:
            nt, fs, t0, t1, ntr = DASutils.read_PASSCAL_SEGY_headers(f)
            if nt == 0:
                # nSamples exceeds uint16 max (OptaSense >65535 samples) — recover from file size
                nt, ntr = _segy_ns_ntr_from_filesize(f, fs if fs > 0 else 2500.0)
                t1 = _segy_end_time_from_samples(t0, nt, fs if fs > 0 else 2500.0)
            rows.append(dict(file=f, startTime=t0, endTime=t1,
                             fs=fs, nTraces=ntr, nSamples=nt))
        except Exception as e:
            print(f"  header error {os.path.basename(f)}: {e}")
    return pd.DataFrame(rows)


def select_files(manifest, date_list, hour_mode):
    hrs = hour_set(hour_mode)
    manifest = manifest.copy()
    manifest["date"] = manifest["startTime"].dt.strftime("%Y-%m-%d")
    manifest["hour"] = manifest["startTime"].dt.hour
    sel = manifest[manifest["date"].isin(date_list) & manifest["hour"].isin(hrs)]
    return sel.sort_values("startTime").reset_index(drop=True)


def load_segy(row, ch_start, ch_end):
    """
    Read one SEG-Y file, return float64 array [nch, nsamples] for channels
    [ch_start, ch_end). Applies detrend only — caller handles diff/decimate.
    """
    ntr = int(row["nTraces"])
    ns  = int(row["nSamples"])
    raw = DASutils.read_PASSCAL_segy(row["file"], ntr, ns)   # (nTraces, nSamples)
    raw = raw.astype(np.float64)

    # Trim to requested channel range
    actual_end = min(ch_end, raw.shape[0])
    if ch_start >= raw.shape[0]:
        raise ValueError(f"ch_start={ch_start} >= nTraces={raw.shape[0]}")
    raw = raw[ch_start:actual_end, :]

    # Linear detrend per channel
    raw -= raw.mean(axis=1, keepdims=True)
    t = np.arange(raw.shape[1], dtype=np.float64)
    t -= t.mean()
    for i in range(raw.shape[0]):
        raw[i] -= t * (np.dot(raw[i], t) / np.dot(t, t))
    return raw


def main():
    if not SEGY_DIR:
        print("Set SEGY_DIR to the directory containing 2017 SEG-Y files.")
        sys.exit(1)

    os.makedirs(OUT_DIR, exist_ok=True)
    target_dates = parse_target_dates()
    date_tag = target_dates[0] if len(target_dates) == 1 else \
               f"{target_dates[0]}_to_{target_dates[-1]}_{len(target_dates)}d"

    print("=== SAFOD DAS Lellouch CC — 2017 SEG-Y variant ===")
    print(f"  SEGY_DIR    : {SEGY_DIR}")
    print(f"  Dates       : {', '.join(target_dates)}")
    print(f"  Hours       : {HOUR_MODE}")
    print(f"  Channels    : [{CH_START}, {CH_END})")
    print(f"  Source ch   : {SOURCE_CH}")
    print(f"  Diff        : {DO_DIFF}")
    print(f"  Decimate    : {DECIMATE_Q}x")
    print(f"  Output      : {OUT_DIR}")

    manifest = build_manifest(SEGY_DIR)
    if manifest.empty:
        print("Empty manifest — no valid SEG-Y files. Aborting.")
        sys.exit(1)

    sel = select_files(manifest, target_dates, HOUR_MODE)
    if sel.empty:
        print(f"No files match dates={target_dates} hours={HOUR_MODE}. Aborting.")
        sys.exit(1)
    print(f"Selected {len(sel)} file(s) for processing.")

    # Diagnostic: show first file dimensions
    row0 = sel.iloc[0]
    print(f"  First file: nSamples={int(row0['nSamples'])}, nTraces={int(row0['nTraces'])}, "
          f"start={row0['startTime']}, end={row0['endTime']}")

    # Determine working fs from first file
    raw_fs  = float(sel.iloc[0]["fs"])
    work_fs = raw_fs / DECIMATE_Q if DECIMATE_Q > 1 else raw_fs
    dt      = 1.0 / work_fs
    print(f"  Raw fs={raw_fs} Hz  →  working fs={work_fs:.2f} Hz  dt={dt:.6f} s")

    # Minimum samples needed for scipy_decimate (padlen=27 at raw fs)
    MIN_SAMPLES_RAW = 50  # well above padlen=27

    # actual channel count may be less than CH_END - CH_START when the array
    # has fewer traces than CH_END (e.g. 118-channel OptaSense file with CH_END=800)
    actual_nch   = min(CH_END, int(sel.iloc[0]["nTraces"])) - CH_START
    isource      = SOURCE_CH - CH_START
    win_npts     = int(WINDOW_SEC / dt)
    step_npts    = int(win_npts * (1 - OVERLAP))
    max_lag_npts = int(MAX_LAG / dt)
    n_lag        = 2 * max_lag_npts + 1
    print(f"  Actual channels: {actual_nch}  (requested {CH_END - CH_START})")

    # ── bad-channel estimation from first few files ───────────────────────────
    rms_accum = np.zeros(actual_nch)
    rms_count = 0
    for _, row in sel.iloc[:min(5, len(sel))].iterrows():
        try:
            X = load_segy(row, CH_START, CH_END)
            print(f"  RMS scan file shape: {X.shape}  ({row['file'].split('/')[-1]})")
            if X.shape[1] < MIN_SAMPLES_RAW:
                print(f"  RMS scan skip: only {X.shape[1]} samples (need >{MIN_SAMPLES_RAW})")
                continue
            if DO_DIFF:
                X = DASutils.preprocess_diff(X, 1.0 / raw_fs)
            if DECIMATE_Q > 1:
                X = scipy_decimate(X, DECIMATE_Q, axis=1, zero_phase=True)
            rms_accum += np.sqrt(np.mean(X ** 2, axis=1))
            rms_count += 1
        except Exception as e:
            print(f"  RMS scan skip: {e}")

    if rms_count == 0:
        print("Could not estimate RMS. Aborting.")
        sys.exit(3)

    rms_per_ch = rms_accum / rms_count
    med = np.median(rms_per_ch)
    mad = np.median(np.abs(rms_per_ch - med)) + 1e-30
    bad_mask = rms_per_ch > med + RMS_MAD_K * 1.4826 * mad
    bad_chs  = np.where(bad_mask)[0] + CH_START
    print(f"Bad channels flagged: {bad_mask.sum()} — {bad_chs.tolist()}")
    if bad_mask[isource]:
        print("ERROR: source channel is bad. Choose a different SANITY_SOURCE_CH.")
        sys.exit(4)

    # ── main CC loop ──────────────────────────────────────────────────────────
    cc_sum  = np.zeros((actual_nch, n_lag), dtype=np.float64)
    n_stack = 0
    n_used  = 0
    n_skip  = 0

    for _, row in tqdm.tqdm(sel.iterrows(), total=len(sel)):
        try:
            X = load_segy(row, CH_START, CH_END)

            if DO_DIFF:
                X = DASutils.preprocess_diff(X, 1.0 / raw_fs)

            if DECIMATE_Q > 1:
                X = scipy_decimate(X, DECIMATE_Q, axis=1, zero_phase=True)

            if bad_mask.any():
                X[bad_mask, :] = 0.0

            npts = X.shape[1]
            if npts < win_npts:
                n_skip += 1
                continue

            for i0 in range(0, npts - win_npts + 1, step_npts):
                Xw = X[:, i0:i0 + win_npts].copy()
                Xw -= Xw.mean(axis=1, keepdims=True)

                # Cosine taper 10% per side
                ntap = int(0.1 * win_npts)
                taper = np.ones(win_npts)
                taper[:ntap]  = 0.5 * (1 - np.cos(np.pi * np.arange(ntap) / ntap))
                taper[-ntap:] = 0.5 * (1 - np.cos(np.pi * np.arange(ntap)[::-1] / ntap))
                Xw *= taper[None, :]

                Xw = DASutils.bandpass2D_c(Xw, FMIN, FMAX, dt, zerophase=True)
                Xw = temporal_normalization(Xw, work_fs, window_time=TN_WINDOW_S)

                cc_sum += computeCC(Xw, dt, MAX_LAG, isource=isource,
                                    ch_buffer_in=actual_nch, whitening_params=None)
                n_stack += 1

            n_used += 1

        except Exception as e:
            print(f"  Skip {os.path.basename(row['file'])}: {e}")
            n_skip += 1

    if n_stack == 0:
        print("No CC windows computed. Aborting.")
        sys.exit(2)

    cc_mean  = cc_sum / n_stack
    lags     = np.arange(-max_lag_npts, max_lag_npts + 1) * dt
    channels = np.arange(CH_START, CH_START + actual_nch)

    out_name = (f"sanity_cc_segy_{date_tag}_hours{HOUR_MODE}"
                f"_ch{CH_START}-{CH_END}_src{SOURCE_CH}.npz")
    out_path = os.path.join(OUT_DIR, out_name)
    np.savez(
        out_path,
        cc=cc_mean, lags=lags, channels=channels,
        date=date_tag, target_dates=np.array(target_dates),
        hour_mode=HOUR_MODE, ch_start=CH_START, ch_end=CH_END,
        isource_idx=isource, source_channel=SOURCE_CH,
        fs=work_fs, dt=dt, raw_fs=raw_fs, decimate_q=DECIMATE_Q,
        do_diff=DO_DIFF,
        fmin=FMIN, fmax=FMAX, window_sec=WINDOW_SEC, overlap=OVERLAP,
        max_lag=MAX_LAG, tn_window_s=TN_WINDOW_S,
        n_stack=n_stack, n_files_used=n_used, n_files_skipped=n_skip,
        bad_channels=bad_chs, bad_channel_threshold_mad=RMS_MAD_K,
        rms_per_channel=rms_per_ch,
    )
    print(f"\nSaved {out_path}")
    print(f"  n_stack={n_stack}  files_used={n_used}  skipped={n_skip}")
    print("Next: run sanity_plot.py against this npz.")


if __name__ == "__main__":
    main()
