"""
Canonical Lellouch et al. 2019 Figure 7c reproduction cross-correlation.

Faithful to the paper as far as possible given the differences in the dataset:
  - Bandpass:                    5-20 Hz
  - Window:                      30 s, 50% overlap
  - Time normalization:          running absolute mean (Bensen 2007), window 0.1 s
                                 (= half the longest period in 5-20 Hz band)
  - Spectral whitening:          OFF by default — Lellouch does not describe a
                                 whitening step (only running-AM). Optional knob
                                 via SANITY_WHITEN=true if you want to test it.
  - Geometry (a):                virtual source at top of usable channel range by default
                                 receivers along array
  - Stack across windows:        per-channel sum across all 30s windows in selected day(s)

Not done in this script (intentionally — handled in sanity_plot.py):
  - Adjacent-channel pre-shifted stacking (R±10 with 3200 m/s shift). That is a
    post-processing operation on this script's output and does not require re-running
    the heavy CC.

Differences from the paper that we accept:
  - Sampling rate: 500 Hz (vs Lellouch's 2500 Hz). Nyquist 250 Hz, so the 5-20 Hz
    CC band is fully usable.
  - Gauge length: 16.3 m (vs Lellouch's 10 m). In the 5-20 Hz band, the
    finite-gauge sinc response is close to flat for both 3200 m/s body waves
    and 1500 m/s guided/tube waves, so this is not expected to be a first-order
    obstacle to reproduction.

Inputs:
  SAFOD_CSV       : CSV manifest path (default oak path)
  SANITY_DATE     : YYYY-MM-DD (default 2024-10-23)
  SANITY_DATES    : optional comma-separated date list or inclusive range
                    YYYY-MM-DD:YYYY-MM-DD. Overrides SANITY_DATE.
  SANITY_HOURS    : all, day, or night (default all). Lellouch-style reproduction
                    should use all continuous files unless testing source timing.
  SANITY_CH_START : top channel of usable in-well array (default 150)
  SANITY_CH_END   : bottom channel exclusive (default 800)
  SANITY_SOURCE_CH: virtual source channel (default SANITY_CH_START)
  SANITY_OUT      : output directory for the daily npz

Output:
  <SANITY_OUT>/sanity_cc_<DATE_TAG>_ch<START>-<END>_src<SOURCE>.npz
  with cc, lags, channels, and parameters.
"""
import os
import sys
import numpy as np
import pandas as pd
import tqdm

import DASutils
from cc_tools import temporal_normalization, computeCC

# ---------------- USER PARAMETERS ----------------
CSV = os.environ.get(
    "SAFOD_CSV",
    "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv",
)
DATA_ROOT_OLD = "/oak/stanford/groups/ettore88/data/SAFODAS1-harddrive-transfer"
DATA_ROOT_NEW = "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer"

TARGET_DATE = os.environ.get("SANITY_DATE", "2024-10-23")
TARGET_DATES_ENV = os.environ.get("SANITY_DATES", "").strip()
HOUR_MODE   = os.environ.get("SANITY_HOURS", "all").strip().lower()
CH_START    = int(os.environ.get("SANITY_CH_START", "150"))
CH_END      = int(os.environ.get("SANITY_CH_END",   "800"))
SOURCE_CH   = int(os.environ.get("SANITY_SOURCE_CH", str(CH_START)))

# Lellouch-faithful CC parameters
FMIN, FMAX  = 5.0, 20.0
WINDOW_SEC  = 30.0
OVERLAP     = 0.5
MAX_LAG     = 1.0          # seconds; Lellouch shows ±0.3 s but 1.0 gives margin
TN_WINDOW_S = 0.1          # running-AM window = half the longest period (1/5 = 0.2 s)

# Spectral whitening: NOT in the Lellouch paper. They only describe running-AM
# (Bensen 2007). Default OFF for strict reproduction. Set USE_WHITEN=true via
# env var if you want to test with phase-only whitening as an additional knob.
USE_WHITEN    = os.environ.get("SANITY_WHITEN", "false").lower() == "true"
WHITEN_WIN_HZ = 0.0        # 0 -> phase-only whitening when USE_WHITEN is True

# Daytime UTC band for filtering files (10 AM - 8 PM PDT = UTC-7 in October)
# October 2024 is PDT (DST), so UTC offset is -7. Day window in UTC is 17:00 - 03:00 next day.
DAY_UTC_HOURS = set(list(range(17, 24)) + list(range(0, 3)))
ALL_UTC_HOURS = set(range(24))
NIGHT_UTC_HOURS = ALL_UTC_HOURS - DAY_UTC_HOURS

# Bad-channel masking. Channels whose median-of-window RMS is more than this
# many MAD-multiples above the per-day median are zeroed out before CC, so they
# don't contaminate the per-window running-AM normalization. The 2024-10-23
# preflight RMS plot showed a clear ~30x spike near channel 310; this catches
# that and any similar ones automatically.
RMS_OUTLIER_MAD_K = float(os.environ.get("SANITY_BADCH_K", "8.0"))

OUT_DIR = os.environ.get(
    "SANITY_OUT",
    "/oak/stanford/groups/ettore88/nberrios/sanity_v1",
)
os.makedirs(OUT_DIR, exist_ok=True)
# -------------------------------------------------


def normalize_path(p):
    p = str(p)
    if os.path.exists(p):
        return p
    if p.startswith(DATA_ROOT_OLD):
        alt = p.replace(DATA_ROOT_OLD, DATA_ROOT_NEW, 1)
        if os.path.exists(alt):
            return alt
    return p


def parse_target_dates():
    """Return selected UTC date strings."""
    if not TARGET_DATES_ENV:
        return [TARGET_DATE]

    if ":" in TARGET_DATES_ENV:
        start, end = [x.strip() for x in TARGET_DATES_ENV.split(":", 1)]
        dates = pd.date_range(start=start, end=end, freq="D")
        return [d.strftime("%Y-%m-%d") for d in dates]

    return [x.strip() for x in TARGET_DATES_ENV.split(",") if x.strip()]


def make_date_tag(date_list):
    if len(date_list) == 1:
        return date_list[0]
    return f"{date_list[0]}_to_{date_list[-1]}_{len(date_list)}d"


def selected_hour_set(hour_mode):
    """Return selected UTC hours for all/day/night modes."""
    if hour_mode == "all":
        return ALL_UTC_HOURS
    if hour_mode == "day":
        return DAY_UTC_HOURS
    if hour_mode == "night":
        return NIGHT_UTC_HOURS
    raise ValueError("SANITY_HOURS must be one of: all, day, night")


def select_files(csv_path, date_list, hour_mode):
    """Continuous-only files for the target date(s) and hour mode."""
    db = pd.read_csv(csv_path, sep=r"\s+").drop_duplicates()
    db = db[db["nSamples"] == 30000].reset_index(drop=True)  # continuous-only — drop ~20s event triggers
    db["startTime_dt"] = pd.to_datetime(db["startTime"], errors="coerce", utc=True)
    db = db.dropna(subset=["startTime_dt"]).reset_index(drop=True)
    db["date"]      = db["startTime_dt"].dt.strftime("%Y-%m-%d")
    db["hour"]      = db["startTime_dt"].dt.hour
    db["file_norm"] = db["file"].map(normalize_path)
    hours = selected_hour_set(hour_mode)
    sel = db[
        (db["date"].isin(date_list))
        & (db["hour"].isin(hours))
    ].copy()
    sel = sel[sel["file_norm"].map(os.path.exists)].reset_index(drop=True)
    sel = sel.sort_values("startTime_dt").reset_index(drop=True)
    return sel


def main():
    target_dates = parse_target_dates()
    date_tag = make_date_tag(target_dates)
    try:
        selected_hour_set(HOUR_MODE)
    except ValueError as e:
        print(str(e))
        sys.exit(6)
    if not (CH_START <= SOURCE_CH < CH_END):
        print(f"SANITY_SOURCE_CH={SOURCE_CH} must be in [{CH_START}, {CH_END}). Aborting.")
        sys.exit(5)

    print("=== SAFOD DAS canonical Lellouch CC pipeline ===")
    print(f"Dates        : {', '.join(target_dates)}")
    print(f"Date tag     : {date_tag}")
    print(f"Hour mode    : {HOUR_MODE}")
    print(f"Channels     : [{CH_START}, {CH_END})")
    print(f"Source ch    : {SOURCE_CH}")
    print(f"Bandpass     : {FMIN}-{FMAX} Hz")
    print(f"Window       : {WINDOW_SEC}s, overlap {OVERLAP}")
    print(f"Running-AM   : {TN_WINDOW_S} s")
    print(f"Whitening    : {'phase-only' if USE_WHITEN else 'off'}")
    print(f"Output       : {OUT_DIR}")

    files_df = select_files(CSV, target_dates, HOUR_MODE)
    files = files_df["file_norm"].tolist()
    if len(files) == 0:
        print(f"No {HOUR_MODE}-hours continuous files for {target_dates}. Aborting.")
        sys.exit(1)
    print(f"Selected {len(files)} {HOUR_MODE}-hours continuous files for CC.")
    print("Files by date:")
    for d, n in files_df.groupby("date").size().items():
        print(f"  {d}: {n}")

    # Read first file to lock down fs/dt
    DAS0, info0 = DASutils.readFile_HDF(
        [files[0]], 0.05, 24.0, verbose=0,
        preproc=True, diff=True, taper=False,
        desampling=True, nChbuffer=900, system="OptaSense",
    )
    fs = float(info0["fs"])
    dt = 1.0 / fs
    print(f"fs={fs} Hz   dt={dt:.6f} s")

    nch          = CH_END - CH_START
    isource      = SOURCE_CH - CH_START
    win_npts     = int(WINDOW_SEC / dt)
    step_npts    = int(win_npts * (1 - OVERLAP))
    max_lag_npts = int(MAX_LAG / dt)
    n_lag        = 2 * max_lag_npts + 1
    print(f"nch={nch}   win_npts={win_npts}   step_npts={step_npts}   n_lag={n_lag}")

    whitening_params = (fs, WHITEN_WIN_HZ, FMIN, FMAX) if USE_WHITEN else None

    cc_sum  = np.zeros((nch, n_lag), dtype=np.float64)
    n_stack = 0
    n_files_used = 0
    n_files_skipped = 0

    # Pass 1: estimate per-channel RMS over the first few files to identify
    # bad channels. We only need a rough estimate, so a handful of files suffices.
    n_rms_files = min(5, len(files))
    rms_accum = np.zeros(nch, dtype=np.float64)
    rms_count = 0
    for f in files[:n_rms_files]:
        try:
            DAS, _ = DASutils.readFile_HDF(
                [f], 0.05, 24.0, verbose=0,
                preproc=True, diff=True, taper=False,
                desampling=True, nChbuffer=900, system="OptaSense",
            )
            Xrms = DAS[CH_START:CH_END, :].astype(np.float64, copy=False)
            rms_accum += np.sqrt(np.mean(Xrms ** 2, axis=1))
            rms_count += 1
        except Exception as e:
            print(f"  RMS pass: skipping {f}: {e}")
    if rms_count == 0:
        print("Could not estimate RMS for bad-channel masking. Aborting.")
        sys.exit(3)
    rms_per_ch = rms_accum / rms_count
    med  = np.median(rms_per_ch)
    mad  = np.median(np.abs(rms_per_ch - med)) + 1e-30
    bad_mask = rms_per_ch > med + RMS_OUTLIER_MAD_K * 1.4826 * mad
    bad_channels = np.where(bad_mask)[0] + CH_START
    print(f"Bad-channel mask: {bad_mask.sum()} channels flagged "
          f"(threshold {RMS_OUTLIER_MAD_K} MAD above median).")
    if bad_mask.sum() > 0:
        print(f"  Bad channels (absolute idx): {bad_channels.tolist()}")
    if bad_mask[isource]:
        print("  WARNING: source channel flagged as bad. "
              "CC will be unusable — pick a different SANITY_SOURCE_CH or lower SANITY_BADCH_K.")
        sys.exit(4)

    for f in tqdm.tqdm(files):
        try:
            DAS, info = DASutils.readFile_HDF(
                [f], 0.05, 24.0, verbose=0,
                preproc=True, diff=True, taper=False,
                desampling=True, nChbuffer=900, system="OptaSense",
            )
            X = DAS[CH_START:CH_END, :].astype(np.float64, copy=False)
            # Zero out bad channels — they will not contribute to CC since
            # CC of zero with anything is zero. They remain in the array shape
            # so the channel index alignment downstream stays simple.
            if bad_mask.any():
                X[bad_mask, :] = 0.0
            npts = X.shape[1]
            if npts < win_npts:
                n_files_skipped += 1
                continue

            for i0 in range(0, npts - win_npts + 1, step_npts):
                Xw = X[:, i0:i0 + win_npts].copy()

                # 1. Detrend per channel (linear)
                Xw -= Xw.mean(axis=1, keepdims=True)
                # Cosine taper (10% per side) to reduce edge artifacts before bandpass
                ntap = int(0.1 * win_npts)
                taper = np.ones(win_npts)
                taper[:ntap]  = 0.5 * (1 - np.cos(np.pi * np.arange(ntap) / ntap))
                taper[-ntap:] = 0.5 * (1 - np.cos(np.pi * np.arange(ntap)[::-1] / ntap))
                Xw *= taper[None, :]

                # 2. Bandpass 5-20 Hz
                Xw = DASutils.bandpass2D_c(Xw, FMIN, FMAX, dt, zerophase=True)

                # 3. Running absolute mean normalization (Bensen 2007), nwin = TN_WINDOW_S * fs
                Xw = temporal_normalization(Xw, fs, window_time=TN_WINDOW_S)

                # 4. Cross-correlation with phase-only whitening
                cc_i = computeCC(
                    Xw, dt, MAX_LAG,
                    isource=isource,
                    ch_buffer_in=nch,
                    whitening_params=whitening_params,
                )
                cc_sum += cc_i
                n_stack += 1

            n_files_used += 1

        except Exception as e:
            print(f"  Skipping file due to error: {e}")
            n_files_skipped += 1
            continue

    if n_stack == 0:
        print("No CC windows computed. Aborting.")
        sys.exit(2)

    cc_mean = cc_sum / n_stack
    lags    = np.arange(-max_lag_npts, max_lag_npts + 1) * dt
    channels = np.arange(CH_START, CH_END)

    out_name = f"sanity_cc_{date_tag}_hours{HOUR_MODE}_ch{CH_START}-{CH_END}_src{SOURCE_CH}.npz"
    out_path = os.path.join(OUT_DIR, out_name)
    np.savez(
        out_path,
        cc=cc_mean,
        lags=lags,
        channels=channels,
        date=date_tag,
        target_dates=np.array(target_dates),
        hour_mode=HOUR_MODE,
        selected_hours=np.array(sorted(selected_hour_set(HOUR_MODE))),
        ch_start=CH_START, ch_end=CH_END,
        isource_idx=isource, source_channel=SOURCE_CH,
        fs=fs, dt=dt,
        fmin=FMIN, fmax=FMAX,
        window_sec=WINDOW_SEC, overlap=OVERLAP,
        max_lag=MAX_LAG,
        tn_window_s=TN_WINDOW_S,
        use_whiten=USE_WHITEN, whiten_win_hz=WHITEN_WIN_HZ,
        n_stack=n_stack,
        n_files_used=n_files_used,
        n_files_skipped=n_files_skipped,
        bad_channels=bad_channels,                    # absolute channel indices
        bad_channel_threshold_mad=RMS_OUTLIER_MAD_K,
        rms_per_channel=rms_per_ch,
    )
    print(f"Saved {out_path}  (n_stack={n_stack}, files_used={n_files_used}, skipped={n_files_skipped})")
    print("Done. Next: run sanity_plot.py against this npz.")


if __name__ == "__main__":
    main()
