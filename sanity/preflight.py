"""
Pre-flight diagnostic for the SAFOD DAS sanity reproduction of Lellouch et al. 2019 Fig 7c.

Goal: BEFORE running cross-correlation on a full day of data, confirm:
  1. What sample rate the data is actually delivered at after readFile_HDF.
  2. Whether body-wave-velocity energy (~3000 m/s apparent vertical velocity)
     exists in the wavefield, or whether tube/casing modes (~1500 m/s) and
     surface modes (<1000 m/s) dominate.
  3. Where the uphole/in-well boundary is (RMS-per-channel transition).
  4. The frequency content of the wavefield in the CC band of interest (5-20 Hz).

If body-wave energy is not present in the F-K plot, no preprocessing trick
will conjure it through cross-correlation. This script answers that question
in ~minutes instead of after another week of CC tuning.

Inputs:
  - SAFOD CSV manifest (default: oak Stanford path, override via SAFOD_CSV env var)
  - Target date and UTC hour range (default: 2024-10-23, 20:00 UTC = 1pm PDT)

Outputs (PNG and one .npz with the F-K matrix), under sanity/preflight_out/:
  - rms_per_channel.png : log-RMS per channel, with ch_start guide
  - spectrogram.png     : spectrogram of one mid-array channel
  - timedist_one_file.png : raw strain rate, time vs channel, one 60s file
  - fk_one_file.png     : 2D F-K spectrum with reference velocity lines
  - preflight_arrays.npz : the F-K matrix and supporting arrays for re-plotting
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import DASutils

# ---------------- USER PARAMETERS ----------------
CSV = os.environ.get(
    "SAFOD_CSV",
    "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv",
)
DATA_ROOT_OLD = "/oak/stanford/groups/ettore88/data/SAFODAS1-harddrive-transfer"
DATA_ROOT_NEW = "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer"

TARGET_DATE   = os.environ.get("SANITY_DATE", "2024-10-23")  # zero EQs in catalog week
HOUR_UTC_LO   = int(os.environ.get("SANITY_HOUR_LO", "20"))  # 20:00 UTC = 1 PM PDT
HOUR_UTC_HI   = int(os.environ.get("SANITY_HOUR_HI", "21"))  # one hour
N_FILES_MAX   = int(os.environ.get("SANITY_NFILES", "60"))   # 60 x 60s = 1 hour

CH_START = int(os.environ.get("SANITY_CH_START", "150"))
CH_END   = int(os.environ.get("SANITY_CH_END",   "800"))

# Bandpass on read — wide, encompasses the 5-20 Hz CC band
READ_FMIN = 0.05
READ_FMAX = 24.0

OUT_DIR = os.environ.get(
    "SANITY_PREFLIGHT_OUT",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "preflight_out"),
)
os.makedirs(OUT_DIR, exist_ok=True)

# Reference apparent velocities to overlay on F-K (m/s along fiber)
REF_VELS = {
    "P body (3200 m/s)":  3200.0,
    "S body (1600 m/s)":  1600.0,
    "Tube (~1500 m/s)":   1500.0,
    "Surface (~500 m/s)": 500.0,
}
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


def select_files(csv_path, date_str, hour_lo, hour_hi, n_max):
    """Pick continuous-only (nSamples == 30000) daytime files for one date/hour band."""
    db = pd.read_csv(csv_path, sep=r"\s+").drop_duplicates()
    db = db[db["nSamples"] == 30000].reset_index(drop=True)  # continuous-only
    db["startTime_dt"] = pd.to_datetime(db["startTime"], errors="coerce", utc=True)
    db = db.dropna(subset=["startTime_dt"]).reset_index(drop=True)
    db["date"]    = db["startTime_dt"].dt.strftime("%Y-%m-%d")
    db["hour"]    = db["startTime_dt"].dt.hour
    db["file_norm"] = db["file"].map(normalize_path)

    sel = db[
        (db["date"] == date_str)
        & (db["hour"] >= hour_lo)
        & (db["hour"] < hour_hi)
    ].copy()
    sel = sel[sel["file_norm"].map(os.path.exists)].reset_index(drop=True)
    if len(sel) == 0:
        raise SystemExit(
            f"No continuous (nSamples==30000) files found for {date_str} "
            f"in UTC hours [{hour_lo}, {hour_hi}). Check date/hour or CSV."
        )
    return sel.head(n_max)


def read_concatenated(files):
    """Read each file with the same flags as stack_daily.py and concatenate in time."""
    parts = []
    info_first = None
    for f in files:
        DAS, info = DASutils.readFile_HDF(
            [f], READ_FMIN, READ_FMAX, verbose=0,
            preproc=True, diff=True, taper=False,
            desampling=True, nChbuffer=900, system="OptaSense",
        )
        if info_first is None:
            info_first = info
        parts.append(DAS[CH_START:CH_END, :].astype(np.float32, copy=False))
    X = np.concatenate(parts, axis=1)
    return X, info_first


def plot_rms(X, fs, out_png):
    rms = np.sqrt(np.mean(X.astype(np.float64) ** 2, axis=1))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogy(np.arange(CH_START, CH_END), rms, lw=0.6)
    ax.set_xlabel("Channel index")
    ax.set_ylabel("RMS (log scale)")
    ax.set_title(
        f"RMS per channel — {TARGET_DATE} UTC {HOUR_UTC_LO:02d}:00–{HOUR_UTC_HI:02d}:00, "
        f"fs={fs:.0f} Hz"
    )
    ax.axvline(CH_START, color="r", ls="--", lw=0.8, label=f"ch_start={CH_START}")
    ax.legend(loc="best")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_spectrogram(X, fs, out_png):
    """Spectrogram of a mid-array channel."""
    from scipy.signal import spectrogram
    nch = X.shape[0]
    ich = nch // 2
    f, t, Sxx = spectrogram(
        X[ich, :].astype(np.float64), fs=fs,
        nperseg=int(8 * fs), noverlap=int(4 * fs),
        scaling="density",
    )
    fig, ax = plt.subplots(figsize=(9, 4))
    Sdb = 10.0 * np.log10(np.maximum(Sxx, 1e-30))
    pcm = ax.pcolormesh(t, f, Sdb, shading="auto", cmap="magma",
                        vmin=np.percentile(Sdb, 5), vmax=np.percentile(Sdb, 99))
    ax.set_ylim(0, min(50, fs / 2))
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title(
        f"Spectrogram — channel {CH_START + ich} (mid-array), "
        f"{TARGET_DATE} {HOUR_UTC_LO:02d}–{HOUR_UTC_HI:02d} UTC"
    )
    ax.axhspan(5, 20, color="cyan", alpha=0.15, label="CC band 5–20 Hz")
    ax.legend(loc="upper right")
    fig.colorbar(pcm, ax=ax, label="Power (dB)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_timedist_one_file(X_one, fs, out_png):
    """Plot one 60-s file (or whatever fits) raw strain rate as time vs channel."""
    nch, npts = X_one.shape
    t = np.arange(npts) / fs
    chans = np.arange(CH_START, CH_END)
    vlim = np.percentile(np.abs(X_one), 98)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.imshow(
        X_one, aspect="auto", origin="upper", cmap="RdBu_r",
        vmin=-vlim, vmax=vlim,
        extent=[t[0], t[-1], chans[-1], chans[0]],
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Channel")
    ax.set_title(
        f"Raw strain rate, one 60-s file — {TARGET_DATE} {HOUR_UTC_LO:02d} UTC"
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def compute_fk(X, fs, dz):
    """2D FFT of the time-channel matrix. Returns (f, k, |F-K|^2)."""
    nch, npts = X.shape
    # Window in space and time to reduce edge artifacts
    win_t = np.hanning(npts)[None, :]
    win_x = np.hanning(nch)[:, None]
    Xw = X * win_t * win_x

    # 2D FFT: time -> f (along axis 1), channel -> k (along axis 0)
    F = np.fft.rfft(Xw, axis=1)            # shape (nch, nf)
    FK = np.fft.fftshift(np.fft.fft(F, axis=0), axes=0)  # shape (nch, nf)
    FK_pow = np.abs(FK) ** 2

    f = np.fft.rfftfreq(npts, d=1.0 / fs)              # >= 0
    k = np.fft.fftshift(np.fft.fftfreq(nch, d=dz))     # cycles / m
    return f, k, FK_pow


def plot_fk(f, k, FK, out_png):
    """Log-power F-K with reference apparent-velocity lines.

    Display has k on the x-axis and f on the y-axis, so a line through the origin
    of slope v_app = f / k shows the apparent along-fiber velocity. Body-wave
    energy appears as steep (near-vertical) ridges; tube/surface modes as shallow
    (near-horizontal) ridges.
    """
    f_max_disp = 25.0
    fmask = f <= f_max_disp
    f_disp = f[fmask]
    FK_disp = FK[:, fmask]
    P = 10.0 * np.log10(np.maximum(FK_disp, 1e-30))

    fig, ax = plt.subplots(figsize=(8, 6))
    pcm = ax.pcolormesh(
        k, f_disp, P.T, shading="auto", cmap="magma",
        vmin=np.percentile(P, 70), vmax=np.percentile(P, 99.5),
    )
    for label, v in REF_VELS.items():
        f_line = np.linspace(0.1, f_max_disp, 200)
        k_line = f_line / v
        ax.plot( k_line, f_line, "--", lw=1.0, label=label)
        ax.plot(-k_line, f_line, "--", lw=1.0, color=ax.lines[-1].get_color())
    ax.set_xlabel("k (cycles / m)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_xlim(-0.05, 0.05)
    ax.set_ylim(0, f_max_disp)
    ax.axhspan(5, 20, color="cyan", alpha=0.15)
    ax.set_title(
        f"F-K spectrum — {TARGET_DATE} {HOUR_UTC_LO:02d}–{HOUR_UTC_HI:02d} UTC, "
        f"channels {CH_START}–{CH_END}\n"
        "Body-wave energy = lines near vertical (steep slope = high apparent velocity)"
    )
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(pcm, ax=ax, label="Power (dB)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    print(f"=== SAFOD DAS pre-flight diagnostic ===")
    print(f"Date:    {TARGET_DATE}")
    print(f"UTC hrs: {HOUR_UTC_LO:02d}:00 to {HOUR_UTC_HI:02d}:00")
    print(f"Channels: [{CH_START}, {CH_END})")
    print(f"Output:  {OUT_DIR}")

    files_df = select_files(CSV, TARGET_DATE, HOUR_UTC_LO, HOUR_UTC_HI, N_FILES_MAX)
    files = files_df["file_norm"].tolist()
    print(f"Selected {len(files)} continuous files.")

    # Read all selected files concatenated for stable spectrogram & RMS,
    # but keep one separately for the time-distance plot.
    DAS0, info0 = DASutils.readFile_HDF(
        [files[0]], READ_FMIN, READ_FMAX, verbose=0,
        preproc=True, diff=True, taper=False,
        desampling=True, nChbuffer=900, system="OptaSense",
    )
    fs0 = float(info0["fs"])
    print(f"After readFile_HDF: fs = {fs0:.4f} Hz, dt = {1.0/fs0:.6f} s")
    print(f"Native acquisition was 10000 Hz × 1/Desample(20) = 500 Hz stored.")
    print(f"If fs above is not 500, desampling is doing something extra — investigate.")
    X_one = DAS0[CH_START:CH_END, :].astype(np.float32, copy=False)

    print("Concatenating all selected files...")
    X, info = read_concatenated(files)
    fs = float(info["fs"])
    nch, npts = X.shape
    dz = float(info.get("dCh", 1.0)) if "dCh" in info else 1.0
    # Fallback: read dCh from CSV if not in info
    if "dCh" not in info:
        dz = float(files_df["dCh"].iloc[0])
    print(f"Concatenated shape: nch={nch}, npts={npts}, fs={fs}, dz={dz}")

    # 1. RMS per channel
    plot_rms(X, fs, os.path.join(OUT_DIR, "rms_per_channel.png"))
    print("Wrote rms_per_channel.png")

    # 2. Spectrogram (mid-array channel)
    plot_spectrogram(X, fs, os.path.join(OUT_DIR, "spectrogram.png"))
    print("Wrote spectrogram.png")

    # 3. Time-distance plot of one file
    plot_timedist_one_file(X_one, fs, os.path.join(OUT_DIR, "timedist_one_file.png"))
    print("Wrote timedist_one_file.png")

    # 4. F-K analysis on concatenated data
    print("Computing F-K spectrum...")
    f, k, FK = compute_fk(X.astype(np.float64), fs, dz)
    plot_fk(f, k, FK, os.path.join(OUT_DIR, "fk_one_file.png"))
    print("Wrote fk_one_file.png")

    # Save arrays for re-plotting
    np.savez(
        os.path.join(OUT_DIR, "preflight_arrays.npz"),
        date=TARGET_DATE, hour_lo=HOUR_UTC_LO, hour_hi=HOUR_UTC_HI,
        ch_start=CH_START, ch_end=CH_END,
        fs=fs, dz=dz,
        f=f, k=k, FK_pow=FK,
    )
    print("Wrote preflight_arrays.npz")
    print("Done.")


if __name__ == "__main__":
    main()
