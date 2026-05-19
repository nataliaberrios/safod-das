"""
F-K diagnostic comparing one 60-s file with a M >= 2 earthquake vs one
60-s file from a known-quiet period on the SAFOD 2024-2025 DAS data.

Purpose:
  Just visualize what the raw wavefield looks like in (f, k) space during
  an earthquake vs during ambient-noise-only conditions. No filtering yet —
  this is the diagnostic that tells us whether filtering would be useful
  and what velocity range to filter.

Earthquake selection:
  Query NCEDC for events near Parkfield (35.9747 N, 120.5522 W) within a
  30-km radius in 2024, M >= 2 (raise floor to ensure a clear signature).
  Match against SAFOD_2024_2025.csv to find a continuous 60-s file whose
  [startTime, endTime] window contains the event's origin time. First match
  wins.

Quiet selection:
  2024-10-23 around 20:00 UTC. The week 2024-10-22 to 2024-10-28 had zero
  M >= 1 events near Parkfield in the prior notebook query.

Outputs (under sanity/fk_compare_out/):
  - fk_compare_eq_quiet.png : 2x2 panel (time-distance + F-K for each case)
  - fk_compare_arrays.npz   : F-K matrices and metadata for re-plotting
  - text summary printed to stdout with the selected event details

Run on Sherlock with the same env vars as preflight (DASutils on PYTHONPATH).
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import DASutils
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

# ---------------- USER PARAMETERS ----------------
CSV = os.environ.get(
    "SAFOD_CSV",
    "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv",
)
DATA_ROOT_OLD = "/oak/stanford/groups/ettore88/data/SAFODAS1-harddrive-transfer"
DATA_ROOT_NEW = "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer"

# Channel range used by sanity_cc and preflight — gives a clean in-well slice.
CH_START = int(os.environ.get("SANITY_CH_START", "150"))
CH_END   = int(os.environ.get("SANITY_CH_END", "800"))

# Earthquake catalog parameters (mirrors notebooks/13_stacks_n_seis.ipynb).
CATALOG       = "NCEDC"
PARKFIELD_LAT = 35.9747
PARKFIELD_LON = -120.5522
RADIUS_KM     = 30.0
MIN_MAG       = 2.0
EQ_YEAR_START = "2024-01-01"
EQ_YEAR_END   = "2024-12-31"

# Quiet-window fallback (used if we can't find an EQ that matches the CSV).
QUIET_DATE  = "2024-10-23"
QUIET_HOUR  = 20

# Bandpass on read — same as preflight/sanity_cc.
READ_FMIN = 0.05
READ_FMAX = 24.0

OUT_DIR = os.environ.get(
    "FK_COMPARE_OUT",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "fk_compare_out"),
)
os.makedirs(OUT_DIR, exist_ok=True)

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


def load_manifest(csv_path):
    db = pd.read_csv(csv_path, sep=r"\s+").drop_duplicates()
    db = db[db["nSamples"] == 30000].reset_index(drop=True)
    db["startTime_dt"] = pd.to_datetime(db["startTime"], errors="coerce", utc=True)
    db["endTime_dt"]   = pd.to_datetime(db["endTime"],   errors="coerce", utc=True)
    db = db.dropna(subset=["startTime_dt", "endTime_dt"]).reset_index(drop=True)
    db["file_norm"] = db["file"].map(normalize_path)
    db = db[db["file_norm"].map(os.path.exists)].reset_index(drop=True)
    return db


def find_eq_with_file(db):
    """Query NCEDC for M>=MIN_MAG events near Parkfield in 2024.
       Return (event_dict, csv_row) for the first event whose origin time
       lies inside a continuous 60-s file window."""
    print(f"Querying {CATALOG} for events near Parkfield {EQ_YEAR_START}..{EQ_YEAR_END}, "
          f"M>={MIN_MAG}, within {RADIUS_KM} km...")
    client = Client(CATALOG)
    catalog = client.get_events(
        starttime=pd.Timestamp(EQ_YEAR_START).to_pydatetime(),
        endtime=pd.Timestamp(EQ_YEAR_END).to_pydatetime(),
        minmagnitude=MIN_MAG,
        latitude=PARKFIELD_LAT,
        longitude=PARKFIELD_LON,
        maxradius=RADIUS_KM / 111.32,
        includearrivals=False,
    )
    print(f"  catalog returned {len(catalog)} events")
    # Iterate in descending magnitude order so we prefer the clearest signal.
    events = []
    for ev in catalog:
        origin = ev.preferred_origin() or ev.origins[0]
        mag    = ev.preferred_magnitude() or ev.magnitudes[0]
        dist_m, _, _ = gps2dist_azimuth(
            PARKFIELD_LAT, PARKFIELD_LON, origin.latitude, origin.longitude,
        )
        events.append({
            "time_utc":  pd.Timestamp(origin.time.datetime, tz="UTC"),
            "lat":       origin.latitude,
            "lon":       origin.longitude,
            "depth_km":  origin.depth / 1000.0 if origin.depth else None,
            "mag":       mag.mag,
            "mag_type":  mag.magnitude_type,
            "dist_km":   dist_m / 1000.0,
        })
    events.sort(key=lambda e: -e["mag"])
    for ev in events:
        if ev["dist_km"] > RADIUS_KM:
            continue
        match = db[
            (db["startTime_dt"] <= ev["time_utc"])
            & (db["endTime_dt"] >= ev["time_utc"])
        ]
        if len(match) > 0:
            return ev, match.iloc[0]
    return None, None


def find_quiet_file(db):
    sel = db[
        (db["startTime_dt"].dt.strftime("%Y-%m-%d") == QUIET_DATE)
        & (db["startTime_dt"].dt.hour == QUIET_HOUR)
    ]
    if len(sel) == 0:
        raise SystemExit(
            f"No continuous file found for {QUIET_DATE} UTC {QUIET_HOUR:02d}:00."
        )
    return sel.iloc[0]


def read_file(path):
    DAS, info = DASutils.readFile_HDF(
        [path], READ_FMIN, READ_FMAX, verbose=0,
        preproc=True, diff=True, taper=False,
        desampling=True, nChbuffer=900, system="OptaSense",
    )
    return DAS[CH_START:CH_END, :].astype(np.float64, copy=False), info


def compute_fk(X, fs, dz):
    nch, npts = X.shape
    win_t = np.hanning(npts)[None, :]
    win_x = np.hanning(nch)[:, None]
    Xw = X * win_t * win_x
    F  = np.fft.rfft(Xw, axis=1)
    FK = np.fft.fftshift(np.fft.fft(F, axis=0), axes=0)
    f  = np.fft.rfftfreq(npts, d=1.0 / fs)
    k  = np.fft.fftshift(np.fft.fftfreq(nch, d=dz))
    return f, k, np.abs(FK) ** 2


def plot_timedist(ax, X, fs, title):
    nch, npts = X.shape
    t = np.arange(npts) / fs
    vlim = np.percentile(np.abs(X), 98)
    ax.imshow(
        X, aspect="auto", origin="upper", cmap="RdBu_r",
        vmin=-vlim, vmax=vlim,
        extent=[t[0], t[-1], CH_END - 1, CH_START],
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Channel")
    ax.set_title(title, fontsize=10)


def plot_fk(ax, f, k, FK, title):
    f_max_disp = 25.0
    fmask = f <= f_max_disp
    f_d   = f[fmask]
    FK_d  = FK[:, fmask]
    P = 10.0 * np.log10(np.maximum(FK_d, 1e-30))
    pcm = ax.pcolormesh(
        k, f_d, P.T, shading="auto", cmap="magma",
        vmin=np.percentile(P, 70), vmax=np.percentile(P, 99.5),
    )
    for label, v in REF_VELS.items():
        f_line = np.linspace(0.1, f_max_disp, 200)
        k_line = f_line / v
        ax.plot( k_line, f_line, "--", lw=0.8, label=label)
        ax.plot(-k_line, f_line, "--", lw=0.8, color=ax.lines[-1].get_color())
    ax.set_xlim(-0.05, 0.05)
    ax.set_ylim(0, f_max_disp)
    ax.axhspan(5, 20, color="cyan", alpha=0.15)
    ax.set_xlabel("k (cycles / m)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title(title, fontsize=10)
    ax.legend(loc="upper right", fontsize=7)
    return pcm


def main():
    print("=== F-K compare: earthquake-window vs quiet-window ===")
    print(f"Channels: [{CH_START}, {CH_END})")
    print(f"Output dir: {OUT_DIR}")

    db = load_manifest(CSV)
    print(f"Loaded {len(db)} continuous (nSamples=30000) files from manifest.")

    ev, eq_row = find_eq_with_file(db)
    if ev is None:
        print("No earthquake matched a continuous file in 2024 with the current "
              "settings. Lower MIN_MAG or widen the year range.")
        sys.exit(2)
    print()
    print(f"Selected earthquake:  {ev['time_utc'].isoformat()}")
    print(f"  M = {ev['mag']:.2f} ({ev['mag_type']})")
    print(f"  depth = {ev['depth_km']:.2f} km, distance from SAFOD = {ev['dist_km']:.2f} km")
    print(f"  file = {eq_row['file_norm']}")
    print(f"  file window = {eq_row['startTime_dt']} to {eq_row['endTime_dt']}")

    quiet_row = find_quiet_file(db)
    print()
    print(f"Selected quiet file:  {quiet_row['startTime_dt'].isoformat()}")
    print(f"  file = {quiet_row['file_norm']}")

    # Read both
    X_eq,    info_eq    = read_file(eq_row["file_norm"])
    X_quiet, info_quiet = read_file(quiet_row["file_norm"])
    fs_eq    = float(info_eq["fs"])
    fs_quiet = float(info_quiet["fs"])
    dz       = float(eq_row["dCh"])
    print(f"fs eq={fs_eq:.2f} Hz, quiet={fs_quiet:.2f} Hz, dz={dz:.3f} m")

    # F-K
    f_e, k_e, FK_e = compute_fk(X_eq,    fs_eq,    dz)
    f_q, k_q, FK_q = compute_fk(X_quiet, fs_quiet, dz)

    # Plot 2x2
    fig, axs = plt.subplots(2, 2, figsize=(15, 11))
    plot_timedist(axs[0, 0], X_eq,    fs_eq,
                  f"EQ: M{ev['mag']:.1f} at {ev['time_utc'].strftime('%Y-%m-%d %H:%M:%S')} UTC "
                  f"(dist {ev['dist_km']:.1f} km, depth {ev['depth_km']:.1f} km)")
    plot_timedist(axs[0, 1], X_quiet, fs_quiet,
                  f"Quiet: {QUIET_DATE} UTC {QUIET_HOUR:02d}:XX (0 M>=1 events that week)")
    pcm_e = plot_fk(axs[1, 0], f_e, k_e, FK_e, "F-K, EQ window")
    pcm_q = plot_fk(axs[1, 1], f_q, k_q, FK_q, "F-K, quiet window")
    fig.colorbar(pcm_e, ax=axs[1, 0], label="Power (dB)")
    fig.colorbar(pcm_q, ax=axs[1, 1], label="Power (dB)")
    fig.suptitle("EQ vs quiet wavefield, SAFOD 2024", fontsize=12)
    fig.tight_layout()
    out_png = os.path.join(OUT_DIR, "fk_compare_eq_quiet.png")
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"\nWrote {out_png}")

    # Save arrays
    np.savez(
        os.path.join(OUT_DIR, "fk_compare_arrays.npz"),
        ev_time=ev["time_utc"].isoformat(), ev_mag=ev["mag"],
        ev_dist_km=ev["dist_km"], ev_depth_km=ev["depth_km"],
        ch_start=CH_START, ch_end=CH_END,
        fs_eq=fs_eq, fs_quiet=fs_quiet, dz=dz,
        f_eq=f_e, k_eq=k_e, FK_eq=FK_e,
        f_quiet=f_q, k_quiet=k_q, FK_quiet=FK_q,
    )
    print("Wrote fk_compare_arrays.npz. Done.")


if __name__ == "__main__":
    main()
