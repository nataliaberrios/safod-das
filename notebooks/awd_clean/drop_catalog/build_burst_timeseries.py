"""Extract one burst as a CONTINUOUS Nano time series, for the drop-marker figure.

Why this exists. Every other figure in this directory windows each drop and
re-zeroes it on its own pick, which is right for stacking but wrong for the
question "do the drops line up with the marks?" -- each panel then has exactly
one mark in it. This pulls the raw record straight through, unwindowed, so a
burst of ~20 drops appears as ~20 impulses on one continuous trace and the
delivered pick times can be drawn on top as vertical lines.

Reads the Nano .pb files spanning the burst, beams over the same 81-439 m
aperture at the same 2,975 m/s used by nano_hierarchical_repeatability.py, and
saves the 1-D trace plus its absolute UTC time axis. Output is small (~a few MB),
so the notebook plots it without touching $OAK.

    sbatch awd_clean/drop_catalog/burst_timeseries_job.sh
"""
import csv
import datetime as dt
import os
import sys
from pathlib import Path

import numpy as np
from scipy.signal import butter, hilbert, sosfiltfilt

# Three DAS-utilities checkouts are live on this system and they are NOT
# interchangeable: only the ettore88 safod_das_git copy defines
# readFile_protobuf. Pick by capability rather than by first-found, and insert
# in reverse so the preferred path ends up at sys.path[0].
for _p in reversed([
        "/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python",
        "/home/groups/edunham/nberrios/safod_das/DAS-utilities/python"]):
    if (Path(_p) / "DASutils.py").exists() and _p not in sys.path:
        sys.path.insert(0, _p)
from DASutils import readFile_protobuf  # noqa: E402
print("DASutils from:", sys.modules["DASutils"].__file__)

HERE = Path(__file__).resolve().parent
AWD = HERE.parent
NANO_DIR = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/Nano")
MANIFEST = AWD / "awd_manifest.csv"
BEAM_NPZ = AWD / "nano_hierarchical_repeatability.npz"
OUT = HERE / "burst_timeseries.npz"

BURST = int(os.environ.get("BURST", "0"))
BAND = (30.0, 60.0)          # same band as the detection beam
VELOCITY = 2975.0            # same moveout
PAD_S = 25.0                 # context either side of the burst


def pb_start(name):
    """Start time from a Nano protobuf filename: parts[1]=date, parts[2]=HH.MM.SS."""
    p = Path(name).name.split("_")
    return dt.datetime.strptime(f"{p[1]} {p[2]}", "%Y-%m-%d %H.%M.%S")


def main():
    rows = [r for r in csv.DictReader(MANIFEST.open()) if int(r["burst_id"]) == BURST]
    if not rows:
        raise SystemExit(f"no drops for burst {BURST}")
    drops = [dt.datetime.fromisoformat(r["utc_time"]).replace(tzinfo=None) for r in rows]
    files = sorted({r["nano_file"] for r in rows}, key=pb_start)
    print(f"burst {BURST}: {len(drops)} drops, "
          f"{drops[0]:%Y-%m-%d %H:%M:%S} -> {drops[-1]:%H:%M:%S} UTC")
    print(f"spanning {len(files)} Nano file(s): {files}")

    z = np.load(BEAM_NPZ, allow_pickle=True)
    ch, dist = z["aperture_channel"], z["aperture_distance_m"]
    print(f"aperture: {len(ch)} channels, {dist.min():.0f}-{dist.max():.0f} m")

    # Read each file once, beam it, then concatenate on the absolute time axis.
    segs, t0s, fs = [], [], None
    for name in files:
        path = NANO_DIR / name
        if not path.exists():
            raise SystemExit(f"missing {path}")
        DAS, info = readFile_protobuf([str(path)], fmin=BAND[0], fmax=BAND[1],
                                      desampling=False, verbose=False)
        fs = float(info["fs"])
        sub = DAS[ch, :]
        # shift each channel back along the 2975 m/s moveout, then average
        beam = np.zeros(sub.shape[1], dtype=np.float64)
        for k in range(sub.shape[0]):
            s = int(round(dist[k] / VELOCITY * fs))
            beam += np.roll(sub[k], -s)
        beam /= sub.shape[0]
        segs.append(beam.astype(np.float32))
        t0s.append(pb_start(name))
        print(f"  {name}: {sub.shape[1]/fs:.0f} s at {fs:.0f} Hz, starts {t0s[-1]:%H:%M:%S}")

    # Build one continuous trace on an absolute-seconds axis.
    origin = t0s[0]
    total = int(round(((t0s[-1] - origin).total_seconds()) * fs)) + len(segs[-1])
    trace = np.full(total, np.nan, dtype=np.float32)
    for seg, t0 in zip(segs, t0s):
        i = int(round((t0 - origin).total_seconds() * fs))
        trace[i:i + len(seg)] = seg
    tsec = np.arange(total) / fs

    # Crop to the burst plus padding.
    lo = (drops[0] - origin).total_seconds() - PAD_S
    hi = (drops[-1] - origin).total_seconds() + PAD_S
    m = (tsec >= lo) & (tsec <= hi)
    trace, tsec = trace[m], tsec[m]
    env = np.abs(hilbert(np.nan_to_num(trace)))

    np.savez_compressed(
        OUT,
        trace=trace, envelope=env.astype(np.float32), t_seconds=tsec.astype(np.float64),
        origin_utc=str(origin), fs=fs, band_hz=np.array(BAND), velocity_mps=VELOCITY,
        burst_id=BURST, n_channels=len(ch),
        aperture_m=np.array([dist.min(), dist.max()]),
        drop_seconds=np.array([(d - origin).total_seconds() for d in drops]),
        drop_utc=np.array([d.isoformat() for d in drops]),
        source_files=np.array(files),
    )
    print(f"\nwrote {OUT.name}: {len(trace)/fs:.0f} s at {fs:.0f} Hz, "
          f"{len(drops)} drop marks, {OUT.stat().st_size/1e6:.1f} MB")


if __name__ == "__main__":
    main()
