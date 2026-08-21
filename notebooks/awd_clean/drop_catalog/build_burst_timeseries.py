"""Extract several bursts as CONTINUOUS Nano time series, for the drop-marker figure.

Why this exists. Every other figure in this directory windows each drop and
re-zeroes it on its own pick, which is right for stacking but wrong for the
question "do the drops line up with the marks?" -- each panel then has exactly
one mark in it. This pulls the raw record straight through, unwindowed, so a
burst of ~20 drops appears as ~20 impulses on one continuous trace and the
delivered pick times can be drawn on top as vertical lines.

Bursts are sampled across the whole 24 h rather than taking one, because
detection is strongly time-varying: burst 30 has 0 of 20 drops detected on Nano
while burst 48 has 16 of 20. One burst is not representative of the survey.

Reads the Nano .pb files spanning each burst, beams over the same 81-439 m
aperture at the same 2,975 m/s used by nano_hierarchical_repeatability.py, and
saves the 1-D trace plus its time axis. Output stays small enough to commit, so
the notebook plots it without touching $OAK.

    sbatch awd_clean/drop_catalog/burst_timeseries_job.sh
    BURSTS=0,5,10 sbatch awd_clean/drop_catalog/burst_timeseries_job.sh
"""
import csv
import datetime as dt
import os
import sys
from pathlib import Path

import numpy as np

# Three DAS-utilities checkouts are live on this system and they are NOT
# interchangeable: only the ettore88 safod_das_git copy defines
# readFile_protobuf. Insert in reverse so the preferred path ends up first.
for _p in reversed([
        "/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python",
        "/home/groups/edunham/nberrios/safod_das/DAS-utilities/python"]):
    if (Path(_p) / "DASutils.py").exists() and _p not in sys.path:
        sys.path.insert(0, _p)
from DASutils import readFile_protobuf  # noqa: E402

HERE = Path(__file__).resolve().parent
AWD = HERE.parent
NANO_DIR = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/Nano")
MANIFEST = AWD / "awd_manifest.csv"
BEAM_NPZ = AWD / "nano_hierarchical_repeatability.npz"
OUT = HERE / "burst_timeseries.npz"

BURSTS = [int(x) for x in os.environ.get("BURSTS", "0,10,20,30,40,48").split(",")]
BAND = (30.0, 60.0)          # same band as the detection beam
VELOCITY = 2975.0            # same moveout
PAD_S = 20.0                 # context either side of each burst


def pb_start(name):
    """Start time from a Nano protobuf filename: parts[1]=date, parts[2]=HH.MM.SS."""
    p = Path(name).name.split("_")
    return dt.datetime.strptime(f"{p[1]} {p[2]}", "%Y-%m-%d %H.%M.%S")


def beam_one(path, ch, dist):
    """Read one .pb file and collapse the aperture onto a single 2975 m/s beam."""
    DAS, info = readFile_protobuf([str(path)], fmin=BAND[0], fmax=BAND[1],
                                  desampling=False, verbose=False)
    fs = float(info["fs"])
    sub = DAS[ch, :]
    beam = np.zeros(sub.shape[1], dtype=np.float64)
    for k in range(sub.shape[0]):
        beam += np.roll(sub[k], -int(round(dist[k] / VELOCITY * fs)))
    return (beam / sub.shape[0]).astype(np.float32), fs


def extract(burst, rows_by_burst, ch, dist):
    rows = rows_by_burst[burst]
    drops = [dt.datetime.fromisoformat(r["utc_time"]).replace(tzinfo=None) for r in rows]
    files = sorted({r["nano_file"] for r in rows}, key=pb_start)
    print(f"\nburst {burst}: {len(drops)} drops, "
          f"{drops[0]:%Y-%m-%d %H:%M:%S}-{drops[-1]:%H:%M:%S} UTC, {len(files)} file(s)")

    segs, t0s, fs = [], [], None
    for name in files:
        path = NANO_DIR / name
        if not path.exists():
            print(f"  MISSING {name} -- skipping burst {burst}")
            return None
        seg, fs = beam_one(path, ch, dist)
        segs.append(seg)
        t0s.append(pb_start(name))
        print(f"  {name}: {len(seg)/fs:.0f} s at {fs:.0f} Hz")

    origin = t0s[0]
    total = int(round((t0s[-1] - origin).total_seconds() * fs)) + len(segs[-1])
    trace = np.full(total, np.nan, dtype=np.float32)
    for seg, t0 in zip(segs, t0s):
        i = int(round((t0 - origin).total_seconds() * fs))
        trace[i:i + len(seg)] = seg
    tsec = np.arange(total) / fs

    lo = (drops[0] - origin).total_seconds() - PAD_S
    hi = (drops[-1] - origin).total_seconds() + PAD_S
    m = (tsec >= lo) & (tsec <= hi)
    # store time relative to the burst's FIRST drop, so panels share an x axis
    return {
        "trace": trace[m],
        "t_rel": (tsec[m] - (drops[0] - origin).total_seconds()).astype(np.float32),
        "drop_rel": np.array([(d - drops[0]).total_seconds() for d in drops]),
        "start_utc": drops[0].isoformat(),
        "n_drops": len(drops),
        "fs": fs,
    }


def main():
    rows_by_burst = {}
    for r in csv.DictReader(MANIFEST.open()):
        rows_by_burst.setdefault(int(r["burst_id"]), []).append(r)

    z = np.load(BEAM_NPZ, allow_pickle=True)
    ch, dist = z["aperture_channel"], z["aperture_distance_m"]
    print(f"aperture: {len(ch)} channels, {dist.min():.0f}-{dist.max():.0f} m")
    print(f"bursts requested: {BURSTS}")

    out, kept = {}, []
    for b in BURSTS:
        if b not in rows_by_burst:
            print(f"\nburst {b}: not in the manifest -- skipped")
            continue
        d = extract(b, rows_by_burst, ch, dist)
        if d is None:
            continue
        kept.append(b)
        for k, v in d.items():
            out[f"b{b}_{k}"] = v

    if not kept:
        raise SystemExit("no bursts extracted")
    out["bursts"] = np.array(kept)
    out["band_hz"] = np.array(BAND)
    out["velocity_mps"] = VELOCITY
    out["aperture_m"] = np.array([dist.min(), dist.max()])
    out["n_channels"] = len(ch)
    np.savez_compressed(OUT, **out)
    print(f"\nwrote {OUT.name}: {len(kept)} bursts {kept}, "
          f"{OUT.stat().st_size/1e6:.1f} MB")


if __name__ == "__main__":
    main()
