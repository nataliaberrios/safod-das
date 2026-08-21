"""Extract epochs as single-channel Nano traces for the GPS pick-QC figure.

This reproduces the "Identifying drops" figure from
``notebooks/JULY24_fml_coupling.ipynb`` cell 17 (Step 3c QC), which is the
figure of record for the picking technique. Parameters are taken from that
notebook and must not drift from it:

    FMIN,FMAX = 1,100   broadband; NOT the 30-60 Hz detection beam
    FILE_DUR_S = 300    Nano .pb files are 5 minutes

**One parameter is deliberately NOT inherited.** That notebook uses
``CH_DETECT = 50``, commented "ch 50 ~ 63 m, best-coupled depth".
``nano_find_wellhead.txt`` puts fibre entry at channel **73** (92 m along
fibre), so channel 50 is **in the air** and its 0.995 neighbour correlation is
surface coupling, not good coupling. ``AUDIT_2026-08-20.md`` 1.4 voids Nano
results built on air channels. The channel is chosen from the data here instead,
restricted to the cemented range.

The single channel is the point. Beaming over 81-439 m averages the drops into
the background and the impulses stop being legible; one well-coupled channel
shows each drop as a clean spike sitting on its GPS mark.

Each panel uses the ONE .pb file holding the median drop of its epoch, and the
x axis is time within that file -- so panels legitimately span different ranges,
and an epoch straddling a file boundary shows fewer marks than it has drops.
That is what the "(N in file)" half of each title reports.

    sbatch awd_clean/drop_catalog/burst_timeseries_job.sh
    BURSTS=0,5,10 sbatch awd_clean/drop_catalog/burst_timeseries_job.sh
"""
import csv
import datetime as dt
import os
import sys
from pathlib import Path

import numpy as np

# Three DAS-utilities checkouts are live and they are NOT interchangeable: only
# the ettore88 safod_das_git copy defines readFile_protobuf. Insert in reverse
# so the preferred path ends up first.
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
OUT = HERE / "burst_timeseries.npz"

# --- band and file duration from JULY24_fml_coupling.ipynb cell 17 ---
FMIN, FMAX = 1.0, 100.0
FILE_DUR_S = 300.0
PAD_S = 5.0                       # the notebook's xlim padding
# Same six epochs the original figure shows.
BURSTS = [int(x) for x in os.environ.get("BURSTS", "0,9,19,28,38,47").split(",")]

# CHANNEL. The original figure used CH_DETECT = 50, annotated there as
# "ch 50 ~ 63 m, best-coupled depth". That is WRONG and the figure must not be
# reproduced with it: `nano_find_wellhead.txt` puts fibre entry at channel 73
# (92 m along fibre), so channels 0-72 are in the AIR. Channel 50's high RMS and
# 0.995 neighbour correlation are surface/air coupling, not good coupling.
# `AUDIT_2026-08-20.md` §1.4 voids Nano results built on air channels.
#
# So the channel is chosen from the data instead of asserted: scan the cemented
# range for the channel whose drop impulses stand highest above their own
# background, and report it. Set CH_DETECT in the environment to override.
CH_MIN = 73                       # wellhead; first channel in the ground
CH_MAX = 347                      # top of the aperture used by the detection beam
CH_OVERRIDE = os.environ.get("CH_DETECT")


def pb_start(name):
    """Start time from a Nano .pb filename: parts[1]=date, parts[2]=HH.MM.SS."""
    p = Path(name).name.split("_")
    return dt.datetime.strptime(f"{p[1]} {p[2]}", "%Y-%m-%d %H.%M.%S")


def pick_channel(DAS, fs, rel):
    """Channel in [CH_MIN, CH_MAX] whose drop impulses stand highest above noise.

    Score is median peak |amplitude| in a 0.3 s window after each drop, divided
    by the channel's own median |amplitude| away from drops. Per-channel
    normalisation means a loud but uninformative channel cannot win.
    """
    n = DAS.shape[1]
    on = np.zeros(n, dtype=bool)
    for r in rel:
        i = int(r * fs)
        on[max(0, i):min(n, i + int(0.3 * fs))] = True
    hi = min(CH_MAX, DAS.shape[0] - 1)
    best, best_score = None, -np.inf
    for c in range(CH_MIN, hi + 1):
        x = np.abs(np.asarray(DAS[c, :], dtype=np.float64))
        bg = np.median(x[~on])
        if not np.isfinite(bg) or bg <= 0:
            continue
        score = np.median([x[max(0, int(r*fs)):min(n, int(r*fs)+int(0.3*fs))].max()
                           for r in rel]) / bg
        if score > best_score:
            best, best_score = c, score
    return best, best_score


def main():
    rows_by_burst = {}
    for r in csv.DictReader(MANIFEST.open()):
        rows_by_burst.setdefault(int(r["burst_id"]), []).append(r)
    print(f"epochs requested: {BURSTS}")
    print(f"band {FMIN:.0f}-{FMAX:.0f} Hz (JULY24 cell 17)")
    how = f"forced to {CH_OVERRIDE}" if CH_OVERRIDE else f"chosen from data in [{CH_MIN}, {CH_MAX}]"
    print(f"channel: {how}")
    print("  NOT channel 50 -- nano_find_wellhead.txt puts fibre entry at 73, "
          "so 0-72 are in the air")
    ch_detect = int(CH_OVERRIDE) if CH_OVERRIDE else None

    out, kept = {}, []
    for b in BURSTS:
        if b not in rows_by_burst:
            print(f"\nepoch {b}: not in the manifest -- skipped")
            continue
        rows = rows_by_burst[b]
        drops = [dt.datetime.fromisoformat(r["utc_time"]).replace(tzinfo=None)
                 for r in rows]

        # The file holding the MEDIAN drop, as the original does -- not every
        # file the epoch touches.
        median_drop = drops[len(drops) // 2]
        name = rows[len(rows) // 2]["nano_file"]
        path = NANO_DIR / name
        if not path.exists():
            print(f"\nepoch {b}: missing {name} -- skipped")
            continue
        t_file = pb_start(name)

        DAS, info = readFile_protobuf([str(path)], fmin=FMIN, fmax=FMAX,
                                      desampling=False, verbose=False)
        fs = float(info["fs"])
        if ch_detect is None:                      # choose once, on the first epoch
            rel0 = np.array([(d - t_file).total_seconds() for d in drops
                             if t_file <= d < t_file + dt.timedelta(seconds=FILE_DUR_S)])
            ch_detect, sc = pick_channel(DAS, fs, rel0)
            print(f"  chose channel {ch_detect} (drop/background {sc:.1f}x); "
                  f"reused for every epoch")
        trace = np.asarray(DAS[ch_detect, :], dtype=np.float32)
        t = np.arange(len(trace)) / fs

        # drops that actually fall inside THIS file
        rel = np.array([(d - t_file).total_seconds() for d in drops
                        if t_file <= d < t_file + dt.timedelta(seconds=FILE_DUR_S)])
        if rel.size == 0:
            print(f"\nepoch {b}: no drops inside {name} -- skipped")
            continue

        # crop to the plotted window so the stored arrays stay small
        m = (t >= rel.min() - PAD_S - 1) & (t <= rel.max() + PAD_S + 1)
        print(f"\nepoch {b:2d}: {len(drops)} GPS drops ({rel.size} in file), "
              f"{median_drop:%Y-%m-%d %H:%M:%S} UTC, {name}")
        print(f"           file window {rel.min():.1f}-{rel.max():.1f} s, "
              f"{m.sum()/fs:.0f} s kept at {fs:.0f} Hz")

        kept.append(b)
        out[f"b{b}_trace"] = trace[m]
        out[f"b{b}_t_file"] = t[m].astype(np.float32)
        out[f"b{b}_drop_in_file"] = rel
        out[f"b{b}_n_gps"] = len(drops)
        out[f"b{b}_start_utc"] = drops[0].isoformat()
        out[f"b{b}_file"] = name

    if not kept:
        raise SystemExit("no epochs extracted")
    out["bursts"] = np.array(kept)
    out["ch_detect"] = ch_detect
    out["band_hz"] = np.array([FMIN, FMAX])
    out["pad_s"] = PAD_S
    np.savez_compressed(OUT, **out)
    print(f"\nwrote {OUT.name}: {len(kept)} epochs {kept}, "
          f"{OUT.stat().st_size/1e6:.1f} MB")


if __name__ == "__main__":
    main()
