#!/usr/bin/env python3
"""Resumable exact-resolution unfiltered Lellouch-style seasonal stacks."""
from __future__ import annotations
import json, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "ambient_transfer"
MANIFEST = OUT / "seasonal_day_selection.json"
LOG = OUT / "seasonal_unfiltered_progress.log"
CHUNK = 10

def main():
    days = json.loads(MANIFEST.read_text())["days"]
    py = sys.executable
    with LOG.open("a") as lf:
        lf.write(f"RUN days={len(days)} chunk={CHUNK} exact_resolution=True norm_seconds=5.0\n"); lf.flush()
        for item in days:
            date = item["date"]; n = int(item["nfiles"])
            for start in range(0, n, CHUNK):
                count = min(CHUNK, n-start)
                stem = OUT / f"transfer_{date}_start{start}_n{count}.json"
                if stem.exists():
                    lf.write(f"SKIP {date} {start} {count} existing\n"); lf.flush(); continue
                cmd = [py, str(ROOT / "ambient_transfer_test.py"), "--date", date,
                       "--start", str(start), "--nfiles", str(count),
                       "--norm-seconds", "5"]
                t=time.time(); lf.write(f"START {date} {start} {count}\n"); lf.flush()
                r=subprocess.run(cmd, cwd=ROOT.parent, capture_output=True, text=True)
                lf.write(r.stdout); lf.write(r.stderr)
                lf.write(f"END {date} {start} rc={r.returncode} elapsed={time.time()-t:.1f}\n"); lf.flush()
                if r.returncode != 0: raise SystemExit(f"chunk failed {date} {start}")
    print("seasonal unfiltered chunks complete")

if __name__ == "__main__": main()
