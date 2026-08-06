#!/usr/bin/env python3
"""Resumable one-day driver for the corrected seasonal signed-lag workflow."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
LEGACY = ROOT / "ambient_transfer"
OUT = LEGACY / "signed_lag_v2"
MANIFEST = LEGACY / "seasonal_day_selection.json"
CHUNK = 10


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--day-index", type=int, required=True)
    args = parser.parse_args()
    days = json.loads(MANIFEST.read_text())["days"]
    item = days[args.day_index]
    date = item["date"]
    nfiles = int(item["nfiles"])
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / f"seasonal_signed_fk_v2_{date}.log"
    with log.open("a") as handle:
        handle.write(f"RUN date={date} files={nfiles} chunk={CHUNK}\n")
        handle.flush()
        for start in range(0, nfiles, CHUNK):
            count = min(CHUNK, nfiles - start)
            stem = OUT / f"signed_fk_v2_{date}_start{start}_n{count}"
            if stem.with_suffix(".json").exists() and stem.with_suffix(".npz").exists():
                handle.write(f"SKIP start={start} n={count}\n")
                handle.flush()
                continue
            command = [
                sys.executable,
                str(ROOT / "ambient_signed_fk_v2.py"),
                "--date", date,
                "--start", str(start),
                "--nfiles", str(count),
                "--output-dir", str(OUT),
            ]
            began = time.time()
            handle.write(f"START start={start} n={count}\n")
            handle.flush()
            result = subprocess.run(command, cwd=ROOT.parent, capture_output=True, text=True)
            handle.write(result.stdout)
            handle.write(result.stderr)
            handle.write(
                f"END start={start} rc={result.returncode} elapsed={time.time()-began:.1f}\n"
            )
            handle.flush()
            if result.returncode != 0:
                raise SystemExit(f"chunk failed: {date} start={start} rc={result.returncode}")
    print(f"complete {date}")


if __name__ == "__main__":
    main()
