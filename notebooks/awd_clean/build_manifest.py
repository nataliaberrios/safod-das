"""Build the authoritative SAFOD AWD drop/burst/file-availability manifest."""

from __future__ import annotations

import csv
import glob
import os
import re
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GPS_FILE = ROOT / "p26.cc9.txt"
NANO_DIR = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/Nano")
DEEP_DIR = Path(
    "/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/"
    "01_--_recording_2026-06-15T230629Z_--_active_source"
)
OUT = Path(__file__).with_name("awd_manifest.csv")

NANO_DURATION_S = 300.0
DEEP_DURATION_S = 62.0
BURST_GAP_S = 60.0
MIN_BURST_DROPS = 5


@dataclass(frozen=True)
class TimedFile:
    start: datetime
    path: str
    duration_s: float

    def contains(self, time: datetime) -> bool:
        return self.start <= time < self.start + timedelta(seconds=self.duration_s + 2)


def nano_time(path: str) -> datetime | None:
    parts = os.path.basename(path).split("_")
    try:
        value = f"{parts[1]}T{parts[2].replace('.', ':')}Z"
        return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    except (IndexError, ValueError):
        return None


def deep_time(path: str) -> datetime | None:
    match = re.search(r"_(\d{4}-\d{2}-\d{2}T\d{6}Z)\.h5$", path)
    if not match:
        return None
    return datetime.strptime(match.group(1), "%Y-%m-%dT%H%M%SZ").replace(
        tzinfo=timezone.utc
    )


def timed_files(pattern: str, parser, duration_s: float) -> list[TimedFile]:
    result = []
    for path in glob.glob(pattern):
        start = parser(path)
        if start is not None:
            result.append(TimedFile(start, path, duration_s))
    return sorted(result, key=lambda item: item.start)


def containing(time: datetime, files: list[TimedFile]) -> TimedFile | None:
    # File counts are modest; explicit selection makes overlap handling visible.
    hits = [item for item in files if item.contains(time)]
    return max(hits, key=lambda item: item.start) if hits else None


def gps_rows() -> list[tuple[datetime, float, float]]:
    rows = []
    with GPS_FILE.open() as stream:
        next(stream)
        for line in stream:
            fields = line.split()
            if len(fields) < 2:
                continue
            time = datetime.fromisoformat(fields[1]).replace(tzinfo=timezone.utc)
            max_cc = float(fields[2]) if len(fields) > 2 else float("nan")
            diff_prev = float(fields[3]) if len(fields) > 3 else float("nan")
            rows.append((time, max_cc, diff_prev))
    return sorted(rows)


def main() -> None:
    nano = timed_files(str(NANO_DIR / "*.pb"), nano_time, NANO_DURATION_S)
    deep = timed_files(str(DEEP_DIR / "*.h5"), deep_time, DEEP_DURATION_S)
    if not nano:
        raise RuntimeError(f"No Nano files found in {NANO_DIR}")

    # The AWD acquisition window is defined by Nano coverage. Retain every
    # qualifying burst; do not force an expected 48-epoch count.
    survey_start = nano[0].start
    survey_end = nano[-1].start + timedelta(seconds=NANO_DURATION_S)
    gps = [row for row in gps_rows() if survey_start <= row[0] <= survey_end]

    grouped: list[list[tuple[datetime, float, float]]] = []
    for row in gps:
        if not grouped or (row[0] - grouped[-1][-1][0]).total_seconds() > BURST_GAP_S:
            grouped.append([row])
        else:
            grouped[-1].append(row)
    bursts = [group for group in grouped if len(group) >= MIN_BURST_DROPS]

    fieldnames = [
        "burst_id", "drop_id", "utc_time", "seconds_from_burst_start",
        "gps_max_cc", "gps_diff_prev_s", "nano_available", "nano_file",
        "nano_offset_s", "deep_available", "deep_file", "deep_offset_s",
        "paired_available",
    ]
    with OUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for burst_id, burst in enumerate(bursts):
            beginning = burst[0][0]
            for drop_id, (time, max_cc, diff_prev) in enumerate(burst):
                nf = containing(time, nano)
                df = containing(time, deep)
                writer.writerow({
                    "burst_id": burst_id,
                    "drop_id": drop_id,
                    "utc_time": time.isoformat(),
                    "seconds_from_burst_start": f"{(time - beginning).total_seconds():.6f}",
                    "gps_max_cc": f"{max_cc:.6f}",
                    "gps_diff_prev_s": f"{diff_prev:.6f}",
                    "nano_available": int(nf is not None),
                    "nano_file": os.path.basename(nf.path) if nf else "",
                    "nano_offset_s": f"{(time - nf.start).total_seconds():.6f}" if nf else "",
                    "deep_available": int(df is not None),
                    "deep_file": os.path.basename(df.path) if df else "",
                    "deep_offset_s": f"{(time - df.start).total_seconds():.6f}" if df else "",
                    "paired_available": int(nf is not None and df is not None),
                })

    total = sum(map(len, bursts))
    print(f"Nano files: {len(nano)}; Deep files: {len(deep)}")
    print(f"GPS drops in Nano coverage: {len(gps)}")
    print(f"Qualifying AWD bursts: {len(bursts)}; drops: {total}")
    print(f"Drops per burst: min={min(map(len, bursts))}, max={max(map(len, bursts))}")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
