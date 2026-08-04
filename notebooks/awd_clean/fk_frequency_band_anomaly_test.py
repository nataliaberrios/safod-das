"""Frequency-band validation of the conditional SAFOD moveout anomaly test.

The same-fiber channel-0/depth convention used in v37 is retained.  A fixed,
deterministic set of blocks is sampled from each selected complete day.  Each
block is processed independently in several signed negative F-K bands, then
the early/late velocities and descriptive breakpoint are summarized by date.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import theilslopes

from ambient_transfer_test import corrected_path, load_segment, normalized_corr_pairs, preprocess


ROOT = Path(__file__).resolve().parent
CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")
MANIFEST = ROOT / "ambient_transfer" / "seasonal_day_selection.json"
OUT = ROOT / "ambient_transfer" / "frequency_band_anomaly_test.json"
BANDS = ((3.0, 8.0), (5.0, 12.0), (8.0, 20.0), (15.0, 30.0))


def fk_negative(x, fs, dx, fmin, fmax):
    x = x[::2, ::2]
    fs, dx = fs / 2.0, dx * 2.0
    f = np.fft.fftfreq(x.shape[1], 1.0 / fs)
    k = np.fft.fftfreq(x.shape[0], dx)
    K, F = np.meshgrid(k, f, indexing="ij")
    af, ak = np.abs(F), np.abs(K)
    v = af / np.maximum(ak, 1e-12)
    mask = (af >= fmin) & (af <= fmax) & (v >= 2500.0) & (v <= 4500.0) & (ak > 0) & ((F * K) < 0)
    return np.fft.ifft2(np.fft.fft2(x) * mask).real, fs, dx


def ridge_metrics(top, lags, dist):
    lag = lags[np.argmax(top, axis=1)]
    q = (dist <= 650.0) & np.isfinite(lag)
    d, t = dist[q], lag[q]

    def velocity(lo, hi):
        qq = (d >= lo) & (d <= hi)
        if qq.sum() < 4:
            return np.nan
        return float(1.0 / theilslopes(t[qq], d[qq]).slope)

    one = np.polyfit(d, t, 1)
    sse1 = np.sum((t - np.polyval(one, d)) ** 2)
    bps = np.arange(250.0, 551.0, 25.0)
    improvements = []
    for bp in bps:
        left, right = d <= bp, d > bp
        if left.sum() < 4 or right.sum() < 4:
            improvements.append(np.nan)
            continue
        p1, p2 = np.polyfit(d[left], t[left], 1), np.polyfit(d[right], t[right], 1)
        sse = np.sum((t[left] - np.polyval(p1, d[left])) ** 2) + np.sum((t[right] - np.polyval(p2, d[right])) ** 2)
        improvements.append(1.0 - sse / sse1)
    return {
        "early_velocity_m_s": velocity(50.0, 350.0),
        "late_velocity_m_s": velocity(450.0, 650.0),
        "best_breakpoint_m": float(bps[int(np.nanargmax(improvements))]),
        "breakpoint_sse_reduction": float(np.nanmax(improvements)),
    }


def main(args):
    manifest = json.loads(MANIFEST.read_text())
    dates = [d["date"] for d in manifest["days"]]
    db = pd.read_csv(CSV, sep=r"\s+")
    db = db[db.nSamples > 0].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    blocks_per_day = args.blocks_per_day
    nfiles = args.nfiles
    results = []
    for date in dates:
        day = db[db.t.dt.strftime("%Y-%m-%d") == date].sort_values("t")
        starts = np.linspace(0, max(0, len(day) - nfiles), blocks_per_day, dtype=int)
        for start in sorted(set(starts.tolist())):
            stacks = {band: None for band in BANDS}
            used = 0
            lags = dist = None
            for row in day.iloc[start : start + nfiles].itertuples(index=False):
                path = corrected_path(row.file)
                if not path.exists():
                    continue
                raw, fs, dx = load_segment(path)
                x = preprocess(raw, fs, fmin=3.0, fmax=30.0, norm_seconds=5.0)
                for band in BANDS:
                    y, fy, dxy = fk_negative(x, fs, dx, *band)
                    targets = [int(round(50.0 * j / dxy)) for j in range(1, 15) if int(round(50.0 * j / dxy)) < y.shape[0]]
                    lag, cc = normalized_corr_pairs(y, [(0, c) for c in targets], fy)
                    if stacks[band] is None:
                        stacks[band], lags, dist = cc, lag, np.asarray(targets) * dxy
                    else:
                        stacks[band] += cc
                used += 1
            if used == 0:
                continue
            for band, stack in stacks.items():
                metric = ridge_metrics(stack / used, lags, dist)
                results.append({"date": date, "start": int(start), "nfiles": int(used), "band_hz": list(band), **metric})
    summary = {"assumption": "Same cemented main-hole fiber; current channel 0 and 1.020952-m spacing inherit Lellouch position/depth convention.", "bands_hz": [list(b) for b in BANDS], "blocks_per_day": blocks_per_day, "nfiles_per_block": nfiles, "results": results}
    OUT.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--blocks-per-day", type=int, default=4)
    parser.add_argument("--nfiles", type=int, default=10)
    main(parser.parse_args())
