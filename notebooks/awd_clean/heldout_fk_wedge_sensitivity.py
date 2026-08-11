#!/usr/bin/env python3
"""Held-out sensitivity test for the ambient signed F-K wedge.

The candidate velocity wedges are evaluated on independent training and test
dates. The training set selects one negative signed wedge by its aggregate
3.2-km/s score; that frozen wedge is then evaluated on held-out dates with a
receiver-permutation null. This addresses selection optimism without claiming
that any candidate wedge is an independent velocity estimate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ambient_transfer_test import CSV, OUT, corrected_path, load_segment, preprocess, normalized_corr_pairs


CANDIDATES = [(2000.0, 3000.0), (2500.0, 3500.0), (2500.0, 4500.0),
              (3000.0, 4000.0), (3500.0, 5000.0)]
TRAIN_DATES = ["2024-12-20", "2025-02-24", "2024-05-11", "2024-10-28"]
TEST_DATES = ["2025-03-04", "2024-06-17", "2024-06-26", "2024-11-30"]


def mask(nz, nt, fs, dx, vmin, vmax, sign=-1):
    f = np.fft.fftfreq(nt, 1.0 / fs)
    k = np.fft.fftfreq(nz, dx)
    K, F = np.meshgrid(k, f, indexing="ij")
    af, ak = np.abs(F), np.abs(K)
    v = af / np.maximum(ak, 1e-12)
    out = (af >= 5.0) & (af <= 20.0) & (v >= vmin) & (v <= vmax) & (ak > 0)
    return out & ((F * K < 0) if sign < 0 else (F * K > 0))


def block_stacks(day, starts, nfiles):
    stacks = {w: None for w in CANDIDATES}
    used = 0
    lags = None
    distances = None
    for start in starts:
        rows = day.iloc[start:start + nfiles]
        for row in rows.itertuples(index=False):
            path = corrected_path(row.file)
            if not path.exists():
                continue
            x, fs, dx = load_segment(path)
            x = preprocess(x, fs, norm_seconds=5.0)
            y = x[::2, ::2]
            fs2, dx2 = fs / 2.0, dx * 2.0
            targets = [int(round(50.0 * j / dx2)) for j in range(1, 15)]
            targets = [c for c in targets if c < y.shape[0]]
            pairs = [(0, c) for c in targets]
            lags, corr = normalized_corr_pairs(y, pairs, fs2)
            distances = np.asarray(targets, dtype=float) * dx2
            for wedge in CANDIDATES:
                filtered = np.fft.ifft2(np.fft.fft2(y) * mask(y.shape[0], y.shape[1], fs2, dx2, *wedge)).real
                _, cc = normalized_corr_pairs(filtered, pairs, fs2)
                stacks[wedge] = cc if stacks[wedge] is None else stacks[wedge] + cc
            used += 1
    if used == 0:
        raise RuntimeError("No usable files")
    return {w: stacks[w] / used for w in CANDIDATES}, lags, distances, used


def velocity_scores(top, lags, distances):
    velocities = np.linspace(1500.0, 5500.0, 161)
    scores = np.array([
        np.nanmedian([row[np.argmin(np.abs(lags - distance / velocity))]
                      for row, distance in zip(top, distances)])
        for velocity in velocities
    ])
    return velocities, scores


def summarize(top, lags, distances, null_seed):
    velocities, scores = velocity_scores(top, lags, distances)
    rng = np.random.default_rng(null_seed)
    null = np.array([
        np.nanmax(velocity_scores(top, lags, rng.permutation(distances))[1])
        for _ in range(500)
    ])
    peak = int(np.nanargmax(scores))
    return {
        "peak_velocity_m_s": float(velocities[peak]),
        "peak_score": float(scores[peak]),
        "score_3200": float(scores[np.argmin(np.abs(velocities - 3200.0))]),
        "null95": float(np.quantile(null, 0.95)),
        "p_peak": float((1.0 + np.sum(null >= scores[peak])) / (len(null) + 1.0)),
        "velocities_m_s": velocities,
        "scores": scores,
        "null": null,
    }


def day_blocks(db, date, nblocks):
    day = db[db.t.dt.strftime("%Y-%m-%d") == date].sort_values("t").reset_index(drop=True)
    if len(day) < 10:
        return day, [0]
    starts = np.linspace(0, max(0, len(day) - 10), nblocks).round().astype(int).tolist()
    return day, sorted(set(starts))


def run(args):
    db = pd.read_csv(CSV, sep=r"\s+")
    db = db[db.nSamples > 0].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    manifest = json.loads((OUT / "seasonal_day_selection.json").read_text())
    nfiles = {item["date"]: int(item["nfiles"]) for item in manifest["days"]}

    aggregate = {}
    block_counts = {}
    for group, dates in (("train", TRAIN_DATES), ("heldout", TEST_DATES)):
        all_stacks = {w: [] for w in CANDIDATES}
        lags = distances = None
        nused = 0
        for date in dates:
            day, starts = day_blocks(db, date, args.blocks_per_day)
            stacks, lags, distances, used = block_stacks(day, starts, args.block_files)
            nused += used
            block_counts[date] = {"starts": starts, "used_files": used, "available": len(day)}
            for wedge in CANDIDATES:
                all_stacks[wedge].append(stacks[wedge] * used)
        aggregate[group] = {
            wedge: sum(all_stacks[wedge]) / nused for wedge in CANDIDATES
        }
        aggregate[group]["lags"] = lags
        aggregate[group]["distances"] = distances
        aggregate[group]["used_files"] = nused

    train_summary = {w: summarize(aggregate["train"][w], aggregate["train"]["lags"], aggregate["train"]["distances"], 20260805 + i)
                     for i, w in enumerate(CANDIDATES)}
    selected = max(CANDIDATES, key=lambda w: train_summary[w]["peak_score"])
    heldout_summary = {w: summarize(aggregate["heldout"][w], aggregate["heldout"]["lags"], aggregate["heldout"]["distances"], 20260905 + i)
                       for i, w in enumerate(CANDIDATES)}

    OUT.mkdir(exist_ok=True)
    stem = "heldout_fk_wedge_sensitivity"
    report = {
        "train_dates": TRAIN_DATES,
        "heldout_dates": TEST_DATES,
        "blocks_per_day": args.blocks_per_day,
        "block_files": args.block_files,
        "candidate_velocity_wedges_m_s": [list(w) for w in CANDIDATES],
        "selected_wedge_m_s": list(selected),
        "block_counts": block_counts,
        "train_used_files": aggregate["train"]["used_files"],
        "heldout_used_files": aggregate["heldout"]["used_files"],
        "train": {str(w): {k: v for k, v in s.items() if k not in {"velocities_m_s", "scores", "null"}}
                  for w, s in train_summary.items()},
        "heldout": {str(w): {k: v for k, v in s.items() if k not in {"velocities_m_s", "scores", "null"}}
                    for w, s in heldout_summary.items()},
    }
    (OUT / f"{stem}.json").write_text(json.dumps(report, indent=2))

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, group, summaries in zip(axes[0], ["train", "heldout"], [train_summary, heldout_summary]):
        for wedge, summary in summaries.items():
            label = f"{wedge[0]/1000:.1f}–{wedge[1]/1000:.1f} km/s"
            ax.plot(summary["velocities_m_s"] / 1000.0, summary["scores"], label=label,
                    lw=2.2 if wedge == selected else 1.0)
        ax.axvline(3.2, color="k", ls=":", lw=0.8)
        ax.set_title(f"{group}: {'selected' if group == 'train' else 'frozen held-out'} wedges")
        ax.set_xlabel("Trial velocity (km/s)")
        ax.set_ylabel("Median normalized correlation")
        ax.legend(frameon=False, fontsize=7)
    for ax, group, summaries in zip(axes[1], ["train", "heldout"], [train_summary, heldout_summary]):
        wedges = [f"{w[0]/1000:.1f}–{w[1]/1000:.1f}" for w in CANDIDATES]
        scores = [summaries[w]["score_3200"] for w in CANDIDATES]
        pvals = [summaries[w]["p_peak"] for w in CANDIDATES]
        xx = np.arange(len(wedges))
        ax.bar(xx, scores, color=["tab:red" if w == selected else "0.6" for w in CANDIDATES])
        ax.set_xticks(xx, wedges, rotation=35, ha="right")
        ax.set_title(f"{group}: score at 3.2 km/s")
        ax.set_ylabel("Score")
        for x, score, p in zip(xx, scores, pvals):
            ax.text(x, score, f"p={p:.3f}", ha="center", va="bottom", fontsize=7)
    fig.suptitle("Held-out signed F-K wedge sensitivity")
    fig.savefig(OUT / f"{stem}.png", dpi=250)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--blocks-per-day", type=int, default=4)
    parser.add_argument("--block-files", type=int, default=10)
    run(parser.parse_args())

