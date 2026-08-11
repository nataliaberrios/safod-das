#!/usr/bin/env python3
"""Resumable array-task driver for pre-F-K full-pipeline null realizations."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ambient_fk_full_pipeline_null_v2 import run


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--date", required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--nfiles", type=int, default=30)
    parser.add_argument("--null-index-offset", type=int, default=0)
    parser.add_argument("--nulls-per-task", type=int, default=1)
    parser.add_argument("--null-methods", default="channel_permutation,circular_time_shift")
    parser.add_argument("--seed", type=int, default=20260805)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.task_id < 0 or args.nulls_per_task < 1:
        raise ValueError("task id must be non-negative and nulls per task positive")

    null_start = args.null_index_offset + args.task_id * args.nulls_per_task
    stem = (
        f"fk_full_pipeline_null_v2_{args.date}_start{args.start}_n{args.nfiles}_"
        f"r{null_start}-{null_start + args.nulls_per_task - 1}"
    )
    json_path = args.output_dir / f"{stem}.json"
    npz_path = args.output_dir / f"{stem}.npz"
    if json_path.exists() and npz_path.exists():
        report = json.loads(json_path.read_text())
        expected = {
            "workflow_version": "ambient_fk_full_pipeline_null_v2",
            "date": args.date,
            "start": args.start,
            "requested_files": args.nfiles,
            "null_start": null_start,
            "null_realizations": args.nulls_per_task,
            "random_seed": args.seed,
        }
        mismatch = {
            key: (report.get(key), value)
            for key, value in expected.items()
            if report.get(key) != value
        }
        if mismatch:
            raise RuntimeError(f"existing output metadata mismatch: {mismatch}")
        print(f"SKIP complete task {args.task_id}: {json_path}")
        return
    if json_path.exists() != npz_path.exists():
        raise RuntimeError(
            f"partial output for task {args.task_id}; inspect before rerunning: "
            f"json={json_path.exists()} npz={npz_path.exists()}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    run(argparse.Namespace(
        date=args.date,
        start=args.start,
        nfiles=args.nfiles,
        null_start=null_start,
        nulls=args.nulls_per_task,
        seed=args.seed,
        null_methods=args.null_methods,
        output_dir=args.output_dir,
        stem=stem,
        overwrite=False,
    ))


if __name__ == "__main__":
    main()
