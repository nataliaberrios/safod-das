#!/usr/bin/env python3
"""Completeness gate plus aggregation for distributed pre-F-K null batches."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--glob", default="fk_full_pipeline_null_v2_*.json")
    parser.add_argument("--expected-start", type=int, required=True)
    parser.add_argument("--expected-count", type=int, required=True)
    parser.add_argument("--output-stem", default="fk_full_pipeline_null_v2_aggregate")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.expected_count < 1:
        raise ValueError("--expected-count must be positive")

    paths = sorted(args.input_dir.glob(args.glob))
    paths = [path for path in paths if path.stem != args.output_stem]
    if not paths:
        raise SystemExit("no null JSON batches found")
    reference_report = None
    reference_npz = None
    observed_reference = None
    ids = []
    for path in paths:
        paired = path.with_suffix(".npz")
        if not paired.exists():
            raise FileNotFoundError(f"partial batch, missing {paired}")
        report = json.loads(path.read_text())
        product = np.load(paired)
        if report.get("workflow_version") != "ambient_fk_full_pipeline_null_v2":
            raise ValueError(f"wrong workflow: {path}")
        count = int(report["null_realizations"])
        batch_ids = np.arange(int(report["null_start"]), int(report["null_start"]) + count)
        if product["null_realization_ids"].shape[0] != count:
            raise ValueError(f"JSON/NPZ null-count mismatch: {path}")
        ids.extend(map(int, batch_ids))
        identity = (
            report["date"], report["start"], report["requested_files"],
            report["used_files"], report["random_seed"], tuple(report["null_methods"]),
        )
        if reference_report is None:
            reference_report = identity
            reference_npz = product
            observed_reference = report["observed"]
        else:
            if identity != reference_report:
                raise ValueError(f"incompatible batch metadata: {path}")
            if report["observed"] != observed_reference:
                raise ValueError(f"observed metrics differ: {path}")
            if not np.allclose(product["lags_s"], reference_npz["lags_s"]):
                raise ValueError(f"lag axis differs: {path}")
            if not np.allclose(
                product["receiver_offsets_m"], reference_npz["receiver_offsets_m"]
            ):
                raise ValueError(f"receiver geometry differs: {path}")

    ids_array = np.asarray(sorted(ids), dtype=int)
    if np.unique(ids_array).size != ids_array.size:
        raise ValueError("duplicate null realization IDs")
    expected = np.arange(args.expected_start, args.expected_start + args.expected_count)
    if not np.array_equal(ids_array, expected):
        missing = np.setdiff1d(expected, ids_array).tolist()
        unexpected = np.setdiff1d(ids_array, expected).tolist()
        raise ValueError(f"incomplete null array: missing={missing}, unexpected={unexpected}")

    command = [
        sys.executable,
        str(Path(__file__).resolve().parent / "aggregate_ambient_fk_full_pipeline_null_v2.py"),
        "--input-dir", str(args.input_dir),
        "--glob", args.glob,
        "--output-stem", args.output_stem,
    ]
    if args.overwrite:
        command.append("--overwrite")
    subprocess.run(command, check=True)
    print(
        f"COMPLETE: aggregated {args.expected_count} distinct null IDs; "
        "the observed result was equality-checked across tasks and retained once."
    )


if __name__ == "__main__":
    main()
