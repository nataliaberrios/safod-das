#!/usr/bin/env python3
"""Aggregate disjoint pre-F-K null batches into one empirical distribution."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ambient_fk_full_pipeline_null_v2 import VELOCITIES, empirical_p


def same_value(first: float, second: float) -> bool:
    return bool(np.isclose(first, second, rtol=1e-10, atol=1e-12))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--glob", default="fk_full_pipeline_null_v2_*.json")
    parser.add_argument("--output-stem", default="fk_full_pipeline_null_v2_aggregate")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    json_paths = sorted(args.input_dir.glob(args.glob))
    json_paths = [path for path in json_paths if path.stem != args.output_stem]
    if not json_paths:
        raise SystemExit(f"no inputs matched {args.input_dir / args.glob}")

    records = []
    for json_path in json_paths:
        npz_path = json_path.with_suffix(".npz")
        if not npz_path.exists():
            raise FileNotFoundError(f"missing paired NPZ: {npz_path}")
        report = json.loads(json_path.read_text())
        if report.get("workflow_version") != "ambient_fk_full_pipeline_null_v2":
            raise ValueError(f"unexpected workflow in {json_path}")
        records.append((report, np.load(npz_path)))

    reference = records[0][0]
    methods = tuple(reference["null_methods"])
    identity_keys = ("date", "start", "requested_files", "used_files", "random_seed")
    seen_ids = set()
    ordered = []
    for report, product in records:
        for key in identity_keys:
            if report[key] != reference[key]:
                raise ValueError(f"incompatible {key}: {report[key]} != {reference[key]}")
        if tuple(report["null_methods"]) != methods:
            raise ValueError("all batches must contain the same null methods")
        for mode in ("negative", "positive"):
            for metric in (
                "peak_absolute_score_in_wedge",
                "peak_velocity_m_s_in_wedge",
                "signed_score_3200",
            ):
                if not same_value(
                    report["observed"][mode][metric],
                    reference["observed"][mode][metric],
                ):
                    raise ValueError(f"observed result mismatch for {mode}/{metric}")
        null_ids = np.arange(
            int(report["null_start"]),
            int(report["null_start"]) + int(report["null_realizations"]),
        )
        overlap = seen_ids.intersection(map(int, null_ids))
        if overlap:
            raise ValueError(f"duplicate null realization ids: {sorted(overlap)}")
        seen_ids.update(map(int, null_ids))
        ordered.append((int(null_ids[0]), null_ids, report, product))
    ordered.sort(key=lambda item: item[0])

    null_ids = np.concatenate([item[1] for item in ordered])
    index_3200 = int(np.argmin(abs(VELOCITIES - 3200.0)))
    aggregate_report = {
        "workflow_version": "ambient_fk_full_pipeline_null_v2_aggregate",
        "source_batches": [path.name for path in json_paths],
        "date": reference["date"],
        "start": reference["start"],
        "requested_files": reference["requested_files"],
        "used_files": reference["used_files"],
        "random_seed": reference["random_seed"],
        "null_realizations": int(null_ids.size),
        "null_realization_ids": null_ids.tolist(),
        "null_methods": list(methods),
        "observed": reference["observed"],
        "observed_familywise_score": reference["observed_familywise_score"],
        "null_results": {},
        "supported_claim": reference["supported_claim"],
        "limitations": reference["limitations"],
    }
    arrays = {
        "null_realization_ids": null_ids,
        "velocities_m_s": VELOCITIES,
        "lags_s": ordered[0][3]["lags_s"],
        "receiver_offsets_m": ordered[0][3]["receiver_offsets_m"],
    }
    for mode in ("negative", "positive"):
        arrays[f"observed_{mode}_top"] = ordered[0][3][f"observed_{mode}_top"]
        arrays[f"observed_{mode}_scores"] = ordered[0][3][f"observed_{mode}_scores"]

    for method in methods:
        method_report = {"branches": {}}
        for mode in ("negative", "positive"):
            score_curves = np.concatenate([
                item[3][f"null_{method}_{mode}_scores"] for item in ordered
            ])
            statistics = np.concatenate([
                item[3][f"null_{method}_{mode}_selection_statistics"]
                for item in ordered
            ])
            observed_peak = float(reference["observed"][mode]["peak_absolute_score_in_wedge"])
            observed_3200 = float(reference["observed"][mode]["absolute_score_3200"])
            fixed_3200 = np.abs(score_curves[:, index_3200])
            method_report["branches"][mode] = {
                "observed_peak_absolute_score_in_wedge": observed_peak,
                "observed_peak_velocity_m_s_in_wedge": float(
                    reference["observed"][mode]["peak_velocity_m_s_in_wedge"]
                ),
                "null95_peak_score_in_wedge": float(np.quantile(statistics, 0.95)),
                "p_peak_in_wedge": empirical_p(observed_peak, statistics),
                "null95_absolute_score_3200": float(np.quantile(fixed_3200, 0.95)),
                "p_absolute_score_3200": empirical_p(observed_3200, fixed_3200),
            }
            arrays[f"null_{method}_{mode}_scores"] = score_curves
            arrays[f"null_{method}_{mode}_selection_statistics"] = statistics
        familywise = np.concatenate([
            item[3][f"null_{method}_familywise_statistics"] for item in ordered
        ])
        observed_familywise = float(reference["observed_familywise_score"])
        method_report["familywise"] = {
            "observed_maximum_across_branches_and_wedge": observed_familywise,
            "null95": float(np.quantile(familywise, 0.95)),
            "p": empirical_p(observed_familywise, familywise),
        }
        arrays[f"null_{method}_familywise_statistics"] = familywise
        aggregate_report["null_results"][method] = method_report

    json_out = args.input_dir / f"{args.output_stem}.json"
    npz_out = args.input_dir / f"{args.output_stem}.npz"
    if not args.overwrite and (json_out.exists() or npz_out.exists()):
        raise FileExistsError("aggregate output exists; pass --overwrite or change --output-stem")
    np.savez_compressed(npz_out, **arrays)
    json_out.write_text(json.dumps(aggregate_report, indent=2))
    print(json.dumps(aggregate_report, indent=2))


if __name__ == "__main__":
    main()
