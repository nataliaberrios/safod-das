#!/usr/bin/env python3
"""Aggregate the frozen F-K mask audit into development/held-out products.

The primary comparison is generalization from 2024-12-20 (development) to
seven independently selected days (held out).  Reported velocity peaks remain
conditional on each F-K mask; no mask is selected from these results.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ambient_fk_mask_sensitivity_v2 import (
    BRANCH_LAG_SIGN,
    MASK_SPECS,
)


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "ambient_transfer" / "seasonal_day_selection.json"
OUT = ROOT / "ambient_transfer" / "fk_mask_sensitivity_v2"
DEVELOPMENT_DATE = "2024-12-20"
VELOCITIES = np.linspace(1200.0, 6000.0, 193)
REFERENCE_VELOCITY = 3200.0
CHUNK = 10


def velocity_scores(top, lags, distance, lag_sign, distances_override=None):
    used_distance = distance if distances_override is None else distances_override
    return np.array(
        [
            np.nanmedian(
                [
                    row[
                        np.argmin(
                            np.abs(lags - lag_sign * offset / velocity)
                        )
                    ]
                    for row, offset in zip(top, used_distance)
                ]
            )
            for velocity in VELOCITIES
        ]
    )


def summarize(top, lags, distance, lag_sign, permutations, seed):
    """Return conditional ridge scores and a receiver-order null.

    This null tests ordered moveout *after* the selected mask.  It is not a
    full-pipeline F-K selection null and is labelled accordingly downstream.
    """
    scores = velocity_scores(top, lags, distance, lag_sign)
    peak_index = int(np.nanargmax(np.abs(scores)))
    reference_index = int(np.argmin(np.abs(VELOCITIES - REFERENCE_VELOCITY)))
    rng = np.random.default_rng(seed)
    fixed_null = np.empty(permutations, dtype=float)
    peak_null = np.empty(permutations, dtype=float)
    for index in range(permutations):
        permuted = rng.permutation(distance)
        null_scores = velocity_scores(
            top, lags, distance, lag_sign, distances_override=permuted
        )
        fixed_null[index] = abs(null_scores[reference_index])
        peak_null[index] = np.nanmax(np.abs(null_scores))
    return {
        "scores": scores,
        "fixed_null": fixed_null,
        "peak_null": peak_null,
        "conditional_peak_velocity_m_s": float(VELOCITIES[peak_index]),
        "conditional_peak_signed_score": float(scores[peak_index]),
        "conditional_peak_absolute_score": float(abs(scores[peak_index])),
        "signed_score_3200": float(scores[reference_index]),
        "absolute_score_3200": float(abs(scores[reference_index])),
        "conditional_null95_fixed_3200": float(np.quantile(fixed_null, 0.95)),
        "conditional_p_fixed_3200": float(
            (1 + np.sum(fixed_null >= abs(scores[reference_index])))
            / (permutations + 1)
        ),
        "conditional_null95_peak_scan": float(np.quantile(peak_null, 0.95)),
        "conditional_p_peak_scan": float(
            (1 + np.sum(peak_null >= abs(scores[peak_index])))
            / (permutations + 1)
        ),
    }


def load_day(date, nfiles, allow_partial):
    weighted_sums = {
        (mask_name, branch): None
        for mask_name in MASK_SPECS
        for branch in BRANCH_LAG_SIGN
    }
    used_counts = {mask_name: 0 for mask_name in MASK_SPECS}
    found_chunks = 0
    expected_chunks = len(range(0, nfiles, CHUNK))
    lags = distance = None
    missing = []
    for start in range(0, nfiles, CHUNK):
        count = min(CHUNK, nfiles - start)
        stem = OUT / f"fk_mask_v2_{date}_start{start}_n{count}"
        json_path = stem.with_suffix(".json")
        npz_path = stem.with_suffix(".npz")
        if not json_path.exists() or not npz_path.exists():
            missing.append({"start": start, "nfiles": count})
            continue
        report = json.loads(json_path.read_text())
        product = np.load(npz_path, allow_pickle=False)
        if report.get("workflow_version") != "ambient_fk_mask_sensitivity_v2":
            raise RuntimeError(f"Unexpected workflow version in {json_path}")
        if lags is None:
            lags = np.asarray(product["lags"])
            distance = np.asarray(product["distance"])
        elif not (
            np.array_equal(lags, product["lags"])
            and np.allclose(distance, product["distance"])
        ):
            raise RuntimeError(f"Axis mismatch in {npz_path}")
        for mask_name in MASK_SPECS:
            used = int(report["mask_products"][mask_name]["used_files"])
            for branch in BRANCH_LAG_SIGN:
                top = np.asarray(product[f"{mask_name}__{branch}_top"])
                key = (mask_name, branch)
                contribution = top * used
                weighted_sums[key] = (
                    contribution
                    if weighted_sums[key] is None
                    else weighted_sums[key] + contribution
                )
            used_counts[mask_name] += used
        found_chunks += 1

    if missing and not allow_partial:
        raise RuntimeError(
            f"{date}: {len(missing)}/{expected_chunks} chunks missing; "
            "rerun with --allow-partial only for diagnostics"
        )
    if found_chunks == 0:
        return None
    tops = {
        mask_name: {
            branch: weighted_sums[(mask_name, branch)] / used_counts[mask_name]
            for branch in BRANCH_LAG_SIGN
        }
        for mask_name in MASK_SPECS
    }
    return {
        "date": date,
        "split_role": "development" if date == DEVELOPMENT_DATE else "held_out",
        "tops": tops,
        "lags": lags,
        "distance": distance,
        "used_files_by_mask": used_counts,
        "found_chunks": found_chunks,
        "expected_chunks": expected_chunks,
        "missing_chunks": missing,
    }


def equal_day_stack(day_products, dates, mask_name, branch):
    selected = [item for item in day_products if item["date"] in dates]
    return np.mean(
        [item["tops"][mask_name][branch] for item in selected], axis=0
    )


def clean_metrics(result):
    return {
        key: value
        for key, value in result.items()
        if key not in {"scores", "fixed_null", "peak_null"}
    }


def plot_summary(group_results, day_results, heldout_dates):
    colors = {
        "production_2p5_4p5": "#111111",
        "narrow_2p8_3p8": "#0072B2",
        "broad_2p0_5p5": "#D55E00",
        "direction_only": "#009E73",
    }
    labels = {
        "production_2p5_4p5": "2.5--4.5 km/s production",
        "narrow_2p8_3p8": "2.8--3.8 km/s narrow",
        "broad_2p0_5p5": "2.0--5.5 km/s broad",
        "direction_only": "direction only",
    }
    fig, axes = plt.subplots(
        2, 3, figsize=(15, 8.5), constrained_layout=True, sharex="col"
    )
    for row, branch in enumerate(("negative", "positive")):
        for column, group in enumerate(("development", "held_out")):
            ax = axes[row, column]
            for mask_name in MASK_SPECS:
                result = group_results[group][mask_name][branch]
                ax.plot(
                    VELOCITIES / 1000.0,
                    np.abs(result["scores"]),
                    color=colors[mask_name],
                    lw=2.2 if mask_name == "production_2p5_4p5" else 1.5,
                    label=labels[mask_name],
                )
            ax.axvline(REFERENCE_VELOCITY / 1000.0, color="0.5", ls=":")
            ax.set_title(f"{group.replace('_', ' ').title()} — {branch} branch")
            ax.set_ylabel("Absolute median correlation")
            ax.grid(alpha=0.2)
            if row == 1:
                ax.set_xlabel("Conditional trial apparent velocity (km/s)")
            if row == 0 and column == 0:
                ax.legend(frameon=False, fontsize=8)

        ax = axes[row, 2]
        x = np.arange(len(heldout_dates))
        for mask_name in MASK_SPECS:
            values = [
                day_results[date][mask_name][branch]["absolute_score_3200"]
                for date in heldout_dates
            ]
            ax.plot(
                x,
                values,
                marker="o",
                ms=4,
                color=colors[mask_name],
                label=labels[mask_name],
            )
        ax.set_xticks(x, [date[5:] for date in heldout_dates], rotation=45)
        ax.set_title(f"Held-out days at 3.2 km/s — {branch}")
        ax.set_ylabel("Absolute median correlation")
        ax.grid(alpha=0.2)
        if row == 1:
            ax.set_xlabel("Held-out date (MM-DD)")

    fig.suptitle(
        "Frozen F-K mask sensitivity: post-filter ridge diagnostics (not an unbiased velocity inversion)",
        fontsize=14,
    )
    fig.savefig(OUT / "ambient_fk_mask_sensitivity_v2.png", dpi=350)
    fig.savefig(OUT / "ambient_fk_mask_sensitivity_v2.pdf")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--permutations", type=int, default=1000)
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    manifest_dates = [item["date"] for item in manifest["days"]]
    if DEVELOPMENT_DATE not in manifest_dates:
        raise RuntimeError(f"Development date {DEVELOPMENT_DATE} not in manifest")
    heldout_dates = [date for date in manifest_dates if date != DEVELOPMENT_DATE]
    if len(heldout_dates) != 7:
        raise RuntimeError(
            f"Expected seven held-out dates; manifest contains {len(heldout_dates)}"
        )

    day_products = []
    for item in manifest["days"]:
        product = load_day(
            item["date"], int(item["nfiles"]), args.allow_partial
        )
        if product is not None:
            day_products.append(product)
    available_dates = [item["date"] for item in day_products]
    if DEVELOPMENT_DATE not in available_dates:
        raise RuntimeError("Development day has no available audit chunks")
    available_heldout = [
        date for date in heldout_dates if date in available_dates
    ]
    if not args.allow_partial and available_heldout != heldout_dates:
        raise RuntimeError("Not all seven held-out days are available")
    if not available_heldout:
        raise RuntimeError("No held-out day is available")

    lags = day_products[0]["lags"]
    distance = day_products[0]["distance"]
    group_dates = {
        "development": [DEVELOPMENT_DATE],
        "held_out": available_heldout,
    }
    group_results = {group: {} for group in group_dates}
    arrays = {
        "lags": lags,
        "distance": distance,
        "velocities_m_s": VELOCITIES,
    }
    seed = 202608050
    for group, dates in group_dates.items():
        for mask_index, mask_name in enumerate(MASK_SPECS):
            group_results[group][mask_name] = {}
            for branch_index, (branch, lag_sign) in enumerate(
                BRANCH_LAG_SIGN.items()
            ):
                top = equal_day_stack(
                    day_products, dates, mask_name, branch
                )
                result = summarize(
                    top,
                    lags,
                    distance,
                    lag_sign,
                    args.permutations,
                    seed + 1000 * (group == "held_out") + 10 * mask_index + branch_index,
                )
                group_results[group][mask_name][branch] = result
                prefix = f"{group}__{mask_name}__{branch}"
                arrays[f"{prefix}__top"] = top
                arrays[f"{prefix}__scores"] = result["scores"]
                arrays[f"{prefix}__fixed_null"] = result["fixed_null"]
                arrays[f"{prefix}__peak_null"] = result["peak_null"]

    day_results = {}
    for day_index, item in enumerate(day_products):
        date = item["date"]
        day_results[date] = {}
        for mask_index, mask_name in enumerate(MASK_SPECS):
            day_results[date][mask_name] = {}
            for branch_index, (branch, lag_sign) in enumerate(
                BRANCH_LAG_SIGN.items()
            ):
                result = summarize(
                    item["tops"][mask_name][branch],
                    lags,
                    distance,
                    lag_sign,
                    max(100, min(args.permutations, 250)),
                    seed + 10000 + 100 * day_index + 10 * mask_index + branch_index,
                )
                day_results[date][mask_name][branch] = clean_metrics(result)

    generalization = {}
    for mask_name in MASK_SPECS:
        generalization[mask_name] = {}
        for branch in BRANCH_LAG_SIGN:
            development_score = group_results["development"][mask_name][branch][
                "absolute_score_3200"
            ]
            heldout_score = group_results["held_out"][mask_name][branch][
                "absolute_score_3200"
            ]
            heldout_values = np.array(
                [
                    day_results[date][mask_name][branch]["absolute_score_3200"]
                    for date in available_heldout
                ]
            )
            generalization[mask_name][branch] = {
                "development_absolute_score_3200": development_score,
                "heldout_equal_day_absolute_score_3200": heldout_score,
                "heldout_day_median_absolute_score_3200": float(
                    np.median(heldout_values)
                ),
                "heldout_day_min_absolute_score_3200": float(
                    np.min(heldout_values)
                ),
                "heldout_day_max_absolute_score_3200": float(
                    np.max(heldout_values)
                ),
                "heldout_to_development_ratio_3200": float(
                    heldout_score / development_score
                ) if development_score > 0 else None,
            }

    report = {
        "workflow_version": "ambient_fk_mask_sensitivity_v2_aggregate",
        "development_date": DEVELOPMENT_DATE,
        "heldout_dates": heldout_dates,
        "available_heldout_dates": available_heldout,
        "split_is_frozen": True,
        "manifest": str(MANIFEST),
        "manifest_selection_criteria": manifest.get("criteria"),
        "predeclared_masks": MASK_SPECS,
        "reference_velocity_m_s": REFERENCE_VELOCITY,
        "velocity_statement": "Velocity peaks are post-filter conditional ridge diagnostics; they are not independent or unbiased estimates of formation velocity.",
        "null_statement": "Receiver-order permutations test ordered moveout after F-K selection; they are not full-pipeline mask-selection nulls.",
        "aggregation": "Each date is stacked internally by file count; held-out dates then receive equal weight.",
        "permutations": args.permutations,
        "partial": args.allow_partial,
        "day_completeness": {
            item["date"]: {
                "found_chunks": item["found_chunks"],
                "expected_chunks": item["expected_chunks"],
                "missing_chunks": item["missing_chunks"],
                "used_files_by_mask": item["used_files_by_mask"],
            }
            for item in day_products
        },
        "group_results": {
            group: {
                mask_name: {
                    branch: clean_metrics(result)
                    for branch, result in branches.items()
                }
                for mask_name, branches in masks.items()
            }
            for group, masks in group_results.items()
        },
        "day_results": day_results,
        "generalization_at_3200_m_s": generalization,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "ambient_fk_mask_sensitivity_v2.json").write_text(
        json.dumps(report, indent=2)
    )
    np.savez_compressed(
        OUT / "ambient_fk_mask_sensitivity_v2.npz", **arrays
    )
    plot_summary(group_results, day_results, available_heldout)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
