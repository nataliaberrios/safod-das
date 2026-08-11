#!/usr/bin/env python3
"""Audit the completed corrected seasonal signed-lag F-K aggregate.

This script is intentionally independent of the aggregation job's reported
branch metrics.  It reloads the saved correlation sections, recomputes the
physical- and opposite-lag velocity scans, verifies the branch sign convention
and complete seasonal file count, and writes a compact promotion report.  It
does not decide wave type, Green's-function convergence, or formation Vp.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "ambient_transfer" / "signed_lag_v2"
MANIFEST = ROOT / "ambient_transfer" / "seasonal_day_selection.json"
JSON_PATH = OUT / "seasonal_signed_fk_v2_aggregate.json"
NPZ_PATH = OUT / "seasonal_signed_fk_v2_aggregate.npz"
FIGURE_PATH = OUT / "seasonal_signed_fk_v2_aggregate.png"
AUDIT_PATH = OUT / "seasonal_signed_fk_v2_completion_audit.json"
VELOCITIES = np.linspace(1200.0, 6000.0, 193)
BRANCH_SIGN = {"negative": 1.0, "positive": -1.0}


def velocity_scores(
    top: np.ndarray,
    lags: np.ndarray,
    distance: np.ndarray,
    lag_sign: float,
) -> np.ndarray:
    """Recompute median moveout scores without importing production code."""
    return np.asarray(
        [
            np.nanmedian(
                [
                    row[np.argmin(np.abs(lags - lag_sign * offset / velocity))]
                    for row, offset in zip(top, distance)
                ]
            )
            for velocity in VELOCITIES
        ]
    )


def close(a: float, b: float, label: str, atol: float = 1e-10) -> None:
    if not np.isclose(a, b, rtol=1e-8, atol=atol):
        raise ValueError(f"{label} mismatch: recomputed={a}, reported={b}")


def main() -> None:
    for path in (JSON_PATH, NPZ_PATH, FIGURE_PATH):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    report = json.loads(JSON_PATH.read_text())
    if report.get("workflow_version") != "signed_lag_v2":
        raise ValueError("unexpected aggregate workflow version")
    if report.get("branch_lag_sign") != BRANCH_SIGN:
        raise ValueError(f"unexpected branch signs: {report.get('branch_lag_sign')}")

    manifest = json.loads(MANIFEST.read_text())["days"]
    expected_days = [item["date"] for item in manifest]
    expected_files = int(sum(int(item["nfiles"]) for item in manifest))
    if report.get("days") != expected_days:
        raise ValueError("aggregate day list differs from frozen manifest")
    if int(report.get("weighted_files", -1)) != expected_files:
        raise ValueError(
            f"aggregate uses {report.get('weighted_files')} of {expected_files} files"
        )

    product = np.load(NPZ_PATH, allow_pickle=False)
    lags = np.asarray(product["lags"])
    distance = np.asarray(product["distance"])
    saved_velocities = np.asarray(product["velocities_m_s"])
    if not np.allclose(saved_velocities, VELOCITIES):
        raise ValueError("saved velocity grid differs from audit grid")
    if not np.all(np.diff(lags) > 0) or not np.all(np.diff(distance) > 0):
        raise ValueError("lag or receiver-offset coordinate is not increasing")

    index_3200 = int(np.argmin(abs(VELOCITIES - 3200.0)))
    audited: dict[str, dict[str, float]] = {}
    for branch, lag_sign in BRANCH_SIGN.items():
        top = np.asarray(product[f"{branch}_top"])
        if top.shape != (distance.size, lags.size):
            raise ValueError(f"unexpected {branch} top shape: {top.shape}")
        physical = velocity_scores(top, lags, distance, lag_sign)
        leakage = velocity_scores(top, lags, distance, -lag_sign)
        peak = int(np.nanargmax(np.abs(physical)))
        saved = report["branches"][branch]
        close(float(VELOCITIES[peak]), float(saved["peak_velocity_m_s"]), f"{branch} peak velocity")
        close(float(physical[peak]), float(saved["peak_signed_score"]), f"{branch} peak score")
        close(float(physical[index_3200]), float(saved["score_3200"]), f"{branch} score 3200")
        close(
            float(leakage[index_3200]),
            float(saved["opposite_lag_leakage_3200"]),
            f"{branch} leakage 3200",
        )
        audited[branch] = {
            "physical_lag_sign": float(lag_sign),
            "peak_velocity_m_s": float(VELOCITIES[peak]),
            "peak_signed_score": float(physical[peak]),
            "absolute_score_3200": float(abs(physical[index_3200])),
            "absolute_opposite_lag_leakage_3200": float(abs(leakage[index_3200])),
            "physical_to_leakage_ratio_3200": float(
                abs(physical[index_3200])
                / (abs(leakage[index_3200]) + np.finfo(float).tiny)
            ),
            "null95": float(saved["null95"]),
            "p_peak": float(saved["p_peak"]),
        }
    product.close()

    negative = audited["negative"]["absolute_score_3200"]
    positive = audited["positive"]["absolute_score_3200"]
    output = {
        "workflow_version": "signed_lag_v2_completion_audit",
        "source_workflow_version": report["workflow_version"],
        "expected_and_used_files": expected_files,
        "days": expected_days,
        "coordinate_convention": (
            "channel 0 is the top virtual source; F*K<0 is evaluated at positive "
            "lag and F*K>0 at negative lag"
        ),
        "branches": audited,
        "positive_to_negative_absolute_score_ratio_3200": float(
            positive / (negative + np.finfo(float).tiny)
        ),
        "signed_branch_asymmetry_3200": float(
            (negative - positive) / (negative + positive + np.finfo(float).tiny)
        ),
        "checks": {
            "frozen_manifest_complete": True,
            "branch_lag_signs_correct": True,
            "saved_metrics_reproduced_from_sections": True,
            "aggregate_figure_present": True,
        },
        "interpretive_boundary": (
            "The audit validates computation and signed-lag geometry. Statistical "
            "significance supports ordered moveout after F-K selection, but does not "
            "alone establish Green's-function convergence, mode identity, causality, "
            "or formation Vp."
        ),
    }
    AUDIT_PATH.write_text(json.dumps(output, indent=2))
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
