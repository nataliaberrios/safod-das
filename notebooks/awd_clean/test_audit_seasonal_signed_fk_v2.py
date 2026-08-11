#!/usr/bin/env python3
"""Synthetic integration test for the signed-lag completion auditor."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np

import audit_seasonal_signed_fk_v2 as audit


def branch_report(
    top: np.ndarray,
    lags: np.ndarray,
    distance: np.ndarray,
    sign: float,
) -> dict[str, float]:
    physical = audit.velocity_scores(top, lags, distance, sign)
    leakage = audit.velocity_scores(top, lags, distance, -sign)
    peak = int(np.nanargmax(np.abs(physical)))
    index_3200 = int(np.argmin(abs(audit.VELOCITIES - 3200.0)))
    return {
        "peak_velocity_m_s": float(audit.VELOCITIES[peak]),
        "peak_signed_score": float(physical[peak]),
        "score_3200": float(physical[index_3200]),
        "opposite_lag_leakage_3200": float(leakage[index_3200]),
        "null95": 0.1,
        "p_peak": 0.002,
    }


def test_audit() -> None:
    temporary = Path(tempfile.mkdtemp(prefix="safod_signed_audit_", dir="/tmp"))
    manifest = temporary / "manifest.json"
    report_path = temporary / "aggregate.json"
    product_path = temporary / "aggregate.npz"
    figure_path = temporary / "aggregate.png"
    output_path = temporary / "audit.json"

    lags = np.linspace(-0.35, 0.35, 701)
    distance = np.arange(50.0, 701.0, 50.0)
    width = 0.003
    velocity = 3200.0
    negative_top = np.stack(
        [np.exp(-0.5 * ((lags - offset / velocity) / width) ** 2) for offset in distance]
    )
    positive_top = np.stack(
        [np.exp(-0.5 * ((lags + offset / velocity) / width) ** 2) for offset in distance]
    )
    report = {
        "workflow_version": "signed_lag_v2",
        "days": ["synthetic-day"],
        "weighted_files": 10,
        "branch_lag_sign": audit.BRANCH_SIGN,
        "branches": {
            "negative": branch_report(negative_top, lags, distance, 1.0),
            "positive": branch_report(positive_top, lags, distance, -1.0),
        },
    }
    manifest.write_text(
        json.dumps({"days": [{"date": "synthetic-day", "nfiles": 10}]})
    )
    report_path.write_text(json.dumps(report))
    np.savez_compressed(
        product_path,
        lags=lags,
        distance=distance,
        velocities_m_s=audit.VELOCITIES,
        negative_top=negative_top,
        positive_top=positive_top,
    )
    figure_path.write_bytes(b"synthetic figure placeholder")

    original = {
        name: getattr(audit, name)
        for name in (
            "MANIFEST",
            "JSON_PATH",
            "NPZ_PATH",
            "FIGURE_PATH",
            "AUDIT_PATH",
        )
    }
    try:
        audit.MANIFEST = manifest
        audit.JSON_PATH = report_path
        audit.NPZ_PATH = product_path
        audit.FIGURE_PATH = figure_path
        audit.AUDIT_PATH = output_path
        audit.main()
        result = json.loads(output_path.read_text())
        assert all(result["checks"].values())
        assert np.isclose(
            result["positive_to_negative_absolute_score_ratio_3200"], 1.0
        )
        assert abs(result["signed_branch_asymmetry_3200"]) < 1e-12

        bad = json.loads(report_path.read_text())
        bad["branches"]["positive"]["score_3200"] += 0.1
        report_path.write_text(json.dumps(bad))
        try:
            audit.main()
        except ValueError as error:
            assert "positive score 3200 mismatch" in str(error)
        else:
            raise AssertionError("auditor accepted inconsistent positive score")
    finally:
        for name, value in original.items():
            setattr(audit, name, value)


if __name__ == "__main__":
    test_audit()
    print("PASS test_audit")
