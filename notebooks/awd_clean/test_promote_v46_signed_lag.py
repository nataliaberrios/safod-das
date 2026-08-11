#!/usr/bin/env python3
"""Temporary-copy integration test for the v46 signed-lag promoter."""
from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import promote_v46_signed_lag as promote


def test_promotion() -> None:
    temporary = Path(tempfile.mkdtemp(prefix="safod_v46_promotion_", dir="/tmp"))
    notebook = temporary / "AWD_results_dashboard.ipynb"
    tex = temporary / "AWD_advisor_figure_guide.tex"
    tex_v46 = temporary / "AWD_advisor_figure_guide_v46.tex"
    product = temporary / "signed_lag_v2"
    product.mkdir()
    shutil.copy2(promote.NOTEBOOK, notebook)
    shutil.copy2(promote.TEX, tex)

    aggregate = {
        "workflow_version": "signed_lag_v2",
        "weighted_files": 11523,
        "branch_lag_sign": {"negative": 1.0, "positive": -1.0},
    }
    branch_common = {
        "absolute_score_3200": 0.20,
        "physical_to_leakage_ratio_3200": 8.0,
        "null95": 0.05,
        "p_peak": 0.002,
    }
    audit = {
        "workflow_version": "signed_lag_v2_completion_audit",
        "expected_and_used_files": 11523,
        "days": [f"date-{index}" for index in range(8)],
        "checks": {
            "frozen_manifest_complete": True,
            "branch_lag_signs_correct": True,
            "saved_metrics_reproduced_from_sections": True,
            "aggregate_figure_present": True,
        },
        "branches": {
            "negative": {
                **branch_common,
                "peak_velocity_m_s": 3075.0,
            },
            "positive": {
                **branch_common,
                "peak_velocity_m_s": 3125.0,
                "absolute_score_3200": 0.18,
                "physical_to_leakage_ratio_3200": 7.0,
            },
        },
        "positive_to_negative_absolute_score_ratio_3200": 0.9,
        "signed_branch_asymmetry_3200": 0.0526315789,
    }
    aggregate_path = product / "seasonal_signed_fk_v2_aggregate.json"
    audit_path = product / "seasonal_signed_fk_v2_completion_audit.json"
    figure_path = product / "seasonal_signed_fk_v2_aggregate.png"
    aggregate_path.write_text(json.dumps(aggregate))
    audit_path.write_text(json.dumps(audit))
    figure_path.write_bytes(b"synthetic figure placeholder")

    original = {
        name: getattr(promote, name)
        for name in ("NOTEBOOK", "TEX", "TEX_V46", "AGGREGATE", "AUDIT", "FIGURE")
    }
    try:
        promote.NOTEBOOK = notebook
        promote.TEX = tex
        promote.TEX_V46 = tex_v46
        promote.AGGREGATE = aggregate_path
        promote.AUDIT = audit_path
        promote.FIGURE = figure_path
        promote.main()

        updated = json.loads(notebook.read_text())
        assert updated["metadata"]["awd_dashboard"]["version"] == "v46"
        assert len(updated["cells"]) == 193
        combined = "\n".join(
            "".join(cell.get("source", [])) for cell in updated["cells"]
        )
        normalized = " ".join(combined.split())
        assert "positive/negative absolute-score ratio at 3.2 km/s is 0.900" in normalized
        assert "independently validated Green's functions" in normalized
        assert "# v46" in combined

        tex_text = tex.read_text()
        assert "version v46" in tex_text
        assert "Corrected eight-day signed-lag comparison (v46)" in tex_text
        assert "0.900" in tex_text
        assert tex_v46.read_text() == tex_text

        try:
            promote.main()
        except RuntimeError as error:
            assert "already v46" in str(error)
        else:
            raise AssertionError("second promotion did not refuse duplicate v46")
    finally:
        for name, value in original.items():
            setattr(promote, name, value)


if __name__ == "__main__":
    test_promotion()
    print("PASS test_promotion")
