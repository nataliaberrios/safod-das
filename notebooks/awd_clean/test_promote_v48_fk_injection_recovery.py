#!/usr/bin/env python3
"""Temporary-copy integration test for the v48 injection promoter."""
from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import promote_v48_fk_injection_recovery as promote


def synthetic_summary() -> dict:
    scenarios = []
    for velocity in (1800.0, 2750.0, 3200.0, 4000.0):
        for direction in (1, -1):
            in_wedge = velocity >= 2500.0
            points = []
            for amplitude in promote.EXPECTED_AMPLITUDES:
                recovered = in_wedge and amplitude >= 0.1
                points.append({
                    "amplitude_ratio": amplitude,
                    "prefilter_uplift_ci95": [0.01 if recovered else -0.01, 0.03],
                    "postfilter_score_uplift_ci95": [0.01 if recovered else -0.01, 0.03],
                    "independent_recovery": recovered,
                    "postfilter_recovery": recovered,
                })
            scenarios.append({
                "velocity_m_s": velocity,
                "direction": direction,
                "inside_production_wedge": in_wedge,
                "minimum_independent_recovery_ratio": 0.1 if in_wedge else None,
                "minimum_postfilter_recovery_ratio": 0.1 if in_wedge else None,
                "amplitudes": points,
            })
    return {
        "workflow_version": "ambient_fk_injection_recovery_v2_aggregate",
        "used_files": 300,
        "scenarios": scenarios,
    }


def test_promotion() -> None:
    temporary = Path(tempfile.mkdtemp(prefix="safod_v48_promotion_", dir="/tmp"))
    notebook = temporary / "AWD_results_dashboard.ipynb"
    tex = temporary / "AWD_advisor_figure_guide.tex"
    tex_v48 = temporary / "AWD_advisor_figure_guide_v48.tex"
    product = temporary / "fk_injection_recovery_v2_n300"
    product.mkdir()
    source_notebook = json.loads(promote.NOTEBOOK.read_text())
    assert source_notebook["metadata"]["awd_dashboard"]["version"] == "v48"
    assert "# v48" in "".join(source_notebook["cells"][-4]["source"])
    source_notebook["cells"] = source_notebook["cells"][:-4]
    source_notebook["metadata"]["awd_dashboard"]["version"] = "v47"
    source_notebook["metadata"]["awd_dashboard"]["last_status_sync"] = "2026-08-11"
    notebook.write_text(json.dumps(source_notebook, indent=1) + "\n")

    source_tex = promote.TEX.read_text()
    section = r"\section{Ambient F--K injection--recovery calibration (v48)}"
    assert source_tex.count(section) == 1 and source_tex.count(r"\end{document}") == 1
    source_tex = source_tex[:source_tex.index(section)] + "\\end{document}\n"
    source_tex = source_tex.replace("version v48", "version v47", 1)
    tex.write_text(source_tex)
    aggregate = product / "ambient_fk_injection_recovery_v2_aggregate.json"
    audit = product / "ambient_fk_injection_recovery_v2_completion_audit.json"
    figure = product / "ambient_fk_injection_recovery_v2_aggregate.png"
    aggregate.write_text(json.dumps(synthetic_summary()))
    audit.write_text(json.dumps({
        "workflow_version": "ambient_fk_injection_recovery_v2_completion_audit",
        "expected_files": 300,
        "all_checks_pass": True,
        "checks": {"synthetic_integration_check": True},
    }))
    figure.write_bytes(b"synthetic figure placeholder")
    original = {
        name: getattr(promote, name)
        for name in ("NOTEBOOK", "TEX", "TEX_V48", "OUT", "AGGREGATE", "AUDIT", "FIGURE")
    }
    try:
        promote.NOTEBOOK = notebook
        promote.TEX = tex
        promote.TEX_V48 = tex_v48
        promote.OUT = product
        promote.AGGREGATE = aggregate
        promote.AUDIT = audit
        promote.FIGURE = figure
        promote.main()
        updated = json.loads(notebook.read_text())
        assert updated["metadata"]["awd_dashboard"]["version"] == "v48"
        assert len(updated["cells"]) == 197
        combined = "\n".join("".join(cell.get("source", [])) for cell in updated["cells"])
        assert "# v48" in combined
        assert "Real-noise plane-wave injection" in combined
        assert "The 1.8 km/s off-wedge control never passed" in combined
        tex_text = tex.read_text()
        assert "version v48" in tex_text
        assert "injection--recovery calibration (v48)" in tex_text
        assert tex_v48.read_text() == tex_text
        try:
            promote.main()
        except RuntimeError as error:
            assert "already v48" in str(error)
        else:
            raise AssertionError("second promotion did not refuse duplicate v48")
    finally:
        for name, value in original.items():
            setattr(promote, name, value)


if __name__ == "__main__":
    test_promotion()
    print("PASS test_promotion")
