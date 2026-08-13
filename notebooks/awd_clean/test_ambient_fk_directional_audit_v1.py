#!/usr/bin/env python3
"""Regression checks for the frozen ambient directional F-K audit."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPORT = (
    ROOT
    / "ambient_transfer"
    / "fk_directional_audit_v1"
    / "ambient_fk_directional_audit_v1.json"
)


def main() -> None:
    report = json.loads(REPORT.read_text())
    assert report["workflow_version"] == "ambient_fk_directional_audit_v1"

    decision = report["decision"]
    assert decision["signed_fk_implementation_validated_by_one_way_synthetics"]
    assert not decision[
        "standard_direction_only_filter_shows_independent_3p2_km_s_result"
    ]
    assert decision["production_velocity_fan_result_is_operator_conditioned"]

    held = report["mask_comparisons"]["held_out"]
    assert held["production_2p5_4p5"]["mirror_correlation"] > 0.99
    assert held["direction_only"]["mirror_correlation"] < 0.5
    assert abs(held["direction_only"]["negative_score_3200"]) < 0.02
    assert abs(held["direction_only"]["positive_score_3200"]) < 0.02

    balance = report["real_pre_filter_direction_balance"]
    assert 0.95 < balance["positive_to_negative_total_power_ratio"] < 1.05
    assert not balance["target_enrichment_passes_both_familywise_nulls"]

    legacy = report["legacy_figure45_mirror_check"]
    assert legacy["mirror_correlation"] > 0.99
    assert legacy["mirror_relative_rms"] < 0.12
    print("ambient directional F-K audit checks passed")


if __name__ == "__main__":
    main()
