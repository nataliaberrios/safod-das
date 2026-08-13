"""Focused regression tests for the v50 F-K defensibility workflow."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ambient_lellouch2019_reproduction_v1 import (
    REFERENCE_VELOCITY_M_S,
    align_and_stack,
    geometry,
)


ROOT = Path(__file__).resolve().parent


def test_lellouch_geometry_reaches_800_m() -> None:
    dx = 1.0209523439407349
    targets, required = geometry(dx, 900)
    assert len(targets) == 16
    assert np.isclose(targets[0] * dx, 50.0, atol=0.6)
    assert np.isclose(targets[-1] * dx, 800.0, atol=0.6)
    assert required[0] == 0
    assert required[-1] < 900


def test_lellouch_alignment_preserves_fixed_moveout() -> None:
    dx = 1.0
    fs = 1000.0
    lags = np.arange(-0.35, 0.351, 1.0 / fs)
    targets, required = geometry(dx, 900)
    correlations = np.empty((required.size, lags.size), dtype=float)
    for index, channel in enumerate(required):
        arrival = channel * dx / REFERENCE_VELOCITY_M_S
        correlations[index] = np.exp(-0.5 * ((lags - arrival) / 0.003) ** 2)
    _, _, aligned = align_and_stack(
        correlations, lags, required, targets, dx
    )
    picked = lags[np.argmax(aligned, axis=1)]
    expected = targets * dx / REFERENCE_VELOCITY_M_S
    assert np.max(np.abs(picked - expected)) <= 1.0 / fs


def test_consolidated_fk_decision_and_frequency_boundary() -> None:
    report = json.loads(
        (
            ROOT
            / "ambient_transfer"
            / "fk_validation_evidence_v1"
            / "fk_validation_evidence_v1.json"
        ).read_text()
    )
    assert report["decision"]["operator_implementation_validated"] is True
    assert (
        report["decision"][
            "both_coordinate_directions_recovered_on_real_background"
        ]
        is True
    )
    assert report["decision"]["real_ambient_input_is_nearly_direction_balanced"] is True
    active = report["awd_empirical_direction"][
        "expected_to_opposite_fixed_tube_energy_ratio"
    ]
    assert active["5_20_hz"] < 2.0
    assert active["25_60_hz"] > 100.0


def test_active_report_does_not_overstate_low_frequency_direction() -> None:
    report = json.loads(
        (
            ROOT
            / "ambient_transfer"
            / "awd_production_operator_validation_v1"
            / "awd_fk_production_operator_validation_v1.json"
        ).read_text()
    )
    decision = report["decision"]
    assert decision["signed_fan_selects_known_downgoing_direction_at_5_20_hz"] is False
    assert decision["signed_fan_selects_known_downgoing_direction_at_25_60_hz"] is True
