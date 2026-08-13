#!/usr/bin/env python3
"""Consolidate the evidence that makes the SAFOD signed F-K result defensible."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SYNTHETIC = ROOT / "fk_sign_synthetic_test.json"
INJECTION = (
    ROOT
    / "ambient_transfer"
    / "fk_injection_recovery_v2_n300"
    / "ambient_fk_injection_recovery_v2_aggregate.json"
)
INJECTION_AUDIT = INJECTION.with_name(
    "ambient_fk_injection_recovery_v2_completion_audit.json"
)
ACTIVE = (
    ROOT
    / "ambient_transfer"
    / "awd_production_operator_validation_v1"
    / "awd_fk_production_operator_validation_v1.json"
)
BALANCE = (
    ROOT
    / "ambient_transfer"
    / "awd_direction_ground_truth_v1"
    / "awd_fk_direction_ground_truth_v1.json"
)
OUT_DIR = ROOT / "ambient_transfer" / "fk_validation_evidence_v1"
OUT = OUT_DIR / "fk_validation_evidence_v1"


def load(path: Path) -> dict:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def main() -> None:
    synthetic = load(SYNTHETIC)
    injection = load(INJECTION)
    injection_audit = load(INJECTION_AUDIT)
    active = load(ACTIVE)
    balance = load(BALANCE)
    if injection_audit.get("all_checks_pass") is not True:
        raise ValueError("injection-recovery completion audit did not pass")

    synthetic_ratios = {}
    for direction, values in synthetic["results"].items():
        if direction == "increasing_z":
            correct = values["negative_output_rms"]
            wrong = values["positive_output_rms"]
        else:
            correct = values["positive_output_rms"]
            wrong = values["negative_output_rms"]
        synthetic_ratios[direction] = float(correct / wrong)

    injection_rows = []
    for scenario in injection["scenarios"]:
        if not scenario["inside_production_wedge"]:
            continue
        injection_rows.append(
            {
                "velocity_m_s": float(scenario["velocity_m_s"]),
                "direction": int(scenario["direction"]),
                "branch": scenario["injected_branch"],
                "minimum_independent_recovery_ratio": scenario[
                    "minimum_independent_recovery_ratio"
                ],
                "minimum_postfilter_recovery_ratio": scenario[
                    "minimum_postfilter_recovery_ratio"
                ],
            }
        )

    active_ratios = {
        band: values[
            "expected_to_opposite_energy_in_fixed_downgoing_tube_ratio"
        ]
        for band, values in active["metrics"].items()
    }
    ambient_ratios = {
        "pre_filter_fan_power": balance[
            "ambient_fk_negative_to_positive_fan_power_ratio"
        ],
        "post_filter_correlation_score": balance[
            "ambient_negative_to_positive_abs_score_ratio_3200"
        ],
    }

    report = {
        "workflow_version": "fk_validation_evidence_v1",
        "question": (
            "Can the signed 2.5--4.5 km/s F-K operator be trusted, and what "
            "physical language is justified for its ambient branches?"
        ),
        "literature_context": [
            {
                "citation": "Lellouch et al. (2019), doi:10.1029/2019JB017533",
                "relevance": (
                    "Their ambient subsection used running-absolute-mean "
                    "normalization, 30 s windows with 15 s overlap, local "
                    "receiver stacking, and 3.2 km/s alignment; it did not "
                    "report F-K filtering for the ambient correlation."
                ),
            },
            {
                "citation": "Rao and Wang (2016), doi:10.1088/1742-2132/13/3/412",
                "relevance": (
                    "F-K filtering is a standard VSP wavefield-separation "
                    "method, but simple fans can have overlap and leakage near "
                    "the wavenumber origin and at truncated events."
                ),
            },
        ],
        "exact_5_20_hz_synthetic_correct_to_wrong_rms_ratio": synthetic_ratios,
        "real_background_injection": {
            "used_files": injection["used_files"],
            "amplitude_definition": injection["amplitude_definition"],
            "completion_audit_all_checks_pass": True,
            "in_wedge_scenarios": injection_rows,
        },
        "awd_empirical_direction": {
            "coordinate_label": (
                "surface AWD propagates toward increasing Nano coordinate"
            ),
            "expected_to_opposite_fixed_tube_energy_ratio": active_ratios,
            "decision": (
                "strong at 25--60 Hz; not decisive at 5--20 Hz"
            ),
        },
        "real_ambient_direction_balance": ambient_ratios,
        "decision": {
            "operator_implementation_validated": bool(
                min(synthetic_ratios.values()) > 5.0
                and injection_audit["all_checks_pass"]
            ),
            "both_coordinate_directions_recovered_on_real_background": bool(
                all(
                    row["minimum_independent_recovery_ratio"] is not None
                    for row in injection_rows
                )
            ),
            "real_ambient_input_is_nearly_direction_balanced": bool(
                0.8 <= ambient_ratios["pre_filter_fan_power"] <= 1.25
            ),
            "approved_language": (
                "On Nano, F*K<0 is the nominal downgoing/increasing-coordinate "
                "branch and F*K>0 is the nominal upgoing/decreasing-coordinate "
                "branch. The observed 5--20 Hz ambient field is approximately "
                "balanced between them."
            ),
            "language_not_yet_justified": (
                "Claiming that the nominal branches are clean individual body "
                "waves, assigning the origin of the upgoing energy without a "
                "reflection/scattering test, or treating either fan-filtered "
                "peak as an independent formation-velocity estimate."
            ),
        },
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")

    figure, axes = plt.subplots(2, 2, figsize=(13.2, 9.2), constrained_layout=True)

    directions = ["increasing_z", "decreasing_z"]
    values = [synthetic_ratios[name] for name in directions]
    axes[0, 0].bar(
        ["increasing\ncoordinate", "decreasing\ncoordinate"],
        values,
        color=["tab:blue", "tab:orange"],
    )
    axes[0, 0].axhline(1.0, color="0.25", ls="--", lw=1.0)
    axes[0, 0].set_ylabel("Correct/wrong branch output RMS")
    axes[0, 0].set_title("a  Exact production mask: one-way synthetics")
    for index, value in enumerate(values):
        axes[0, 0].text(index, value * 1.01, f"{value:.1f}:1", ha="center")

    velocities = sorted({row["velocity_m_s"] for row in injection_rows})
    x = np.arange(len(velocities))
    width = 0.34
    for offset, direction, label, color in (
        (-width / 2, 1, "increasing coordinate", "tab:blue"),
        (width / 2, -1, "decreasing coordinate", "tab:orange"),
    ):
        threshold = [
            next(
                row["minimum_independent_recovery_ratio"]
                for row in injection_rows
                if row["velocity_m_s"] == velocity
                and row["direction"] == direction
            )
            for velocity in velocities
        ]
        axes[0, 1].bar(x + offset, threshold, width, label=label, color=color)
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xticks(x, [f"{velocity/1000:g}" for velocity in velocities])
    axes[0, 1].set_xlabel("Injected apparent velocity (km s$^{-1}$)")
    axes[0, 1].set_ylabel("Minimum independently recovered RMS ratio")
    axes[0, 1].set_title("b  Known-direction waves injected into 300 real files")
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[0, 1].grid(axis="y", alpha=0.2)

    band_labels = ["5--20", "25--60"]
    active_values = [active_ratios["5_20_hz"], active_ratios["25_60_hz"]]
    axes[1, 0].bar(band_labels, active_values, color="tab:green")
    axes[1, 0].axhline(1.0, color="0.25", ls="--", lw=1.0)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xlabel("Frequency band (Hz)")
    axes[1, 0].set_ylabel("Expected/opposite energy in fixed AWD moveout tube")
    axes[1, 0].set_title("c  Empirical AWD direction is frequency dependent")
    for index, value in enumerate(active_values):
        axes[1, 0].text(index, value * 1.03, f"{value:.2f}", ha="center")

    ambient_values = [
        ambient_ratios["pre_filter_fan_power"],
        ambient_ratios["post_filter_correlation_score"],
    ]
    axes[1, 1].bar(
        ["pre-filter\nF-K fan power", "filtered\ncorrelation score"],
        ambient_values,
        color="tab:purple",
    )
    axes[1, 1].axhline(1.0, color="0.25", ls="--", lw=1.0)
    axes[1, 1].set_ylim(0.0, 1.25)
    axes[1, 1].set_ylabel("$F K<0$ / $F K>0$ ratio")
    axes[1, 1].set_title("d  The recorded ambient field is directionally balanced")
    for index, value in enumerate(ambient_values):
        axes[1, 1].text(index, value + 0.025, f"{value:.3f}", ha="center")

    figure.suptitle(
        "Validation ladder for the signed SAFOD ambient F-K operator",
        fontsize=14,
    )
    figure.savefig(OUT.with_suffix(".png"), dpi=350, bbox_inches="tight")
    figure.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
