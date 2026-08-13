#!/usr/bin/env python3
"""Calibrate ambient signed F-K directions with the known surface AWD source.

The surface AWD is a ground-truth direction experiment: on the cemented Nano
fiber, whose channel coordinate increases away from the top of the well, the
first coherent source arrival must propagate toward increasing coordinate.
This script compares that known-direction active wavefield with the signed
power balance of the 2024--2025 ambient archive.  It uses frozen products and
does not refilter or modify raw waveforms.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
ACTIVE = ROOT / "fk_dispersion.npz"
AMBIENT = (
    ROOT
    / "ambient_transfer"
    / "fk_prefilter_energy_v1_n300_r20"
    / "fk_prefilter_energy_v1_aggregate.npz"
)
SIGNED = (
    ROOT
    / "ambient_transfer"
    / "signed_lag_v2"
    / "seasonal_signed_fk_v2_aggregate.npz"
)
SIGNED_JSON = SIGNED.with_suffix(".json")
OUT_DIR = ROOT / "ambient_transfer" / "awd_direction_ground_truth_v1"
OUT = OUT_DIR / "awd_fk_direction_ground_truth_v1"

BANDS_HZ = ((5.0, 20.0), (20.0, 60.0), (25.0, 60.0))
FAN_M_S = (2500.0, 4500.0)


def active_ratios(
    frequency: np.ndarray,
    slowness: np.ndarray,
    power: np.ndarray,
) -> dict[str, dict[str, float]]:
    """Return expected/opposite signed-power ratios for fixed bands."""
    result = {}
    fan = (
        (np.abs(slowness) >= 1.0 / FAN_M_S[1])
        & (np.abs(slowness) <= 1.0 / FAN_M_S[0])
    )
    for fmin, fmax in BANDS_HZ:
        selected_frequency = (frequency >= fmin) & (frequency <= fmax)
        expected = np.ix_(selected_frequency, slowness > 0.0)
        opposite = np.ix_(selected_frequency, slowness < 0.0)
        expected_fan = np.ix_(selected_frequency, fan & (slowness > 0.0))
        opposite_fan = np.ix_(selected_frequency, fan & (slowness < 0.0))
        result[f"{fmin:g}_{fmax:g}_hz"] = {
            "all_slowness_expected_to_opposite_power_ratio": float(
                power[expected].sum() / power[opposite].sum()
            ),
            "fan_2p5_4p5_expected_to_opposite_power_ratio": float(
                power[expected_fan].sum() / power[opposite_fan].sum()
            ),
        }
    return result


def main() -> None:
    for path in (ACTIVE, AMBIENT, SIGNED, SIGNED_JSON):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    with np.load(ACTIVE, allow_pickle=False) as product:
        frequency = np.asarray(product["frequency"])
        slowness = np.asarray(product["slowness"])
        active_power = np.asarray(product["nano_80_440_m_power"])
        active_window = np.asarray(product["time_window"])
    with np.load(AMBIENT, allow_pickle=False) as product:
        ambient_frequency = np.asarray(product["frequency_hz"])
        ambient_wavenumber = np.asarray(product["wavenumber_cycles_m"])
        ambient_power = np.asarray(product["observed_power"])
        negative_target = np.asarray(product["negative_target_support"]).astype(bool)
        positive_target = np.asarray(product["positive_target_support"]).astype(bool)
        used_files = int(np.asarray(product["used_files"]).size)
    with np.load(SIGNED, allow_pickle=False) as product:
        signed_velocities = np.asarray(product["velocities_m_s"])
        index_3200 = int(np.argmin(np.abs(signed_velocities - 3200.0)))
        negative_score = float(
            np.asarray(product["negative_physical_scores"])[index_3200]
        )
        positive_score = float(
            np.asarray(product["positive_physical_scores"])[index_3200]
        )
    signed_report = json.loads(SIGNED_JSON.read_text())
    signed_files = int(signed_report["weighted_files"])

    ratios = active_ratios(frequency, slowness, active_power)
    ambient_fan_ratio = float(
        ambient_power[negative_target].sum() / ambient_power[positive_target].sum()
    )
    ambient_total_negative = float(
        ambient_power[ambient_wavenumber < 0.0].sum()
    )
    ambient_total_positive = float(
        ambient_power[ambient_wavenumber > 0.0].sum()
    )
    ambient_total_ratio = ambient_total_negative / ambient_total_positive
    ambient_correlation_ratio = abs(negative_score) / abs(positive_score)

    report = {
        "workflow_version": "awd_fk_direction_ground_truth_v1",
        "physical_question": (
            "Does the signed F-K convention select the known downgoing AWD "
            "arrival, and what does that imply for the strong opposite ambient branch?"
        ),
        "coordinate_assumption": (
            "Nano channel coordinate increases from the top of the cemented "
            "fiber into the well, following the Lellouch/main-hole convention."
        ),
        "expected_active_direction": {
            "phase_beam_slowness_sign": "positive",
            "fft_mask_sign": "F*K < 0",
            "physical_label": "toward increasing fiber coordinate (downgoing for Nano)",
        },
        "active_input": str(ACTIVE.relative_to(ROOT)),
        "active_aperture_m": [80.0, 440.0],
        "active_time_window_relative_to_catalog_reference_s": active_window.tolist(),
        "active_expected_to_opposite_power_ratios": ratios,
        "ambient_input": str(AMBIENT.relative_to(ROOT)),
        "ambient_prefilter_files": used_files,
        "ambient_fk_negative_to_positive_total_power_ratio": ambient_total_ratio,
        "ambient_fk_negative_to_positive_fan_power_ratio": ambient_fan_ratio,
        "ambient_signed_correlation_files": signed_files,
        "ambient_negative_to_positive_abs_score_ratio_3200": ambient_correlation_ratio,
        "decision": {
            "production_sign_convention_empirically_validated": bool(
                ratios["25_60_hz"][
                    "fan_2p5_4p5_expected_to_opposite_power_ratio"
                ]
                > 10.0
            ),
            "ambient_input_is_directionally_balanced_in_target_fan": bool(
                0.8 <= ambient_fan_ratio <= 1.25
            ),
            "strong_opposite_ambient_branch_is_explained_by_sign_bug": False,
            "interpretation": (
                "The same sign convention strongly selects the known downgoing "
                "AWD arrival, whereas the ambient field contains comparable signed "
                "power. The opposite ambient branch is therefore a property of the "
                "recorded ambient wavefield and correlation geometry, not evidence "
                "that the signed F-K code simply mirrors one branch."
            ),
        },
        "interpretive_boundary": (
            "This calibration validates coordinate direction. It does not by "
            "itself determine whether the decreasing-coordinate ambient energy is "
            "generated below the aperture, reflected, scattered, or produced by a "
            "different source distribution than in the 2017 experiment."
        ),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.4), constrained_layout=True)

    active_db = 10.0 * np.log10(
        active_power / np.nanmax(active_power) + np.finfo(float).tiny
    )
    active_image = axes[0, 0].pcolormesh(
        slowness * 1e3,
        frequency,
        active_db,
        shading="auto",
        cmap="magma",
        vmin=-25.0,
        vmax=0.0,
    )
    for sign in (-1.0, 1.0):
        axes[0, 0].axvspan(
            sign / FAN_M_S[1] * 1e3,
            sign / FAN_M_S[0] * 1e3,
            color="cyan" if sign > 0 else "white",
            alpha=0.12,
        )
    axes[0, 0].axvline(0.0, color="white", lw=0.8)
    axes[0, 0].set_xlim(-0.8, 0.8)
    axes[0, 0].set_ylim(5.0, 65.0)
    axes[0, 0].set_xlabel("Signed slowness (ms m$^{-1}$)")
    axes[0, 0].set_ylabel("Frequency (Hz)")
    axes[0, 0].set_title("a  AWD ground truth: expected direction is positive slowness")
    fig.colorbar(active_image, ax=axes[0, 0], label="Relative phase-beam power (dB)")

    labels = ["5--20", "20--60", "25--60"]
    all_ratios = [
        ratios[key]["all_slowness_expected_to_opposite_power_ratio"]
        for key in ("5_20_hz", "20_60_hz", "25_60_hz")
    ]
    fan_ratios = [
        ratios[key]["fan_2p5_4p5_expected_to_opposite_power_ratio"]
        for key in ("5_20_hz", "20_60_hz", "25_60_hz")
    ]
    x = np.arange(3)
    width = 0.36
    axes[0, 1].bar(x - width / 2, all_ratios, width, label="all slownesses")
    axes[0, 1].bar(x + width / 2, fan_ratios, width, label="2.5--4.5 km/s")
    axes[0, 1].axhline(1.0, color="0.25", ls="--", lw=1.0)
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xticks(x, labels)
    axes[0, 1].set_xlabel("Frequency band (Hz)")
    axes[0, 1].set_ylabel("Expected/opposite signed power")
    axes[0, 1].set_title("b  Known surface source selects one direction")
    axes[0, 1].legend(frameon=False)
    axes[0, 1].grid(axis="y", alpha=0.2)

    order = np.argsort(ambient_wavenumber)
    ambient_db = 10.0 * np.log10(
        ambient_power[order] / np.nanmax(ambient_power) + np.finfo(float).tiny
    )
    ambient_image = axes[1, 0].pcolormesh(
        ambient_frequency,
        ambient_wavenumber[order] * 1000.0,
        ambient_db,
        shading="auto",
        cmap="magma",
        vmin=-35.0,
        vmax=0.0,
    )
    for velocity, linestyle in ((FAN_M_S[0], "--"), (FAN_M_S[1], ":")):
        k = ambient_frequency / velocity * 1000.0
        axes[1, 0].plot(ambient_frequency, k, "w", ls=linestyle, lw=0.8)
        axes[1, 0].plot(ambient_frequency, -k, "w", ls=linestyle, lw=0.8)
    axes[1, 0].axhline(0.0, color="white", lw=0.8)
    axes[1, 0].set_ylim(-9.0, 9.0)
    axes[1, 0].set_xlabel("Frequency (Hz)")
    axes[1, 0].set_ylabel("Wavenumber (cycles km$^{-1}$)")
    axes[1, 0].set_title("c  Ambient input: comparable signed F-K power")
    fig.colorbar(ambient_image, ax=axes[1, 0], label="Relative power (dB)")

    values = [
        fan_ratios[0],
        fan_ratios[2],
        ambient_fan_ratio,
        ambient_correlation_ratio,
    ]
    labels = [
        "AWD\n5--20 Hz",
        "AWD\n25--60 Hz",
        "Ambient pre-filter\nfan power",
        "Ambient filtered\nscore",
    ]
    colors = ["tab:blue", "tab:blue", "tab:orange", "tab:orange"]
    bars = axes[1, 1].bar(np.arange(4), values, color=colors)
    axes[1, 1].axhline(1.0, color="0.25", ls="--", lw=1.0)
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xticks(np.arange(4), labels)
    axes[1, 1].set_ylabel("Expected/opposite or $F K<0$/$F K>0$ ratio")
    axes[1, 1].set_title("d  Active asymmetry versus ambient balance")
    axes[1, 1].grid(axis="y", alpha=0.2)
    for bar, value in zip(bars, values):
        axes[1, 1].text(
            bar.get_x() + bar.get_width() / 2,
            value * 1.12,
            f"{value:.2f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    fig.suptitle(
        "Signed F-K direction calibration: the AWD validates the sign; the ambient field is bidirectional",
        fontsize=14,
    )
    fig.savefig(OUT.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
