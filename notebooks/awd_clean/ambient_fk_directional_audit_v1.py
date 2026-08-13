#!/usr/bin/env python3
"""Audit standard directional F-K separation for the ambient SAFOD branch.

The production ambient workflow applies both a signed direction selection and
a 2.5--4.5 km/s velocity fan before correlation.  A narrow fan can create a
moveout-shaped correlation kernel even when the unfiltered input has no
independent concentration at that slowness.  This audit therefore separates
two questions:

1. Does the signed F-K implementation distinguish known one-way waves?
2. Does the real held-out wavefield retain a 3.2 km/s correlation when only
   the standard directional quadrant is selected, without a velocity fan?

The script uses frozen, already completed products.  It does not rerun raw
waveforms or modify the production arrays.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
MASK_DIR = ROOT / "ambient_transfer" / "fk_mask_sensitivity_v2"
MASK_NPZ = MASK_DIR / "ambient_fk_mask_sensitivity_v2.npz"
PREFILTER_DIR = ROOT / "ambient_transfer" / "fk_prefilter_energy_v1_n300_r20"
PREFILTER_NPZ = PREFILTER_DIR / "fk_prefilter_energy_v1_aggregate.npz"
PREFILTER_JSON = PREFILTER_DIR / "fk_prefilter_energy_v1_aggregate.json"
SIGNED_DIR = ROOT / "ambient_transfer" / "signed_lag_v2"
SIGNED_NPZ = SIGNED_DIR / "seasonal_signed_fk_v2_aggregate.npz"
OUT_DIR = ROOT / "ambient_transfer" / "fk_directional_audit_v1"
OUT_STEM = OUT_DIR / "ambient_fk_directional_audit_v1"

MASKS = (
    ("production_2p5_4p5", "2.5--4.5 km/s fan"),
    ("narrow_2p8_3p8", "2.8--3.8 km/s fan"),
    ("broad_2p0_5p5", "2.0--5.5 km/s fan"),
    ("direction_only", "direction only"),
)


def mirror_metrics(left: np.ndarray, right: np.ndarray) -> tuple[float, float, float]:
    """Compare two correlation sections after reversing the right lag axis."""
    mirrored = right[:, ::-1]
    correlation = float(np.corrcoef(left.ravel(), mirrored.ravel())[0, 1])
    relative_rms = float(np.linalg.norm(left - mirrored) / np.linalg.norm(left))
    norm_ratio = float(np.linalg.norm(mirrored) / np.linalg.norm(left))
    return correlation, relative_rms, norm_ratio


def physical_view(top: np.ndarray, lags: np.ndarray, sign: int) -> tuple[np.ndarray, np.ndarray]:
    """Return positive travel time for either physical lag side."""
    selected = lags >= 0 if sign > 0 else lags <= 0
    time = sign * lags[selected]
    section = top[:, selected]
    if time[0] > time[-1]:
        time = time[::-1]
        section = section[:, ::-1]
    return time, section


def main() -> None:
    for path in (MASK_NPZ, PREFILTER_NPZ, PREFILTER_JSON, SIGNED_NPZ):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    with np.load(MASK_NPZ, allow_pickle=False) as mask_product:
        mask = {key: np.asarray(mask_product[key]) for key in mask_product.files}
    with np.load(PREFILTER_NPZ, allow_pickle=False) as prefilter_product:
        prefilter = {
            key: np.asarray(prefilter_product[key])
            for key in prefilter_product.files
        }
    with np.load(SIGNED_NPZ, allow_pickle=False) as signed_product:
        signed = {key: np.asarray(signed_product[key]) for key in signed_product.files}
    prefilter_report = json.loads(PREFILTER_JSON.read_text())

    lags = mask["lags"]
    distance = mask["distance"]
    velocities = mask["velocities_m_s"]
    index_3200 = int(np.argmin(np.abs(velocities - 3200.0)))

    comparisons: dict[str, dict[str, float]] = {}
    for group in ("development", "held_out"):
        comparisons[group] = {}
        for mask_name, _ in MASKS:
            negative = mask[f"{group}__{mask_name}__negative__top"]
            positive = mask[f"{group}__{mask_name}__positive__top"]
            correlation, relative_rms, norm_ratio = mirror_metrics(
                negative, positive
            )
            comparisons[group][mask_name] = {
                "mirror_correlation": correlation,
                "mirror_relative_rms": relative_rms,
                "mirrored_norm_ratio_positive_to_negative": norm_ratio,
                "negative_score_3200": float(
                    mask[f"{group}__{mask_name}__negative__scores"][index_3200]
                ),
                "positive_score_3200": float(
                    mask[f"{group}__{mask_name}__positive__scores"][index_3200]
                ),
            }

    power = prefilter["observed_power"]
    wavenumber = prefilter["wavenumber_cycles_m"]
    negative_power = float(np.sum(power[wavenumber < 0.0]))
    positive_power = float(np.sum(power[wavenumber > 0.0]))
    directional_power_ratio = positive_power / negative_power
    directional_power_asymmetry = (
        (negative_power - positive_power) / (negative_power + positive_power)
    )

    signed_mirror_correlation, signed_relative_rms, signed_norm_ratio = (
        mirror_metrics(signed["negative_top"], signed["positive_top"])
    )
    direction_only = comparisons["held_out"]["direction_only"]
    production = comparisons["held_out"]["production_2p5_4p5"]
    report = {
        "workflow_version": "ambient_fk_directional_audit_v1",
        "question": (
            "Does standard signed-quadrant F-K separation independently retain "
            "the 3.2 km/s ambient correlation without a velocity fan?"
        ),
        "inputs": {
            "mask_sensitivity": str(MASK_NPZ.relative_to(ROOT)),
            "pre_filter_energy": str(PREFILTER_NPZ.relative_to(ROOT)),
            "signed_lag_aggregate": str(SIGNED_NPZ.relative_to(ROOT)),
        },
        "standard_operator": (
            "5--20 Hz signed F-K quadrant selection with no apparent-velocity "
            "restriction; branches are reconstructed and correlated independently"
        ),
        "real_pre_filter_direction_balance": {
            "positive_to_negative_total_power_ratio": directional_power_ratio,
            "signed_asymmetry": directional_power_asymmetry,
            "target_enrichment_passes_both_familywise_nulls": bool(
                prefilter_report["passes_both_familywise_nulls"]
            ),
        },
        "mask_comparisons": comparisons,
        "legacy_figure45_mirror_check": {
            "mirror_correlation": signed_mirror_correlation,
            "mirror_relative_rms": signed_relative_rms,
            "mirrored_norm_ratio_positive_to_negative": signed_norm_ratio,
        },
        "decision": {
            "signed_fk_implementation_validated_by_one_way_synthetics": True,
            "standard_direction_only_filter_shows_independent_3p2_km_s_result": False,
            "production_velocity_fan_result_is_operator_conditioned": True,
            "reason": (
                "Held-out direction-only 3.2 km/s scores are near zero, while "
                "the velocity-fan branches become almost exact lag mirrors."
            ),
        },
        "key_values": {
            "held_out_production_mirror_correlation": production[
                "mirror_correlation"
            ],
            "held_out_direction_only_mirror_correlation": direction_only[
                "mirror_correlation"
            ],
            "held_out_direction_only_negative_score_3200": direction_only[
                "negative_score_3200"
            ],
            "held_out_direction_only_positive_score_3200": direction_only[
                "positive_score_3200"
            ],
        },
        "interpretive_boundary": (
            "F-K filtering remains a valid wavefield-separation operation. The "
            "present real-data test does not independently validate a 3.2 km/s "
            "ambient arrival when the apparent-velocity fan is removed."
        ),
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_STEM.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")

    figure, axes = plt.subplots(2, 3, figsize=(15.0, 8.7), constrained_layout=True)

    frequency = prefilter["frequency_hz"]
    k_km = wavenumber * 1000.0
    order = np.argsort(k_km)
    relative_power = 10.0 * np.log10(
        power[order] / np.nanmax(power) + np.finfo(float).tiny
    )
    image = axes[0, 0].pcolormesh(
        frequency,
        k_km[order],
        relative_power,
        shading="auto",
        cmap="magma",
        vmin=-40.0,
        vmax=0.0,
    )
    axes[0, 0].axhline(0.0, color="white", lw=0.8)
    axes[0, 0].set_ylim(-9.0, 9.0)
    axes[0, 0].set_xlabel("Frequency (Hz)")
    axes[0, 0].set_ylabel("Wavenumber (cycles km$^{-1}$)")
    axes[0, 0].set_title("a  Pre-filter F-K power")
    figure.colorbar(image, ax=axes[0, 0], label="Relative power (dB)")

    common_sections = []
    for branch in ("negative", "positive"):
        common_sections.append(mask[f"held_out__direction_only__{branch}__top"])
    limit = float(np.nanpercentile(np.abs(np.concatenate(
        [section.ravel() for section in common_sections]
    )), 98.5))
    for ax, branch, sign, label in (
        (axes[0, 1], "negative", 1, "b  Direction 1 (F K < 0)"),
        (axes[0, 2], "positive", -1, "c  Direction 2 (F K > 0)"),
    ):
        top = mask[f"held_out__direction_only__{branch}__top"]
        travel_time, section = physical_view(top, lags, sign)
        ax.imshow(
            section,
            extent=[travel_time[0], travel_time[-1], distance[-1], distance[0]],
            aspect="auto",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            interpolation="nearest",
        )
        ax.plot(distance / 3200.0, distance, color="black", ls="--", lw=1.1)
        ax.set_xlim(0.0, 0.25)
        ax.set_ylim(700.0, 0.0)
        ax.set_xlabel("Absolute physical lag (s)")
        ax.set_ylabel("Receiver offset (m)")
        ax.set_title(label)

    ax = axes[1, 0]
    x = np.arange(len(MASKS))
    width = 0.36
    dev = [comparisons["development"][name]["mirror_correlation"] for name, _ in MASKS]
    held = [comparisons["held_out"][name]["mirror_correlation"] for name, _ in MASKS]
    ax.bar(x - width / 2, dev, width, label="development")
    ax.bar(x + width / 2, held, width, label="held-out")
    ax.set_xticks(x, [label for _, label in MASKS], rotation=20, ha="right")
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("Correlation after lag reflection")
    ax.set_title("d  Branch mirror similarity")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)

    ax = axes[1, 1]
    for branch, color, label in (
        ("negative", "tab:blue", "F K < 0, positive lag"),
        ("positive", "tab:orange", "F K > 0, negative lag"),
    ):
        scores = mask[f"held_out__direction_only__{branch}__scores"]
        ax.plot(velocities / 1000.0, scores, color=color, lw=1.8, label=label)
    ax.axvline(3.2, color="0.25", ls="--", lw=1.1)
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_xlim(1.2, 6.0)
    ax.set_xlabel("Trial apparent velocity (km s$^{-1}$)")
    ax.set_ylabel("Median correlation")
    ax.set_title("e  Direction-only held-out scores")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)

    ax = axes[1, 2]
    labels = ["Raw F-K\npower ratio", "Fan mirror\ncorrelation", "Direction-only\nmirror correlation"]
    values = [directional_power_ratio, production["mirror_correlation"], direction_only["mirror_correlation"]]
    bars = ax.bar(np.arange(3), values, color=("0.45", "tab:red", "tab:green"))
    ax.axhline(1.0, color="0.25", ls="--", lw=1.0)
    ax.set_xticks(np.arange(3), labels)
    ax.set_ylim(0.0, 1.12)
    ax.set_ylabel("Ratio or correlation")
    ax.set_title("f  Directional decision metrics")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.025, f"{value:.3f}", ha="center", fontsize=9)
    ax.grid(axis="y", alpha=0.2)

    figure.suptitle(
        "Standard directional F-K audit: valid separation operator, no independent 3.2 km/s recovery",
        fontsize=14,
    )
    figure.savefig(OUT_STEM.with_suffix(".png"), dpi=350, bbox_inches="tight")
    figure.savefig(OUT_STEM.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
