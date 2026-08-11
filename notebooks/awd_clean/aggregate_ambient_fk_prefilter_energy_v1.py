#!/usr/bin/env python3
"""Aggregate and plot distributed pre-filter F-K energy null realizations."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


METHODS = ("channel_permutation", "circular_time_shift")
BRANCHES = ("negative", "positive")


def empirical_upper_p(observed: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= observed)) / (null.size + 1))


def require_close(reference: np.ndarray, value: np.ndarray, label: str) -> None:
    if reference.shape != value.shape or not np.allclose(
        reference, value, rtol=5e-5, atol=1e-8, equal_nan=True
    ):
        raise ValueError(f"distributed products disagree for {label}")


def load_batches(input_dir: Path, expected_nulls: int) -> tuple[list[dict], list[dict]]:
    reports = []
    products = []
    for json_path in sorted(input_dir.glob("fk_prefilter_energy_v1_*.json")):
        report = json.loads(json_path.read_text())
        if report.get("workflow_version") != "ambient_fk_prefilter_energy_v1":
            continue
        if int(report.get("null_realizations", -1)) != 1:
            raise ValueError(f"expected one realization per batch: {json_path}")
        npz_path = json_path.with_suffix(".npz")
        if not npz_path.is_file():
            raise FileNotFoundError(npz_path)
        product = np.load(npz_path, allow_pickle=False)
        reports.append(report)
        products.append({key: np.asarray(product[key]) for key in product.files})
        product.close()
    ids = [int(report["null_start"]) for report in reports]
    if len(reports) != expected_nulls or sorted(ids) != list(range(expected_nulls)):
        raise ValueError(
            f"expected null IDs 0..{expected_nulls - 1}; found {sorted(ids)}"
        )
    order = np.argsort(ids)
    return [reports[i] for i in order], [products[i] for i in order]


def validate_common(reports: list[dict], products: list[dict]) -> None:
    first_report = reports[0]
    exact_keys = (
        "date",
        "start",
        "requested_files",
        "used_files",
        "random_seed",
        "passband_hz",
        "target_velocity_m_s",
        "reference_velocity_m_s",
        "statistic",
        "preprocessing",
    )
    for report in reports[1:]:
        for key in exact_keys:
            if report.get(key) != first_report.get(key):
                raise ValueError(f"distributed metadata disagree for {key}")
    first = products[0]
    array_keys = (
        "frequency_hz",
        "wavenumber_cycles_m",
        "negative_target_support",
        "positive_target_support",
        "negative_reference_support",
        "positive_reference_support",
    )
    for product in products[1:]:
        for key in array_keys:
            require_close(first[key], product[key], key)
        if not np.array_equal(first["used_files"], product["used_files"]):
            raise ValueError("distributed products used different input files")


def aggregate(reports: list[dict], products: list[dict]) -> tuple[dict, dict[str, np.ndarray]]:
    first_report = reports[0]
    first = products[0]
    observed_per_file = {
        branch: np.mean(
            np.stack([product[f"observed_{branch}_per_file"] for product in products]),
            axis=0,
        )
        for branch in BRANCHES
    }
    observed_means_by_product = {
        branch: np.asarray(
            [np.mean(product[f"observed_{branch}_per_file"]) for product in products]
        )
        for branch in BRANCHES
    }
    maximum_observed_mean_spread = max(
        float(np.ptp(values)) for values in observed_means_by_product.values()
    )
    if maximum_observed_mean_spread > 1e-3:
        raise ValueError(
            "distributed observed means differ by more than 1e-3: "
            f"{maximum_observed_mean_spread}"
        )
    observed = {
        branch: float(np.mean(observed_per_file[branch])) for branch in BRANCHES
    }
    observed_familywise = max(observed.values())
    observed_power_stack = np.stack([product["observed_power"] for product in products])
    observed_power = np.mean(observed_power_stack, axis=0)
    maximum_observed_power_relative_spread = float(
        np.max(np.ptp(observed_power_stack, axis=0))
        / (np.max(observed_power) + np.finfo(float).tiny)
    )
    arrays: dict[str, np.ndarray] = {
        "null_realization_ids": np.arange(len(reports), dtype=int),
        "frequency_hz": first["frequency_hz"],
        "wavenumber_cycles_m": first["wavenumber_cycles_m"],
        "observed_power": observed_power,
        "negative_target_support": first["negative_target_support"],
        "positive_target_support": first["positive_target_support"],
        "negative_reference_support": first["negative_reference_support"],
        "positive_reference_support": first["positive_reference_support"],
        "observed_negative_per_file": observed_per_file["negative"],
        "observed_positive_per_file": observed_per_file["positive"],
        "used_files": first["used_files"],
    }
    report: dict[str, object] = {
        "workflow_version": "ambient_fk_prefilter_energy_v1_aggregate",
        "date": first_report["date"],
        "start": first_report["start"],
        "requested_files": first_report["requested_files"],
        "used_files": first_report["used_files"],
        "null_realizations": len(reports),
        "null_realization_ids": list(range(len(reports))),
        "random_seed": first_report["random_seed"],
        "passband_hz": first_report["passband_hz"],
        "target_velocity_m_s": first_report["target_velocity_m_s"],
        "reference_velocity_m_s": first_report["reference_velocity_m_s"],
        "statistic": first_report["statistic"],
        "preprocessing": first_report["preprocessing"],
        "observed": {
            branch: {
                "mean_log_enrichment": observed[branch],
                "geometric_enrichment_ratio": float(np.exp(observed[branch])),
            }
            for branch in BRANCHES
        },
        "observed_familywise": observed_familywise,
        "distributed_reproducibility": {
            "maximum_observed_mean_log_enrichment_spread": maximum_observed_mean_spread,
            "maximum_observed_power_relative_spread": maximum_observed_power_relative_spread,
            "input_file_lists_identical": True,
            "handling": "average cross-node observed products after bounded-spread checks",
        },
        "null_results": {},
        "decision_rule": (
            "The raw input supports target-slowness enrichment only if both separately "
            "reported familywise null probabilities are <= 0.05; passing does not "
            "establish Green's-function convergence or formation velocity."
        ),
        "limitations": first_report["limitations"],
    }
    for method in METHODS:
        method_report: dict[str, object] = {"branches": {}}
        familywise = np.full(len(products), -np.inf, dtype=float)
        for branch in BRANCHES:
            values = np.asarray([
                product[f"null_{method}_{branch}_mean"][0]
                for product in products
            ])
            arrays[f"null_{method}_{branch}"] = values
            familywise = np.maximum(familywise, values)
            method_report["branches"][branch] = {
                "observed_mean_log_enrichment": observed[branch],
                "observed_geometric_enrichment_ratio": float(np.exp(observed[branch])),
                "null_median": float(np.median(values)),
                "null95": float(np.quantile(values, 0.95)),
                "p_upper": empirical_upper_p(observed[branch], values),
            }
        arrays[f"null_{method}_familywise"] = familywise
        method_report["familywise"] = {
            "observed_maximum_across_branches": observed_familywise,
            "null_median": float(np.median(familywise)),
            "null95": float(np.quantile(familywise, 0.95)),
            "p_upper": empirical_upper_p(observed_familywise, familywise),
        }
        report["null_results"][method] = method_report
    report["passes_both_familywise_nulls"] = all(
        float(report["null_results"][method]["familywise"]["p_upper"]) <= 0.05
        for method in METHODS
    )
    return report, arrays


def plot(report: dict, arrays: dict[str, np.ndarray], output_stem: Path) -> None:
    frequency = arrays["frequency_hz"]
    wavenumber = arrays["wavenumber_cycles_m"] * 1000.0
    power = 10.0 * np.log10(
        arrays["observed_power"] / np.nanmax(arrays["observed_power"])
        + np.finfo(float).tiny
    )
    spatial_order = np.argsort(wavenumber)
    wavenumber = wavenumber[spatial_order]
    power = power[spatial_order]
    figure, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), constrained_layout=True)
    image = axes[0, 0].pcolormesh(
        frequency,
        wavenumber,
        power,
        shading="auto",
        cmap="magma",
        vmin=-40.0,
        vmax=0.0,
    )
    for velocity, style in ((2500.0, "--"), (4500.0, ":")):
        axes[0, 0].plot(frequency, frequency / velocity * 1000.0, "w" + style, lw=1.4)
        axes[0, 0].plot(frequency, -frequency / velocity * 1000.0, "w" + style, lw=1.4)
    axes[0, 0].set_ylim(-9.0, 9.0)
    axes[0, 0].set_xlabel("Frequency (Hz)")
    axes[0, 0].set_ylabel("Wavenumber (cycles km$^{-1}$)")
    axes[0, 0].set_title("Mean pre-filter F-K power")
    figure.colorbar(image, ax=axes[0, 0], label="Relative power (dB)")

    observed = {
        branch: float(report["observed"][branch]["mean_log_enrichment"])
        for branch in BRANCHES
    }
    for ax, method in zip((axes[0, 1], axes[1, 0]), METHODS):
        colors = {"negative": "tab:blue", "positive": "tab:orange"}
        for branch in BRANCHES:
            values = arrays[f"null_{method}_{branch}"]
            ax.hist(values, bins="auto", alpha=0.45, color=colors[branch], label=f"{branch} null")
            ax.axvline(observed[branch], color=colors[branch], lw=2.0)
        family = report["null_results"][method]["familywise"]
        ax.axvline(family["null95"], color="0.2", ls=":", lw=1.8, label="familywise null 95%")
        label = "Channel permutation" if method == "channel_permutation" else "Independent time shift"
        ax.set_title(f"{label}: familywise p={family['p_upper']:.3f}")
        ax.set_xlabel("Mean log target/reference enrichment")
        ax.set_ylabel("Null count")
        ax.legend(fontsize=8)

    ax = axes[1, 1]
    x = np.arange(2)
    values = [np.exp(observed[branch]) for branch in BRANCHES]
    ax.bar(x, values, color=("tab:blue", "tab:orange"))
    ax.axhline(1.0, color="0.2", lw=1.2)
    ax.set_xticks(x, ("F×K < 0", "F×K > 0"))
    ax.set_ylabel("Geometric target/reference power ratio")
    ax.set_title(
        "Raw target-band enrichment\n"
        + ("passes both stated nulls" if report["passes_both_familywise_nulls"] else "does not pass both stated nulls")
    )
    figure.suptitle(
        f"Pre-filter F-K energy test: {report['used_files']} one-minute files",
        fontsize=14,
    )
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_stem.with_suffix(".png"), dpi=320)
    figure.savefig(output_stem.with_suffix(".pdf"))
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--expected-nulls", type=int, default=20)
    parser.add_argument("--output-stem", type=Path, required=True)
    args = parser.parse_args()
    reports, products = load_batches(args.input_dir, args.expected_nulls)
    validate_common(reports, products)
    report, arrays = aggregate(reports, products)
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    args.output_stem.with_suffix(".json").write_text(json.dumps(report, indent=2))
    np.savez_compressed(args.output_stem.with_suffix(".npz"), **arrays)
    plot(report, arrays, args.output_stem.with_name(args.output_stem.name + "_figure"))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
