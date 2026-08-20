#!/usr/bin/env python3
"""Aggregate exact-stack ordered products against their matched INPUT-LEVEL nulls.

The decision statistic is the maximum causal envelope moveout score over the
declared scan.  The maximum is retaken in every null realization, which accounts
for the fact that the velocity was selected after inspecting Figure 7c.

WHY THE INPUT-LEVEL NULL IS THE ONE THAT COUNTS.  The receiver-order permutation
that `exact_stack` computes internally (`receiver_order_null_maxima`) permutes
the finished gather, so it tests whether the moveout is specific to the true
receiver order.  It cannot test whether the OPERATOR manufactures moveout,
because the operator has already run.  The nulls summarized here are built
before the operator runs -- `white_noise` replaces the data, `channel_permutation`
scrambles channel identity pre-correlation -- so the operator is inside the null.
That distinction is exactly what the ambient F-K fan failed, and the reason it is
not reported as a Lellouch reproduction.

REWRITTEN 2026-08-19.  The previous version was keyed to a schema that no longer
exists and could not run:
    arrays  fixed_simple, fixed_coordinates_m, fixed_simple_velocity_scores
            -> now r_plus_minus_10_correlation, offsets_m, causal_moveout_scores
    keys    fk_mode, used_files, used_windows, null_realization
            -> fk_mode is gone (option removed, fan QC-rejected); the rest are
               files, windows_30_s_15_s_step, null_realization
    null    circular_time_shift was never valid for exact_stack; it exists only
            in ambient_fk_full_pipeline_null_v2.py
The comparison axis is now the common-mode branch rather than the F-K mode, in
step with ambient_lellouch2019_matched_{ordered,nulls}.sbatch.

Requires the null products those two launchers write.  It reports what is
missing and exits 1 rather than raising, so an incomplete sweep is legible.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
FULL = ROOT / "ambient_transfer" / "lellouch2019_exact_stack"
NULLS = ROOT / "ambient_transfer" / "lellouch2019_exact_stack_nulls"
OUT = ROOT / "ambient_transfer" / "lellouch2019_exact_stack_summary"

# (label, common_mode_removal) -- the axis the launchers now sweep
BRANCHES = (("config0_no_removal", False), ("config3_median_removal", True))
NULL_METHODS = ("white_noise", "channel_permutation")
EXPECTED_REALIZATIONS = 20
V_REF = 3200.0
COLORS = {"white_noise": "tab:green", "channel_permutation": "tab:blue"}


def records(directory: Path) -> list[dict]:
    output = []
    if not directory.is_dir():
        return output
    for path in sorted(directory.glob("lellouch_exact_*.json")):
        report = json.loads(path.read_text())
        report["_json"] = path
        report["_npz"] = path.with_suffix(".npz")
        if report["_npz"].is_file():
            output.append(report)
    return output


def empirical_p(observed: float, null: np.ndarray) -> float:
    """One-sided, with the +1 correction so p can never be reported as 0."""
    return float((1 + np.sum(null >= observed)) / (null.size + 1))


def field(item: dict, *names, default=None):
    """First present key among `names`; tolerates the older report schema."""
    for name in names:
        if name in item:
            return item[name]
    if default is not None:
        return default
    raise KeyError("none of %s in %s" % (names, item.get("_json")))


def main() -> int:
    available = records(FULL) + records(NULLS)
    if not available:
        print("no exact-stack products found under:\n  %s\n  %s" % (FULL, NULLS))
        print("run ambient_lellouch2019_matched_ordered.sbatch and "
              "ambient_lellouch2019_matched_nulls.sbatch first")
        return 1

    problems: list[str] = []
    summary = {
        "workflow_version": "summarize_lellouch2019_exact_stack_v2",
        "decision_statistic": (
            "maximum causal envelope moveout score over the declared velocity "
            "scan; velocity selection repeated independently in every null"
        ),
        "null_class": (
            "input-level surrogates built BEFORE the operator runs, so the "
            "operator is inside the null and cannot inflate the p-value"
        ),
        "branches": {},
    }
    loaded: dict[tuple, object] = {}

    for label, removal in BRANCHES:
        def same_branch(item):
            return (
                bool(field(item, "common_mode_removal", default=False)) == removal
                and field(item, "spectral_mode", default="cross_correlation")
                == "cross_correlation"
            )

        ordered_records = [
            item for item in available
            if same_branch(item)
            and field(item, "null_method", default="ordered") == "ordered"
        ]
        if len(ordered_records) != 1:
            problems.append(
                "%s: expected exactly one ordered product, found %d"
                % (label, len(ordered_records))
            )
            continue

        ordered_report = ordered_records[0]
        ordered = np.load(ordered_report["_npz"])
        scores = ordered["causal_moveout_scores"]
        grid = ordered["velocity_grid_m_s"]
        observed_stat = float(np.max(scores))
        observed_i = int(np.argmax(scores))
        ordered_windows = field(
            ordered_report, "windows_30_s_15_s_step", "used_windows"
        )

        branch_report = {
            "ordered_product": str(ordered_report["_npz"]),
            "ordered_windows": ordered_windows,
            "ordered_selected_velocity_m_s": float(grid[observed_i]),
            "ordered_maximum_score": observed_stat,
            # The pedestal diagnostic. A statistic that rises monotonically with
            # trial velocity is measuring proximity to the zero-lag lobe, not
            # moveout, and its p-value is not interpretable as a detection.
            "pedestal_corr_velocity_score": float(np.corrcoef(grid, scores)[0, 1]),
            "selected_at_scan_edge": observed_i in (0, len(scores) - 1),
            "nulls": {},
        }
        loaded[(label, "ordered")] = ordered

        for method in NULL_METHODS:
            candidates = [
                item for item in available
                if same_branch(item)
                and field(item, "null_method", default="ordered") == method
            ]
            candidates.sort(key=lambda it: field(it, "null_realization", default=0))
            if len(candidates) != EXPECTED_REALIZATIONS:
                problems.append(
                    "%s/%s: expected %d realizations, found %d"
                    % (label, method, EXPECTED_REALIZATIONS, len(candidates))
                )
                if not candidates:
                    continue
            statistics, selected, curves = [], [], []
            for item in candidates:
                item_windows = field(
                    item, "windows_30_s_15_s_step", "used_windows"
                )
                if item_windows != ordered_windows:
                    # A null on fewer windows is a weaker stack, not a null.
                    problems.append(
                        "%s/%s: window mismatch in %s (%s vs ordered %s)"
                        % (label, method, item["_json"].name,
                           item_windows, ordered_windows)
                    )
                    continue
                product = np.load(item["_npz"])
                curve = product["causal_moveout_scores"]
                index = int(np.argmax(curve))
                statistics.append(float(curve[index]))
                selected.append(float(product["velocity_grid_m_s"][index]))
                curves.append(curve)
            if not statistics:
                continue
            statistics = np.asarray(statistics)
            branch_report["nulls"][method] = {
                "realizations": len(statistics),
                "null95_maximum_score": float(np.quantile(statistics, 0.95)),
                "empirical_p": empirical_p(observed_stat, statistics),
                "selected_velocity_m_s_median": float(np.median(selected)),
            }
            loaded[(label, method)] = {
                "statistics": statistics,
                "curves": np.asarray(curves),
            }
        summary["branches"][label] = branch_report

    summary["problems"] = problems
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "lellouch2019_exact_stack_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    complete = [lab for lab, _ in BRANCHES
                if lab in summary["branches"]
                and len(summary["branches"][lab]["nulls"]) == len(NULL_METHODS)]
    if complete:
        n = len(complete)
        fig, axes = plt.subplots(3, n, figsize=(6.3 * n, 12), squeeze=False,
                                 constrained_layout=True)
        for column, label in enumerate(complete):
            ordered = loaded[(label, "ordered")]
            rep = summary["branches"][label]
            lags = ordered["lags_s"]
            distance = ordered["offsets_m"]
            section = ordered["r_plus_minus_10_correlation"]
            limit = float(np.percentile(np.abs(section), 99.0))
            axes[0, column].imshow(
                section,
                extent=[lags[0], lags[-1], distance[-1], distance[0]],
                aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit,
                interpolation="nearest",
            )
            axes[0, column].plot(distance / V_REF, distance, "k--", lw=1.2)
            axes[0, column].set(
                xlim=(0.0, 0.35), xlabel="Correlation lag (s)",
                ylabel="Receiver offset (m)",
                title="%s: ordered R+-10 stack" % label,
            )

            velocity = ordered["velocity_grid_m_s"] / 1000.0
            axes[1, column].plot(velocity, ordered["causal_moveout_scores"],
                                 color="black", lw=2, label="ordered")
            for method in NULL_METHODS:
                if (label, method) not in loaded:
                    continue
                curves = loaded[(label, method)]["curves"]
                lo, hi = np.quantile(curves, [0.05, 0.95], axis=0)
                axes[1, column].fill_between(velocity, lo, hi,
                                             color=COLORS[method], alpha=0.16,
                                             label=method)
            axes[1, column].axvline(V_REF / 1000.0, color="0.35", ls=":")
            axes[1, column].set(
                xlabel="Trial apparent velocity (km/s)",
                ylabel="Envelope moveout score",
                title="Ordered curve vs matched-null 5-95%%\npedestal corr(v,score) = %+.3f"
                      % rep["pedestal_corr_velocity_score"],
            )
            axes[1, column].legend(frameon=False, fontsize=8)

            observed = rep["ordered_maximum_score"]
            top = max([observed] + [float(np.max(loaded[(label, m)]["statistics"]))
                                    for m in NULL_METHODS if (label, m) in loaded])
            bins = np.linspace(0.0, top * 1.05, 40)
            for method in NULL_METHODS:
                if (label, method) not in loaded:
                    continue
                axes[2, column].hist(
                    loaded[(label, method)]["statistics"], bins=bins,
                    color=COLORS[method], alpha=0.55, label="%s (p=%.3f)"
                    % (method, rep["nulls"][method]["empirical_p"]),
                )
            axes[2, column].axvline(observed, color="black", lw=2, label="observed")
            axes[2, column].set(xlabel="Maximum moveout score",
                                ylabel="Realizations",
                                title="Matched-null distributions")
            axes[2, column].legend(frameon=False, fontsize=8)
        stem = OUT / "lellouch2019_exact_stack_summary"
        fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
        print("wrote %s.{png,json}" % stem)
    else:
        print("no branch has a complete null sweep; wrote JSON only")

    for label in summary["branches"]:
        rep = summary["branches"][label]
        print("\n=== %s ===" % label)
        print("  ordered max %.4f at %.0f m/s | pedestal corr(v,score) %+.3f | edge %s"
              % (rep["ordered_maximum_score"], rep["ordered_selected_velocity_m_s"],
                 rep["pedestal_corr_velocity_score"], rep["selected_at_scan_edge"]))
        for method, nrep in rep["nulls"].items():
            print("  %-22s n=%2d  null95 %.4f  p = %.4f"
                  % (method, nrep["realizations"],
                     nrep["null95_maximum_score"], nrep["empirical_p"]))
        if abs(rep["pedestal_corr_velocity_score"]) >= 0.5:
            print("  NOTE: statistic is pedestal-dominated; its p-value is not a")
            print("        detection. See FIG7C_MULTIDAY_RESULT.md.")

    if problems:
        print("\nincomplete or inconsistent (%d):" % len(problems))
        for line in problems:
            print("  " + line)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
