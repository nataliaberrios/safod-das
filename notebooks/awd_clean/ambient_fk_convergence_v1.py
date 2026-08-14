#!/usr/bin/env python3
"""Quantify ambient-stack convergence before and after the production F-K fan.

The analysis enumerates every combination of the eight completed seasonal
days.  It compares each partial stack with the eight-day reference and reports
the 3.2 km/s moveout score.  The filtered result is explicitly conditional on
the 2.5--4.5 km/s fan; convergence of that product is not treated as evidence
for Green's-function convergence when pre-filter surrogate nulls fail.
"""
from __future__ import annotations

from itertools import combinations
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
TRANSFER = ROOT / "ambient_transfer"
SIGNED = TRANSFER / "signed_lag_v2"
OUT = TRANSFER / "fk_qc_notebook_v2"
DATES = (
    "2024-05-11", "2024-06-17", "2024-06-26", "2024-10-28",
    "2024-11-30", "2024-12-20", "2025-02-24", "2025-03-04",
)


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a.ravel(), b.ravel())[0, 1])


def score(top: np.ndarray, lags: np.ndarray, distance: np.ndarray,
          velocity: float, lag_sign: float = 1.0) -> float:
    values = [
        row[np.argmin(np.abs(lags - lag_sign * offset / velocity))]
        for row, offset in zip(top, distance)
    ]
    return float(np.nanmedian(values))


def load_products():
    unfiltered = []
    filtered = []
    weights = []
    for date in DATES:
        un = np.load(TRANSFER / f"transfer_seasonal_{date}.npz")
        fk = np.load(SIGNED / f"seasonal_signed_fk_v2_{date}.npz")
        meta = json.loads(
            (SIGNED / f"seasonal_signed_fk_v2_{date}.json").read_text()
        )
        unfiltered.append(un["top_stack"])
        filtered.append(fk["negative_top"])
        weights.append(meta["used_files"])
    return (
        np.stack(unfiltered), np.stack(filtered), np.asarray(weights, float),
        un["lags"], un["distances"], fk["lags"], fk["distance"],
    )


def combination_metrics(tops, weights, lags, distance, final):
    rows = []
    for n_days in range(1, len(tops) + 1):
        for chosen in combinations(range(len(tops)), n_days):
            index = np.asarray(chosen)
            partial = np.average(tops[index], axis=0, weights=weights[index])
            rows.append((
                n_days,
                correlation(partial, final),
                score(partial, lags, distance, 3200.0),
            ))
    return np.asarray(rows, float)


def quantiles(rows, column):
    return np.asarray([
        np.quantile(rows[rows[:, 0] == n, column], [0.1, 0.5, 0.9])
        for n in range(1, 9)
    ])


def main():
    un, fk, weights, un_lags, un_dist, fk_lags, fk_dist = load_products()
    un_final = np.average(un, axis=0, weights=weights)
    fk_final = np.average(fk, axis=0, weights=weights)
    un_rows = combination_metrics(un, weights, un_lags, un_dist, un_final)
    fk_rows = combination_metrics(fk, weights, fk_lags, fk_dist, fk_final)

    five_un = np.load(TRANSFER / "interim_2024-12-20_n300_unfiltered.npz")
    five_fk = np.load(TRANSFER / "interim_2024-12-20_n300_fk_negative.npz")
    five_hour = {
        "unfiltered_section_correlation": correlation(five_un["top"], un_final),
        "filtered_section_correlation": correlation(five_fk["top"], fk_final),
        "unfiltered_score_3200": score(
            five_un["top"], five_un["lags"], five_un["dist"], 3200.0
        ),
        "filtered_score_3200": score(
            five_fk["top"], five_fk["lags"], five_fk["dist"], 3200.0
        ),
    }

    day_axis = np.arange(1, 9)
    un_corr = quantiles(un_rows, 1)
    fk_corr = quantiles(fk_rows, 1)
    un_score = quantiles(un_rows, 2)
    fk_score = quantiles(fk_rows, 2)

    # A purely descriptive convergence rule: at least 90% of day combinations
    # correlate >=0.90 with the final section.  This rule is evaluated for both
    # products but is not allowed to override the independent surrogate tests.
    def first_stable(q):
        passing = q[:, 0] >= 0.90
        for i in range(len(passing)):
            if np.all(passing[i:]):
                return int(i + 1)
        return None

    report = {
        "workflow_version": "ambient_fk_convergence_v1",
        "dates": list(DATES),
        "five_hour": five_hour,
        "descriptive_rule": (
            "first n for which the 10th percentile section correlation with "
            "the final eight-day stack is >=0.90 and remains so"
        ),
        "unfiltered_first_stable_days": first_stable(un_corr),
        "production_fan_first_stable_days": first_stable(fk_corr),
        "decision": (
            "The production-fan section is descriptively stable much sooner "
            "than the unfiltered section, but this is not physical signal "
            "convergence because the same fan fails pre-filter surrogate nulls. "
            "The unfiltered 3.2 km/s observable is not independently recovered "
            "by the eight-day aggregate."
        ),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUT / "ambient_fk_convergence_v1.npz",
        unfiltered_rows=un_rows,
        filtered_rows=fk_rows,
        unfiltered_corr_quantiles=un_corr,
        filtered_corr_quantiles=fk_corr,
        unfiltered_score_quantiles=un_score,
        filtered_score_quantiles=fk_score,
    )
    (OUT / "ambient_fk_convergence_v1.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5), constrained_layout=True)
    for q, color, label in (
        (un_corr, "0.25", "unfiltered"),
        (fk_corr, "tab:blue", "2.5–4.5 km/s signed fan"),
    ):
        axes[0].plot(day_axis, q[:, 1], "o-", color=color, label=label)
        axes[0].fill_between(day_axis, q[:, 0], q[:, 2], color=color, alpha=0.18)
    axes[0].axhline(0.90, color="tab:red", ls="--", lw=1.0)
    axes[0].scatter([5 / 24], [five_hour["unfiltered_section_correlation"]],
                    color="0.25", marker="s")
    axes[0].scatter([5 / 24], [five_hour["filtered_section_correlation"]],
                    color="tab:blue", marker="s")
    axes[0].set(xlabel="Stacked recording duration (days)",
                ylabel="Correlation with eight-day section",
                title="Section stability across all day combinations")
    axes[0].set_xlim(0, 8.25)
    axes[0].legend(frameon=False)

    for q, color, label in (
        (un_score, "0.25", "unfiltered"),
        (fk_score, "tab:blue", "2.5–4.5 km/s signed fan"),
    ):
        axes[1].plot(day_axis, q[:, 1], "o-", color=color, label=label)
        axes[1].fill_between(day_axis, q[:, 0], q[:, 2], color=color, alpha=0.18)
    axes[1].scatter([5 / 24], [five_hour["unfiltered_score_3200"]],
                    color="0.25", marker="s")
    axes[1].scatter([5 / 24], [five_hour["filtered_score_3200"]],
                    color="tab:blue", marker="s")
    axes[1].axhline(0.0, color="0.6", lw=0.8)
    axes[1].set(xlabel="Stacked recording duration (days)",
                ylabel="Median correlation on 3.2 km/s trajectory",
                title="Moveout score across all day combinations")
    axes[1].set_xlim(0, 8.25)
    axes[1].legend(frameon=False)
    fig.suptitle(
        "Ambient-stack convergence: selected-fan stability is conditional",
        fontsize=13,
    )
    fig.savefig(OUT / "ambient_fk_convergence_v1.png", dpi=300)
    fig.savefig(OUT / "ambient_fk_convergence_v1.pdf")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
