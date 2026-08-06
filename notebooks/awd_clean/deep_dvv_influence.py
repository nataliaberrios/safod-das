"""Single-aperture and single-burst influence on the Deep sensitivity result.

Pre-registration failure criterion: "strong dependence on one aperture or one
burst".  Nothing in the main script measured that, so this diagnostic closes the
gap.  It was written before any recovered value was inspected.

Two jackknives, both on the primary pass of the held-out population:

* leave-one-aperture-out -- refit each trial with each aperture dropped in turn,
  and report the largest resulting shift in the recovered fractional change.
* leave-one-burst-out    -- drop each burst from the population and report the
  largest resulting shift in the per-level median and in the empirical 95%
  zero-injection threshold.

The influence is reported as a ratio to that leg's own empirical null threshold,
because that is the scale at which a shift would actually change a conclusion.
No pass/fail cutoff is invented here: the pre-registration names the criterion
but does not quantify it, so the numbers are reported and the comparison is made
explicit rather than being decided by a threshold chosen after the fact.

Reads the blinded gathers and the sealed truth table. Changes no band,
trajectory, aperture geometry, or quality threshold.
"""

from __future__ import annotations

import csv

import numpy as np

import deep_dvv_injection_recovery as D


POPULATION = "heldout"
BAND = D.PRIMARY_BAND
OUT_CSV = D.HERE / "deep_dvv_influence.csv"
OUT_TXT = D.HERE / "deep_dvv_influence.txt"


def _aperture_delays(target, reference, q, fs, config):
    """Per-aperture delay and correlation, identical to the primary estimator."""
    use = (q >= config["recovery_q_s"][0]) & (q <= config["recovery_q_s"][1])
    max_lag = config["max_lag_s"]
    n = target.shape[0]
    delay = np.full(n, np.nan)
    correlation = np.full(n, np.nan)
    for index in range(n):
        delay[index], correlation[index] = D._normalized_correlation_lag(
            reference[index][use], target[index][use], fs, max_lag
        )
    good = (
        np.isfinite(delay)
        & np.isfinite(correlation)
        & (correlation >= D.MIN_APERTURE_CORRELATION)
        & (np.abs(delay) <= D.BOUNDARY_FRACTION * max_lag)
    )
    return delay, correlation, good


def _fit(travel, delay, correlation, mask):
    if np.sum(mask) < D.MIN_APERTURES:
        return np.nan
    weight = np.maximum(correlation[mask], 0.0) ** 2
    beta, _, _ = D._robust_line(travel[mask], delay[mask], weight)
    return -float(beta[1])


def main() -> None:
    # Composite key, for the same reason the main script needs one: trial ids
    # alone are not unique across populations.
    truth = D._assert_unique_truth(D._read_csv(D.TRUTH))
    config = D.BAND_CONFIG[BAND]
    rows = []
    per_leg = {}

    for leg in D.LEGS:
        path = D._blind_path(leg, BAND, POPULATION)
        if not path.exists():
            print(f"missing {path}, skipping {leg}")
            continue
        with np.load(path) as blind:
            trial_ids = blind["trial_id"].astype(str)
            trials = blind["trial_gather"]
            references = blind["reference_loo"]
            target_index = blind["target_index"]
            q = blind["q_s"]
            coordinate = blind["coordinate_m"]
            fs = float(blind["fs"])
            intercept = float(blind["intercept_s"])
            slowness = float(blind["slowness_s_per_m"])
            epochs = blind["epochs"]

            layout = D._aperture_layout(coordinate, D.APERTURE_M)
            centers = np.asarray([c for _, c in layout])
            travel = intercept + slowness * centers

            cache: dict[int, np.ndarray] = {}
            records = []
            for itrial, trial_id in enumerate(trial_ids):
                meta = truth[(POPULATION, leg, BAND, trial_id)]
                if meta["trial_kind"] != "dvv":
                    continue
                index = int(target_index[itrial])
                if index not in cache:
                    cache[index] = D._beams(references[index], q, coordinate, layout, None)
                target = D._beams(trials[itrial], q, coordinate, layout, None)
                delay, correlation, good = _aperture_delays(
                    target, cache[index], q, fs, config
                )
                full = _fit(travel, delay, correlation, good)
                worst = 0.0
                worst_center = np.nan
                for j in np.flatnonzero(good):
                    dropped = good.copy()
                    dropped[j] = False
                    value = _fit(travel, delay, correlation, dropped)
                    if np.isfinite(value) and np.isfinite(full):
                        if abs(value - full) > worst:
                            worst = abs(value - full)
                            worst_center = float(centers[j])
                records.append({
                    "leg": leg,
                    "trial_id": trial_id,
                    "epoch": int(epochs[index]),
                    "injected_dvv": float(meta["injected_value"]),
                    "estimated_dvv": full,
                    "max_leave_one_aperture_shift": worst,
                    "most_influential_aperture_center_m": worst_center,
                    "n_apertures": int(np.sum(good)),
                })
            per_leg[leg] = records
            print(f"influence {leg:9s} {len(records)} dvv trials, {len(layout)} apertures")

    for records in per_leg.values():
        rows.extend(records)
    if not rows:
        raise RuntimeError("No blinded held-out gathers found; run --stage inject first")
    with OUT_CSV.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Deep guided mode: single-aperture and single-burst influence",
        "=" * 70,
        f"Population: {POPULATION}; band: {BAND.replace('_', '-')} Hz; primary pass.",
        "",
    ]
    for leg, records in per_leg.items():
        injected = np.asarray([r["injected_dvv"] for r in records])
        estimated = np.asarray([r["estimated_dvv"] for r in records])
        epoch = np.asarray([r["epoch"] for r in records])
        aperture_shift = np.asarray([r["max_leave_one_aperture_shift"] for r in records])
        yields = np.asarray([r["n_apertures"] for r in records])

        null = estimated[(injected == 0) & np.isfinite(estimated)]
        if null.size < 5:
            lines.append(f"{leg}: too few finite null trials to characterise influence")
            continue
        centered = null - np.median(null)
        threshold = float(np.quantile(np.abs(centered), 0.95))

        finite = np.isfinite(aperture_shift) & (aperture_shift > 0)
        median_influence = float(np.median(aperture_shift[finite]))
        p95_influence = float(np.quantile(aperture_shift[finite], 0.95))

        # Leave-one-burst-out on the per-level median and on the null threshold.
        worst_level_shift, worst_level, worst_burst = 0.0, np.nan, np.nan
        for level in np.unique(injected):
            selected = (injected == level) & np.isfinite(estimated)
            base = float(np.median(estimated[selected]))
            for burst in np.unique(epoch):
                keep = selected & (epoch != burst)
                if np.sum(keep) < 5:
                    continue
                shift = abs(float(np.median(estimated[keep])) - base)
                if shift > worst_level_shift:
                    worst_level_shift, worst_level, worst_burst = shift, level, burst
        worst_threshold_shift, worst_threshold_burst = 0.0, np.nan
        null_epoch = epoch[(injected == 0) & np.isfinite(estimated)]
        for burst in np.unique(null_epoch):
            keep = null_epoch != burst
            if np.sum(keep) < 5:
                continue
            subset = null[keep]
            value = float(np.quantile(np.abs(subset - np.median(subset)), 0.95))
            if abs(value - threshold) > worst_threshold_shift:
                worst_threshold_shift, worst_threshold_burst = abs(value - threshold), burst

        lines += [
            f"{leg}",
            "-" * 70,
            f"  empirical 95% null threshold        {threshold:.4e}",
            f"  aperture yield  median {np.median(yields):.1f}, "
            f"min {yields.min()}, "
            f"below MIN_APERTURES={D.MIN_APERTURES}: "
            f"{float(np.mean(yields < D.MIN_APERTURES)):.3f} of trials",
            f"  leave-one-aperture-out shift        "
            f"median {median_influence:.3e} ({median_influence / threshold:.2f}x threshold), "
            f"p95 {p95_influence:.3e} ({p95_influence / threshold:.2f}x threshold)",
            f"  leave-one-burst-out, level median   "
            f"worst {worst_level_shift:.3e} ({worst_level_shift / threshold:.2f}x threshold) "
            f"at injected {worst_level:.1e}, epoch {worst_burst:.0f}",
            f"  leave-one-burst-out, null threshold "
            f"worst {worst_threshold_shift:.3e} "
            f"({worst_threshold_shift / threshold:.2f}x threshold), epoch {worst_threshold_burst:.0f}",
            "",
        ]
    lines += [
        "Reading these numbers: a ratio well below 1 means no single aperture or",
        "burst can move the result by as much as the null threshold, so the estimate",
        "is not carried by one measurement. A ratio near or above 1 satisfies the",
        "pre-registered 'strong dependence on one aperture or one burst' failure",
        "criterion. The pre-registration did not fix a numerical cutoff, so none is",
        "asserted here.",
        "",
    ]
    report = "\n".join(lines) + "\n"
    OUT_TXT.write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
