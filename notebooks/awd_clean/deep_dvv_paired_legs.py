"""Paired outbound/return combination for the Deep guided-mode observable.

A reader will ask why two legs sampling the same source were analysed
separately instead of combined. This answers that with a measurement.

Outbound is the benchmark, not the thing being replaced. The combined estimator
has to beat it. Per the frozen addendum in DEEP_DVV_PREREGISTRATION.md, three
combinations are compared, all at the level of the per-burst estimate and never
at the level of raw traces:

    equal       (eps_out + eps_ret) / 2
    invvar      inverse-variance weighted, legs treated as independent
    joint       (1' C^-1 eps) / (1' C^-1 1), covariance-aware

sigma and C are estimated leave-one-burst-out from the zero-injection trials, so
the weights applied to a burst never come from the trial being evaluated. They
derive only from the null and are therefore identical at every injection level.

The 23 outbound and 23 return estimates are 23 *paired* observations of the same
23 source bursts. They are never concatenated into 46 independent ones.

Every threshold, bias and reliable-detection number is produced by the same
``_analyse_group`` used for the single-leg results, so the comparison rows in
the output table are computed identically.
"""

from __future__ import annotations

import argparse
import csv

import numpy as np

import deep_dvv_injection_recovery as D


BAND = D.PRIMARY_BAND
PASS = "primary"

# Set by main().  "heldout" is the primary, strictly held-out population; the
# "allbursts" population is trajectory-contaminated for a sensitivity claim but
# is legitimate for a tidal amplitude fit, where the extra bursts buy precision
# and trajectory selection is time-independent.
POPULATION = "heldout"
OUT_CSV = OUT_TABLE = OUT_TXT = None


def _set_population(population: str) -> None:
    global POPULATION, OUT_CSV, OUT_TABLE, OUT_TXT
    POPULATION = population
    suffix = "" if population == "heldout" else f"_{population}"
    OUT_CSV = D.HERE / f"deep_dvv_paired_legs{suffix}.csv"
    OUT_TABLE = D.HERE / f"deep_dvv_paired_comparison{suffix}.csv"
    OUT_TXT = D.HERE / f"deep_dvv_paired_legs{suffix}.txt"


def _load_pairs() -> tuple[dict, dict]:
    """Return {(epoch, level): {leg: eps}} for dvv trials and for shift controls."""
    truth = D._assert_unique_truth(D._read_csv(D.TRUTH))
    dvv: dict = {}
    controls: dict = {}
    for row in D._read_csv(D.RECOVERY):
        if (
            row["population"] != POPULATION
            or row["band"] != BAND
            or row["pass"] != PASS
        ):
            continue
        meta = truth[D._truth_key(row)]
        key = (int(meta["target_epoch"]), float(meta["injected_value"]))
        target = dvv if meta["trial_kind"] == "dvv" else controls
        if meta["trial_kind"] not in ("dvv",):
            key = (int(meta["target_epoch"]), meta["trial_kind"])
        target.setdefault(key, {})[row["leg"]] = float(row["estimated_dvv"])
    return dvv, controls


def _covariance(pairs: np.ndarray) -> np.ndarray:
    """2x2 covariance of (outbound, return) null errors."""
    centered = pairs - np.median(pairs, axis=0, keepdims=True)
    return np.cov(centered, rowvar=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--population", choices=("heldout", "allbursts"),
                        default="heldout")
    _set_population(parser.parse_args().population)
    print(f"population: {POPULATION}")
    dvv, controls = _load_pairs()

    epochs = sorted({epoch for epoch, _ in dvv})
    levels = sorted({level for _, level in dvv})
    print(f"{len(epochs)} paired bursts, {len(levels)} injection levels")

    # Null pairs, one per burst, used only for the leave-one-burst-out weights.
    null_pairs = {}
    for epoch in epochs:
        entry = dvv.get((epoch, 0.0), {})
        if "outbound" in entry and "return" in entry:
            null_pairs[epoch] = (entry["outbound"], entry["return"])
    if len(null_pairs) < 10:
        raise RuntimeError("Too few paired zero-injection trials for weighting")

    full = np.asarray([null_pairs[e] for e in epochs if e in null_pairs])
    full_cov = _covariance(full)
    correlation = float(
        full_cov[0, 1] / np.sqrt(full_cov[0, 0] * full_cov[1, 1])
    )

    weights = {}
    for epoch in epochs:
        others = np.asarray([v for e, v in null_pairs.items() if e != epoch])
        cov = _covariance(others)
        variance = np.diag(cov)
        inv_var = 1.0 / variance
        w_invvar = inv_var / inv_var.sum()
        try:
            precision = np.linalg.inv(cov)
            ones = np.ones(2)
            w_joint = (precision @ ones) / float(ones @ precision @ ones)
        except np.linalg.LinAlgError:
            w_joint = w_invvar
        weights[epoch] = {
            "equal": np.array([0.5, 0.5]),
            "invvar": w_invvar,
            "joint": w_joint,
            "sigma_out": float(np.sqrt(variance[0])),
            "sigma_ret": float(np.sqrt(variance[1])),
            "rho": float(cov[0, 1] / np.sqrt(variance[0] * variance[1])),
        }

    estimators = ("outbound", "return", "equal", "invvar", "joint")
    injected_by: dict[str, list] = {name: [] for name in estimators}
    estimated_by: dict[str, list] = {name: [] for name in estimators}
    rows = []
    dropped = 0
    for epoch in epochs:
        w = weights[epoch]
        for level in levels:
            entry = dvv.get((epoch, level), {})
            if "outbound" not in entry or "return" not in entry:
                dropped += 1
                continue
            vector = np.array([entry["outbound"], entry["return"]])
            if not np.all(np.isfinite(vector)):
                dropped += 1
                continue
            values = {
                "outbound": vector[0],
                "return": vector[1],
                "equal": float(w["equal"] @ vector),
                "invvar": float(w["invvar"] @ vector),
                "joint": float(w["joint"] @ vector),
            }
            for name in estimators:
                injected_by[name].append(level)
                estimated_by[name].append(values[name])
            rows.append({
                "epoch": epoch,
                "injected_dvv": level,
                **{f"eps_{name}": values[name] for name in estimators},
                "w_invvar_out": w["invvar"][0],
                "w_joint_out": w["joint"][0],
                "sigma_out": w["sigma_out"],
                "sigma_ret": w["sigma_ret"],
                "rho_loo": w["rho"],
            })
    if dropped:
        print(f"dropped {dropped} incomplete or non-finite pairs")

    with OUT_CSV.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    analyses = {}
    for name in estimators:
        analyses[name] = D._analyse_group(
            np.asarray(injected_by[name]), np.asarray(estimated_by[name])
        )

    nano = np.load(D.HERE / "nano_dvv_injection_recovery.npz", allow_pickle=True)
    nano_threshold = float(nano["empirical_95pct_threshold_dvv"])
    outbound_threshold = analyses["outbound"]["threshold"]

    labels = {
        "outbound": "Deep outbound (benchmark)",
        "return": "Deep return",
        "equal": "Deep paired, equal weight",
        "invvar": "Deep paired, inverse-variance",
        "joint": "Deep paired, covariance-aware",
    }
    table = []
    for name in estimators:
        analysis = analyses[name]
        zero_bias = next(
            r["bias"] for r in analysis["summary_rows"] if r["injected_dvv"] == 0
        )
        scatter = next(
            r["robust_scatter_1p4826_mad"]
            for r in analysis["summary_rows"]
            if r["injected_dvv"] == 0
        )
        table.append({
            "observable": labels[name],
            "null_threshold_dvv": f"{analysis['threshold']:.6g}",
            "reliable_detection": f"{analysis['reliable_detection']:.6g}",
            "reliable_sign_only": f"{analysis['reliable_sign']:.6g}",
            "null_scatter_mad": f"{scatter:.6g}",
            "zero_injection_bias": f"{zero_bias:.6g}",
            "false_positive_rate": f"{analysis['false_positive_rate']:.4f}",
            "vs_outbound_threshold": f"{outbound_threshold / analysis['threshold']:.3f}",
            "vs_nano_threshold": f"{nano_threshold / analysis['threshold']:.3f}",
        })
    with OUT_TABLE.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)

    best = min(estimators, key=lambda n: analyses[n]["threshold"])
    ratio = outbound_threshold / analyses[best]["threshold"]
    if best == "outbound" or ratio < 1.05:
        verdict = (
            "Combining the legs does not improve on the outbound branch alone."
        )
    elif best in ("equal", "invvar", "joint"):
        verdict = (
            f"The {labels[best]} estimator improves on outbound by "
            f"{ratio:.2f}x in null threshold."
        )
    worse = [
        n for n in ("equal", "invvar", "joint")
        if analyses[n]["threshold"] > outbound_threshold * 1.05
    ]

    mean_w = np.mean([weights[e]["joint"][0] for e in epochs])
    lines = [
        "Deep guided mode: paired outbound/return combination",
        "=" * 70,
        f"Population {POPULATION}; band {BAND.replace('_', '-')} Hz; pass {PASS}.",
        f"{len(epochs)} paired bursts x {len(levels)} levels; "
        f"weights leave-one-burst-out from zero-injection pairs.",
        "",
        f"Between-leg null error correlation: rho = {correlation:+.3f}",
        f"Mean covariance-aware weight on outbound: {mean_w:.3f} "
        f"(0.5 would be equal weight)",
        "",
        f"{'observable':32s} {'threshold':>11s} {'reliable':>10s} "
        f"{'scatter':>10s} {'bias':>11s} {'vs out':>7s}",
        "-" * 88,
    ]
    for row in table:
        lines.append(
            f"{row['observable']:32s} {float(row['null_threshold_dvv']):11.4e} "
            f"{float(row['reliable_detection']):10.1e} "
            f"{float(row['null_scatter_mad']):10.4e} "
            f"{float(row['zero_injection_bias']):+11.3e} "
            f"{float(row['vs_outbound_threshold']):7.3f}"
        )
    lines += [
        "",
        f"Nano reference threshold {nano_threshold:.4e}, reliable detection "
        f"{float(nano['symmetric_95pct_detection_limit_dvv']):.1e}",
        "",
        "Verdict",
        "-" * 70,
        f"  {verdict}",
    ]
    if worse:
        lines.append(
            f"  Degraded relative to outbound: {', '.join(labels[n] for n in worse)}."
        )
    lines += [
        "",
        "Reading rho: the two legs record the same source impacts, so their",
        "errors are correlated and the second leg carries less independent",
        "information than its own scatter alone would suggest. That correlation",
        "bounds how much any combination can gain, whatever the weighting.",
        "",
        "Outbound was the pre-registered benchmark. A paired estimator that does",
        "not beat it does not replace it.",
        "",
    ]
    report = "\n".join(lines) + "\n"
    OUT_TXT.write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
