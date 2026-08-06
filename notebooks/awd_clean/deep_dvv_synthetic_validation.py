"""Per-leg synthetic geometry and injection validation for the Deep observable.

Pre-registration criterion "SNR-limited 1": before a weaker return-leg result
may be attributed to signal-to-noise, the return geometry and injection must be
shown to work on noiseless synthetic data.  A leg that fails here has a
trajectory, geometry, or estimator problem, not an SNR problem.

Two independent paths are compared at each injected level:

* ``physical``  -- the guided mode is re-synthesised with its speed genuinely
  scaled by ``1 - eps``, then run through the whole pipeline.
* ``pipeline``  -- the unperturbed mode is perturbed by the injection routine
  used in the real analysis.

If the two agree, the injection is a faithful proxy for a medium change rather
than a shift the estimator is guaranteed to find.  Both must also return the
injected value with the correct sign.

Reads only the frozen trajectory.  Touches no recovered value.
"""

from __future__ import annotations

import csv
import json

import numpy as np

import deep_dvv_injection_recovery as D


TEST_LEVELS = (-5e-3, -1e-3, -1e-4, 0.0, 1e-4, 1e-3, 5e-3)
CENTER_FREQUENCY_HZ = {"15_30": 22.0, "3_15": 9.0, "60_120": 90.0}
TOLERANCE = 0.05

OUT_CSV = D.HERE / "deep_dvv_synthetic_validation.csv"
OUT_TXT = D.HERE / "deep_dvv_synthetic_validation.txt"


def _ricker(t: np.ndarray, f0: float) -> np.ndarray:
    a = (np.pi * f0 * t) ** 2
    return (1 - 2 * a) * np.exp(-a)


def _synthetic(coordinate, sample_time, intercept, slowness, f0, eps):
    """A clean guided mode whose speed is faster by the fraction ``eps``."""
    travel = (intercept + slowness * coordinate) * (1.0 - eps)
    return np.stack([_ricker(sample_time - t, f0) for t in travel])


def main() -> None:
    frozen = json.loads(D.FROZEN.read_text())
    fs = float(frozen["fs_hz"])
    dx = float(frozen["dx_deep_m"])
    sample_time = np.arange(3500, dtype=float) / fs - D.PRE_S
    coordinate = np.arange(
        D.COORD_RANGE_M[0], D.COORD_RANGE_M[1], dx * D.CHANNEL_STRIDE
    )
    layout = D._aperture_layout(coordinate, D.APERTURE_M)
    centers = np.asarray([c for _, c in layout])

    rows = []
    for leg in D.LEGS:
        for band_tag, config in D.BAND_CONFIG.items():
            trajectory = frozen["trajectories"][f"{leg}|{band_tag}"]
            intercept = trajectory["intercept_s"]
            slowness = trajectory["slowness_s_per_m"]
            f0 = CENTER_FREQUENCY_HZ[band_tag]
            q0, q1 = config["extract_q_s"]
            q = np.arange(int(round((q1 - q0) * fs)), dtype=float) / fs + q0
            travel = intercept + slowness * coordinate
            travel_centers = intercept + slowness * centers

            clean = _synthetic(coordinate, sample_time, intercept, slowness, f0, 0.0)
            reference = D._align(clean, sample_time, coordinate, intercept, slowness, q)
            reference_beams = D._beams(reference, q, coordinate, layout, None)

            for eps in TEST_LEVELS:
                physical = D._align(
                    _synthetic(coordinate, sample_time, intercept, slowness, f0, eps),
                    sample_time, coordinate, intercept, slowness, q,
                )
                pipeline = D._shift_gather(reference, q, -eps * travel)
                for path, gather in (("physical", physical), ("pipeline", pipeline)):
                    result = D._estimate(
                        D._beams(gather, q, coordinate, layout, None),
                        reference_beams, layout, travel_centers, q, fs, config,
                    )
                    estimated = result["estimated_dvv"]
                    rows.append({
                        "leg": leg,
                        "band": band_tag,
                        "role": config["role"],
                        "path": path,
                        "injected_dvv": eps,
                        "estimated_dvv": estimated,
                        "ratio": (estimated / eps) if eps else np.nan,
                        "n_apertures": result["n_apertures"],
                        "median_aperture_correlation": result["median_aperture_correlation"],
                    })

    with OUT_CSV.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Deep guided mode: per-leg synthetic geometry and injection validation",
        "=" * 70,
        "Pre-registration criterion 'SNR-limited 1'. Noiseless synthetic mode on the",
        "real Deep geometry, using each leg's own frozen trajectory.",
        "",
    ]
    verdicts = {}
    for leg in D.LEGS:
        for band_tag in D.BAND_CONFIG:
            subset = [r for r in rows if r["leg"] == leg and r["band"] == band_tag]
            nonzero = [r for r in subset if r["injected_dvv"] != 0]
            zero = [r for r in subset if r["injected_dvv"] == 0]
            ratios = np.asarray([r["ratio"] for r in nonzero], dtype=float)
            signs_ok = all(
                np.sign(r["estimated_dvv"]) == np.sign(r["injected_dvv"])
                for r in nonzero
                if np.isfinite(r["estimated_dvv"])
            )
            worst = float(np.nanmax(np.abs(ratios - 1.0))) if ratios.size else np.nan
            zero_worst = float(
                np.nanmax([abs(r["estimated_dvv"]) for r in zero])
            ) if zero else np.nan
            agreement = float(
                np.nanmax([
                    abs(a["estimated_dvv"] - b["estimated_dvv"])
                    for a in subset if a["path"] == "physical"
                    for b in subset if b["path"] == "pipeline"
                    and a["injected_dvv"] == b["injected_dvv"]
                ])
            )
            passed = bool(
                signs_ok and np.isfinite(worst) and worst < TOLERANCE
            )
            verdicts[(leg, band_tag)] = passed
            lines.append(
                f"{leg:9s} {band_tag:6s}  "
                f"{'PASS' if passed else 'FAIL'}  "
                f"worst scale error {worst:.4f}, "
                f"zero-injection |eps| {zero_worst:.2e}, "
                f"physical-vs-pipeline {agreement:.2e}, "
                f"signs {'all correct' if signs_ok else 'INCORRECT'}"
            )
    lines += [
        "",
        "A leg that fails here cannot have a weaker real result attributed to",
        "signal-to-noise; the pre-registration classifies that as a trajectory,",
        "geometry, or estimator failure instead.",
        "",
    ]
    report = "\n".join(lines) + "\n"
    OUT_TXT.write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
