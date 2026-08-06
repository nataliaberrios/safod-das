"""Solid-Earth-tide regression on the Deep guided-mode per-burst estimates.

Paper 1 reports a pre-registered tidal null on the Nano observable: amplitude
3.4e-4 +/- 2.3e-4, 95% upper limit 7.98e-4, or 1.6x the De Fazio et al. (1973)
scale of 5e-4. The Deep observable has a lower per-burst scatter, so a
scatter-ratio projection suggests its upper limit could fall near or below the
De Fazio scale. This script replaces that projection with a fit.

What this is, and is not
------------------------
This is the same *style* of regression applied to a different observable. It is
NOT a re-run of Paper 1's registered procedure: that fit used the depth-median
dv/v chain with common-mode removal over 46 epochs, whereas the primary
population here is the 23 held-out bursts of the Deep injection-recovery test.
The two upper limits are therefore comparable in scale but not identical in
construction, and the comparison should be stated that way.

The measurement
---------------
The zero-injection trials already are a dv/v time series: for burst b, the
recovered eps is the apparent-velocity change of that burst against a
leave-one-out reference built from the others. No new waveform processing is
needed.

One correction is applied. A leave-one-out reference makes

    eps_b = N/(N-1) * (x_b - mean(x)),

so the raw estimates are inflated by N/(N-1) relative to the per-burst deviation
from the mean. With N = 23 that is 4.3%, and it is divided out before fitting.

The tide model is the compact degree-2 formulation from ``safod_tides.ipynb``,
reproduced here so the fit does not depend on a notebook outside the repository.
It yields a scalar areal-strain proxy, used only as a *shape* -- normalised to
unit amplitude -- so the fitted coefficient carries the units of dv/v. No attempt
is made to convert strain to a signed velocity change through a constitutive law.

Degeneracy warning
------------------
The survey is 24 h long, so a diurnal tidal component is strongly degenerate with
instrumental drift, and this experiment has known drift (compaction and plate
settlement). The primary model therefore includes a linear trend. The no-trend
model is reported alongside because the difference between them measures how much
of the apparent tidal amplitude is really drift.

Inference is by surrogate null: the tide predictor is re-evaluated at randomly
shifted times and the fit repeated, giving an empirical p-value that does not
assume Gaussian independent errors.
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timedelta, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import deep_dvv_injection_recovery as D


# SAFOD site and physical constants, from safod_tides.ipynb.
LATITUDE_DEG = 35.9746
LONGITUDE_DEG = -120.5516
G = 6.674e-11
EARTH_RADIUS_M = 6.371e6
SURFACE_GRAVITY_M_S2 = 9.81
MOON_MASS_KG = 7.342e22
SUN_MASS_KG = 1.989e30
ASTRONOMICAL_UNIT_M = 1.495978707e11
H2 = 0.603741
L2 = 0.084010
AREAL_LOVE_FACTOR = 2.0 * H2 - 6.0 * L2

# Benchmarks.
DE_FAZIO_DVV = 5.0e-4
NIU_DVV = 5.76e-5  # 2.4e-7 /Pa * 240 Pa
PAPER1_NANO_UL = 7.98e-4

N_SURROGATE = 999
SEED = 20260806

ESTIMATORS = ("outbound", "return", "equal", "invvar", "joint")
LABELS = {
    "outbound": "Deep outbound",
    "return": "Deep return",
    "equal": "Deep paired, equal weight",
    "invvar": "Deep paired, inverse-variance",
    "joint": "Deep paired, covariance-aware",
}

OUT_CSV = D.HERE / "deep_dvv_tidal_fit.csv"
OUT_TXT = D.HERE / "deep_dvv_tidal_fit.txt"
OUT_PNG = D.HERE / "deep_dvv_tidal_fit.png"


# --------------------------------------------------------------------------
# Tide model (degree-2 potential -> areal-strain proxy)
# --------------------------------------------------------------------------
def julian_date(times) -> np.ndarray:
    return np.array([t.timestamp() for t in times], dtype=float) / 86400.0 + 2440587.5


def sun_position(jd):
    d = jd - 2451545.0
    mean_longitude = np.deg2rad((280.460 + 0.9856474 * d) % 360.0)
    mean_anomaly = np.deg2rad((357.528 + 0.9856003 * d) % 360.0)
    ecliptic_longitude = (
        mean_longitude
        + np.deg2rad(1.915) * np.sin(mean_anomaly)
        + np.deg2rad(0.020) * np.sin(2.0 * mean_anomaly)
    )
    distance_au = (
        1.00014 - 0.01671 * np.cos(mean_anomaly) - 0.00014 * np.cos(2.0 * mean_anomaly)
    )
    obliquity = np.deg2rad(23.4393 - 3.563e-7 * d)
    ra = np.arctan2(
        np.cos(obliquity) * np.sin(ecliptic_longitude), np.cos(ecliptic_longitude)
    )
    dec = np.arcsin(np.sin(obliquity) * np.sin(ecliptic_longitude))
    return ra, dec, distance_au * ASTRONOMICAL_UNIT_M


def moon_position(jd):
    d = jd - 2451543.5
    node = np.deg2rad((125.1228 - 0.0529538083 * d) % 360.0)
    inclination = np.deg2rad(5.1454)
    perigee = np.deg2rad((318.0634 + 0.1643573223 * d) % 360.0)
    semimajor_axis = 60.2666
    eccentricity = 0.054900
    mean_anomaly = np.deg2rad((115.3654 + 13.0649929509 * d) % 360.0)
    eccentric_anomaly = mean_anomaly.copy()
    for _ in range(3):
        eccentric_anomaly -= (
            eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly) - mean_anomaly
        ) / (1.0 - eccentricity * np.cos(eccentric_anomaly))
    x_orbit = semimajor_axis * (np.cos(eccentric_anomaly) - eccentricity)
    y_orbit = semimajor_axis * np.sqrt(1.0 - eccentricity**2) * np.sin(eccentric_anomaly)
    radius = np.hypot(x_orbit, y_orbit)
    argument = np.arctan2(y_orbit, x_orbit) + perigee
    x_ecl = radius * (
        np.cos(node) * np.cos(argument)
        - np.sin(node) * np.sin(argument) * np.cos(inclination)
    )
    y_ecl = radius * (
        np.sin(node) * np.cos(argument)
        + np.cos(node) * np.sin(argument) * np.cos(inclination)
    )
    z_ecl = radius * np.sin(argument) * np.sin(inclination)
    obliquity = np.deg2rad(23.4393 - 3.563e-7 * (jd - 2451545.0))
    y_eq = y_ecl * np.cos(obliquity) - z_ecl * np.sin(obliquity)
    z_eq = y_ecl * np.sin(obliquity) + z_ecl * np.cos(obliquity)
    ra = np.arctan2(y_eq, x_ecl)
    dec = np.arctan2(z_eq, np.hypot(x_ecl, y_eq))
    return ra, dec, radius * EARTH_RADIUS_M


def cos_zenith_angle(jd, ra, dec):
    d = jd - 2451545.0
    lst = np.deg2rad((280.46061837 + 360.98564736629 * d + LONGITUDE_DEG) % 360.0)
    hour_angle = lst - ra
    latitude = np.deg2rad(LATITUDE_DEG)
    return np.sin(latitude) * np.sin(dec) + np.cos(latitude) * np.cos(dec) * np.cos(
        hour_angle
    )


def areal_strain(times) -> np.ndarray:
    jd = julian_date(times)
    sun_ra, sun_dec, sun_distance = sun_position(jd)
    moon_ra, moon_dec, moon_distance = moon_position(jd)

    def potential(mass, distance, cosz):
        return G * mass / distance**3 * EARTH_RADIUS_M**2 * 0.5 * (3.0 * cosz**2 - 1.0)

    total = potential(
        SUN_MASS_KG, sun_distance, cos_zenith_angle(jd, sun_ra, sun_dec)
    ) + potential(MOON_MASS_KG, moon_distance, cos_zenith_angle(jd, moon_ra, moon_dec))
    return AREAL_LOVE_FACTOR / (SURFACE_GRAVITY_M_S2 * EARTH_RADIUS_M) * total


# --------------------------------------------------------------------------
# Regression
# --------------------------------------------------------------------------
def _fit(design: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ coefficients
    dof = max(1, y.size - design.shape[1])
    sigma2 = float(residual @ residual) / dof
    covariance = sigma2 * np.linalg.pinv(design.T @ design)
    return coefficients, np.sqrt(np.maximum(np.diag(covariance), 0.0))


def _design(shape: np.ndarray, hours: np.ndarray, with_trend: bool) -> np.ndarray:
    columns = [shape, np.ones_like(shape)]
    if with_trend:
        columns.append(hours - hours.mean())
    return np.column_stack(columns)


def main() -> None:
    rng = np.random.default_rng(SEED)

    # Burst times for the held-out epochs.
    with np.load(D.STACKS) as data:
        begtimes = data["begtimes_str"].astype(str)
    epoch_time = {}
    for index, text in enumerate(begtimes):
        if not text or text == "":
            continue
        epoch_time[index] = datetime.fromisoformat(text).astimezone(timezone.utc)

    rows = list(csv.DictReader(open(D.HERE / "deep_dvv_paired_legs.csv")))
    null_rows = [r for r in rows if float(r["injected_dvv"]) == 0.0]
    epochs = [int(r["epoch"]) for r in null_rows]
    missing = [e for e in epochs if e not in epoch_time]
    if missing:
        raise RuntimeError(f"No timestamp for epochs {missing}")

    times = [epoch_time[e] for e in epochs]
    order = np.argsort([t.timestamp() for t in times])
    epochs = [epochs[i] for i in order]
    times = [times[i] for i in order]
    null_rows = [null_rows[i] for i in order]

    n = len(epochs)
    loo_correction = (n - 1) / n  # undo the leave-one-out inflation
    hours = np.array(
        [(t - times[0]).total_seconds() / 3600.0 for t in times], dtype=float
    )
    strain = areal_strain(times)
    shape = strain - strain.mean()
    scale = np.max(np.abs(shape))
    shape = shape / scale

    print(f"{n} held-out bursts spanning {hours[-1]:.2f} h")
    print(
        f"tide proxy over the burst times: {np.ptp(strain) * 1e9:.1f} nanostrain "
        f"peak to peak"
    )

    results = []
    for name in ESTIMATORS:
        y = np.array([float(r[f"eps_{name}"]) for r in null_rows]) * loo_correction
        for with_trend in (True, False):
            design = _design(shape, hours, with_trend)
            coefficients, errors = _fit(design, y)
            amplitude = float(coefficients[0])
            sigma = float(errors[0])

            surrogate = np.empty(N_SURROGATE)
            for k in range(N_SURROGATE):
                offset = float(rng.uniform(1.0, 23.0))
                shifted_times = [t + timedelta(hours=offset) for t in times]
                shifted = areal_strain(shifted_times)
                shifted = shifted - shifted.mean()
                shifted = shifted / np.max(np.abs(shifted))
                surrogate[k] = _fit(_design(shifted, hours, with_trend), y)[0][0]
            p_value = float(
                (1 + np.sum(np.abs(surrogate) >= abs(amplitude))) / (N_SURROGATE + 1)
            )
            upper_limit = abs(amplitude) + 1.96 * sigma
            surrogate_limit = float(np.quantile(np.abs(surrogate), 0.95))

            results.append(
                {
                    "estimator": name,
                    "label": LABELS[name],
                    "model": "with_linear_trend" if with_trend else "no_trend",
                    "n_bursts": n,
                    "amplitude_dvv": amplitude,
                    "sigma_dvv": sigma,
                    "p_surrogate": p_value,
                    "upper_limit_95_dvv": upper_limit,
                    "surrogate_95_dvv": surrogate_limit,
                    "vs_de_fazio": upper_limit / DE_FAZIO_DVV,
                    "vs_niu": upper_limit / NIU_DVV,
                    "vs_paper1_nano_ul": upper_limit / PAPER1_NANO_UL,
                }
            )
            print(
                f"  {name:9s} {'trend' if with_trend else 'notrend':8s} "
                f"A={amplitude:+.3e} +/- {sigma:.3e}  p={p_value:.3f}  "
                f"UL={upper_limit:.3e} ({upper_limit / DE_FAZIO_DVV:.2f}x De Fazio)"
            )

    with OUT_CSV.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    primary = [r for r in results if r["model"] == "with_linear_trend"]
    best = min(primary, key=lambda r: r["upper_limit_95_dvv"])

    # ---------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.6), constrained_layout=True)
    ax = axes[0]
    dense = [times[0] + timedelta(minutes=10 * i) for i in range(int(hours[-1] * 6) + 1)]
    dense_strain = areal_strain(dense)
    ax.plot(
        [(t - times[0]).total_seconds() / 3600.0 for t in dense],
        (dense_strain - dense_strain.mean()) * 1e9,
        color="0.4", lw=1.2, label="degree-2 areal-strain proxy",
    )
    ax.plot(hours, (strain - strain.mean()) * 1e9, "o", color="#b2182b", ms=5,
            label="burst times")
    ax.set(xlabel="hours from first held-out burst", ylabel="nanostrain",
           title="A  Tidal forcing over the survey")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    y = np.array([float(r["eps_outbound"] if False else r["eps_joint"]) for r in null_rows]) * loo_correction
    fit = [r for r in primary if r["estimator"] == "joint"][0]
    ax.plot(hours, y, "o", color="#2166ac", ms=5, label="Deep paired, covariance-aware")
    ax.plot(hours, fit["amplitude_dvv"] * shape, "-", color="#b2182b",
            label=f"fitted tidal term, A={fit['amplitude_dvv']:+.2e}")
    ax.axhline(0, color="0.3", lw=0.8)
    ax.set(xlabel="hours from first held-out burst",
           ylabel="apparent velocity change", title="B  Observation and tidal fit")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[2]
    names = [r["label"] for r in primary]
    limits = [r["upper_limit_95_dvv"] for r in primary]
    ax.barh(range(len(names)), limits, color="#4393c3")
    ax.axvline(DE_FAZIO_DVV, color="#b2182b", ls="--", lw=1.4,
               label=f"De Fazio scale {DE_FAZIO_DVV:.0e}")
    ax.axvline(NIU_DVV, color="#7b3294", ls=":", lw=1.4,
               label=f"Niu SAFOD scale {NIU_DVV:.2e}")
    ax.axvline(PAPER1_NANO_UL, color="0.35", ls="-.", lw=1.2,
               label=f"Paper 1 Nano UL {PAPER1_NANO_UL:.2e}")
    ax.set_xscale("log")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set(xlabel="95% upper limit on tidal dv/v amplitude",
           title="C  Upper limits against benchmarks")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.suptitle(
        "Deep guided mode: solid-Earth-tide regression on 23 held-out bursts\n"
        "upper limits on an apparent guided-mode velocity change, not a formation "
        "velocity change", fontsize=11,
    )
    fig.savefig(OUT_PNG, dpi=250)
    plt.close(fig)

    # ---------------------------------------------------------------- report
    lines = [
        "Deep guided mode: solid-Earth-tide regression",
        "=" * 74,
        f"{n} held-out bursts over {hours[-1]:.2f} h; "
        f"tide proxy {np.ptp(strain) * 1e9:.1f} nanostrain peak to peak.",
        f"Leave-one-out inflation removed (factor {loo_correction:.4f}).",
        f"Surrogate null: {N_SURROGATE} random time shifts of the tide predictor.",
        "",
        "Primary model, tidal term plus constant plus linear trend",
        "-" * 74,
        f"{'observable':32s} {'amplitude':>12s} {'sigma':>11s} {'p':>7s} "
        f"{'95% UL':>11s} {'/DeFazio':>9s}",
    ]
    for row in primary:
        lines.append(
            f"{row['label']:32s} {row['amplitude_dvv']:+12.3e} "
            f"{row['sigma_dvv']:11.3e} {row['p_surrogate']:7.3f} "
            f"{row['upper_limit_95_dvv']:11.3e} {row['vs_de_fazio']:9.2f}"
        )
    lines += [
        "",
        "Sensitivity to the trend term (24 h span makes diurnal tide and drift",
        "degenerate; a large change here means the tidal term is absorbing drift)",
        "-" * 74,
    ]
    for name in ESTIMATORS:
        a = [r for r in results if r["estimator"] == name and r["model"] == "with_linear_trend"][0]
        b = [r for r in results if r["estimator"] == name and r["model"] == "no_trend"][0]
        lines.append(
            f"  {LABELS[name]:32s} UL {a['upper_limit_95_dvv']:.3e} with trend, "
            f"{b['upper_limit_95_dvv']:.3e} without "
            f"({b['upper_limit_95_dvv'] / a['upper_limit_95_dvv']:.2f}x)"
        )
    lines += [
        "",
        "Benchmarks",
        "-" * 74,
        f"  De Fazio et al. (1973) scale      {DE_FAZIO_DVV:.2e}",
        f"  Niu et al. (2008) SAFOD scale     {NIU_DVV:.2e}",
        f"  Paper 1 Nano 95% upper limit      {PAPER1_NANO_UL:.2e} "
        f"({PAPER1_NANO_UL / DE_FAZIO_DVV:.2f}x De Fazio)",
        "",
        f"  Best Deep upper limit: {best['label']}, {best['upper_limit_95_dvv']:.3e} "
        f"= {best['vs_de_fazio']:.2f}x De Fazio, {best['vs_niu']:.1f}x Niu",
        "",
        "Reading this",
        "-" * 74,
        "  Every surrogate p-value should be read as a null result unless it is",
        "  small: no tidal detection is claimed here. The quantity of interest is",
        "  the upper limit, which states what the experiment could have excluded.",
        "  This is the same style of regression as Paper 1's registered tidal null,",
        "  not a re-run of it: that fit used the depth-median dv/v chain with",
        "  common-mode removal over 46 epochs. The limits are comparable in scale",
        "  but not identical in construction.",
        "  The limit applies to an apparent guided-mode velocity change and is not",
        "  a formation Vp or Vs limit.",
        "",
    ]
    report = "\n".join(lines) + "\n"
    OUT_TXT.write_text(report)
    print()
    print(report, end="")


if __name__ == "__main__":
    main()
