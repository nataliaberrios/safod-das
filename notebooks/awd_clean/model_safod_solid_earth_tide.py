"""Model the first-pass solid-Earth tide at SAFOD for the June 2026 AWD test.

This is a transparent, low-precision body-tide calculation.  It computes
approximate Sun and Moon ephemerides, the degree-2 tidal potential, and the
scalar areal-strain proxy documented in ``solid_earth_tide_model.tex``.

It is a forcing/template model, not a complete downhole DAS strain model:
the latter would calculate the full strain tensor and project it onto the
local fiber direction.

Run from this directory with::

    python model_safod_solid_earth_tide.py

Outputs are written beside this script:
``safod_solid_earth_tide_2026-06-17.png``
and ``safod_solid_earth_tide_2026-06-17.csv``.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent

# SAFOD site used in the methods note.
LATITUDE_DEG = 35.9746
LONGITUDE_DEG = -120.5516

# Exact window represented in awd_earth_tide_starter.ipynb / the methods note.
START_UTC = datetime(2026, 6, 16, 23, 48, tzinfo=timezone.utc)
END_UTC = datetime(2026, 6, 17, 23, 44, tzinfo=timezone.utc)
SAMPLE_MINUTES = 1

# Physical constants and degree-2 Love numbers from the methods note.
G = 6.674e-11
EARTH_RADIUS_M = 6.371e6
SURFACE_GRAVITY_M_S2 = 9.81
MOON_MASS_KG = 7.342e22
SUN_MASS_KG = 1.989e30
MOON_DISTANCE_UNIT_M = EARTH_RADIUS_M
H2 = 0.6032
L2 = 0.0839
AREAL_LOVE_FACTOR = 2.0 * H2 - 6.0 * L2


def datetime_range(start: datetime, end: datetime, step_minutes: int) -> list[datetime]:
    """Return inclusive UTC datetimes at a fixed minute spacing."""
    count = int((end - start).total_seconds() // (60.0 * step_minutes))
    return [start + timedelta(minutes=step_minutes * i) for i in range(count + 1)]


def julian_date(times: list[datetime]) -> np.ndarray:
    """Convert timezone-aware UTC datetimes to Julian Date."""
    unix_seconds = np.array([time.timestamp() for time in times], dtype=float)
    return unix_seconds / 86400.0 + 2440587.5


def sun_position(jd: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return approximate Sun right ascension, declination, and distance.

    Angles are radians and distance is metres.  This short-period solar
    ephemeris is adequate for a one-day diagnostic tide plot.
    """
    d = jd - 2451545.0
    mean_longitude = np.deg2rad((280.460 + 0.9856474 * d) % 360.0)
    mean_anomaly = np.deg2rad((357.528 + 0.9856003 * d) % 360.0)
    ecliptic_longitude = (
        mean_longitude
        + np.deg2rad(1.915) * np.sin(mean_anomaly)
        + np.deg2rad(0.020) * np.sin(2.0 * mean_anomaly)
    )
    distance_au = (
        1.00014
        - 0.01671 * np.cos(mean_anomaly)
        - 0.00014 * np.cos(2.0 * mean_anomaly)
    )
    obliquity = np.deg2rad(23.4393 - 3.563e-7 * d)
    right_ascension = np.arctan2(
        np.cos(obliquity) * np.sin(ecliptic_longitude),
        np.cos(ecliptic_longitude),
    )
    declination = np.arcsin(
        np.sin(obliquity) * np.sin(ecliptic_longitude)
    )
    astronomical_unit_m = 1.495978707e11
    return right_ascension, declination, distance_au * astronomical_unit_m


def moon_position(jd: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an approximate lunar position using a truncated orbital model.

    The equations follow the compact orbital-element construction used in the
    methods note: one Newton-style correction for eccentric anomaly followed
    by rotations from the lunar orbital plane to equatorial coordinates.
    Angles are radians and distance is metres.
    """
    # Days from 2000 January 0.0, the epoch used by this compact lunar model.
    d = jd - 2451543.5
    ascending_node = np.deg2rad((125.1228 - 0.0529538083 * d) % 360.0)
    inclination = np.deg2rad(5.1454)
    argument_perigee = np.deg2rad((318.0634 + 0.1643573223 * d) % 360.0)
    semimajor_axis = 60.2666
    eccentricity = 0.054900
    mean_anomaly = np.deg2rad((115.3654 + 13.0649929509 * d) % 360.0)

    # First-order eccentric-anomaly correction with the standard e*cos(M)
    # refinement.  This is more than sufficient for the plotted one-day phase.
    eccentric_anomaly = mean_anomaly + eccentricity * np.sin(mean_anomaly) * (
        1.0 + eccentricity * np.cos(mean_anomaly)
    )
    x_orbit = semimajor_axis * (np.cos(eccentric_anomaly) - eccentricity)
    y_orbit = semimajor_axis * np.sqrt(1.0 - eccentricity**2) * np.sin(
        eccentric_anomaly
    )
    orbital_radius = np.hypot(x_orbit, y_orbit)
    true_anomaly = np.arctan2(y_orbit, x_orbit)
    argument = true_anomaly + argument_perigee

    x_ecliptic = orbital_radius * (
        np.cos(ascending_node) * np.cos(argument)
        - np.sin(ascending_node) * np.sin(argument) * np.cos(inclination)
    )
    y_ecliptic = orbital_radius * (
        np.sin(ascending_node) * np.cos(argument)
        + np.cos(ascending_node) * np.sin(argument) * np.cos(inclination)
    )
    z_ecliptic = orbital_radius * np.sin(argument) * np.sin(inclination)

    obliquity = np.deg2rad(23.4393 - 3.563e-7 * (jd - 2451545.0))
    x_equatorial = x_ecliptic
    y_equatorial = y_ecliptic * np.cos(obliquity) - z_ecliptic * np.sin(obliquity)
    z_equatorial = y_ecliptic * np.sin(obliquity) + z_ecliptic * np.cos(obliquity)
    right_ascension = np.arctan2(y_equatorial, x_equatorial)
    declination = np.arctan2(
        z_equatorial, np.hypot(x_equatorial, y_equatorial)
    )
    distance_m = orbital_radius * MOON_DISTANCE_UNIT_M
    return right_ascension, declination, distance_m


def local_zenith_cosine(
    jd: np.ndarray,
    right_ascension: np.ndarray,
    declination: np.ndarray,
) -> np.ndarray:
    """Return cos(theta), where theta is body-to-site zenith angle."""
    d = jd - 2451545.0
    local_sidereal_time = np.deg2rad(
        (280.46061837 + 360.98564736629 * d + LONGITUDE_DEG) % 360.0
    )
    hour_angle = local_sidereal_time - right_ascension
    latitude = np.deg2rad(LATITUDE_DEG)
    return (
        np.sin(latitude) * np.sin(declination)
        + np.cos(latitude) * np.cos(declination) * np.cos(hour_angle)
    )


def degree_two_potential(
    mass_kg: float,
    distance_m: np.ndarray,
    cos_zenith: np.ndarray,
) -> np.ndarray:
    """Calculate the degree-2 tidal potential in m^2/s^2."""
    legendre_p2 = 0.5 * (3.0 * cos_zenith**2 - 1.0)
    return G * mass_kg / distance_m**3 * EARTH_RADIUS_M**2 * legendre_p2


def model_tide(times: list[datetime]) -> dict[str, np.ndarray]:
    """Calculate Sun, Moon, and total scalar tidal strain time series."""
    jd = julian_date(times)
    sun_ra, sun_dec, sun_distance = sun_position(jd)
    moon_ra, moon_dec, moon_distance = moon_position(jd)
    sun_cos_zenith = local_zenith_cosine(jd, sun_ra, sun_dec)
    moon_cos_zenith = local_zenith_cosine(jd, moon_ra, moon_dec)

    sun_potential = degree_two_potential(
        SUN_MASS_KG, sun_distance, sun_cos_zenith
    )
    moon_potential = degree_two_potential(
        MOON_MASS_KG, moon_distance, moon_cos_zenith
    )
    forcing = (sun_potential + moon_potential) / (
        SURFACE_GRAVITY_M_S2 * EARTH_RADIUS_M
    )
    strain_scale = AREAL_LOVE_FACTOR * forcing
    return {
        "jd": jd,
        "sun_potential": sun_potential,
        "moon_potential": moon_potential,
        "sun_strain": AREAL_LOVE_FACTOR * sun_potential / (SURFACE_GRAVITY_M_S2 * EARTH_RADIUS_M),
        "moon_strain": AREAL_LOVE_FACTOR * moon_potential / (SURFACE_GRAVITY_M_S2 * EARTH_RADIUS_M),
        "total_strain": strain_scale,
        "sun_distance_m": sun_distance,
        "moon_distance_m": moon_distance,
        "sun_zenith_deg": np.rad2deg(np.arccos(np.clip(sun_cos_zenith, -1.0, 1.0))),
        "moon_zenith_deg": np.rad2deg(np.arccos(np.clip(moon_cos_zenith, -1.0, 1.0))),
    }


def save_csv(times: list[datetime], result: dict[str, np.ndarray], path: Path) -> None:
    """Write the modeled quantities for later template use or inspection."""
    columns = [
        np.array([time.isoformat() for time in times], dtype=object),
        result["total_strain"],
        result["sun_strain"],
        result["moon_strain"],
        result["sun_zenith_deg"],
        result["moon_zenith_deg"],
    ]
    header = "utc,total_areal_strain,sun_areal_strain,moon_areal_strain,sun_zenith_deg,moon_zenith_deg"
    np.savetxt(path, np.column_stack(columns), delimiter=",", header=header, fmt="%s")


def plot_result(times: list[datetime], result: dict[str, np.ndarray], path: Path) -> None:
    """Plot body contributions and the normalized matched-filter template."""
    total = result["total_strain"]
    mean_removed = total - np.mean(total)
    template = mean_removed / np.max(np.abs(mean_removed))
    time_array = np.array(times)
    nano = 1e9

    fig, axes = plt.subplots(
        2, 1, figsize=(11, 7), sharex=True, constrained_layout=True
    )
    axes[0].plot(time_array, result["sun_strain"] * nano, label="Sun", lw=1.8)
    axes[0].plot(time_array, result["moon_strain"] * nano, label="Moon", lw=1.8)
    axes[0].plot(
        time_array,
        total * nano,
        label="Sun + Moon",
        color="black",
        lw=2.2,
    )
    axes[0].axhline(0.0, color="0.5", lw=0.7)
    axes[0].set_ylabel("nanostrain (10^-9 strain)\n[dimensionless strain]")
    axes[0].set_title(
        "First-pass solid-Earth tide at SAFOD\n"
        "35.9746°N, 120.5516°W — June 17, 2026 experiment window\n"
        "This is tidal strain, not seismic δv/v"
    )
    axes[0].legend(frameon=False, ncol=3, loc="upper right")
    axes[0].grid(alpha=0.2)

    axes[1].plot(time_array, template, color="#8c2d04", lw=2.0)
    axes[1].axhline(0.0, color="0.5", lw=0.7)
    axes[1].set_ylabel("normalized template")
    axes[1].set_xlabel("UTC")
    axes[1].set_ylim(-1.1, 1.1)
    axes[1].grid(alpha=0.2)
    axes[1].text(
        0.01,
        0.08,
        "Template = (ε_tide − mean) / max|ε_tide − mean|\n"
        "Use this shape when fitting an amplitude A in δv/v",
        transform=axes[1].transAxes,
        fontsize=9,
        color="#8c2d04",
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )
    axes[1].xaxis.set_major_locator(mdates.HourLocator(interval=3))
    axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:%M"))
    for label in axes[1].get_xticklabels():
        label.set_rotation(30)
        label.set_ha("right")

    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    times = datetime_range(START_UTC, END_UTC, SAMPLE_MINUTES)
    result = model_tide(times)
    png_path = HERE / "safod_solid_earth_tide_2026-06-17.png"
    csv_path = HERE / "safod_solid_earth_tide_2026-06-17.csv"
    plot_result(times, result, png_path)
    save_csv(times, result, csv_path)

    total = result["total_strain"]
    print(f"Window: {times[0].isoformat()} -> {times[-1].isoformat()}")
    print(f"Samples: {len(times)} at {SAMPLE_MINUTES}-minute spacing")
    print(f"Love-number areal factor: {AREAL_LOVE_FACTOR:.4f}")
    print(f"Total scalar strain min: {total.min():.6e}")
    print(f"Total scalar strain max: {total.max():.6e}")
    print(f"Total scalar strain peak-to-peak: {np.ptp(total):.6e}")
    print(f"Total scalar strain half-range: {0.5 * np.ptp(total):.6e}")
    print(f"Saved {png_path}")
    print(f"Saved {csv_path}")


if __name__ == "__main__":
    main()
