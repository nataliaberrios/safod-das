"""Physically phased astronomical tide predictor and exploratory CCA for AWD Nano delays.

This is intentionally a low-precision, deterministic astronomical phase test.  It
uses the actual UTC time of each burst and a standard low-precision Sun/Moon
ephemeris to evaluate the degree-2 tide-generating potential at the SAFOD
latitude/longitude.  It is a scalar potential predictor, not a fault-normal
traction or opening/closing calculation: the local fault orientation and stress
boundary conditions are not supplied by the AWD dataset.
"""
from pathlib import Path
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
LAT_DEG = 35.98       # approximate SAFOD latitude; used only for phase geometry
LON_DEG = -120.55     # approximate SAFOD longitude, east-positive convention below
EARTH_RADIUS_M = 6_371_000.0
G = 6.67430e-11
M_SUN = 1.98847e30
M_MOON = 7.342e22
AU_M = 149_597_870_700.0
EARTH_RADIUS_KM = 6378.14
BLOCK_LEN = 3
N_PERM = 5000
SEED = 20260802


def z(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / v.std(ddof=1)


def julian_date(dt):
    """UTC Julian date for a pandas/DatetimeIndex-like array."""
    # datetime64[ns] is exact enough for this low-precision ephemeris.
    ns = pd.to_datetime(dt, utc=True).astype("int64").to_numpy()
    return 2440587.5 + ns / 86_400_000_000_000.0


def sun_radec(jd):
    """Low-precision apparent geocentric solar RA/Dec and distance (AU)."""
    d = np.asarray(jd, float) - 2451543.5
    g = np.deg2rad((357.529 + 0.98560028 * d) % 360.0)
    q = (280.459 + 0.98564736 * d) % 360.0
    lam = np.deg2rad((q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g)) % 360.0)
    r = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)
    eps = np.deg2rad(23.4393 - 3.563e-7 * d)
    ra = np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam))
    dec = np.arcsin(np.sin(eps) * np.sin(lam))
    return ra, dec, r * AU_M


def moon_radec(jd):
    """Low-precision geocentric lunar RA/Dec and distance.

    This is the compact orbital-element algorithm commonly used for quick
    sky-position calculations.  Its phase accuracy is adequate for a 24-hour
    screening test, but it is not a substitute for a precision tide package.
    """
    d = np.asarray(jd, float) - 2451543.5
    N = np.deg2rad((125.1228 - 0.0529538083 * d) % 360.0)
    inc = np.deg2rad(5.1454)
    w = np.deg2rad((318.0634 + 0.1643573223 * d) % 360.0)
    a = 60.2666
    ecc = 0.054900
    M = np.deg2rad((115.3654 + 13.0649929509 * d) % 360.0)
    # Two Newton steps are sufficient for this eccentricity.
    E = M.copy()
    for _ in range(3):
        E -= (E - ecc * np.sin(E) - M) / (1.0 - ecc * np.cos(E))
    xv = a * (np.cos(E) - ecc)
    yv = a * np.sqrt(1.0 - ecc * ecc) * np.sin(E)
    v = np.arctan2(yv, xv)
    r = np.sqrt(xv * xv + yv * yv)
    lon = v + w
    xh = r * (np.cos(N) * np.cos(lon) - np.sin(N) * np.sin(lon) * np.cos(inc))
    yh = r * (np.sin(N) * np.cos(lon) + np.cos(N) * np.sin(lon) * np.cos(inc))
    zh = r * (np.sin(lon) * np.sin(inc))
    # Ecliptic-to-equatorial rotation.
    eps = np.deg2rad(23.4393 - 3.563e-7 * d)
    xe = xh
    ye = yh * np.cos(eps) - zh * np.sin(eps)
    ze = yh * np.sin(eps) + zh * np.cos(eps)
    ra = np.arctan2(ye, xe)
    dec = np.arctan2(ze, np.sqrt(xe * xe + ye * ye))
    return ra, dec, r * EARTH_RADIUS_KM * 1000.0


def gmst_rad(jd):
    d = np.asarray(jd, float) - 2451545.0
    gmst_deg = (280.46061837 + 360.98564736629 * d + 0.000387933 * (d / 36525.0) ** 2) % 360.0
    return np.deg2rad(gmst_deg)


def degree2_potential(jd, body):
    """Degree-2 tide-generating potential U (m^2/s^2) at SAFOD."""
    if body == "sun":
        ra, dec, dist = sun_radec(jd)
        mass = M_SUN
    elif body == "moon":
        ra, dec, dist = moon_radec(jd)
        mass = M_MOON
    else:
        raise ValueError(body)
    lat = np.deg2rad(LAT_DEG)
    lst = gmst_rad(jd) + np.deg2rad(LON_DEG)
    hour_angle = lst - ra
    cos_zenith = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(hour_angle)
    p2 = 0.5 * (3.0 * cos_zenith * cos_zenith - 1.0)
    coeff = G * mass * EARTH_RADIUS_M**2 / dist**3
    return coeff * p2, np.rad2deg(ra), np.rad2deg(dec), dist


def load():
    m = pd.read_csv(HERE / "awd_manifest.csv", parse_dates=["utc_time"])
    b = m.groupby("burst_id", as_index=False).agg(utc_time=("utc_time", "median"))
    r = pd.read_csv(HERE / "nano_burst_repeatability_hierarchical.csv")
    d = b.merge(r, on="burst_id").sort_values("utc_time").reset_index(drop=True)
    d["hours"] = (d.utc_time - d.utc_time.iloc[0]).dt.total_seconds() / 3600.0
    d["cumulative_drops"] = d.n_drops.cumsum() - d.n_drops.iloc[0]
    d["delay_ms"] = 1000.0 * d.loo_burst_delay_s
    B = np.c_[np.ones(len(d)), z(d.cumulative_drops), z(d.burst_signal_rms)]
    d["delay_source_corrected_ms"] = d.delay_ms - B @ np.linalg.lstsq(B, d.delay_ms, rcond=None)[0]
    jd = julian_date(d.utc_time)
    d["moon_potential_m2s2"], d["moon_ra_deg"], d["moon_dec_deg"], d["moon_distance_m"] = degree2_potential(jd, "moon")
    d["sun_potential_m2s2"], d["sun_ra_deg"], d["sun_dec_deg"], d["sun_distance_m"] = degree2_potential(jd, "sun")
    d["total_potential_m2s2"] = d.moon_potential_m2s2 + d.sun_potential_m2s2
    # A radial degree-2 acceleration proxy is supplied only as a sensitivity
    # coordinate; without a fault orientation it must not be called traction.
    d["radial_tide_accel_nm_s2"] = 2.0 * d.total_potential_m2s2 / EARTH_RADIUS_M * 1e9
    return d


def matrices(d):
    # Body-specific scalar potentials retain the astronomical phase; no fitted
    # sine/cosine phase is introduced.  Standardization is implicit in CCA.
    X = d[["moon_potential_m2s2", "sun_potential_m2s2"]].to_numpy(float)
    Y = np.c_[z(d.delay_source_corrected_ms), z(d.burst_signal_rms)]
    return X, Y


def invsqrt(S):
    w, v = np.linalg.eigh(S)
    return (v * (1.0 / np.sqrt(np.maximum(w, 1e-8)))) @ v.T


def cca(X, Y):
    xm, ym = X.mean(0), Y.mean(0)
    xx, yy = X - xm, Y - ym
    n = len(X)
    sxx, syy, sxy = xx.T @ xx / (n - 1), yy.T @ yy / (n - 1), xx.T @ yy / (n - 1)
    A = invsqrt(sxx) @ sxy @ invsqrt(syy)
    u, s, vt = np.linalg.svd(A, full_matrices=False)
    wx, wy = invsqrt(sxx) @ u[:, 0], invsqrt(syy) @ vt.T[:, 0]
    sx, sy = xx @ wx, yy @ wy
    corr = float(np.corrcoef(sx, sy)[0, 1])
    if corr < 0:
        wy, sy, corr = -wy, -sy, -corr
    return {"corr": corr, "wx": wx, "wy": wy, "xscore": sx, "yscore": sy,
            "singular_values": s, "xmean": xm, "ymean": ym}


def block_order(rng, n):
    starts = np.arange(0, n, BLOCK_LEN)
    out = []
    for j in rng.permutation(len(starts)):
        out.extend((starts[j] + np.arange(BLOCK_LEN)) % n)
    return np.asarray(out[:n])


def permutation_null(X, Y):
    rng = np.random.default_rng(SEED)
    n = len(Y)
    row, block = np.empty(N_PERM), np.empty(N_PERM)
    for k in range(N_PERM):
        row[k] = cca(X, Y[rng.permutation(n)])["corr"]
        block[k] = cca(X, Y[block_order(rng, n)])["corr"]
    return row, block


def blocked_cv(X, Y, nfold=5):
    values = []
    for test in np.array_split(np.arange(len(X)), nfold):
        train = np.setdiff1d(np.arange(len(X)), test)
        fit = cca(X[train], Y[train])
        sx = (X[test] - fit["xmean"]) @ fit["wx"]
        sy = (Y[test] - fit["ymean"]) @ fit["wy"]
        values.append(float(abs(np.corrcoef(sx, sy)[0, 1])))
    return np.asarray(values)


def make_figure(d, fit, row, block, cv, out):
    t = d.hours.to_numpy(float)
    fig, ax = plt.subplots(2, 2, figsize=(11.8, 8.2), constrained_layout=True)
    a = ax[0, 0]
    a.plot(t, z(d.total_potential_m2s2), color="#c0392b", lw=2, label="astronomical total potential")
    a.plot(t, z(d.delay_source_corrected_ms), color="#34495e", lw=1.5, label="source-corrected delay")
    a.axhline(0, color=".3", lw=.7, ls=":"); a.grid(alpha=.2)
    a.set(xlabel="hours since first burst", ylabel="standardized value", title="(a) Physical predictor and response")
    a.legend(frameon=False, fontsize=8)
    a = ax[0, 1]
    a.scatter(fit["xscore"], fit["yscore"], c=t, cmap="plasma", s=40, edgecolor="white", lw=.35)
    a.set(xlabel="astronomical-potential canonical score", ylabel="delay/coupling canonical score",
          title=f"(b) Canonical variates (r={fit['corr']:.3f})"); a.grid(alpha=.2)
    a.text(.04, .94, "color = elapsed time", transform=a.transAxes, fontsize=8, va="top")
    a = ax[1, 0]
    labels = ["Moon U₂", "Sun U₂"]
    x = np.arange(2)
    a.bar(x - .18, fit["wx"], .36, color="#c0392b", label="astronomical predictors")
    a.bar(x + .18, [fit["wy"][0], fit["wy"][1]], .36, color="#34495e", label="responses (delay, RMS)")
    a.set_xticks(x, labels); a.axhline(0, color=".3", lw=.7); a.grid(axis="y", alpha=.2)
    a.set(ylabel="CCA weight", title="(c) Canonical weights"); a.legend(frameon=False, fontsize=8)
    a = ax[1, 1]
    a.hist(row, bins=30, color="#95a5a6", alpha=.45, density=True, label="row permutation")
    a.hist(block, bins=30, color="#5d6d7e", alpha=.40, density=True, label="3-burst block permutation")
    a.axvline(fit["corr"], color="#8e44ad", lw=2, label=f"observed r={fit['corr']:.3f}")
    a.set(xlabel="first canonical correlation", ylabel="density", title="(d) Serial-dependence null tests")
    a.grid(alpha=.2); a.legend(frameon=False, fontsize=7)
    fig.suptitle("Physically phased astronomical tide potential versus source-corrected AWD delay", fontsize=14, fontweight="bold")
    fig.savefig(out, dpi=260); fig.savefig(out.with_suffix(".pdf")); plt.close(fig)


def main():
    d = load(); X, Y = matrices(d); fit = cca(X, Y); row, block = permutation_null(X, Y); cv = blocked_cv(X, Y)
    p_row = float((1 + (row >= fit["corr"]).sum()) / (len(row) + 1))
    p_block = float((1 + (block >= fit["corr"]).sum()) / (len(block) + 1))
    out_csv = HERE / "nano_physical_tide_cca.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["quantity", "value", "units", "interpretation"])
        w.writerows([
            ("n_bursts", len(d), "bursts", "burst-level observations"),
            ("canonical_correlation_1", fit["corr"], "r", "full-sample CCA"),
            ("row_permutation_p", p_row, "probability", "rows independently permuted"),
            ("block_permutation_p", p_block, "probability", "three-burst blocks permuted"),
            ("blocked_cv_corr_mean", float(np.mean(cv)), "r", "five contiguous test folds"),
            ("blocked_cv_corr_median", float(np.median(cv)), "r", "five contiguous test folds"),
            ("blocked_cv_corr_values", ";".join(f"{v:.3f}" for v in cv), "r", "fold-by-fold values"),
            ("x_set", "Moon and Sun degree-2 tidal potential", "predictor", "physically phased from UTC ephemeris"),
            ("y_set", "source-corrected delay and burst signal RMS", "response", "delay corrected with history/coupling baseline"),
            ("safod_latitude", LAT_DEG, "deg", "approximate phase geometry"),
            ("safod_longitude", LON_DEG, "deg", "approximate phase geometry"),
            ("status", "exploratory physical-phase CCA; not traction detection", "status", "fault orientation and compaction measurement unavailable"),
        ])
    txt = HERE / "nano_physical_tide_cca.txt"
    with txt.open("w") as f:
        f.write("SAFOD AWD Nano physically phased astronomical tide CCA\n")
        f.write("Status: exploratory timing/association test; not a fault-traction, opening/closing, or permeability detection.\n\n")
        f.write("Predictors: degree-2 Moon and Sun tide-generating potentials evaluated at the actual UTC burst times using a low-precision ephemeris.\n")
        f.write(f"Approximate geometry: latitude={LAT_DEG:.2f} deg, longitude={LON_DEG:.2f} deg.\n")
        f.write("Response: Nano burst delay after cumulative-drop and burst-RMS source-history correction, plus burst signal RMS.\n")
        f.write("No independent source-ground compaction/settlement time series was found locally; cumulative drops remain a mechanical-history proxy.\n\n")
        f.write(f"Full-sample first canonical correlation: {fit['corr']:.3f}.\n")
        f.write(f"Independent row-permutation p={p_row:.4f}; three-burst block-permutation p={p_block:.4f}.\n")
        f.write(f"Five contiguous-fold test correlations: {', '.join(f'{v:.3f}' for v in cv)}; mean={np.mean(cv):.3f}, median={np.median(cv):.3f}.\n")
        f.write("Interpretation: a physically phased association is screened, but significance must be judged by the serial block null and held-out stability. The scalar potential does not provide a fault-normal sign; an opening/closing interpretation requires surveyed fault orientation, stress convention, and ideally an independent compaction measurement.\n")
    with (HERE / "nano_physical_tide_cca.json").open("w") as f:
        json.dump({"status": "exploratory_physical_phase_cca_not_tidal_traction_detection", "n_bursts": len(d),
                   "canonical_correlation": fit["corr"], "row_permutation_p": p_row, "block_permutation_p": p_block,
                   "blocked_cv_correlations": cv.tolist(), "x_weights": fit["wx"].tolist(), "y_weights": fit["wy"].tolist(),
                   "x_set": "Moon and Sun degree-2 tide-generating potential from UTC ephemeris", "y_set": "source-corrected delay and burst signal RMS",
                   "geometry": {"latitude_deg": LAT_DEG, "longitude_deg": LON_DEG},
                   "limitations": ["No independent ground-level compaction measurement was found.", "Low-precision ephemeris is suitable for phase screening, not precision gravimetry.", "Scalar potential is not fault-normal traction.", "Fault strike/dip and stress convention are unavailable.", "Serial block null and blocked CV control the interpretation."]}, f, indent=2)
    d.to_csv(HERE / "nano_physical_tide_predictors.csv", index=False)
    np.savez(HERE / "nano_physical_tide_cca.npz", hours=d.hours.to_numpy(float), delay_ms=d.delay_ms.to_numpy(float),
             source_corrected_delay_ms=d.delay_source_corrected_ms.to_numpy(float), moon_potential_m2s2=d.moon_potential_m2s2.to_numpy(float),
             sun_potential_m2s2=d.sun_potential_m2s2.to_numpy(float), total_potential_m2s2=d.total_potential_m2s2.to_numpy(float),
             radial_tide_accel_nm_s2=d.radial_tide_accel_nm_s2.to_numpy(float), x_score=fit["xscore"], y_score=fit["yscore"],
             x_weights=fit["wx"], y_weights=fit["wy"], row_null=row, block_null=block, blocked_cv=cv)
    make_figure(d, fit, row, block, cv, HERE / "nano_physical_tide_cca.png")
    print(txt.read_text())


if __name__ == "__main__":
    main()
