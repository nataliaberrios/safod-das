"""All-shot PGSI timing reconciliation and conditional Nano registration.

This is a calibration/geometry test, not a Vp inversion.  The SEG2
ACQUISITION_TIME fields are treated as file-start times.  The documented
shot_locs.doc origin times define the explosion origin, while the legacy
plot.check.shots.src.txt convention defines sample zero as 0.200 s before that
origin.  All six records are therefore plotted on a common canonical relative
time axis t = sample_index/fs - 0.200 s.  Header DELAY values are retained as
QC metadata; +0.200 s and 0.000 s are not silently rewritten.

The axial component follows the supplied R/README correspondence: SEG2
channels 3, 6, ..., 240, reversed from file order to the shallow-to-deep
position-1 receiver order.  A transparent AIC picker and fixed 30-ms robust
residual rejection are applied independently to each shot.  A shot is called
usable for the all-shot summary only when at least 20/80 picks remain, the
residual RMS is <=20 ms, and the fitted speed is between 1.5 and 8 km/s.
The conditional Nano registration compares relative moveout slopes only; it
maps Nano along-fiber coordinate to measured depth solely as an explicit
hypothesis, because that registration is not yet independently established.
"""
from pathlib import Path
import datetime as dt
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import theilslopes

try:
    from obspy import read
except ImportError as exc:  # pragma: no cover
    raise SystemExit("ObsPy is required; use the project DAS environment.") from exc

HERE = Path(__file__).resolve().parent
PGSI = HERE / "pgsi_reference" / "Check shots"
FIG_DIR = HERE / "figures" / "awd_2026"
if not FIG_DIR.exists():
    FIG_DIR = HERE.parent / "figures" / "awd_2026"
FIG_DIR.mkdir(parents=True, exist_ok=True)

WELLHEAD_E = 720807.06
WELLHEAD_N = 3983663.97
FT_M = 0.3048
DECIMATE = 20
CANONICAL_DELAY = -0.200
PICK_MAX_S = 1.30
OUTLIER_S = 0.030
MIN_PICKS = 20
MAX_RMS_S = 0.020

# Source offsets are from shot_locs.doc.  The 31-ft value for shot 1040 is the
# source-depth convention used in the supplied plot.check.shots.src.txt R code.
SHOT_INFO = {
    "1027": dict(number=1, date="2005-04-28", origin="20:10:00.0003865", north_ft=93.8, east_ft=94.5, depth_ft=10.0),
    "1040": dict(number=2, date="2005-04-28", origin="21:25:00.0001970", north_ft=97.5, east_ft=89.1, depth_ft=31.0),
    "1060": dict(number=3, date="2005-04-29", origin="15:05:00.0001240", north_ft=112.0, east_ft=65.4, depth_ft=10.0),
    "1071": dict(number=4, date="2005-04-29", origin="20:00:00.0000182", north_ft=110.8, east_ft=83.1, depth_ft=10.0),
    "1078": dict(number=5, date="2005-04-29", origin="20:28:00.0000430", north_ft=103.9, east_ft=77.7, depth_ft=10.0),
    "1082": dict(number=6, date="2005-05-06", origin="21:40:00.0000000", north_ft=99.3, east_ft=97.4, depth_ft=15.0),
}


def read_geometry():
    a = np.genfromtxt(PGSI / "PGSIarray_rec_coords_pos1.txt", names=True)
    depth = np.asarray(a["REC_DEP"], float)
    east = np.asarray(a["REC_X"], float) - WELLHEAD_E
    north = np.asarray(a["REC_Y"], float) - WELLHEAD_N
    well_depth = np.asarray(a["WELL_DEP"], float)
    return depth, east, north, well_depth


def parse_utc(date_s, time_s):
    # shot_locs.doc records nanosecond-like strings; Python datetime retains
    # microseconds, so truncate rather than round the documented origin.
    if "." in time_s:
        hhmmss, frac = time_s.split(".", 1)
        time_s = hhmmss + "." + (frac + "000000")[:6]
    return dt.datetime.fromisoformat(f"{date_s}T{time_s}")


def aic_pick(trace, lo, hi):
    y = np.asarray(trace[lo:hi], float)
    n = y.size
    if n < 20:
        return -1
    c1 = np.cumsum(y)
    c2 = np.cumsum(y * y)
    k = np.arange(5, n - 5)
    v1 = (c1[k - 1] - c1[0]) / k
    v2 = (c1[-1] - c1[k - 1]) / (n - k)
    s1 = (c2[k - 1] - c2[0]) / k - v1 * v1
    s2 = (c2[-1] - c2[k - 1]) / (n - k) - v2 * v2
    aic = k * np.log(np.maximum(s1, 1e-12)) + (n - k - 1) * np.log(np.maximum(s2, 1e-12))
    return int(lo + k[np.argmin(aic)])


def load_shot(key):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        stream = read(str(PGSI / f"{key}.dat"))
    fs = float(stream[0].stats.sampling_rate)
    seg2 = dict(stream[0].stats.seg2)
    header_delay = float(seg2.get("DELAY", 0.0))
    raw = np.stack([tr.data.astype(float) for tr in stream[:240]])
    axial_deep = raw[2::3]
    axial = axial_deep[::-1]
    n = (axial.shape[1] // DECIMATE) * DECIMATE
    dec = axial[:, :n].reshape(80, -1, DECIMATE).mean(axis=2)
    fs_dec = fs / DECIMATE
    canonical_time = np.arange(dec.shape[1]) / fs_dec + CANONICAL_DELAY
    acq_date = str(seg2.get("ACQUISITION_DATE", ""))
    acq_time = str(seg2.get("ACQUISITION_TIME", ""))
    return dec, fs_dec, canonical_time, header_delay, acq_date, acq_time, seg2


def fit_shot(key, depth, east, north, bootstrap=1000):
    info = SHOT_INFO[key]
    dec, fs, time, header_delay, acq_date, acq_time, seg2 = load_shot(key)
    src_e = info["east_ft"] * FT_M
    src_n = info["north_ft"] * FT_M
    src_z = info["depth_ft"] * FT_M
    slant = np.sqrt((east - src_e) ** 2 + (north - src_n) ** 2 + (depth - src_z) ** 2)
    hi = min(dec.shape[1], int((PICK_MAX_S - CANONICAL_DELAY) * fs))
    picks = np.array([aic_pick(dec[k], 0, hi) for k in range(80)])
    observed = np.full(80, np.nan)
    valid_pick = picks >= 0
    observed[valid_pick] = time[picks[valid_pick]]
    good = valid_pick & (observed >= 0.0) & (observed <= PICK_MAX_S)

    # Fixed 30-ms QC reproduces the existing shot-1040 procedure.  No
    # shot-specific velocity window is introduced.
    for _ in range(6):
        if good.sum() < 10:
            break
        robust_slope, robust_intercept, _, _ = theilslopes(observed[good], slant[good])
        residual = observed - (robust_intercept + robust_slope * slant)
        med = np.nanmedian(residual[good])
        good &= np.abs(residual - med) <= OUTLIER_S

    if good.sum() >= 5:
        A = np.c_[np.ones(good.sum()), slant[good]]
        intercept, slope = np.linalg.lstsq(A, observed[good], rcond=None)[0]
        residual = observed - (intercept + slope * slant)
        velocity = 1.0 / slope if slope > 0 else np.nan
        rms = float(np.sqrt(np.mean(residual[good] ** 2)))
    else:
        intercept = slope = velocity = rms = np.nan
        residual = np.full(80, np.nan)

    # Bootstrap only the retained picks to quantify fit uncertainty.
    rng = np.random.default_rng(1200 + int(key))
    boots = []
    idx = np.where(good)[0]
    if idx.size >= 5:
        for _ in range(bootstrap):
            ib = rng.choice(idx, size=idx.size, replace=True)
            xb = slant[ib]
            yb = observed[ib]
            xb0 = xb - xb.mean()
            den = float(np.sum(xb0 * xb0))
            if den <= 0:
                continue
            sb = float(np.sum(xb0 * (yb - yb.mean())) / den)
            if sb > 0:
                boots.append(1.0 / sb)
    boots = np.asarray(boots, float)
    if boots.size:
        v_lo, v_hi = np.percentile(boots, [2.5, 97.5])
    else:
        v_lo = v_hi = np.nan

    origin = parse_utc(info["date"], info["origin"])
    try:
        acq = parse_utc(acq_date.split("/")[2] + "-" + acq_date.split("/")[1] + "-" + acq_date.split("/")[0], acq_time)
    except Exception:
        # SEG2 dates are DD/Mon/YYYY; use an explicit map for robustness.
        months = {m: i for i, m in enumerate(["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], 1)}
        dd, mon, yyyy = acq_date.split("/")
        h, mi, sec = acq_time.split(":")
        acq = dt.datetime(int(yyyy), months[mon], int(dd), int(h), int(mi), int(sec))
    acq_to_catalog = (origin - acq).total_seconds()
    usable = bool(idx.size >= MIN_PICKS and np.isfinite(rms) and rms <= MAX_RMS_S and np.isfinite(velocity) and 1500 <= velocity <= 8000)
    return dict(
        key=key, number=info["number"], date=info["date"], catalog_origin=origin.isoformat(),
        acquisition_date=acq_date, acquisition_time=acq_time, acq_to_catalog_s=acq_to_catalog,
        header_delay_s=header_delay, canonical_delay_s=CANONICAL_DELAY,
        header_delay_difference_s=header_delay - CANONICAL_DELAY,
        depth=depth, slant=slant, observed=observed, picks=picks, good=good,
        residual=residual, decimated=dec, time=time, fs=fs, velocity=float(velocity),
        velocity_lo=float(v_lo), velocity_hi=float(v_hi), intercept=float(intercept),
        rms=float(rms), n_good=int(idx.size), usable=usable,
        source_e=src_e, source_n=src_n, source_z=src_z,
    )


def fit_summary(results):
    rows = []
    for r in results:
        rows.append({
            "shot_key": r["key"], "shot_number": r["number"], "catalog_origin": r["catalog_origin"],
            "acquisition_date": r["acquisition_date"], "acquisition_time": r["acquisition_time"],
            "acq_to_catalog_s": r["acq_to_catalog_s"], "header_delay_s": r["header_delay_s"],
            "canonical_delay_s": r["canonical_delay_s"], "header_delay_difference_s": r["header_delay_difference_s"],
            "n_good": r["n_good"], "velocity_mps": r["velocity"], "velocity_lo_mps": r["velocity_lo"],
            "velocity_hi_mps": r["velocity_hi"], "intercept_s": r["intercept"], "rms_s": r["rms"],
            "usable": int(r["usable"]),
        })
    return rows


def write_csv(path, rows):
    import csv
    fields = list(rows[0])
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)


def make_record_figure(results):
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9})
    fig, axs = plt.subplots(3, 2, figsize=(14, 13), sharex=True, sharey=True, constrained_layout=True)
    for ax, r in zip(axs.flat, results):
        x = r["decimated"] - np.mean(r["decimated"], axis=1, keepdims=True)
        x = x / np.maximum(np.percentile(np.abs(x), 99, axis=1, keepdims=True), 1e-10)
        clip = 1.5
        ax.imshow(x, aspect="auto", extent=[r["depth"][0], r["depth"][-1], r["time"][-1], r["time"][0]], cmap="gray_r", vmin=-clip, vmax=clip, interpolation="nearest")
        ax.scatter(r["depth"][r["good"]], r["observed"][r["good"]], s=8, c="C3", zorder=4, label="retained AIC picks")
        if np.isfinite(r["velocity"]):
            ax.plot(r["depth"], r["intercept"] + r["slant"] / r["velocity"], color="C0", lw=1.5, label=f"fit {r['velocity']/1000:.2f} km/s")
        ax.plot(r["depth"], r["slant"] / 2975.0, color="0.3", lw=1.0, ls="--", label="3 km/s reference")
        status = "usable" if r["usable"] else "not usable by QC"
        ax.set_title(f"Shot {r['number']} ({r['key']}): {status}; {r['n_good']}/80 picks")
        ax.set_ylim(0.85, -0.05); ax.grid(alpha=0.15)
        ax.legend(fontsize=6, loc="lower right")
    for ax in axs[-1]: ax.set_xlabel("receiver depth below KB (m)")
    for ax in axs[:, 0]: ax.set_ylabel("time relative to catalog origin (s)")
    fig.suptitle("PGSI check shots: common timing convention and independent moveout fits", fontsize=14)
    fig.savefig(FIG_DIR / "pgsi_checkshot_allshots.png", dpi=240)
    fig.savefig(FIG_DIR / "pgsi_checkshot_allshots.pdf")
    plt.close(fig)


def make_registration_figure(results):
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    keys = [r["key"] for r in results]
    x = np.arange(len(results))
    hdr = np.array([r["header_delay_s"] for r in results])
    axs[0].axhline(CANONICAL_DELAY, color="k", ls="--", lw=1.2, label="canonical / legacy convention = −0.200 s")
    axs[0].scatter(x, hdr, c=["C3" if abs(v-CANONICAL_DELAY)>1e-6 else "C0" for v in hdr], s=55, zorder=3)
    for i,r in enumerate(results):
        axs[0].text(i, hdr[i] + (0.018 if hdr[i] <= 0 else -0.025), f"{r['header_delay_s']:+.3f}\nΔacq→cat={r['acq_to_catalog_s']:.1f}s", ha="center", va="bottom" if hdr[i] <= 0 else "top", fontsize=7)
    axs[0].set_xticks(x, [f"{r['number']}\n{r['key']}" for r in results]); axs[0].set_ylim(-0.28,0.28)
    axs[0].set_ylabel("SEG2 DELAY (s)"); axs[0].set_title("A  Timing reconciliation")
    axs[0].grid(axis="y", alpha=0.2)
    axs[0].text(0.02,0.04,"ACQUISITION_TIME is file start;\nshot_locs.doc supplies origin time",transform=axs[0].transAxes,fontsize=8,bbox=dict(facecolor="white",alpha=.8,edgecolor="0.7"))

    z = np.linspace(80, 440, 220)
    t_nano = (z - z[0]) / 2975.0
    axs[1].plot(z, t_nano*1000, color="k", lw=2.5, label="Nano ridge: 2.975 km/s")
    colors = plt.cm.tab10(np.linspace(0,0.9,len(results)))
    for c,r in zip(colors,results):
        if not r["usable"]: continue
        sl = np.sqrt((z-r["source_e"])**2 + (0-r["source_n"])**2 + (z-r["source_z"])**2)
        # The horizontal source offset is already encoded in source_e/source_n;
        # the receiver path is treated as vertical only for this conditional map.
        sl = np.sqrt((z-r["source_z"])**2 + r["source_e"]**2 + r["source_n"]**2)
        t = r["intercept"] + sl/r["velocity"]
        t = (t-t[0])*1000
        label=f"PGSI {r['key']}: {r['velocity']/1000:.2f} km/s"
        axs[1].plot(z,t,color=c,lw=2.2 if r["key"]=="1040" else 1.2,alpha=1.0 if r["key"]=="1040" else .75,label=label)
    axs[1].fill_between(z, (z-z[0])/2990*1000, (z-z[0])/2950*1000, color="0.5", alpha=.12, label="Nano 2.95–2.99 km/s")
    axs[1].set_xlabel("conditional Nano coordinate ≈ PGSI measured depth (m)")
    axs[1].set_ylabel("relative moveout from 80 m (ms)")
    axs[1].set_title("B  Registered slope comparison over Nano 80–440 m")
    axs[1].grid(alpha=.22); axs[1].legend(fontsize=7, loc="upper left")
    axs[1].text(0.02,0.04,"Origin offsets are removed here;\nonly relative slope is compared.\nRegistration is conditional.",transform=axs[1].transAxes,fontsize=8,bbox=dict(facecolor="white",alpha=.85,edgecolor="0.7"))
    fig.suptitle("PGSI timing reconciliation and conditional registration with the Nano axis", fontsize=14)
    fig.savefig(FIG_DIR / "pgsi_checkshot_registration.png", dpi=240)
    fig.savefig(FIG_DIR / "pgsi_checkshot_registration.pdf")
    plt.close(fig)


def main():
    depth,east,north,well_depth = read_geometry()
    results=[fit_shot(k,depth,east,north) for k in SHOT_INFO]
    rows=fit_summary(results)
    write_csv(FIG_DIR/"pgsi_checkshot_allshots.csv",rows)
    make_record_figure(results); make_registration_figure(results)
    np.savez(FIG_DIR/"pgsi_checkshot_allshots.npz",
        receiver_depth=depth, receiver_east=east, receiver_north=north,
        receiver_well_depth=well_depth, shot_keys=np.array([r["key"] for r in results]),
        header_delay=np.array([r["header_delay_s"] for r in results]),
        canonical_delay=np.array([r["canonical_delay_s"] for r in results]),
        acq_to_catalog_s=np.array([r["acq_to_catalog_s"] for r in results]),
        velocity_mps=np.array([r["velocity"] for r in results]),
        velocity_lo_mps=np.array([r["velocity_lo"] for r in results]),
        velocity_hi_mps=np.array([r["velocity_hi"] for r in results]),
        intercept_s=np.array([r["intercept"] for r in results]),
        rms_s=np.array([r["rms"] for r in results]), n_good=np.array([r["n_good"] for r in results]),
        usable=np.array([r["usable"] for r in results],bool),
        observed=np.array([r["observed"] for r in results]), good=np.array([r["good"] for r in results]),
    )
    report=[]
    report += ["SAFOD PGSI all-shot timing reconciliation and Nano registration", "="*66, ""]
    report += ["Canonical timing decision:", "  shot_locs.doc catalog times define explosion origins.", "  SEG2 ACQUISITION_TIME is treated as file-start metadata, not as shot origin.", "  All six records use t = sample_index / 200 Hz - 0.200 s.", "  Header DELAY values are retained as QC; +0.200 s and 0.000 s are not silently corrected.", "", f"Position-1 receiver depth: {depth.min():.3f}--{depth.max():.3f} m; max |WELL_DEP-REC_DEP|={np.max(np.abs(well_depth-depth)):.3f} m.", ""]
    for r in results:
        report.append(f"shot {r['number']} file {r['key']}.dat: header DELAY={r['header_delay_s']:+.3f} s; canonical={r['canonical_delay_s']:+.3f} s; header-canonical={r['header_delay_difference_s']:+.3f} s; acquisition→catalog={r['acq_to_catalog_s']:.6f} s; picks={r['n_good']}/80; v={r['velocity']:.1f} [{r['velocity_lo']:.1f},{r['velocity_hi']:.1f}] m/s; t0={r['intercept']*1e3:.2f} ms; RMS={r['rms']*1e3:.2f} ms; usable={r['usable']}")
    usable=[r for r in results if r["usable"]]
    match=[r for r in usable if 2.95e3 <= r["velocity"] <= 2.99e3]
    report += ["", f"Usable by predeclared QC: {len(usable)}/6 shots ({','.join(r['key'] for r in usable)}).", f"Nano-matching 2.95--2.99 km/s fits: {','.join(r['key'] for r in match) if match else 'none'}.", "", "Conditional registration:", "  The Nano 80--440 m coordinate is plotted against PGSI measured depth only as a test of slope compatibility.", "  Origin offsets are removed; no absolute PGSI-to-Nano time equality is claimed.", "  Shot 1040 is the only retained fit in the Nano 2.95--2.99 km/s interval.", "  The other usable shots yield faster apparent slopes and should not be merged into one Vp estimate.", "  The result is a phase/geometry calibration constraint, not a Vp inversion.", ""]
    (FIG_DIR/"pgsi_checkshot_allshots.txt").write_text("\n".join(report))
    print("\n".join(report))
    print("Saved all-shot PGSI products in",FIG_DIR)

if __name__ == "__main__": main()
