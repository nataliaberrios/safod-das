"""PGSI check-shot geometry and timing audit for the SAFOD AWD project.

This is deliberately an audit, not a Vp inversion. The 2005 PGSI explosions
are independent of the 2026 AWD experiment, but they provide surveyed receiver
coordinates, documented source offsets, and SEG2 timing fields. Those pieces
make the direct-P moveout/depth-registration hypothesis testable before it is
used to interpret the Nano DAS mode.

The script uses the legacy R-code channel correspondence: the axial geophone
at each level is SEG2 channel 3, 6, ..., 240 (deep-to-shallow in the file),
which is reversed to the top-to-bottom coordinate-file order. It reproduces
the supplied R workflow at 200 samples/s and applies a transparent AIC picker
with a 30-ms robust residual rejection. The resulting picks are shown as a
diagnostic only; the figure and report retain the header-delay discrepancies
so that no timing correction is hidden.
"""

from pathlib import Path
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import theilslopes

try:
    from obspy import read
except ImportError as exc:  # pragma: no cover - environment-specific
    raise SystemExit(
        "ObsPy is required. Run with the SAFOD das environment, e.g. "
        "/home/users/nberrios/miniconda3/envs/das/bin/python."
    ) from exc


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
PICK_MAX_S = 1.30
OUTLIER_S = 0.030

# Source locations extracted from shot_locs.doc. The source depth is the
# legacy R-code value for shot 1040 (31 ft below the KB); the other entries use
# the documented hole depths and are used for the geometry inventory only.
SHOT_INFO = {
    "1027": dict(number=1, north_ft=93.8, east_ft=94.5, depth_ft=10.0),
    "1040": dict(number=2, north_ft=97.5, east_ft=89.1, depth_ft=31.0),
    "1060": dict(number=3, north_ft=112.0, east_ft=65.4, depth_ft=10.0),
    "1071": dict(number=4, north_ft=110.8, east_ft=83.1, depth_ft=10.0),
    "1078": dict(number=5, north_ft=103.9, east_ft=77.7, depth_ft=10.0),
    "1082": dict(number=6, north_ft=99.3, east_ft=97.4, depth_ft=15.0),
}


def read_geometry():
    """Return top-to-bottom position-1 receiver geometry in local metres."""
    path = PGSI / "PGSIarray_rec_coords_pos1.txt"
    a = np.genfromtxt(path, names=True)
    depth = np.asarray(a["REC_DEP"], float)
    east = np.asarray(a["REC_X"], float) - WELLHEAD_E
    north = np.asarray(a["REC_Y"], float) - WELLHEAD_N
    return depth, east, north


def aic_pick(trace, lo, hi):
    """Return an AIC change-point index in a decimated trace."""
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
    aic = k * np.log(np.maximum(s1, 1e-12)) + (n - k - 1) * np.log(
        np.maximum(s2, 1e-12)
    )
    return int(lo + k[np.argmin(aic)])


def load_shot(key):
    """Load one SEG2 file and return axial records in shallow-to-deep order."""
    stream = read(str(PGSI / f"{key}.dat"))
    fs = float(stream[0].stats.sampling_rate)
    delay = float(stream[0].stats.seg2.get("DELAY", 0.0))
    # 240 seismic traces occupy channels 1--240; channel 252 is a timing
    # trace. The R orientation file documents the axial pins as 240,237,...,3.
    raw = np.stack([tr.data.astype(float) for tr in stream[:240]])
    axial_deep = raw[2::3]
    axial = axial_deep[::-1]
    n = (axial.shape[1] // DECIMATE) * DECIMATE
    dec = axial[:, :n].reshape(80, -1, DECIMATE).mean(axis=2)
    fs_dec = fs / DECIMATE
    time = np.arange(dec.shape[1]) / fs_dec + delay
    return dec, fs_dec, time, delay, stream


def pick_1040(depth, east, north):
    """Pick the documented shot-1040 direct-looking onset and robustly fit it."""
    dec, fs, time, delay, stream = load_shot("1040")
    info = SHOT_INFO["1040"]
    src_e = info["east_ft"] * FT_M
    src_n = info["north_ft"] * FT_M
    src_z = info["depth_ft"] * FT_M
    slant = np.sqrt((east - src_e) ** 2 + (north - src_n) ** 2 + (depth - src_z) ** 2)
    # Reproduce the supplied R workflow: AIC is evaluated from file start
    # (nominally -0.200 s for shot 1040), then pre-origin picks are rejected
    # explicitly below. Starting at t=0 biases the picker toward later phases.
    lo = 0
    hi = min(dec.shape[1], int((PICK_MAX_S - delay) * fs))
    picks = np.array([aic_pick(dec[k], lo, hi) for k in range(dec.shape[0])])
    observed = np.full(80, np.nan)
    ok = picks >= 0
    observed[ok] = time[picks[ok]]
    good = ok & (observed >= 0.0) & (observed <= PICK_MAX_S)

    # The AIC picks can land on strong later arrivals. Fit a robust line and
    # reject only points more than 30 ms from it; this is a QC flag, not a
    # hand-tuned acceptance of a preferred velocity.
    for _ in range(4):
        if good.sum() < 10:
            break
        slope, intercept, _, _ = theilslopes(observed[good], slant[good])
        residual = observed - (intercept + slope * slant)
        med = np.nanmedian(residual[good])
        good &= np.abs(residual - med) <= OUTLIER_S
    if good.sum() >= 5:
        A = np.c_[np.ones(good.sum()), slant[good]]
        intercept, slope = np.linalg.lstsq(A, observed[good], rcond=None)[0]
        residual = observed - (intercept + slope * slant)
        velocity = 1.0 / slope
        rms = float(np.sqrt(np.mean(residual[good] ** 2)))
    else:
        intercept = slope = velocity = rms = np.nan
        residual = np.full(80, np.nan)

    return dict(
        decimated=dec,
        fs=fs,
        time=time,
        delay=delay,
        slant=slant,
        observed=observed,
        good=good,
        residual=residual,
        intercept=float(intercept),
        velocity=float(velocity),
        rms=rms,
        stream=stream,
    )


def main():
    depth, east, north = read_geometry()
    pick = pick_1040(depth, east, north)

    # Header inventory: this is part of the result. The supplied R note
    # describes -0.200 s, but two files carry different SEG2 DELAY fields.
    headers = []
    for key, info in SHOT_INFO.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            st = read(str(PGSI / f"{key}.dat"), headonly=True)
        tr = st[0]
        delay = float(tr.stats.seg2.get("DELAY", 0.0))
        headers.append((key, info["number"], delay, tr.stats.seg2.get("ACQUISITION_TIME", "")))

    # Figure layout: geometry; timing inventory; corrected 1040 gather with
    # candidate direct-P curves; and observed-pick residuals.
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9})
    fig = plt.figure(figsize=(13.5, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.82, 1.18])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    # A: plan-view source/receiver geometry.
    ax0.plot(0, 0, "k+", ms=12, mew=2, label="SAFOD wellhead")
    ax0.plot(east, north, color="0.35", lw=1.2, label="PGSI receiver path")
    for key, info in SHOT_INFO.items():
        se = info["east_ft"] * FT_M
        sn = info["north_ft"] * FT_M
        ax0.plot(se, sn, "o", ms=5, label=f"shot {info['number']} ({key})")
        ax0.text(se + 1.5, sn + 1.5, str(info["number"]), fontsize=8)
    ax0.set(xlabel="east of wellhead (m)", ylabel="north of wellhead (m)",
            title="A  Surveyed source/receiver geometry")
    ax0.set_aspect("equal", adjustable="box")
    ax0.grid(alpha=0.25)
    ax0.legend(fontsize=6.5, ncol=2, loc="upper left")

    # B: explicit SEG2 delay inventory.
    keys = [h[0] for h in headers]
    delays = np.array([h[2] for h in headers])
    cols = ["C0" if abs(x + 0.2) < 1e-6 else "C3" for x in delays]
    ax1.axhline(-0.2, color="k", ls="--", lw=1, label="documented/R note: −0.200 s")
    ax1.scatter(np.arange(6), delays, c=cols, s=45, zorder=3)
    for i, (key, number, delay, acq) in enumerate(headers):
        ax1.text(i, delay + (0.012 if delay <= 0 else -0.018), f"{key}\n{delay:+.3f}",
                 ha="center", va="bottom" if delay <= 0 else "top", fontsize=7)
    ax1.set_xticks(np.arange(6), [f"shot {h[1]}" for h in headers])
    ax1.set_ylim(-0.27, 0.27)
    ax1.set_ylabel("SEG2 DELAY (s)")
    ax1.set_title("B  Timing metadata must be audited before joint picking")
    ax1.grid(axis="y", alpha=0.25)
    ax1.legend(fontsize=7, loc="upper right")

    # C: gather. Normalize each trace only for visibility; keep the actual
    # AIC picks and curves in the physical time coordinate.
    gather = pick["decimated"]
    norm = gather / np.maximum(np.percentile(np.abs(gather), 99, axis=1, keepdims=True), 1e-12)
    clip = np.percentile(np.abs(norm), 98)
    extent = [depth[0], depth[-1], pick["time"][-1], pick["time"][0]]
    ax2.imshow(norm, aspect="auto", cmap="gray_r", extent=extent,
               vmin=-clip, vmax=clip, interpolation="nearest")
    ax2.scatter(depth[pick["good"]], pick["observed"][pick["good"]], s=8,
                c="C3", label="AIC picks retained by 30-ms QC", zorder=4)
    for v, color in [(2500, "C0"), (3000, "C1"), (4000, "C2"), (5000, "C4")]:
        ax2.plot(depth, pick["slant"] / v, color=color, lw=1.1, label=f"{v/1000:g} km/s")
    ax2.axhline(0, color="k", lw=0.7, alpha=0.5)
    ax2.set(ylim=(0.70, -0.02), xlabel="receiver depth below KB (m)",
            ylabel="time relative to shot origin (s)",
            title="C  Shot 1040: geometry-corrected moveout")
    ax2.grid(alpha=0.18)
    ax2.legend(fontsize=7, ncol=2, loc="lower right")

    # D: residuals and fit, separating origin-time offset from slope.
    good = pick["good"]
    ax3.axhline(0, color="k", lw=0.8)
    ax3.scatter(pick["slant"][good], pick["residual"][good] * 1e3, s=14,
                c="C3", label="retained picks")
    ax3.scatter(pick["slant"][~good], pick["residual"][~good] * 1e3, s=12,
                facecolors="none", edgecolors="0.5", label="rejected picks")
    ax3.set(xlabel="source–receiver slant distance (m)", ylabel="fit residual (ms)",
            title=(f"D  Robust line: v={pick['velocity']/1000:.2f} km/s, "
                   f"t₀={pick['intercept']*1e3:.1f} ms; N={good.sum()}/80"))
    ax3.grid(alpha=0.25)
    ax3.legend(fontsize=7, loc="upper right")

    fig.suptitle("SAFOD PGSI check shots: geometry and timing decision test", fontsize=13)
    for i, ax in enumerate((ax0, ax1, ax2, ax3)):
        ax.text(-0.08, 1.05, "ABCD"[i], transform=ax.transAxes, fontweight="bold", fontsize=11)

    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"pgsi_checkshot_geometry_timing.{ext}", dpi=220 if ext == "png" else None)

    a = np.genfromtxt(PGSI / "PGSIarray_rec_coords_pos1.txt", names=True)
    well_depth = np.asarray(a["WELL_DEP"], float)
    np.savez(
        FIG_DIR / "pgsi_checkshot_geometry_timing.npz",
        depth=depth, receiver_east=east, receiver_north=north,
        shot_keys=np.array(keys), shot_delays=delays,
        shot1040_slant=pick["slant"], shot1040_observed=pick["observed"],
        shot1040_good=pick["good"], shot1040_residual=pick["residual"],
        shot1040_velocity=pick["velocity"], shot1040_intercept=pick["intercept"],
        shot1040_rms=pick["rms"],
    )

    report = [
        "SAFOD PGSI check-shot geometry/timing audit",
        "===========================================",
        "Decision: use the controlled PGSI explosions as the next highest-leverage",
        "test because surveyed receiver geometry and explicit SEG2 timing make a",
        "direct-P moveout/depth-registration hypothesis testable.",
        "",
        "Geometry: position-1 receiver depths %.3f--%.3f m; max |WELL_DEP-REC_DEP| %.3f m."
        % (depth.min(), depth.max(), np.max(np.abs(well_depth - depth))),
        "Source offsets are from shot_locs.doc (UTM NAD27 relative to the wellhead).",
        "The 1040 source depth follows the supplied R code (31 ft below KB); other",
        "hole depths are documented as 10 ft (shot 6: 15 ft) and are not inverted here.",
        "",
        "SEG2 timing inventory (do not silently harmonize):",
    ]
    for key, number, delay, acq in headers:
        report.append(f"  shot {number} file {key}.dat: DELAY={delay:+.3f} s; ACQUISITION_TIME={acq}")
    report += [
        "The legacy R note explicitly states that file zero is 200 ms before",
        "the shot origin (DELAY=-0.200 s). Two supplied files instead carry",
        "+0.200 s and 0.000 s, so a multi-shot timing correction remains required.",
        "",
        "Shot 1040 corrected axial channel test:",
        "  axial mapping: channels 3,6,...,240, reversed to shallow-to-deep",
        f"  retained AIC picks: {pick['good'].sum()}/80",
        f"  robust slope velocity: {pick['velocity']:.1f} m/s",
        f"  fitted origin offset: {pick['intercept']*1e3:.1f} ms",
        f"  residual RMS: {pick['rms']*1e3:.1f} ms",
        f"Interpretation: the corrected gather has a coherent, approximately {pick['velocity']/1000:.2f} km/s moveout slope after geometry correction,",
        f"but the fitted origin offset ({pick['intercept']*1e3:.1f} ms), {pick['good'].sum()}/80 retained picks, and inconsistent delays mean this is",
        "a validation result,",
        "not a standalone Vp inversion or proof that the 2026 Nano mode is direct P.",
        "The next step is to reconcile SEG2 header/shot-origin conventions across",
        "all six files, then compare the registered PGSI moveout with the Nano depth axis.",
    ]
    report_text = "\n".join(report) + "\n"
    (FIG_DIR / "pgsi_checkshot_geometry_timing.txt").write_text(report_text)
    print(report_text)
    print("Saved", FIG_DIR / "pgsi_checkshot_geometry_timing.png")


if __name__ == "__main__":
    main()

