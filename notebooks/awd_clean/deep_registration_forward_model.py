"""Deep-fiber registration and bounded tube-wave forward consistency test.

This is deliberately a registration/consistency analysis, not a velocity or
permeability inversion.  The raw Deep HDF5 metadata provide an interrogator
StartLocusIndex and a data-derived hairpin channel, but no surveyed borehole
depth, fiber path, or casing/fluid properties.  We therefore keep the measured
coordinate explicit and compare the validated slow-mode candidates with a
bounded fluid/compliance envelope.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE
DX_M = 2.0419
N_CHANNELS = 3200
START_LOCUS_INDEX = 1800
TURNAROUND_CH = 1702
FLUID_MIN = 1300.0
FLUID_MAX = 1800.0
# A phenomenological compliance reduction factor.  It is an envelope, not a
# fitted borehole model: actual tube-wave speed requires fluid, casing, radius,
# cement, and formation-compliance parameters that are not in this workspace.
COMPLIANCE_MIN = 0.75
COMPLIANCE_MAX = 1.00


def read_candidates(path: Path):
    rows = []
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            row["velocity_mps"] = float(row["velocity_mps"])
            row["center_m"] = float(row["center_m"])
            row["band_low_hz"] = float(row["band_low_hz"])
            row["band_high_hz"] = float(row["band_high_hz"])
            row["rank_discovery"] = int(row["rank_discovery"])
            row["validation_fixed_power"] = float(row["validation_fixed_power"])
            row["validation_negative_power"] = float(row["validation_negative_power"])
            row["band_split_null_p"] = float(row["band_split_null_p"])
            rows.append(row)
    return rows


def read_model_arrays(model_dir: Path):
    arrays = {}
    for name in ("Vp_model.dat", "Vs_model.dat", "Vpvs_model.dat"):
        arrays[name] = np.loadtxt(model_dir / name)
    return arrays


def model_stats(arrays):
    stats = {}
    for name, array in arrays.items():
        finite = array[np.isfinite(array)]
        stats[name] = {
            "n": int(finite.size),
            "min": float(np.min(finite)),
            "median": float(np.median(finite)),
            "max": float(np.max(finite)),
            "p05": float(np.percentile(finite, 5)),
            "p95": float(np.percentile(finite, 95)),
        }
    return stats


def build_registration():
    channel = np.arange(N_CHANNELS)
    coordinate = channel * DX_M
    start_locus_m = START_LOCUS_INDEX * DX_M
    turnaround_m = TURNAROUND_CH * DX_M
    end_m = (N_CHANNELS - 1) * DX_M
    return channel, coordinate, start_locus_m, turnaround_m, end_m


def make_figure(rows, model_stats_dict, output: Path):
    channel, coordinate, start_locus_m, turnaround_m, end_m = build_registration()
    velocities = np.asarray([row["velocity_mps"] for row in rows])
    rng = np.random.default_rng(20260802)
    xj = np.asarray([row["center_m"] for row in rows]) + rng.normal(0, 18, len(rows))
    colors = {"outbound": "#1565c0", "return": "#d95f02"}

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.4), constrained_layout=True)
    ax = axes[0, 0]
    ax.axhspan(-0.22, 0.22, color="#eaf2f8")
    ax.plot(coordinate, np.zeros_like(coordinate), color="#273746", lw=2.2)
    ax.axvspan(0, turnaround_m, color="#1565c0", alpha=0.13, label="outbound leg")
    ax.axvspan(turnaround_m, end_m, color="#d95f02", alpha=0.13, label="return leg (raw order)")
    ax.axvline(start_locus_m, color="#7d3c98", ls="--", lw=1.7)
    ax.axvline(turnaround_m, color="#17202a", ls=":", lw=2.1)
    ax.text(start_locus_m, 0.11, f"StartLocusIndex=1800\n{start_locus_m:,.1f} m", ha="center", va="bottom", fontsize=9, color="#7d3c98")
    ax.text(turnaround_m, -0.11, f"hairpin channel 1702\n{turnaround_m:,.1f} m", ha="center", va="top", fontsize=9)
    ax.set(yticks=[], xlabel="Interrogator fiber coordinate from Deep channel 0 (m)",
           title="(a) What is measured: coordinate, not surveyed depth")
    ax.legend(loc="upper left", frameon=False, ncol=2, fontsize=8)
    ax.set_xlim(-100, end_m + 100); ax.set_ylim(-0.28, 0.29)

    ax = axes[0, 1]
    s_out = np.linspace(0, turnaround_m, 200)
    s_ret = np.linspace(0, end_m - turnaround_m, 200)
    ax.plot(s_out, np.zeros_like(s_out), color="#1565c0", lw=6, label="outbound: channel 0 → hairpin")
    ax.plot(s_ret, np.ones_like(s_ret), color="#d95f02", lw=6, label="return: hairpin → channel 3199 after reversal")
    ax.set_yticks([0, 1], ["outbound", "return (reversed)"])
    ax.set_xlabel("Distance along provisional leg coordinate (m)")
    ax.set_title("(b) Orientation used by the Deep analyses")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.text(0.02, 0.05, "Both origins and their relation to measured depth remain unknown.", transform=ax.transAxes, fontsize=8.5)

    ax = axes[1, 0]
    cf = np.linspace(FLUID_MIN, FLUID_MAX, 200)
    ax.fill_between(cf, cf * COMPLIANCE_MIN, cf * COMPLIANCE_MAX, color="#76b7b2", alpha=0.34,
                    label="bounded fluid/compliance envelope")
    ax.axhspan(FLUID_MIN, FLUID_MAX, color="#a9cce3", alpha=0.17, label="fluid-speed prior (1.3–1.8 km/s)")
    for leg in ("outbound", "return"):
        leg_rows = [row for row in rows if row["leg"] == leg]
        for band in ((3.0, 15.0), (15.0, 30.0)):
            sub = [row for row in leg_rows if row["band_low_hz"] == band[0]]
            yy = [row["velocity_mps"] for row in sub]
            xx = [row["center_m"] for row in sub]
            ax.scatter(xx, yy, s=36, color=colors[leg], marker="o" if band[0] == 3.0 else "s",
                       alpha=0.82, edgecolor="white", linewidth=0.35,
                       label=f"{leg}, {band[0]:g}–{band[1]:g} Hz" )
    ax.set(xlabel="Provisional leg coordinate of candidate center (m)", ylabel="Candidate phase velocity (m/s)",
           title="(c) Validated slow mode versus bounded tube-wave prior")
    ax.set_ylim(1050, 2050); ax.grid(alpha=0.22)
    handles, labels = ax.get_legend_handles_labels()
    # Keep the legend compact and avoid duplicate envelope labels.
    ax.legend(handles[:6], labels[:6], frameon=False, fontsize=7.5, ncol=2, loc="upper left")
    ax.text(0.02, 0.04, "Compatibility is not identification: no casing/fluid fit is possible here.", transform=ax.transAxes, fontsize=8.3)

    ax = axes[1, 1]
    names = ["Vp (km/s)", "Vs (km/s)", "Vp/Vs"]
    keys = ["Vp_model.dat", "Vs_model.dat", "Vpvs_model.dat"]
    med = [model_stats_dict[key]["median"] for key in keys]
    p05 = [model_stats_dict[key]["p05"] for key in keys]
    p95 = [model_stats_dict[key]["p95"] for key in keys]
    y = np.arange(3)
    ax.errorbar(med, y, xerr=[np.asarray(med)-np.asarray(p05), np.asarray(p95)-np.asarray(med)],
                fmt="o", color="#2c3e50", ecolor="#85929e", capsize=4, lw=2)
    ax.axvspan(1.3, 1.8, color="#76b7b2", alpha=0.18)
    ax.set_yticks(y, names); ax.set_xlabel("Supplied model value (Vp and Vs in km/s)")
    ax.set_title("(d) Formation model: background only, not depth registered")
    ax.grid(axis="x", alpha=0.25)
    ax.text(0.02, 0.04, "No borehole log or fiber-to-depth transform was supplied.", transform=ax.transAxes, fontsize=8.3)

    fig.suptitle("Deep AWD fiber registration and bounded tube-wave consistency test", fontsize=15, fontweight="bold")
    fig.savefig(output, dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def main():
    rows = read_candidates(OUT / "deep_tube_candidates.csv")
    arrays = read_model_arrays(OUT / "vel_model")
    stats = model_stats(arrays)
    channel, coordinate, start_locus_m, turnaround_m, end_m = build_registration()
    velocities = np.asarray([row["velocity_mps"] for row in rows])
    fluid_compat = (velocities >= FLUID_MIN * COMPLIANCE_MIN) & (velocities <= FLUID_MAX * COMPLIANCE_MAX)
    # The stricter direct fluid prior is reported separately from the envelope.
    direct_compat = (velocities >= FLUID_MIN) & (velocities <= FLUID_MAX)

    out_csv = OUT / "deep_registration_forward_model.csv"
    with out_csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["quantity", "value", "units", "interpretation"])
        writer.writerow(["dx", DX_M, "m/channel", "instrument spatial sampling interval"])
        writer.writerow(["n_channels", N_CHANNELS, "channels", "raw Deep array used by stack builder"])
        writer.writerow(["start_locus_index", START_LOCUS_INDEX, "index", "interrogator absolute locus metadata"])
        writer.writerow(["start_locus_coordinate", start_locus_m, "m", "absolute interrogator coordinate; not borehole depth"])
        writer.writerow(["turnaround_channel", TURNAROUND_CH, "index", "data-derived hairpin split used by Deep analyses"])
        writer.writerow(["turnaround_coordinate", turnaround_m, "m", "coordinate from channel 0; not surveyed depth"])
        writer.writerow(["array_end_coordinate", end_m, "m", "(N-1)*dx"])
        writer.writerow(["candidate_count", len(rows), "candidates", "six ranked candidates per leg and band"])
        writer.writerow(["candidate_velocity_median", float(np.median(velocities)), "m/s", "validated slow-mode candidate median"])
        writer.writerow(["candidate_velocity_p05", float(np.percentile(velocities, 5)), "m/s", "candidate distribution"])
        writer.writerow(["candidate_velocity_p95", float(np.percentile(velocities, 95)), "m/s", "candidate distribution"])
        writer.writerow(["direct_fluid_compatibility_fraction", float(np.mean(direct_compat)), "fraction", "within 1.3–1.8 km/s prior"])
        writer.writerow(["bounded_envelope_compatibility_fraction", float(np.mean(fluid_compat)), "fraction", "within 0.75–1.0 times 1.3–1.8 km/s envelope"])
        writer.writerow(["physical_depth_registration", "unresolved", "status", "raw metadata contain no surveyed fiber-to-depth transform"])
        writer.writerow(["borehole_log_comparison", "not available", "status", "no borehole logs found in workspace"])

    with (OUT / "deep_registration_forward_model.txt").open("w") as stream:
        stream.write("SAFOD AWD Deep registration and bounded tube-wave forward test\n")
        stream.write("Status: conditional consistency test; no physical depth registration or inversion.\n\n")
        stream.write(f"Deep metadata: dx={DX_M:.4f} m/channel, N={N_CHANNELS}, StartLocusIndex={START_LOCUS_INDEX} "
                     f"({start_locus_m:.2f} m interrogator coordinate), hairpin channel={TURNAROUND_CH} "
                     f"({turnaround_m:.2f} m from channel 0), end={(end_m):.2f} m.\n")
        stream.write("The StartLocusIndex and hairpin coordinate are not measured depth. The HDF5 metadata provide no surveyed fiber path, orientation, or depth origin.\n")
        stream.write(f"Validated candidate velocity median={np.median(velocities):.1f} m/s; 5th–95th percentile="
                     f"{np.percentile(velocities,5):.1f}–{np.percentile(velocities,95):.1f} m/s.\n")
        stream.write(f"Direct 1.3–1.8 km/s fluid-prior compatibility={np.mean(direct_compat):.3f}; bounded 0.75–1.0 compliance envelope compatibility={np.mean(fluid_compat):.3f}.\n")
        stream.write("Interpretation: the slow mode is numerically compatible with a tube-wave-like speed range, but this result cannot distinguish fluid, casing, formation compliance, or another guided mode.\n")
        stream.write("The supplied Vp/Vs arrays are plotted only as formation-background context; no candidate can be assigned to a depth or fracture.\n")

    with (OUT / "deep_registration_forward_model.json").open("w") as stream:
        json.dump({
            "status": "conditional_consistency_test",
            "physical_depth_registration": "unresolved",
            "borehole_logs_available": False,
            "dx_m": DX_M, "n_channels": N_CHANNELS,
            "start_locus_index": START_LOCUS_INDEX, "start_locus_coordinate_m": start_locus_m,
            "turnaround_channel": TURNAROUND_CH, "turnaround_coordinate_m": turnaround_m,
            "array_end_coordinate_m": end_m,
            "candidate_count": len(rows),
            "candidate_velocity_median_mps": float(np.median(velocities)),
            "candidate_velocity_p05_mps": float(np.percentile(velocities, 5)),
            "candidate_velocity_p95_mps": float(np.percentile(velocities, 95)),
            "direct_fluid_prior_mps": [FLUID_MIN, FLUID_MAX],
            "compliance_factor_envelope": [COMPLIANCE_MIN, COMPLIANCE_MAX],
            "direct_fluid_compatibility_fraction": float(np.mean(direct_compat)),
            "bounded_envelope_compatibility_fraction": float(np.mean(fluid_compat)),
            "velocity_model_stats": stats,
            "model_limit": "No surveyed Deep fiber-to-depth transform or borehole logs found in workspace.",
        }, stream, indent=2)

    np.savez(OUT / "deep_registration_forward_model.npz",
             candidate_velocity_mps=velocities, candidate_center_m=np.asarray([row["center_m"] for row in rows]),
             candidate_leg=np.asarray([row["leg"] for row in rows]),
             candidate_band_low_hz=np.asarray([row["band_low_hz"] for row in rows]),
             coordinate=np.asarray(coordinate), channel=np.asarray(channel),
             start_locus_coordinate_m=start_locus_m, turnaround_coordinate_m=turnaround_m,
             array_end_coordinate_m=end_m, fluid_prior_mps=np.asarray([FLUID_MIN, FLUID_MAX]),
             compliance_factor=np.asarray([COMPLIANCE_MIN, COMPLIANCE_MAX]))
    make_figure(rows, stats, OUT / "deep_registration_forward_model.png")
    print((OUT / "deep_registration_forward_model.txt").read_text())


if __name__ == "__main__":
    main()
