"""Targeted Deep tube-wave observability test at the provisional SAFOD strand depths.

This is a depth-conditional analysis.  The Deep interrogator coordinate is not
surveyed measured depth, so the report labels every depth as *assumed* and
records the transform used here: outbound s=depth and reversed return
s=turnaround-depth.  The purpose is to test whether the existing slow-mode
observable is present in the windows that cover the SDZ/CDZ, not to claim that
the current experiment detected creep.

Discovery and validation use disjoint non-empty stack epochs.  Each targeted
window is scanned in both signed slowness directions and in two frequency
bands.  Channel-order permutations of the validation window provide a
target-specific null at 3192 and 3302 m.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, sosfiltfilt

from fk_dispersion import weighted_stack


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"

PRE_S = 0.5
TIME_WINDOW = (-0.08, 2.85)
TURNAROUND_CH = 1702
DX_M = 2.0419
TURNAROUND_M = TURNAROUND_CH * DX_M

TARGET_MIN_M = 3000.0
TARGET_MAX_M = 3450.0
SDZ_M = 3192.0
CDZ_M = 3302.0
DAMAGE_MIN_M = 3150.0
DAMAGE_MAX_M = 3414.0

# A 200-m aperture localizes the target interval while leaving enough channels
# for phase coherence.  Target centers are sampled every 50 m; strand depths
# are included exactly in the output even when they fall between grid points.
APERTURE_M = 200.0
CENTER_STEP_M = 50.0
BANDS = ((3.0, 15.0), (15.0, 30.0))
P_POS = np.linspace(1.0 / 1800.0, 1.0 / 1300.0, 21)
TIME_DECIMATION = 4
CHANNEL_STRIDE = 2
WINDOW_S = 0.12
N_PERM = 499
SEED = 20260802

OUT_CSV = HERE / "deep_target_scan.csv"
OUT_NPZ = HERE / "deep_target_scan.npz"
OUT_PNG = HERE / "deep_target_scan.png"
OUT_PDF = HERE / "deep_target_scan.pdf"
OUT_TXT = HERE / "deep_target_scan.txt"
OUT_JSON = HERE / "deep_target_scan.json"


def prepare(section: np.ndarray, fs: float, band: tuple[float, float]):
    sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, np.nan_to_num(section), axis=-1)[:, ::TIME_DECIMATION]


def aligned(section: np.ndarray, coordinate: np.ndarray, fs: float, p: float):
    center = coordinate.mean()
    shifts = np.rint(p * (coordinate - center) * fs).astype(int)
    out = np.full_like(section, np.nan)
    for i, shift in enumerate(shifts):
        if shift > 0:
            out[i, :-shift] = section[i, shift:]
        elif shift < 0:
            out[i, -shift:] = section[i, :shift]
        else:
            out[i] = section[i]
    return out


def semblance_series(section, coordinate, fs, p):
    gather = aligned(section, coordinate, fs, p)
    valid = np.isfinite(gather)
    values = np.nan_to_num(gather)
    numerator = np.sum(values, axis=0) ** 2
    denominator = valid.sum(axis=0) * np.sum(values**2, axis=0)
    length = max(3, int(round(WINDOW_S * fs)))
    numerator = uniform_filter1d(numerator, length, mode="constant")
    denominator = uniform_filter1d(denominator, length, mode="constant")
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)


def scan(section, coordinate, fs, time):
    allowed = (time >= 0.0) & (time <= 2.75)
    result = {}
    for direction, p_values in (("positive", P_POS), ("negative", -P_POS)):
        best = (-np.inf, np.nan, np.nan)
        for p in p_values:
            series = semblance_series(section, coordinate, fs, float(p))
            if not np.any(allowed):
                continue
            index = np.flatnonzero(allowed)[np.argmax(series[allowed])]
            if series[index] > best[0]:
                best = (float(series[index]), float(p), float(time[index]))
        result[direction] = best
    return result


def fixed_score(section, coordinate, fs, time, p, target_time):
    series = semblance_series(section, coordinate, fs, p)
    allowed = np.abs(time - target_time) <= 0.06
    return float(np.max(series[allowed])) if np.any(allowed) else np.nan


def trajectory_rms(section, coordinate, fs, p, target_time):
    """RMS in a 60-ms aligned window, normalized by all-window RMS."""
    gather = aligned(section, coordinate, fs, p)
    center = int(np.argmin(np.abs(np.arange(gather.shape[1]) / fs + TIME_WINDOW[0] - target_time)))
    half = max(2, int(round(0.03 * fs)))
    signal = gather[:, max(0, center - half): center + half + 1]
    all_rms = np.sqrt(np.nanmean(gather**2))
    signal_rms = np.sqrt(np.nanmean(signal**2))
    return float(signal_rms / all_rms) if np.isfinite(all_rms) and all_rms > 0 else np.nan


def target_centers():
    grid = np.arange(TARGET_MIN_M, TARGET_MAX_M + 0.1, CENTER_STEP_M)
    return np.unique(np.r_[grid, SDZ_M, CDZ_M])


def window_for_target(n_channels: int, target_s: float, dx: float):
    width = max(16, int(round(APERTURE_M / dx)))
    target_index = int(round(target_s / dx))
    start = target_index - width // 2
    start = max(0, min(start, n_channels - width))
    stop = start + width
    return start, stop, width


def leg_section(stack, leg: str, i0: int, i1: int):
    if leg == "outbound":
        return stack[:TURNAROUND_CH, i0:i1], np.arange(TURNAROUND_CH) * DX_M
    section = stack[TURNAROUND_CH:, i0:i1][::-1]
    return section, np.arange(section.shape[0]) * DX_M


def target_to_leg_coordinate(depth_m: float, leg: str):
    if leg == "outbound":
        return float(depth_m)
    return float(TURNAROUND_M - depth_m)


def target_permutation(section, coordinate, fs, time, p, target_time, rng):
    observed = fixed_score(section, coordinate, fs, time, p, target_time)
    null = np.empty(N_PERM, dtype=float)
    for k in range(N_PERM):
        null[k] = fixed_score(section[rng.permutation(section.shape[0])], coordinate, fs, time, p, target_time)
    p_value = (1.0 + np.sum(null >= observed)) / (N_PERM + 1.0)
    return observed, null, float(p_value)


def make_figure(rows, nulls, full_sections, time):
    fig = plt.figure(figsize=(14.5, 11.0), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=(1.0, 1.0, 1.15))
    colors = {"positive": "#1f77b4", "negative": "#d62728"}

    # Target-specific observability map.
    ax = fig.add_subplot(grid[0, :])
    for leg, ls in (("outbound", "-"), ("return", "--")):
        for direction in ("positive", "negative"):
            q = [r for r in rows if r["leg"] == leg and r["direction"] == direction]
            q.sort(key=lambda r: r["assumed_depth_m"])
            ax.plot([r["assumed_depth_m"] for r in q], [r["validation_semblance"] for r in q],
                    linestyle=ls, marker="o", ms=3.4, lw=1.5, color=colors[direction],
                    alpha=0.8, label=f"{leg}, {direction}")
    ax.axvspan(DAMAGE_MIN_M, DAMAGE_MAX_M, color="#f4b183", alpha=0.25, label="published damage-zone interval")
    ax.axvline(SDZ_M, color="#7f6000", lw=1.4, ls=":")
    ax.axvline(CDZ_M, color="#7f6000", lw=1.4, ls="--")
    ax.text(SDZ_M, 0.98, "SDZ", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=9)
    ax.text(CDZ_M, 0.98, "CDZ", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=9)
    ax.set(xlim=(TARGET_MIN_M, TARGET_MAX_M), ylim=(0, 1.02), ylabel="even-epoch validation semblance",
           title="Targeted Deep tube-wave observability under provisional depth mapping")
    ax.grid(alpha=0.2); ax.legend(frameon=False, ncol=3, fontsize=8, loc="upper left")

    # Record sections at the two proposed strand depths.
    for col, depth in enumerate((SDZ_M, CDZ_M)):
        ax = fig.add_subplot(grid[1, col])
        section, coordinate, p, target_time = full_sections[depth]
        shown = section[::2]
        scale = np.nanpercentile(np.abs(shown), 99.0, axis=1, keepdims=True)
        scale = np.maximum(scale, np.nanmedian(scale[np.isfinite(scale) & (scale > 0)]) * 1e-3)
        shown = shown / scale
        ax.imshow(shown, aspect="auto", origin="upper", cmap="RdBu_r", vmin=-1, vmax=1,
                  extent=[time[0], time[-1], coordinate[-1], coordinate[0]])
        ax.set(xlabel="time relative to stack alignment (s)", ylabel="local provisional coordinate (m)",
               title=f"Validation window near assumed depth {depth:.0f} m")
        ax.axvline(target_time, color="lime", lw=1.2, label="best trajectory time")
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, fontsize=8, loc="upper right")

    # Nulls at the strands, summarized by leg/direction/band.
    ax = fig.add_subplot(grid[2, :])
    labels = []
    positions = []
    for depth in (SDZ_M, CDZ_M):
        for key in sorted(nulls):
            if key[0] == depth:
                labels.append(f"{depth:.0f} m\n{key[1]} {key[2]} {key[3]:g}–{key[4]:g} Hz")
                positions.append((key, len(labels)))
    for key, pos in positions:
        observed, null, p_value = nulls[key]
        ax.hist(null, bins=30, alpha=0.34, density=True, color="#9e9e9e")
        ax.axvline(observed, color="#c00000", lw=1.5)
        ax.text(pos, 0.98, f"p={p_value:.3f}", transform=ax.get_xaxis_transform(), rotation=90,
                va="top", ha="right", fontsize=7)
    ax.set_xlabel("fixed-trajectory validation semblance under channel-order permutation")
    ax.set_ylabel("density")
    ax.set_title("Target-depth channel-order nulls (red lines: observed validation statistic)")
    ax.grid(alpha=0.2)
    fig.suptitle("SAFOD AWD Deep: targeted test of tube-wave observability at the provisional SDZ/CDZ depths",
                 fontsize=15, fontweight="bold")
    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    plt.close(fig)


def main():
    rng = np.random.default_rng(SEED)
    with np.load(STACKS) as data:
        stacks = data["deep_stacks"]
        counts = data["n_common"]
        fs0 = float(data["fs"])
        dx = float(data["dx_deep"])
        epochs = np.flatnonzero(counts > 0)
        discovery_counts = counts.copy()
        validation_counts = counts.copy()
        discovery_counts[epochs[epochs % 2 == 1]] = 0
        validation_counts[epochs[epochs % 2 == 0]] = 0
        discovery_full = weighted_stack(stacks, discovery_counts)
        validation_full = weighted_stack(stacks, validation_counts)
        full = weighted_stack(stacks, counts)

    i0 = int(round((PRE_S + TIME_WINDOW[0]) * fs0))
    i1 = int(round((PRE_S + TIME_WINDOW[1]) * fs0))
    fs = fs0 / TIME_DECIMATION
    time = np.arange(0, i1 - i0, TIME_DECIMATION) / fs0 + TIME_WINDOW[0]
    rows = []
    nulls = {}
    full_sections = {}
    target_grid = target_centers()

    for leg in ("outbound", "return"):
        disc_leg, coord_leg = leg_section(discovery_full, leg, i0, i1)
        val_leg, _ = leg_section(validation_full, leg, i0, i1)
        full_leg, full_coord = leg_section(full, leg, i0, i1)
        for depth in target_grid:
            target_s = target_to_leg_coordinate(depth, leg)
            start, stop, width = window_for_target(disc_leg.shape[0], target_s, dx)
            channels = np.arange(start, stop, CHANNEL_STRIDE)
            coordinate = np.arange(start, stop) * dx
            coordinate = coordinate[::CHANNEL_STRIDE]
            dsection = prepare(disc_leg[start:stop], fs0, (3.0, 15.0))  # replaced below per band
            for band in BANDS:
                dsec = prepare(disc_leg[start:stop], fs0, band)[::CHANNEL_STRIDE]
                vsec = prepare(val_leg[start:stop], fs0, band)[::CHANNEL_STRIDE]
                result = scan(dsec, coordinate, fs, time)
                for direction in ("positive", "negative"):
                    dsem, p, tc = result[direction]
                    vsem = fixed_score(vsec, coordinate, fs, time, p, tc)
                    trms = trajectory_rms(vsec, coordinate, fs, p, tc)
                    row = {
                        "leg": leg,
                        "band_low_hz": band[0],
                        "band_high_hz": band[1],
                        "assumed_depth_m": float(depth),
                        "target_leg_coordinate_m": float(target_s),
                        "window_start_m": float(start * dx),
                        "window_stop_m": float((stop - 1) * dx),
                        "window_center_m": float((start + 0.5 * (width - 1)) * dx),
                        "n_channels": int(len(channels)),
                        "direction": direction,
                        "discovery_slowness_ms_per_m": float(p * 1e3),
                        "discovery_velocity_mps": float(1.0 / p),
                        "discovery_center_time_s": float(tc),
                        "discovery_semblance": float(dsem),
                        "validation_semblance": float(vsem),
                        "validation_trajectory_rms_ratio": float(trms),
                        "target_null_permutation_p": np.nan,
                    }
                    rows.append(row)
                    if depth in (SDZ_M, CDZ_M):
                        observed, null, pv = target_permutation(vsec, coordinate, fs, time, p, tc, rng)
                        key = (float(depth), leg, direction, band[0], band[1])
                        nulls[key] = (observed, null, pv)
                        row["target_null_permutation_p"] = pv

        # Store full-band record sections only for the target windows. The
        # plot is intentionally display-only; all statistics use the reduced
        # arrays above.
        for depth in (SDZ_M, CDZ_M):
            target_s = target_to_leg_coordinate(depth, leg)
            start, stop, _ = window_for_target(full_leg.shape[0], target_s, dx)
            if leg == "outbound":
                # Use the 3–15 Hz validation window and its best positive pick.
                band = BANDS[0]
                dsec = prepare(disc_leg[start:stop], fs0, band)[::CHANNEL_STRIDE]
                result = scan(dsec, np.arange(start, stop, CHANNEL_STRIDE) * dx, fs, time)
                target_time = result["positive"][2]
            else:
                band = BANDS[0]
                dsec = prepare(disc_leg[start:stop], fs0, band)[::CHANNEL_STRIDE]
                result = scan(dsec, np.arange(start, stop, CHANNEL_STRIDE) * dx, fs, time)
                target_time = result["positive"][2]
            sec = prepare(full_leg[start:stop], fs0, band)[::CHANNEL_STRIDE]
            full_sections.setdefault(float(depth), (sec, np.arange(start, stop, CHANNEL_STRIDE) * dx, result["positive"][1], target_time))

    with OUT_CSV.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    np.savez(OUT_NPZ,
             assumed_depth_m=np.asarray(sorted({r["assumed_depth_m"] for r in rows})),
             target_min_m=TARGET_MIN_M, target_max_m=TARGET_MAX_M,
             sdz_m=SDZ_M, cdz_m=CDZ_M, damage_min_m=DAMAGE_MIN_M, damage_max_m=DAMAGE_MAX_M,
             turnaround_m=TURNAROUND_M, aperture_m=APERTURE_M, center_step_m=CENTER_STEP_M,
             bands=np.asarray(BANDS), p_pos=P_POS, n_perm=N_PERM, seed=SEED,
             row_count=len(rows), null_keys=np.asarray([str(k) for k in nulls]))

    make_figure(rows, nulls, full_sections, time)
    null_summary = []
    for key, (observed, null, pv) in sorted(nulls.items()):
        null_summary.append({"depth_m": key[0], "leg": key[1], "direction": key[2],
                             "band_low_hz": key[3], "band_high_hz": key[4],
                             "observed": float(observed), "null_median": float(np.median(null)),
                             "null_p95": float(np.percentile(null, 95)), "p_value": float(pv)})
    report_lines = [
        "SAFOD AWD Deep targeted tube-wave observability scan",
        "Status: target-specific validation; provisional coordinate-to-depth mapping only.",
        f"Assumed mapping: outbound s=depth; reversed return s={TURNAROUND_M:.2f}-depth (m).",
        f"Target range: {TARGET_MIN_M:.0f}–{TARGET_MAX_M:.0f} m; SDZ={SDZ_M:.0f} m; CDZ={CDZ_M:.0f} m; damage interval={DAMAGE_MIN_M:.0f}–{DAMAGE_MAX_M:.0f} m.",
        f"Aperture={APERTURE_M:.0f} m; target center step={CENTER_STEP_M:.0f} m; bands={BANDS}; channel-order permutations={N_PERM}.",
        f"Rows={len(rows)}; non-empty epochs={len(epochs)}; discovery/validation split by epoch parity.",
        "A non-significant target test means the target window was tested and did not validate; it is not evidence that the array lacks a tube wave elsewhere.",
        "The current 24-hour experiment cannot measure the rate of aseismic creep; it establishes a candidate baseline and a repeat-survey sensitivity test.",
        "",
        "Target null summaries:",
    ]
    for item in null_summary:
        report_lines.append(f"{item['depth_m']:.0f} m {item['leg']} {item['direction']} {item['band_low_hz']:.0f}–{item['band_high_hz']:.0f} Hz: observed={item['observed']:.4f}, null median={item['null_median']:.4f}, null p95={item['null_p95']:.4f}, p={item['p_value']:.4f}")
    OUT_TXT.write_text("\n".join(report_lines) + "\n")
    OUT_JSON.write_text(json.dumps({
        "status": "target_specific_validation_provisional_depth",
        "mapping": "outbound s=depth; return s=turnaround-depth after reversal",
        "turnaround_m": TURNAROUND_M, "target_range_m": [TARGET_MIN_M, TARGET_MAX_M],
        "sdz_m": SDZ_M, "cdz_m": CDZ_M, "damage_interval_m": [DAMAGE_MIN_M, DAMAGE_MAX_M],
        "aperture_m": APERTURE_M, "center_step_m": CENTER_STEP_M,
        "bands_hz": [list(b) for b in BANDS], "n_permutation": N_PERM, "seed": SEED,
        "row_count": len(rows), "null_summaries": null_summary,
    }, indent=2) + "\n")
    print(OUT_TXT.read_text())
    print(f"Saved {OUT_CSV.name}, {OUT_NPZ.name}, {OUT_PNG.name}, {OUT_PDF.name}, {OUT_JSON.name}")


if __name__ == "__main__":
    main()
