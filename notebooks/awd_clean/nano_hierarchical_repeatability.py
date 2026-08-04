"""Hierarchical repeatability of the SAFOD Nano AWD mode at the drop and burst scales.

This analysis is deliberately phase-neutral: it measures the already observed
30--60 Hz moveout at 2,975 m/s over 80--440 m without calling that arrival a P
wave or a guided mode.  Each individual AWD observation is referenced to the
literal ``UTC_Date`` value in ``p26.cc9.txt``.  That value is used only as an
alignment timestamp; its physical meaning is not assumed here.

Raw Nano files are read once per required contiguous file segment.  The data
are filtered to 30--60 Hz by the established DAS reader, then shifted along the
fixed 2,975 m/s trajectory and averaged over the aperture.  Only the resulting
short one-dimensional beam waveforms and scalar features are retained.  This
avoids creating an approximately terabyte-scale individual-drop DAS cube.

Outputs
-------
nano_drop_repeatability.csv
    One row per successfully extracted individual drop.  Metrics use a
    leave-one-drop-out template from the same burst.
nano_burst_repeatability_hierarchical.csv
    One row per burst.  Metrics use a leave-one-burst-out reference.
nano_stack_convergence.csv
    Repeated disjoint within-burst substack comparisons for 1, 2, 4, and 8
    drops per substack.
nano_hierarchical_repeatability.npz
    Reduced beam waveforms, templates, coordinates, and configuration.
nano_hierarchical_repeatability.png
    Publication-style six-panel summary.
nano_hierarchical_repeatability.txt
    Human-readable processing and result audit.
"""

from __future__ import annotations

import csv
import glob
import math
import os
import sys
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MANIFEST = HERE / "awd_manifest.csv"
P26 = ROOT / "p26.cc9.txt"
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
NANO_DIR = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/Nano")

UTILITY_DIRS = (
    Path("/home/groups/edunham/nberrios/safod_das/DAS-utilities/python"),
    Path("/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python"),
)
for directory in UTILITY_DIRS:
    if directory.exists() and str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from DASutils import readFile_protobuf  # noqa: E402


BAND_HZ = (30.0, 60.0)
APERTURE_M = (80.0, 440.0)
VELOCITY_MPS = 2975.0
INTERCEPT_S = -0.022
TAU_RANGE_S = (-0.40, 0.20)
SIGNAL_TAU_S = (-0.04, 0.16)
NOISE_TAU_S = (-0.36, -0.16)
MAX_LAG_S = 0.012
RAW_FILE_DURATION_S = 300.0
FILTER_PAD_S = 2.0
CONVERGENCE_N = (1, 2, 4, 8)
CONVERGENCE_REPEATS = 100
RANDOM_SEED = 20260802

DROP_CSV = HERE / "nano_drop_repeatability.csv"
BURST_CSV = HERE / "nano_burst_repeatability_hierarchical.csv"
CONVERGENCE_CSV = HERE / "nano_stack_convergence.csv"
OUT_NPZ = HERE / "nano_hierarchical_repeatability.npz"
OUT_PNG = HERE / "nano_hierarchical_repeatability.png"
OUT_TXT = HERE / "nano_hierarchical_repeatability.txt"


def parse_nano_start(path: str | Path) -> datetime:
    fields = os.path.basename(path).split("_")
    value = f"{fields[1]}T{fields[2].replace('.', ':')}Z"
    return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(
        tzinfo=timezone.utc
    )


def read_rows() -> list[dict[str, str]]:
    with MANIFEST.open(newline="") as stream:
        rows = [row for row in csv.DictReader(stream) if int(row["nano_available"]) == 1]
    if not rows:
        raise RuntimeError(f"No Nano observations found in {MANIFEST}")
    rows.sort(key=lambda row: (int(row["burst_id"]), int(row["drop_id"])))
    return rows


def p26_utc_dates() -> set[datetime]:
    result: set[datetime] = set()
    with P26.open() as stream:
        next(stream)
        for line in stream:
            fields = line.split()
            if len(fields) >= 2:
                value = datetime.fromisoformat(fields[1])
                if value.tzinfo is None:
                    value = value.replace(tzinfo=timezone.utc)
                result.add(value.astimezone(timezone.utc))
    return result


def audit_timestamp_provenance(rows: list[dict[str, str]]) -> None:
    source = p26_utc_dates()
    absent = []
    for row in rows:
        value = datetime.fromisoformat(row["utc_time"]).astimezone(timezone.utc)
        if value not in source:
            absent.append(row["utc_time"])
    if absent:
        raise RuntimeError(
            f"{len(absent)} manifest UTC values are not literal UTC_Date entries in {P26}"
        )


def raw_file_index() -> tuple[list[Path], list[datetime]]:
    paths = sorted(Path(path) for path in glob.glob(str(NANO_DIR / "*.pb")))
    starts = [parse_nano_start(path) for path in paths]
    if not paths:
        raise RuntimeError(f"No Nano protobuf files found in {NANO_DIR}")
    return paths, starts


def segment_for_drop(
    utc_date: datetime,
    paths: list[Path],
    starts: list[datetime],
    min_arrival_s: float,
    max_arrival_s: float,
) -> tuple[tuple[Path, ...], datetime, datetime]:
    requested_start = utc_date + timedelta(
        seconds=min_arrival_s + TAU_RANGE_S[0] - FILTER_PAD_S
    )
    requested_stop = utc_date + timedelta(
        seconds=max_arrival_s + TAU_RANGE_S[1] + FILTER_PAD_S
    )
    selected = tuple(
        path
        for path, start in zip(paths, starts)
        if start < requested_stop
        and start + timedelta(seconds=RAW_FILE_DURATION_S) > requested_start
    )
    if not selected:
        raise RuntimeError(
            f"No raw Nano file segment spans requested window around {utc_date.isoformat()}"
        )
    selected_starts = [parse_nano_start(path) for path in selected]
    coverage_start = min(selected_starts)
    coverage_stop = max(selected_starts) + timedelta(seconds=RAW_FILE_DURATION_S)
    if requested_start < coverage_start or requested_stop > coverage_stop:
        raise RuntimeError(
            f"Raw Nano coverage is incomplete around {utc_date.isoformat()}"
        )
    return selected, requested_start, requested_stop


def fixed_masks(tau_s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    signal = (tau_s >= SIGNAL_TAU_S[0]) & (tau_s <= SIGNAL_TAU_S[1])
    noise = (tau_s >= NOISE_TAU_S[0]) & (tau_s <= NOISE_TAU_S[1])
    if signal.sum() != noise.sum():
        raise RuntimeError("Signal and noise windows must contain equal sample counts")
    return signal, noise


def moveout_beam(
    section: np.ndarray,
    section_start: datetime,
    utc_date: datetime,
    fs: float,
    channels: np.ndarray,
    distance_m: np.ndarray,
    tau_s: np.ndarray,
) -> np.ndarray:
    """Linearly interpolate and mean the fixed-trajectory aperture."""
    reference_s = (utc_date - section_start).total_seconds()
    # The accepted scan defines INTERCEPT_S at the first channel of the
    # 80--440 m aperture, because the scan coordinate was zeroed there.
    arrival_s = INTERCEPT_S + (distance_m - distance_m[0]) / VELOCITY_MPS
    sample = (
        reference_s * fs
        + arrival_s[:, None] * fs
        + tau_s[None, :] * fs
    )
    lower = np.floor(sample).astype(np.int64)
    fraction = sample - lower
    if lower.min() < 0 or lower.max() + 1 >= section.shape[1]:
        raise RuntimeError(
            f"Moveout window around {utc_date.isoformat()} exceeds the read segment"
        )
    values0 = section[channels[:, None], lower]
    values1 = section[channels[:, None], lower + 1]
    aligned = values0 + fraction * (values1 - values0)
    return np.mean(aligned, axis=0, dtype=np.float64).astype(np.float32)


def paired_metric(
    waveform: np.ndarray,
    reference: np.ndarray,
    mask: np.ndarray,
    fs: float,
    max_lag_s: float = MAX_LAG_S,
) -> tuple[float, float, float]:
    """Return signed NCC, delay of waveform relative to reference, and LS gain.

    A positive lag means that the waveform is delayed relative to the reference.
    The gain is the least-squares multiplier that maps the lag-matched reference
    to the waveform; it is not an absolutely calibrated source amplitude.
    """
    x_full = np.asarray(waveform[mask], dtype=np.float64)
    y_full = np.asarray(reference[mask], dtype=np.float64)
    maximum = int(round(max_lag_s * fs))
    correlations = []
    gains = []
    lags = np.arange(-maximum, maximum + 1)
    for lag in lags:
        if lag > 0:
            x, y = x_full[lag:], y_full[:-lag]
        elif lag < 0:
            x, y = x_full[:lag], y_full[-lag:]
        else:
            x, y = x_full, y_full
        xd = x - np.mean(x)
        yd = y - np.mean(y)
        denominator = np.linalg.norm(xd) * np.linalg.norm(yd)
        correlations.append(float(np.dot(xd, yd) / denominator) if denominator else np.nan)
        yy = float(np.dot(y, y))
        gains.append(float(np.dot(x, y) / yy) if yy else np.nan)
    correlations = np.asarray(correlations)
    if not np.any(np.isfinite(correlations)):
        return np.nan, np.nan, np.nan
    index = int(np.nanargmax(correlations))
    fractional = 0.0
    if 0 < index < correlations.size - 1:
        left, center, right = correlations[index - 1:index + 2]
        curvature = left - 2.0 * center + right
        if np.isfinite(curvature) and curvature < 0:
            fractional = float(np.clip(0.5 * (left - right) / curvature, -1.0, 1.0))
    delay_s = (lags[index] + fractional) / fs
    return float(correlations[index]), float(delay_s), float(gains[index])


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values, dtype=np.float64) ** 2)))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for {path.name}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def percentile_band(values: np.ndarray) -> tuple[float, float, float]:
    finite = values[np.isfinite(values)]
    if not finite.size:
        return np.nan, np.nan, np.nan
    return tuple(float(value) for value in np.percentile(finite, [16, 50, 84]))


def one_sided_sign_probability(n_positive: int, n_total: int) -> float:
    """Exact P[X >= n_positive] for X ~ Binomial(n_total, 0.5)."""
    return float(
        sum(math.comb(n_total, value) for value in range(n_positive, n_total + 1))
        / (2 ** n_total)
    )


def plot_results(
    drop_rows: list[dict[str, object]],
    burst_rows: list[dict[str, object]],
    convergence_rows: list[dict[str, object]],
) -> None:
    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "savefig.facecolor": "white",
    })
    blue = "#0072B2"
    orange = "#D55E00"
    gray = "0.62"
    burst = np.asarray([row["burst_id"] for row in drop_rows], dtype=int)
    drop_order = np.asarray([row["drop_id"] for row in drop_rows], dtype=int)
    signal_corr = np.asarray([row["loo_signal_ncc"] for row in drop_rows], float)
    noise_corr = np.asarray([row["loo_noise_ncc"] for row in drop_rows], float)
    lag_ms = 1e3 * np.asarray([row["loo_signal_delay_s"] for row in drop_rows], float)
    rel_amp = np.asarray([row["relative_amplitude"] for row in drop_rows], float)
    snr_db = np.asarray([row["beam_snr_db"] for row in drop_rows], float)

    fig, axes = plt.subplots(3, 2, figsize=(13.0, 11.2), constrained_layout=True)
    ax = axes[0, 0]
    jitter = (drop_order / np.maximum(1, np.asarray([row["n_drops_in_burst"] for row in drop_rows], int) - 1) - 0.5) * 0.58
    ax.scatter(burst + jitter, noise_corr, s=7, color=gray, alpha=0.30, linewidths=0, label="equal-duration pre-reference noise")
    ax.scatter(burst + jitter, signal_corr, s=8, color=blue, alpha=0.48, linewidths=0, label="30–60 Hz mode window")
    for burst_id in np.unique(burst):
        use = burst == burst_id
        ax.plot(burst_id, np.nanmedian(signal_corr[use]), "_", color="k", ms=7, mew=1.0)
    ax.axhline(0, color="0.2", lw=0.7, ls="--")
    ax.set(xlim=(-1, np.max(burst) + 1), ylim=(-0.5, 1.02), xlabel="AWD burst index", ylabel="leave-one-drop-out normalized correlation", title="A  Individual-drop waveform repeatability within bursts")
    ax.legend(frameon=False, loc="lower right")

    ax = axes[0, 1]
    ax.scatter(burst + jitter, lag_ms, s=8, color=blue, alpha=0.48, linewidths=0)
    for burst_id in np.unique(burst):
        use = burst == burst_id
        ax.plot(burst_id, np.nanmedian(lag_ms[use]), "_", color="k", ms=7, mew=1.0)
    ax.axhline(0, color="0.2", lw=0.7, ls="--")
    ax.set(xlim=(-1, np.max(burst) + 1), xlabel="AWD burst index", ylabel="delay relative to within-burst LOO template (ms)", title="B  Drop-to-drop timing stability")

    ax = axes[1, 0]
    ax.scatter(drop_order, rel_amp, s=9, color=orange, alpha=0.26, linewidths=0)
    unique_order = np.unique(drop_order)
    order_counts = np.asarray([np.sum(drop_order == item) for item in unique_order])
    summary_order = unique_order[order_counts >= 10]
    med = [np.nanmedian(rel_amp[drop_order == item]) for item in summary_order]
    p16 = [np.nanpercentile(rel_amp[drop_order == item], 16) for item in summary_order]
    p84 = [np.nanpercentile(rel_amp[drop_order == item], 84) for item in summary_order]
    ax.fill_between(summary_order, p16, p84, color=orange, alpha=0.16, linewidth=0)
    ax.plot(summary_order, med, color=orange, lw=1.7)
    ax.axhline(1, color="0.2", lw=0.7, ls="--")
    ax.set(xlabel="drop order within burst", ylabel="signal RMS / within-burst median RMS", title="C  Relative amplitude through each burst")
    ax.text(0.98, 0.97, "median and 16–84% shown for ≥10 bursts",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)

    ax = axes[1, 1]
    ax.scatter(drop_order, snr_db, s=9, color=blue, alpha=0.26, linewidths=0)
    med = [np.nanmedian(snr_db[drop_order == item]) for item in summary_order]
    p16 = [np.nanpercentile(snr_db[drop_order == item], 16) for item in summary_order]
    p84 = [np.nanpercentile(snr_db[drop_order == item], 84) for item in summary_order]
    ax.fill_between(summary_order, p16, p84, color=blue, alpha=0.16, linewidth=0)
    ax.plot(summary_order, med, color=blue, lw=1.7)
    ax.axhline(0, color="0.2", lw=0.7, ls="--")
    ax.set(xlabel="drop order within burst", ylabel="moveout-beam SNR (dB)", title="D  Frequency-band observability of individual drops")
    ax.text(0.98, 0.97, "median and 16–84% shown for ≥10 bursts",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)

    ax = axes[2, 0]
    burst_id = np.asarray([row["burst_id"] for row in burst_rows], int)
    burst_corr = np.asarray([row["loo_burst_signal_ncc"] for row in burst_rows], float)
    burst_lag = 1e3 * np.asarray([row["loo_burst_delay_s"] for row in burst_rows], float)
    points = ax.scatter(burst_id, burst_corr, c=burst_lag, cmap="coolwarm", vmin=-6, vmax=6, s=32, edgecolors="white", linewidths=0.35)
    ax.axhline(0, color="0.2", lw=0.7, ls="--")
    ax.set(xlim=(-1, np.max(burst_id) + 1), ylim=(-0.1, 1.02), xlabel="AWD burst index", ylabel="leave-one-burst-out normalized correlation", title="E  Burst-template repeatability across the experiment")
    colorbar = fig.colorbar(points, ax=ax, pad=0.01)
    colorbar.set_label("burst-template delay (ms)")

    ax = axes[2, 1]
    n_values = np.asarray(sorted({int(row["n_drops_per_substack"]) for row in convergence_rows}))
    conv_corr = np.asarray([row["independent_substack_ncc"] for row in convergence_rows], float)
    conv_n = np.asarray([row["n_drops_per_substack"] for row in convergence_rows], int)
    low, center, high = [], [], []
    for n_drop in n_values:
        q16, q50, q84 = percentile_band(conv_corr[conv_n == n_drop])
        low.append(q16); center.append(q50); high.append(q84)
    ax.fill_between(n_values, low, high, color=blue, alpha=0.18, linewidth=0)
    ax.plot(n_values, center, marker="o", color=blue, lw=1.8)
    ax.set(xscale="log", xticks=n_values, xticklabels=[str(value) for value in n_values], ylim=(-0.1, 1.02), xlabel="drops in each independent within-burst substack", ylabel="normalized correlation between disjoint substacks", title="F  Waveform convergence with stacking")

    fig.suptitle(
        "SAFOD AWD Nano hierarchical repeatability: 30–60 Hz, 80–440 m, fixed 2,975 m/s moveout",
        fontsize=13,
    )
    fig.savefig(OUT_PNG, dpi=260)
    plt.close(fig)


def main() -> None:
    rows = read_rows()
    audit_timestamp_provenance(rows)
    with np.load(STACKS) as product:
        fs_expected = float(product["fs"])
        dx = float(product["dx_nano"])
        n_channels_expected = int(product["nano_stacks"].shape[1])

    c0 = int(np.ceil(APERTURE_M[0] / dx))
    c1 = min(n_channels_expected, int(np.floor(APERTURE_M[1] / dx)) + 1)
    channels = np.arange(c0, c1, dtype=int)
    distance_m = channels * dx
    min_arrival_s = float(INTERCEPT_S)
    max_arrival_s = float(
        INTERCEPT_S + (distance_m.max() - distance_m.min()) / VELOCITY_MPS
    )
    tau_s = np.arange(
        int(round(TAU_RANGE_S[0] * fs_expected)),
        int(round(TAU_RANGE_S[1] * fs_expected)) + 1,
    ) / fs_expected
    signal_mask, noise_mask = fixed_masks(tau_s)

    paths, starts = raw_file_index()
    grouped: dict[tuple[Path, ...], list[dict[str, str]]] = defaultdict(list)
    bounds: dict[tuple[Path, ...], list[datetime]] = {}
    for row in rows:
        utc_date = datetime.fromisoformat(row["utc_time"]).astimezone(timezone.utc)
        segment, requested_start, requested_stop = segment_for_drop(
            utc_date, paths, starts, min_arrival_s, max_arrival_s
        )
        grouped[segment].append(row)
        if segment not in bounds:
            bounds[segment] = [requested_start, requested_stop]
        else:
            bounds[segment][0] = min(bounds[segment][0], requested_start)
            bounds[segment][1] = max(bounds[segment][1], requested_stop)

    beam_waveforms = []
    valid_rows = []
    for segment_index, (segment, segment_rows) in enumerate(
        sorted(grouped.items(), key=lambda item: parse_nano_start(item[0][0]))
    ):
        read_start, read_stop = bounds[segment]
        print(
            f"segment {segment_index + 1:03d}/{len(grouped):03d}: "
            f"{len(segment_rows):2d} drops, {len(segment)} raw file(s), "
            f"{read_start.isoformat()} to {read_stop.isoformat()}",
            flush=True,
        )
        section, info = readFile_protobuf(
            [str(path) for path in segment],
            fmin=BAND_HZ[0],
            fmax=BAND_HZ[1],
            desampling=False,
            verbose=False,
            minTime=read_start,
            maxTime=read_stop,
            tapering=False,
            order=4,
            zerophase=True,
            median=True,
            detrend=True,
        )
        fs = float(info["fs"])
        if not np.isclose(fs, fs_expected):
            raise RuntimeError(f"Unexpected sampling rate {fs:g} Hz")
        if section.shape[0] != n_channels_expected:
            raise RuntimeError(
                f"Raw segment contains {section.shape[0]} rather than "
                f"{n_channels_expected} Nano channels"
            )
        if not np.isclose(float(info["dx"]), dx):
            raise RuntimeError("Raw and canonical Nano channel spacing disagree")
        section_start = info["begTime"]
        if section_start.tzinfo is None:
            section_start = section_start.replace(tzinfo=timezone.utc)
        for row in segment_rows:
            utc_date = datetime.fromisoformat(row["utc_time"]).astimezone(timezone.utc)
            beam_waveforms.append(
                moveout_beam(
                    section, section_start, utc_date, fs, channels,
                    distance_m, tau_s,
                )
            )
            valid_rows.append(row)
        del section

    order = np.asarray(sorted(
        range(len(valid_rows)),
        key=lambda index: (int(valid_rows[index]["burst_id"]), int(valid_rows[index]["drop_id"])),
    ))
    waveforms = np.asarray(beam_waveforms, dtype=np.float32)[order]
    valid_rows = [valid_rows[index] for index in order]

    by_burst: dict[int, list[int]] = defaultdict(list)
    for index, row in enumerate(valid_rows):
        by_burst[int(row["burst_id"])].append(index)
    burst_ids = np.asarray(sorted(by_burst), dtype=int)
    burst_templates = np.asarray(
        [np.mean(waveforms[by_burst[burst_id]], axis=0) for burst_id in burst_ids],
        dtype=np.float32,
    )

    drop_output: list[dict[str, object]] = []
    for burst_id in burst_ids:
        indices = by_burst[int(burst_id)]
        group = waveforms[indices]
        if len(indices) < 2:
            raise RuntimeError(f"Burst {burst_id} has fewer than two extracted drops")
        total = np.sum(group, axis=0, dtype=np.float64)
        signal_rms = np.asarray([rms(item[signal_mask]) for item in group])
        median_signal_rms = float(np.median(signal_rms))
        for local_index, global_index in enumerate(indices):
            waveform = waveforms[global_index]
            loo = ((total - waveform) / (len(indices) - 1)).astype(np.float32)
            corr, delay_s, gain = paired_metric(waveform, loo, signal_mask, fs_expected)
            noise_corr, _, _ = paired_metric(waveform, loo, noise_mask, fs_expected)
            noise_rms = rms(waveform[noise_mask])
            beam_snr_db = 20.0 * np.log10(
                (signal_rms[local_index] + np.finfo(float).tiny)
                / (noise_rms + np.finfo(float).tiny)
            )
            row = valid_rows[global_index]
            drop_output.append({
                "burst_id": int(burst_id),
                "drop_id": int(row["drop_id"]),
                "utc_date_from_p26": row["utc_time"],
                "nano_file_listed_in_manifest": row["nano_file"],
                "n_drops_in_burst": len(indices),
                "loo_signal_ncc": corr,
                "loo_noise_ncc": noise_corr,
                "loo_signal_delay_s": delay_s,
                "loo_least_squares_gain": gain,
                "signal_rms": signal_rms[local_index],
                "noise_rms": noise_rms,
                "relative_amplitude": signal_rms[local_index] / median_signal_rms,
                "beam_snr_db": beam_snr_db,
            })

    normalized_templates = []
    burst_signal_rms = []
    for template in burst_templates:
        scale = rms(template[signal_mask])
        burst_signal_rms.append(scale)
        normalized_templates.append(template / scale if scale > 0 else template)
    normalized_templates = np.asarray(normalized_templates, dtype=np.float32)
    burst_signal_rms = np.asarray(burst_signal_rms)
    template_sum = np.sum(normalized_templates, axis=0, dtype=np.float64)
    burst_output: list[dict[str, object]] = []
    for index, burst_id in enumerate(burst_ids):
        loo = ((template_sum - normalized_templates[index]) / (len(burst_ids) - 1)).astype(np.float32)
        corr, delay_s, gain = paired_metric(
            normalized_templates[index], loo, signal_mask, fs_expected
        )
        noise_corr, _, _ = paired_metric(
            normalized_templates[index], loo, noise_mask, fs_expected
        )
        burst_output.append({
            "burst_id": int(burst_id),
            "n_drops": len(by_burst[int(burst_id)]),
            "loo_burst_signal_ncc": corr,
            "loo_burst_noise_ncc": noise_corr,
            "loo_burst_delay_s": delay_s,
            "loo_burst_least_squares_gain_after_rms_normalization": gain,
            "burst_signal_rms": burst_signal_rms[index],
            "burst_relative_amplitude": burst_signal_rms[index] / np.median(burst_signal_rms),
        })

    rng = np.random.default_rng(RANDOM_SEED)
    convergence_output: list[dict[str, object]] = []
    for burst_id in burst_ids:
        group = waveforms[by_burst[int(burst_id)]]
        for n_drop in CONVERGENCE_N:
            if 2 * n_drop > group.shape[0]:
                continue
            for trial in range(CONVERGENCE_REPEATS):
                selection = rng.permutation(group.shape[0])[: 2 * n_drop]
                first = np.mean(group[selection[:n_drop]], axis=0)
                second = np.mean(group[selection[n_drop:]], axis=0)
                corr, delay_s, gain = paired_metric(
                    first, second, signal_mask, fs_expected
                )
                convergence_output.append({
                    "burst_id": int(burst_id),
                    "n_drops_per_substack": int(n_drop),
                    "trial": trial,
                    "independent_substack_ncc": corr,
                    "independent_substack_delay_s": delay_s,
                    "independent_substack_least_squares_gain": gain,
                })

    write_csv(DROP_CSV, drop_output)
    write_csv(BURST_CSV, burst_output)
    write_csv(CONVERGENCE_CSV, convergence_output)
    plot_results(drop_output, burst_output, convergence_output)

    np.savez_compressed(
        OUT_NPZ,
        beam_waveforms=waveforms,
        burst_templates=burst_templates,
        normalized_burst_templates=normalized_templates,
        burst_id=np.asarray([int(row["burst_id"]) for row in valid_rows]),
        drop_id=np.asarray([int(row["drop_id"]) for row in valid_rows]),
        utc_date=np.asarray([row["utc_time"] for row in valid_rows]),
        tau_s=tau_s,
        signal_mask=signal_mask,
        noise_mask=noise_mask,
        aperture_channel=channels,
        aperture_distance_m=distance_m,
        fs=fs_expected,
        dx=dx,
        band_hz=BAND_HZ,
        aperture_m=APERTURE_M,
        velocity_mps=VELOCITY_MPS,
        intercept_s=INTERCEPT_S,
        max_lag_s=MAX_LAG_S,
        convergence_n=CONVERGENCE_N,
        convergence_repeats=CONVERGENCE_REPEATS,
        random_seed=RANDOM_SEED,
        alignment_source="p26.cc9.txt UTC_Date",
    )

    drop_corr = np.asarray([row["loo_signal_ncc"] for row in drop_output], float)
    drop_noise_corr = np.asarray([row["loo_noise_ncc"] for row in drop_output], float)
    drop_lag_ms = 1e3 * np.asarray([row["loo_signal_delay_s"] for row in drop_output], float)
    drop_amp = np.asarray([row["relative_amplitude"] for row in drop_output], float)
    drop_snr = np.asarray([row["beam_snr_db"] for row in drop_output], float)
    burst_corr = np.asarray([row["loo_burst_signal_ncc"] for row in burst_output], float)
    burst_lag_ms = 1e3 * np.asarray([row["loo_burst_delay_s"] for row in burst_output], float)
    burst_signal_minus_noise = []
    for burst_id in burst_ids:
        subset = [
            row for row in drop_output if int(row["burst_id"]) == int(burst_id)
        ]
        burst_signal_minus_noise.append(
            np.median([float(row["loo_signal_ncc"]) for row in subset])
            - np.median([float(row["loo_noise_ncc"]) for row in subset])
        )
    n_bursts_signal_greater = int(np.sum(np.asarray(burst_signal_minus_noise) > 0))
    sign_probability = one_sided_sign_probability(
        n_bursts_signal_greater, len(burst_signal_minus_noise)
    )
    report = [
        "SAFOD Nano AWD hierarchical repeatability",
        "=========================================",
        f"Input observations: {len(rows)} Nano-available rows from {MANIFEST.name}",
        f"Extracted individual drops: {len(drop_output)} in {len(burst_output)} bursts",
        "Alignment timestamp: literal p26.cc9.txt UTC_Date; physical source-time meaning is not assumed",
        f"Fixed mode test: {BAND_HZ[0]:g}-{BAND_HZ[1]:g} Hz, {APERTURE_M[0]:g}-{APERTURE_M[1]:g} m, {VELOCITY_MPS:g} m/s, intercept {INTERCEPT_S:+.3f} s at the first aperture channel",
        "Reduced observable: unnormalized equal-channel mean after trajectory alignment; no individual DAS waveform cube is saved",
        f"Signal tau window: {SIGNAL_TAU_S[0]:+.3f} to {SIGNAL_TAU_S[1]:+.3f} s",
        f"Equal-duration noise tau window: {NOISE_TAU_S[0]:+.3f} to {NOISE_TAU_S[1]:+.3f} s",
        f"Cross-correlation lag search: +/-{MAX_LAG_S * 1e3:.1f} ms; positive delay means the observation is later than its LOO reference",
        "",
        f"Within-burst drop signal NCC (16/50/84%): {percentile_band(drop_corr)}",
        f"Within-burst drop noise NCC (16/50/84%): {percentile_band(drop_noise_corr)}",
        f"Bursts whose median signal NCC exceeds median noise NCC: {n_bursts_signal_greater}/{len(burst_signal_minus_noise)}",
        f"Exact one-sided burst-level sign-test probability: {sign_probability:.6g}",
        f"Within-burst drop delay, ms (16/50/84%): {percentile_band(drop_lag_ms)}",
        f"Within-burst relative amplitude (16/50/84%): {percentile_band(drop_amp)}",
        f"Individual-drop moveout-beam SNR, dB (16/50/84%): {percentile_band(drop_snr)}",
        f"Across-burst template NCC (16/50/84%): {percentile_band(burst_corr)}",
        f"Across-burst template delay, ms (16/50/84%): {percentile_band(burst_lag_ms)}",
        "",
        "Independent within-burst substack convergence:",
    ]
    conv_n = np.asarray([row["n_drops_per_substack"] for row in convergence_output], int)
    conv_corr = np.asarray([row["independent_substack_ncc"] for row in convergence_output], float)
    for n_drop in sorted(set(conv_n)):
        report.append(
            f"  n={n_drop:g} drops per substack, NCC 16/50/84%: "
            f"{percentile_band(conv_corr[conv_n == n_drop])}"
        )
    report.extend([
        "",
        "Interpretive boundary: these metrics quantify repeatability of a fixed phase-neutral moveout beam. They do not by themselves identify the mode as direct P, infer formation V_P, or establish absolute source amplitude.",
        f"Saved: {DROP_CSV.name}, {BURST_CSV.name}, {CONVERGENCE_CSV.name}, {OUT_NPZ.name}, {OUT_PNG.name}",
    ])
    OUT_TXT.write_text("\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    main()
