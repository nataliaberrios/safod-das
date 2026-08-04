"""Split-sample validation of the low-frequency Deep DAS tube-wave candidate."""

from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, detrend, hilbert, sosfiltfilt, windows

from fk_dispersion import weighted_stack


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
PRE_S = 0.5
TIME_WINDOW = (-0.08, 2.85)
TURNAROUND_CH = 1702
APERTURE_M = 400.0
STEP_M = 200.0
BANDS = ((3.0, 15.0), (15.0, 30.0))
P_GRID = np.linspace(-1.0e-3, 1.0e-3, 201)
TUBE_POSITIVE = (1 / 1800.0, 1 / 1300.0)
N_CANDIDATE = 6
N_PERM = 499
SEED = 20260802


def phase_spectrum(section, fs, band):
    section = np.nan_to_num(section.astype(float))
    section = detrend(section, axis=-1, type="linear")
    section *= windows.tukey(section.shape[0], 0.15)[:, None]
    section *= windows.tukey(section.shape[1], 0.10)[None, :]
    nfft = int(2 ** np.ceil(np.log2(section.shape[1])))
    spectrum = np.fft.rfft(section, n=nfft, axis=1)
    frequency = np.fft.rfftfreq(nfft, 1 / fs)
    indices = np.flatnonzero((frequency >= band[0]) & (frequency <= band[1]))[::2]
    amplitude = np.abs(spectrum[:, indices])
    valid = amplitude > np.nanmedian(amplitude, axis=0, keepdims=True) * 1e-6
    unit = np.divide(spectrum[:, indices], amplitude, out=np.zeros_like(spectrum[:, indices]), where=valid)
    return frequency[indices], unit


def beam_profile(frequency, unit, coordinate, p_grid=P_GRID):
    x = coordinate - coordinate.mean()
    power = np.zeros(p_grid.size)
    used = 0
    for j, freq in enumerate(frequency):
        valid = np.abs(unit[:, j]) > 0
        if valid.sum() < 8:
            continue
        steering = np.exp(2j * np.pi * freq * np.outer(p_grid, x[valid]))
        power += (np.abs(steering @ unit[valid, j]) / valid.sum()) ** 2
        used += 1
    return power / max(1, used)


def fixed_power(frequency, unit, coordinate, p):
    x = coordinate - coordinate.mean()
    values = []
    for j, freq in enumerate(frequency):
        valid = np.abs(unit[:, j]) > 0
        if valid.sum() < 8:
            continue
        steering = np.exp(2j * np.pi * freq * p * x[valid])
        values.append((np.abs(steering @ unit[valid, j]) / valid.sum()) ** 2)
    return float(np.mean(values)) if values else np.nan


def apertures(n_channels, dx):
    width = max(16, int(round(APERTURE_M / dx)))
    step = max(8, int(round(STEP_M / dx)))
    starts = np.arange(0, n_channels - width + 1, step)
    centers = (starts + 0.5 * (width - 1)) * dx
    return starts, width, centers


def analyze_stack(section, dx, fs):
    starts, width, centers = apertures(section.shape[0], dx)
    power = np.empty((len(BANDS), starts.size, P_GRID.size), np.float32)
    for ib, band in enumerate(BANDS):
        for ia, start in enumerate(starts):
            stop = start + width
            frequency, unit = phase_spectrum(section[start:stop], fs, band)
            coordinate = np.arange(start, stop) * dx
            power[ib, ia] = beam_profile(frequency, unit, coordinate)
    return starts, width, centers, power


def select_candidates(discovery, centers):
    tube = (P_GRID >= TUBE_POSITIVE[0]) & (P_GRID <= TUBE_POSITIVE[1])
    candidates = {}
    for ib, band in enumerate(BANDS):
        local = discovery[ib][:, tube]
        peak_index = np.argmax(local, axis=1)
        peak_power = local[np.arange(local.shape[0]), peak_index]
        peak_p = P_GRID[tube][peak_index]
        selected = np.argsort(peak_power)[::-1][:min(N_CANDIDATE, len(centers))]
        candidates[ib] = (selected, peak_p[selected], peak_power[selected])
    return candidates


def permutation_test(section, starts, width, dx, fs, candidates, rng):
    results = {}
    for ib, band in enumerate(BANDS):
        selected, fixed_p, _ = candidates[ib]
        observed_values = []
        cached = []
        for index, p in zip(selected, fixed_p):
            start = starts[index]
            stop = start + width
            frequency, unit = phase_spectrum(section[start:stop], fs, band)
            coordinate = np.arange(start, stop) * dx
            observed_values.append(fixed_power(frequency, unit, coordinate, p))
            cached.append((frequency, unit, coordinate, p))
        observed = float(np.nanmedian(observed_values))
        null = np.empty(N_PERM)
        for k in range(N_PERM):
            values = []
            for frequency, unit, coordinate, p in cached:
                values.append(fixed_power(frequency, unit[rng.permutation(unit.shape[0])], coordinate, p))
            null[k] = np.nanmedian(values)
        probability = (1 + np.sum(null >= observed)) / (N_PERM + 1)
        results[ib] = (observed, null, probability, np.asarray(observed_values))
    return results


def bandpass(section, fs, band):
    sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, np.nan_to_num(section), axis=-1)


def line_pick(section, coordinate, fs, band, p):
    filtered = bandpass(section, fs, band)
    time = np.arange(filtered.shape[1]) / fs + TIME_WINDOW[0]
    intercepts = np.arange(-0.08, 0.401, 0.004)
    half = max(12, int(round(0.05 * fs)))
    best = (-np.inf, np.nan)
    for intercept in intercepts:
        centers = np.rint((intercept + p * coordinate - TIME_WINDOW[0]) * fs).astype(int)
        if np.any(centers - half < 0) or np.any(centers + half >= filtered.shape[1]):
            continue
        gather = np.stack([filtered[i, c-half:c+half] for i, c in enumerate(centers)])
        denominator = gather.shape[0] * np.sum(gather ** 2)
        score = np.sum(np.sum(gather, axis=0) ** 2) / denominator if denominator else np.nan
        if score > best[0]:
            best = (score, intercept)
    return best


def burst_statistics(stacks, counts, legs, dx, fs, starts_by_leg, width_by_leg, candidates_by_leg):
    rows = []
    for epoch in np.flatnonzero(counts > 0):
        for leg_name, slicer in legs:
            section = slicer(stacks[epoch])
            starts = starts_by_leg[leg_name]
            width = width_by_leg[leg_name]
            for ib, band in enumerate(BANDS):
                selected, fixed_p, _ = candidates_by_leg[leg_name][ib]
                positive = []
                negative = []
                for index, p in zip(selected, fixed_p):
                    start = starts[index]
                    stop = start + width
                    frequency, unit = phase_spectrum(section[start:stop], fs, band)
                    coordinate = np.arange(start, stop) * dx
                    positive.append(fixed_power(frequency, unit, coordinate, p))
                    negative.append(fixed_power(frequency, unit, coordinate, -p))
                rows.append({
                    "epoch": int(epoch), "n_common_drops": int(counts[epoch]),
                    "leg": leg_name, "band_low_hz": band[0], "band_high_hz": band[1],
                    "positive_power": float(np.nanmedian(positive)),
                    "negative_power": float(np.nanmedian(negative)),
                    "positive_minus_negative": float(np.nanmedian(positive) - np.nanmedian(negative)),
                })
    return rows


def main():
    rng = np.random.default_rng(SEED)
    with np.load(STACKS) as data:
        stacks = data["deep_stacks"]
        counts = data["n_common"]
        fs = float(data["fs"])
        dx = float(data["dx_deep"])
        valid_epochs = np.flatnonzero(counts > 0)
        discovery_counts = counts.copy(); discovery_counts[valid_epochs[valid_epochs % 2 == 0]] = 0
        validation_counts = counts.copy(); validation_counts[valid_epochs[valid_epochs % 2 == 1]] = 0
        full = weighted_stack(stacks, counts)
        discovery_full = weighted_stack(stacks, discovery_counts)
        validation_full = weighted_stack(stacks, validation_counts)

        i0 = int(round((PRE_S + TIME_WINDOW[0]) * fs))
        i1 = int(round((PRE_S + TIME_WINDOW[1]) * fs))
        legs = (
            ("outbound", lambda x: x[:TURNAROUND_CH, i0:i1]),
            ("return", lambda x: x[TURNAROUND_CH:, i0:i1][::-1]),
        )
        products = {}
        candidates_by_leg = {}
        starts_by_leg = {}; width_by_leg = {}; centers_by_leg = {}
        null_by_leg = {}
        candidate_rows = []
        line_segments = {}
        for leg_name, slicer in legs:
            discovery = slicer(discovery_full)
            validation = slicer(validation_full)
            starts, width, centers, discovery_power = analyze_stack(discovery, dx, fs)
            _, _, _, validation_power = analyze_stack(validation, dx, fs)
            candidates = select_candidates(discovery_power, centers)
            nulls = permutation_test(validation, starts, width, dx, fs, candidates, rng)
            starts_by_leg[leg_name] = starts; width_by_leg[leg_name] = width
            centers_by_leg[leg_name] = centers; candidates_by_leg[leg_name] = candidates
            null_by_leg[leg_name] = nulls
            products[f"{leg_name}_centers_m"] = centers
            products[f"{leg_name}_discovery_power"] = discovery_power
            products[f"{leg_name}_validation_power"] = validation_power
            for ib, band in enumerate(BANDS):
                selected, fixed_p, discovery_peak = candidates[ib]
                observed, null, probability, validation_values = nulls[ib]
                products[f"{leg_name}_band{ib}_null"] = null
                products[f"{leg_name}_band{ib}_selected"] = selected
                for rank, (index, p, dp, vp) in enumerate(zip(selected, fixed_p, discovery_peak, validation_values), 1):
                    neg = validation_power[ib, index, np.argmin(np.abs(P_GRID + p))]
                    start = starts[index]; stop = start + width
                    coordinate = np.arange(start, stop) * dx
                    score, intercept = line_pick(slicer(full)[start:stop], coordinate, fs, band, p)
                    candidate_rows.append({
                        "leg": leg_name, "band_low_hz": band[0], "band_high_hz": band[1],
                        "rank_discovery": rank, "center_m": centers[index],
                        "start_m": start * dx, "stop_m": (stop - 1) * dx,
                        "discovery_p_ms_per_m": p * 1e3, "velocity_mps": 1 / p,
                        "discovery_power": dp, "validation_fixed_power": vp,
                        "validation_negative_power": float(neg),
                        "validation_pos_neg_ratio": float(vp / (neg + 1e-30)),
                        "time_semblance": score, "surface_intercept_s": intercept,
                        "band_split_null_p": probability,
                    })
                    line_segments[(leg_name, ib, int(index))] = (start * dx, (stop - 1) * dx, p, intercept)

        burst_rows = burst_statistics(
            stacks, counts, legs, dx, fs, starts_by_leg, width_by_leg, candidates_by_leg
        )

    with (HERE / "deep_tube_candidates.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(candidate_rows[0]))
        writer.writeheader(); writer.writerows(candidate_rows)
    with (HERE / "deep_tube_burst_repeatability.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(burst_rows[0]))
        writer.writeheader(); writer.writerows(burst_rows)

    # Figure 1: time-domain branch mapping.
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    time = np.arange(i1 - i0) / fs + TIME_WINDOW[0]
    for row, (leg_name, slicer) in enumerate(legs):
        section = slicer(full)
        distance = np.arange(section.shape[0]) * dx
        for col, band in enumerate(BANDS):
            filtered = bandpass(section, fs, band)
            shown = filtered[::4]
            dshown = distance[::4]
            # Display-only robust trace scaling is required because a shallow
            # high-amplitude interval otherwise controls the full-leg color
            # scale. Quantitative beams and null tests remain unnormalized.
            trace_scale = np.percentile(np.abs(shown), 99.0, axis=1, keepdims=True)
            positive_scale = trace_scale[np.isfinite(trace_scale) & (trace_scale > 0)]
            floor = np.median(positive_scale) * 1e-3 if positive_scale.size else 1.0
            shown = shown / np.maximum(trace_scale, floor)
            clip = 1.0
            ax = axes[row, col]
            ax.imshow(shown, aspect="auto", origin="upper", cmap="RdBu_r",
                      extent=[time[0], time[-1], dshown[-1], dshown[0]],
                      vmin=-clip, vmax=clip)
            selected = candidates_by_leg[leg_name][col][0]
            for index in selected:
                x0, x1, p, intercept = line_segments[(leg_name, col, int(index))]
                ax.plot([intercept + p*x0, intercept + p*x1], [x0, x1],
                        color="lime", lw=1.2, alpha=0.9)
            ax.set_title(f"{leg_name.capitalize()}, {band[0]:.0f}–{band[1]:.0f} Hz")
            ax.set_xlabel("time relative to catalog pick (s)")
            ax.set_ylabel("distance from leg start (m)")
    fig.suptitle(
        "SAFOD AWD Deep: time-domain test of split-sample slow-mode trajectories\n"
        "display only: each trace divided by its 99th-percentile absolute amplitude"
    )
    fig.tight_layout(); fig.savefig(HERE / "deep_tube_record_sections.png", dpi=200)

    # Figure 2: independent discovery and validation maps.
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex="col", sharey="row")
    tube = (P_GRID >= TUBE_POSITIVE[0]) & (P_GRID <= TUBE_POSITIVE[1])
    for row, leg_name in enumerate(("outbound", "return")):
        centers = centers_by_leg[leg_name]
        for col, band in enumerate(BANDS):
            d = products[f"{leg_name}_discovery_power"][col][:, tube]
            v = products[f"{leg_name}_validation_power"][col][:, tube]
            dp = np.max(d, axis=1); vp = np.max(v, axis=1)
            ax = axes[row, col]
            ax.plot(centers, dp, "o-", color="0.45", ms=4, label="odd-burst discovery")
            ax.plot(centers, vp, "o-", color="C0", ms=4, label="even-burst validation")
            selected = candidates_by_leg[leg_name][col][0]
            ax.scatter(centers[selected], dp[selected], facecolors="none", edgecolors="crimson", s=55, label="preselected")
            ax.set_title(f"{leg_name.capitalize()}, {band[0]:.0f}–{band[1]:.0f} Hz")
            ax.set_xlabel("aperture-center distance (m)")
            ax.set_ylabel("maximum tube-range beam power")
            ax.grid(alpha=0.2)
            if row == 0 and col == 0: ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Deep low-frequency ridge: independent burst-split repeatability")
    fig.tight_layout(); fig.savefig(HERE / "deep_tube_repeatability.png", dpi=200)

    # Figure 3: fixed-p spatial permutation nulls plus burst directionality.
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    burst_array = burst_rows
    for row, leg_name in enumerate(("outbound", "return")):
        for col, band in enumerate(BANDS):
            observed, null, probability, _ = null_by_leg[leg_name][col]
            ax = axes[row, col]
            ax.hist(null, bins=30, color="0.72", edgecolor="white")
            ax.axvline(observed, color="crimson", lw=2,
                       label=f"even-burst observed={observed:.3f}\np={probability:.3f}")
            subset = [r for r in burst_array if r["leg"] == leg_name and r["band_low_hz"] == band[0]]
            fraction = np.mean([r["positive_minus_negative"] > 0 for r in subset])
            ax.set_title(f"{leg_name.capitalize()}, {band[0]:.0f}–{band[1]:.0f} Hz\n{fraction:.0%} bursts: +p power > −p power")
            ax.set_xlabel("median fixed-slowness beam power")
            ax.set_ylabel("spatial permutations")
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Deep slow-mode validation against channel-order nulls")
    fig.tight_layout(); fig.savefig(HERE / "deep_tube_null_tests.png", dpi=200)

    np.savez(HERE / "deep_tube_validation.npz", p_grid=P_GRID, bands=np.asarray(BANDS),
             time_window=np.asarray(TIME_WINDOW), seed=SEED, n_perm=N_PERM, **products)
    lines = ["SAFOD AWD Deep tube-wave candidate: split-sample validation",
             f"discovery epochs: {np.count_nonzero(discovery_counts)}; validation epochs: {np.count_nonzero(validation_counts)}",
             f"spatial permutations: {N_PERM}; seed: {SEED}"]
    for leg_name in ("outbound", "return"):
        for ib, band in enumerate(BANDS):
            observed, _, probability, _ = null_by_leg[leg_name][ib]
            subset = [r for r in burst_rows if r["leg"] == leg_name and r["band_low_hz"] == band[0]]
            fraction = np.mean([r["positive_minus_negative"] > 0 for r in subset])
            lines.append(f"{leg_name} {band[0]:.0f}-{band[1]:.0f} Hz: validation fixed-p power={observed:.5f}, permutation p={probability:.4f}, bursts +p>-p={fraction:.3f}")
    report = "\n".join(lines) + "\n"
    (HERE / "deep_tube_validation.txt").write_text(report)
    print(report)
    print("Saved Deep tube-wave validation figures, tables, NPZ, and report")


if __name__ == "__main__":
    main()
