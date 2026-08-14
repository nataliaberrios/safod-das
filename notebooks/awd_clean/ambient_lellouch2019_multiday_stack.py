#!/usr/bin/env python3
"""Coherent multi-day stack of the paper-faithful Figure 7c observable.

WHY THIS IS NOT THE PER-DAY TABLE. `FIG7C_MULTIDAY_RESULT.md` reports five
independent one-day tests and combines their p values with Fisher's method.
Combining p values is *not* the same experiment as stacking the data. A weak
coherent arrival grows with the number of stacked windows while incoherent noise
averages down, so a single 120-hour stack is strictly more sensitive to it than
five separate 24-hour tests are, however their p values are pooled. Fisher can
only tell you whether five noisy tests were jointly surprising; it cannot build
the SNR that a coherent sum builds.

This is therefore the most sensitive test of Figure 7c available from this
archive, and the last one outstanding.

HOW IT STAYS FAITHFUL. Chunk files store summed cross spectra, which are exactly
additive across days: summing them and dividing by the total window count is
identical to having processed one long record. Every downstream operation --
Equation 6 handling, the inverse transform, the 5-20 Hz bandpass on the stacked
correlations, the moveout scan and the receiver-order null -- is imported from
`ambient_lellouch2019_exact_stack` rather than reimplemented, so the arithmetic
is bit-identical to the per-day aggregates it is compared against.

Only configuration 0 (paper baseline) chunks are stacked, and only from days
whose manifest is exactly continuous. Chunk metadata is checked for agreement
across days exactly as the single-day aggregate checks it across hours; date is
not part of that metadata, so a genuine mismatch still raises.

Output: ambient_lellouch2019_multiday_stack.{npz,png,txt}
"""
from __future__ import annotations

import glob
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ambient_lellouch2019_exact_stack import (  # noqa: E402
    VELOCITY_GRID_M_S,
    correlation_from_spectrum,
    moveout_scores,
    prepare_aggregate_spectra,
    receiver_order_null,
)

STEM = HERE / "ambient_lellouch2019_multiday_stack"
DIRS = [
    HERE / "ambient_transfer" / "lellouch2019_exact_stack",
    HERE / "ambient_transfer" / "lellouch2019_exact_stack_days",
]
PATTERN = "chunk_*_src23_ram0p1_cross_correlation_ordered_r0_start*_n0060.npz"
V_REF = 3200.0
NULL_COUNT = 10000
SEED = 20260814


def collect():
    paths = []
    for d in DIRS:
        for p in sorted(glob.glob(str(d / PATTERN))):
            if "_cm_" in Path(p).name:
                continue
            paths.append(Path(p))
    by_date = {}
    for p in paths:
        m = re.search(r"chunk_(\d{4}-\d{2}-\d{2})_", p.name)
        by_date.setdefault(m.group(1), []).append(p)
    return {d: sorted(v) for d, v in by_date.items()}


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    by_date = collect()
    say("Coherent multi-day stack of the paper-faithful Figure 7c observable")
    say("configuration 0 only; no F-K filter; math imported from the single-day operator")
    say("")
    for d in sorted(by_date):
        say("  %s : %d hourly chunks" % (d, len(by_date[d])))
    dates = [d for d in sorted(by_date) if len(by_date[d]) == 24]
    for d in sorted(by_date):
        if len(by_date[d]) != 24:
            say("  skipping %s (incomplete: %d/24)" % (d, len(by_date[d])))
            dates.remove(d) if d in dates else None

    # The archive is NOT homogeneous in acquisition rate: 2024-05-11 was recorded
    # at 5000 Hz (n_fft 524288, ram_samples 501) against 500 Hz on the other days.
    # Days at a different rate cannot be coherently stacked -- their cross spectra
    # live on a different frequency grid -- so group by rate and stack the
    # dominant configuration only, reporting what that excludes.
    rates = {}
    for d in dates:
        with np.load(by_date[d][0], allow_pickle=False) as z:
            rates.setdefault((float(z["fs_hz"]), int(z["n_fft"])), []).append(d)
    dominant = max(rates, key=lambda k: len(rates[k]))
    for key, ds in sorted(rates.items()):
        mark = "STACKED" if key == dominant else "EXCLUDED, different acquisition rate"
        say("  fs %.0f Hz, n_fft %d : %s  -> %s" % (key[0], key[1], ", ".join(ds), mark))
    dates = sorted(rates[dominant])
    if not dates:
        raise SystemExit("no complete days")

    central = neighbor = source_power = None
    n_windows = 0
    core_files = 0
    ref = None
    for d in dates:
        for path in by_date[d]:
            with np.load(path, allow_pickle=False) as z:
                meta = {
                    "fs_hz": float(z["fs_hz"]),
                    "dx_m": float(z["dx_m"]),
                    "n_fft": int(z["n_fft"]),
                    "source_channel": int(z["source_channel"]),
                    "center_channels": z["center_channels"].tolist(),
                    "offsets_m": z["offsets_m"].tolist(),
                    "ram_samples": int(z["ram_samples"]),
                    "spectral_mode": str(z["spectral_mode"]),
                    "null_method": str(z["null_method"]),
                }
                if ref is None:
                    ref = meta
                    central = np.zeros_like(z["central_cross_spectrum_sum"])
                    neighbor = np.zeros_like(z["neighbor_cross_spectrum_sum"])
                elif meta != ref:
                    raise ValueError("Chunk metadata mismatch at %s" % path)
                central += z["central_cross_spectrum_sum"]
                neighbor += z["neighbor_cross_spectrum_sum"]
                n_windows += int(z["n_windows"])
                core_files += int(z["used_core_files"])

    fs = ref["fs_hz"]; n_fft = ref["n_fft"]
    offsets = np.asarray(ref["offsets_m"], dtype=float)
    say("")
    say("stacked %d days, %d files, %d windows (%.1f hours)"
        % (len(dates), core_files, n_windows, core_files / 60.0))

    (central_spec, neighbor_spec, spa, spd, floored) = prepare_aggregate_spectra(
        central, neighbor, source_power, n_windows, n_fft, fs,
        ref["spectral_mode"], 1e-3,
    )
    lags, _ = correlation_from_spectrum(central_spec, 1, n_fft, fs, apply_bandpass=True)
    _, neighbors = correlation_from_spectrum(neighbor_spec, 1, n_fft, fs, apply_bandpass=True)

    causal = moveout_scores(neighbors, lags, offsets, sign=1.0)
    acausal = moveout_scores(neighbors, lags, offsets, sign=-1.0)
    k = int(np.argmax(causal))
    peak = float(causal[k]); vpeak = float(VELOCITY_GRID_M_S[k])
    nulls, p = receiver_order_null(neighbors, lags, offsets, peak, SEED, NULL_COUNT)
    c32 = float(np.interp(V_REF, VELOCITY_GRID_M_S, causal))
    a32 = float(np.interp(V_REF, VELOCITY_GRID_M_S, acausal))

    say("")
    say("--- result ---")
    say("  peak causal score %.4f at %.0f m/s   (paper: a packet near 3200 m/s)" % (peak, vpeak))
    say("  at 3200 m/s: causal %.4f, acausal %.4f, ratio %.2f  (paper needs > 1)"
        % (c32, a32, c32 / a32 if a32 else float("nan")))
    say("  receiver-order familywise null 95th %.4f | observed %.4f | p = %.4f"
        % (float(np.percentile(nulls, 95)), peak, p))
    say("  detectability: peak / null95 = %.2f  (must reach 1.00)"
        % (peak / max(float(np.percentile(nulls, 95)), 1e-30)))
    say("")
    if p < 0.05:
        say("  VERDICT: the 3200 m/s arrival IS recovered by the coherent multi-day stack.")
    else:
        say("  VERDICT: not recovered. Coherently stacking %.0f hours -- the most"
            % (core_files / 60.0))
        say("  sensitive test this archive supports -- still does not clear the null.")

    fig, ax = plt.subplots(1, 2, figsize=(13, 7), constrained_layout=True)
    gn = neighbors / np.abs(neighbors).max(axis=1, keepdims=True)
    for o, row in zip(offsets, gn):
        ax[0].plot(lags, -o + row * 45, "k-", lw=0.8)
    ax[0].plot(offsets / V_REF, -offsets, "--", color="crimson", lw=1.6, label="3200 m/s (paper)")
    ax[0].set_xlim(-0.4, 0.7); ax[0].set_xlabel("Time [s]")
    ax[0].set_ylabel("Depth below wellhead [m]"); ax[0].legend(fontsize=8)
    ax[0].set_title("coherent %d-day stack, %d windows\nno F-K filter" % (len(dates), n_windows))
    ax[1].plot(VELOCITY_GRID_M_S, causal, "k-", label="causal")
    ax[1].plot(VELOCITY_GRID_M_S, acausal, color="gray", lw=0.9, ls="--", label="acausal")
    ax[1].axhline(np.percentile(nulls, 95), color="crimson", ls="--", lw=1, label="null 95th")
    ax[1].axvline(V_REF, color="steelblue", ls=":", label="3200 m/s")
    ax[1].set_xlabel("trial apparent velocity (m/s)"); ax[1].set_ylabel("moveout score")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    ax[1].set_title("scan vs receiver-order null, p = %.4f" % p)
    fig.savefig(str(STEM) + ".png", dpi=200)

    np.savez(str(STEM) + ".npz", lags=lags, section=neighbors, offsets_m=offsets,
             velocity_grid_m_s=VELOCITY_GRID_M_S, causal_moveout_scores=causal,
             acausal_moveout_scores=acausal, receiver_order_null_maxima=nulls,
             p_value=p, n_windows=n_windows, n_files=core_files,
             dates=np.array(dates))
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
