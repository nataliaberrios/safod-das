#!/usr/bin/env python3
"""Pooled multi-day Figure 7c stack, scored by the exact-stack code itself.

WHY THIS IS THE TEST THAT WAS STILL MISSING. On 2024-12-20 the paper baseline
reaches 0.97 of its own 95 % detection threshold — right at the line. Chunk files
store *summed cross spectra*, so days are exactly additive: pooling N comparable
days multiplies the coherent sum by N while incoherent noise grows as sqrt(N),
roughly a sqrt(N) gain. Four to five days is therefore a ~2x SNR increase over the
single day already tested, which is precisely the factor that could carry 0.97
past 1.0. Every earlier "more data" argument in this thread was run on a defective
pipeline; this one runs on the faithful operator.

Lellouch et al. (2019) themselves used one day for Figure 7c and seven days for
uncertainty, so pooling is consistent with the paper's own design rather than an
escalation beyond it.

NO REIMPLEMENTATION. The pooled spectra are handed to the same functions the
single-day aggregate uses — `prepare_aggregate_spectra`, `correlation_from_spectrum`
(with `apply_bandpass=True`, so the 5–20 Hz band lands on the stacked correlations),
`moveout_scores`, and `receiver_order_null`. The only new code is the summation and
the reporting, so a difference from the single-day numbers cannot come from a
different scorer.

Only configuration 0 (paper baseline: wellhead source channel 23, RAM 0.1 s,
cross-correlation spectral mode, ordered) is pooled, and only exactly continuous
days are eligible. No F-K filter anywhere.

Outputs: ambient_lellouch2019_pooled_days.{npz,png,txt}
"""
from __future__ import annotations

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

STEM = HERE / "ambient_lellouch2019_pooled_days"
MAIN = HERE / "ambient_transfer" / "lellouch2019_exact_stack"
DAYS = HERE / "ambient_transfer" / "lellouch2019_exact_stack_days"
V_REF = 3200.0
NULL_COUNT = 10000

# (date, directory) for configuration 0 only. Ordered so the published day is
# first and each cumulative stack is a superset of the one before it.
SOURCES = [
    ("2024-12-20", MAIN),
    ("2024-06-17", DAYS),
    ("2024-06-26", DAYS),
    ("2025-02-24", DAYS),
    ("2024-05-11", DAYS),
]
PATTERN = "chunk_{date}_src23_ram0p1_cross_correlation_ordered_r0_start*_n*.npz"


def load_day(date: str, directory: Path):
    """Sum one day's chunk cross spectra. Returns None if the day is incomplete."""
    files = sorted(directory.glob(PATTERN.format(date=date)))
    files = [f for f in files if "_cm_" not in f.name]
    if len(files) != 24:
        return None
    central = neighbor = None
    n_windows = 0
    core_files = 0
    meta = None
    for path in files:
        with np.load(path, allow_pickle=False) as d:
            this = (float(d["fs_hz"]), int(d["n_fft"]), int(d["source_channel"]),
                    tuple(np.asarray(d["offsets_m"]).tolist()), float(d["ram_seconds"]),
                    str(d["spectral_mode"]), bool(d["common_mode"]))
            if meta is None:
                meta = this
                central = np.zeros_like(d["central_cross_spectrum_sum"])
                neighbor = np.zeros_like(d["neighbor_cross_spectrum_sum"])
            elif this != meta:
                raise ValueError(f"metadata mismatch within {date} at {path.name}")
            central += d["central_cross_spectrum_sum"]
            neighbor += d["neighbor_cross_spectrum_sum"]
            n_windows += int(d["n_windows"])
            core_files += int(d["used_core_files"])
    return central, neighbor, n_windows, core_files, meta


def score(central, neighbor, n_windows, meta, seed):
    fs, n_fft, _src, offsets, _ram, spectral_mode, _cm = meta
    offsets = np.asarray(offsets, dtype=float)
    central_spec, neighbor_spec, _spa, _spd, _ff = prepare_aggregate_spectra(
        central, neighbor, None, n_windows, n_fft, fs, spectral_mode, 1e-3,
    )
    lags, _central_filtered = correlation_from_spectrum(
        central_spec, 1, n_fft, fs, apply_bandpass=True
    )
    _, neighbors = correlation_from_spectrum(
        neighbor_spec, 1, n_fft, fs, apply_bandpass=True
    )
    causal = moveout_scores(neighbors, lags, offsets, sign=1.0)
    acausal = moveout_scores(neighbors, lags, offsets, sign=-1.0)
    k = int(np.argmax(causal))
    peak = float(causal[k])
    nulls, p = receiver_order_null(neighbors, lags, offsets, peak, seed, NULL_COUNT)
    return dict(
        lags=lags, section=neighbors, offsets=offsets, causal=causal, acausal=acausal,
        peak=peak, v_peak=float(VELOCITY_GRID_M_S[k]),
        c32=float(np.interp(V_REF, VELOCITY_GRID_M_S, causal)),
        a32=float(np.interp(V_REF, VELOCITY_GRID_M_S, acausal)),
        null95=float(np.percentile(nulls, 95)), p=float(p), n_windows=n_windows,
    )


def main():
    log = []

    def say(text=""):
        print(text, flush=True)
        log.append(text)

    say("Pooled multi-day Figure 7c stack -- paper baseline (configuration 0)")
    say("scored by the exact-stack module's own functions; no F-K filter")
    say()

    loaded = []
    for date, directory in SOURCES:
        day = load_day(date, directory)
        if day is None:
            say(f"  {date}: incomplete (not 24 chunks) -- excluded")
            continue
        loaded.append((date, day))
        say(f"  {date}: 24 chunks, {day[2]} windows, {day[3]} files")
    if not loaded:
        raise SystemExit("no complete days")

    say()
    say("--- cumulative pooled stack ---")
    say("  days  windows   peak    at(m/s)  c/ac@3200   null95   peak/null95      p")
    central = neighbor = None
    n_windows = 0
    meta = None
    rows = []
    for index, (date, (c, nb, nw, _cf, m)) in enumerate(loaded, start=1):
        if central is None:
            central = np.zeros_like(c)
            neighbor = np.zeros_like(nb)
            meta = m
        elif m != meta:
            say(f"  {date}: metadata differs from the first day -- stopping pool here")
            break
        central += c
        neighbor += nb
        n_windows += nw
        r = score(central, neighbor, n_windows, meta, seed=20260814 + index)
        rows.append((index, date, r))
        ratio = r["c32"] / r["a32"] if r["a32"] else float("nan")
        say("  %4d %8d %7.3f %8.0f %10.2f %8.3f %12.2f %7.4f"
            % (index, r["n_windows"], r["peak"], r["v_peak"], ratio,
               r["null95"], r["peak"] / r["null95"], r["p"]))

    final = rows[-1][2]
    say()
    say("--- verdict ---")
    if final["p"] < 0.05:
        say("  The pooled stack CLEARS alpha = 0.05: p = %.4f, peak %.3f at %.0f m/s."
            % (final["p"], final["peak"], final["v_peak"]))
        say("  Check the peak velocity against the paper's ~3200 m/s before claiming")
        say("  reproduction: a significant peak at the edge of the scan is not Figure 7c.")
    else:
        best = min(r[2]["p"] for r in rows)
        top = max(r[2]["peak"] / r[2]["null95"] for r in rows)
        say("  No cumulative stack clears alpha = 0.05. Best p = %.4f; the statistic"
            % best)
        say("  reaches at most %.2f of its own 95%% threshold." % top)
        say("  Pooling %d days does not carry the 3,200 m/s packet into detection."
            % len(rows))

    fig, ax = plt.subplots(1, 3, figsize=(17, 6), constrained_layout=True)
    section = final["section"]
    lags = final["lags"]
    offsets = final["offsets"]
    norm = section / np.abs(section).max(axis=1, keepdims=True)
    for offset, trace in zip(offsets, norm):
        ax[0].plot(lags, -offset + trace * 45.0, "k-", lw=0.8)
    ax[0].plot(offsets / V_REF, -offsets, "--", color="crimson", lw=1.6,
               label="3200 m/s (paper)")
    ax[0].set_xlim(-0.4, 0.7)
    ax[0].set_xlabel("Time [s]")
    ax[0].set_ylabel("Depth below wellhead [m]")
    ax[0].legend(fontsize=8)
    ax[0].set_title("pooled %d days, %d windows\nno F-K filter" % (len(rows), final["n_windows"]))

    ax[1].plot(VELOCITY_GRID_M_S, final["causal"], "k-", label="causal")
    ax[1].plot(VELOCITY_GRID_M_S, final["acausal"], color="grey", lw=1, ls="--",
               label="acausal")
    ax[1].axhline(final["null95"], color="crimson", ls="--", lw=1, label="null 95th")
    ax[1].axvline(V_REF, color="steelblue", ls=":", label="3200 m/s")
    ax[1].set_xlabel("trial apparent velocity (m/s)")
    ax[1].set_ylabel("median normalized envelope")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)
    ax[1].set_title("pooled moveout scan")

    ax[2].plot([r[0] for r in rows], [r[2]["peak"] / r[2]["null95"] for r in rows],
               "o-", color="k")
    ax[2].axhline(1.0, color="crimson", ls="--", lw=1.5, label="detection threshold")
    ax[2].set_ylim(0, 1.3)
    ax[2].set_xlabel("days pooled")
    ax[2].set_ylabel("peak / null 95th percentile")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3)
    ax[2].set_title("does pooling reach detection?")
    fig.savefig(str(STEM) + ".png", dpi=200)

    np.savez(str(STEM) + ".npz",
             dates=np.array([r[1] for r in rows]),
             n_days=np.array([r[0] for r in rows]),
             peak=np.array([r[2]["peak"] for r in rows]),
             v_peak=np.array([r[2]["v_peak"] for r in rows]),
             null95=np.array([r[2]["null95"] for r in rows]),
             p=np.array([r[2]["p"] for r in rows]),
             n_windows=np.array([r[2]["n_windows"] for r in rows]),
             lags=lags, section=section, offsets=offsets,
             causal=final["causal"], acausal=final["acausal"],
             velocity_grid_m_s=VELOCITY_GRID_M_S)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say()
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
