"""Plot the semblance surface that chose the Deep trajectory.

The freeze stage of deep_dvv_injection_recovery.py scans start time against
slowness, takes the peak, writes two numbers to deep_dvv_frozen_trajectory.json
and throws the surface away. So the single most-asked question about this
analysis — where did 1544.6 m/s come from? — has no figure answering it.

This recomputes that surface using the analysis's own functions and constants,
plots it, and asserts the peak reproduces the frozen JSON. It changes nothing:
the frozen trajectory is read, not re-selected.

    python make_semblance_figure.py    ->  figures/fig03b_semblance_scan.png
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
sys.path.insert(0, str(ANALYSIS))

import deep_dvv_injection_recovery as D   # noqa: E402
from fk_dispersion import weighted_stack  # noqa: E402

OUT = HERE / "figures" / "fig03b_semblance_scan.png"
BAND = D.PRIMARY_BAND
ACCENT = "#eb6834"
INK, MUTED = "#0b0b0b", "#52514e"

plt.rcParams.update({
    "font.size": 13, "axes.titlesize": 15, "axes.labelsize": 13.5,
    "xtick.labelsize": 12, "ytick.labelsize": 12,
    "text.color": INK, "axes.labelcolor": INK,
})


def surface(section, coordinate, fs, config):
    """Recompute the freeze-stage scan, keeping the whole surface."""
    filtered = D._rms_normalize(D._bandpass(section, fs, config["band"]))
    sample_time = np.arange(filtered.shape[-1], dtype=float) / fs - D.PRE_S
    intercepts = np.arange(
        D.INTERCEPT_GRID_S[0],
        D.INTERCEPT_GRID_S[1] + 0.5 * D.INTERCEPT_STEP_S,
        D.INTERCEPT_STEP_S,
    )
    pad = config["semblance_half_s"] + 2.0 / fs
    q = np.arange(
        int(round(((intercepts[-1] + pad) - (intercepts[0] - pad)) * fs)) + 1,
        dtype=float,
    ) / fs + (intercepts[0] - pad)
    slowness = np.linspace(D.SLOWNESS_RANGE[0], D.SLOWNESS_RANGE[1], D.N_SLOWNESS)

    grid = np.empty((slowness.size, intercepts.size), dtype=float)
    for i, p in enumerate(slowness):
        aligned = D._align(filtered, sample_time, coordinate, 0.0, p, q)
        grid[i] = D._semblance_profile(
            aligned, q, fs, config["semblance_half_s"], intercepts
        )
    return grid, slowness, intercepts


def main() -> None:
    frozen = json.loads((ANALYSIS / "deep_dvv_frozen_trajectory.json").read_text())
    config = D.BAND_CONFIG[BAND]

    with np.load(D.STACKS) as data:
        counts = np.asarray(data["n_common"], dtype=int)
        fs = float(data["fs"])
        dx = float(data["dx_deep"])
        valid = np.flatnonzero(counts > 0)
        discovery_counts = counts.copy()
        discovery_counts[valid[valid % 2 == 0]] = 0     # keep epoch % 2 == 1
        discovery = weighted_stack(data["deep_stacks"], discovery_counts)
        n_channels = discovery.shape[0]

    channels = D._leg_channels(n_channels, dx)

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 6.0), constrained_layout=True)
    for ax, leg in zip(axes, D.LEGS):
        absolute, coordinate = channels[leg]
        grid, slowness, intercepts = surface(discovery[absolute], coordinate, fs, config)
        speed = 1.0 / slowness

        peak = np.unravel_index(np.argmax(grid), grid.shape)
        found_speed, found_t0 = speed[peak[0]], intercepts[peak[1]]
        expect = frozen["trajectories"][f"{leg}|{BAND}"]
        assert abs(found_speed - expect["velocity_mps"]) < 1.0, (
            f"{leg}: recomputed {found_speed:.1f} m/s vs frozen "
            f"{expect['velocity_mps']:.1f} — the scan no longer reproduces the freeze"
        )
        print(f"{leg:9s} peak {found_speed:7.1f} m/s at t0 {found_t0:+.3f} s  "
              f"(frozen {expect['velocity_mps']:.1f}, {expect['intercept_s']:+.3f}) "
              f"semblance {grid[peak]:.4f}")

        im = ax.pcolormesh(intercepts, speed, grid, cmap="Blues", shading="auto")
        ax.plot(found_t0, found_speed, "o", ms=16, mfc="none", mec=ACCENT, mew=3.0)
        ax.annotate(f"{found_speed:.0f} m s$^{{-1}}$\n{found_t0:+.3f} s",
                    xy=(found_t0, found_speed),
                    xytext=(found_t0 + 0.075, found_speed + 95),
                    color=ACCENT, fontsize=13, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=ACCENT, lw=1.6))
        ax.set(xlabel="start time $t_0$ (s)",
               ylabel="apparent speed (m s$^{-1}$)" if leg == "outbound" else "",
               title=f"{leg.capitalize()} leg")
        plt.colorbar(im, ax=ax, label="semblance")

    fig.suptitle(
        f"Where the Deep trajectory comes from: semblance over start time and speed, "
        f"{BAND.replace('_', '–')} Hz\n"
        "discovery half only (23 odd-indexed bursts) — the peak is the frozen "
        "trajectory, circled",
        fontsize=14)
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=200)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(HERE)}")


if __name__ == "__main__":
    main()
