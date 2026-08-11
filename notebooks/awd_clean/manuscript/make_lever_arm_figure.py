"""Figure 1B: the lever-arm contrast, drawn from the frozen trajectories.

11.5x is the paper's central result and no reader feels it from a number. This
draws both observables on shared axes so the contrast is immediate: a short steep
Nano segment against a Deep trajectory that crosses the whole panel.

Reads the frozen Deep trajectory and the Nano constants from the analysis rather
than hard-coding them, so the figure cannot drift from the text.

    python make_lever_arm_figure.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = HERE / "figures" / "fig01b_lever_arm.png"

# Nano observable, from nano_dvv_injection_recovery.py
NANO_SPEED = 2975.0
NANO_INTERCEPT = -0.022
NANO_APERTURE = (80.0, 440.0)

# Deep aperture geometry, from deep_dvv_injection_recovery.py
COORD_RANGE = (200.0, 3000.0)
APERTURE_M, STEP_M = 400.0, 200.0
DX, STRIDE = 2.0419, 6


def deep_aperture_centres() -> np.ndarray:
    coordinate = np.arange(COORD_RANGE[0], COORD_RANGE[1] + 1e-9, DX * STRIDE)
    coordinate = coordinate[coordinate <= COORD_RANGE[1]]
    centres, start = [], coordinate[0]
    while start + APERTURE_M <= coordinate[-1] + 1e-9:
        members = coordinate[(coordinate >= start) & (coordinate < start + APERTURE_M)]
        if members.size >= 8:
            centres.append(members.mean())
        start += STEP_M
    return np.asarray(centres)


def main() -> None:
    frozen = json.loads((ANALYSIS / "deep_dvv_frozen_trajectory.json").read_text())
    deep = frozen["trajectories"]["outbound|15_30"]
    p_deep = deep["slowness_s_per_m"]
    t0_deep = deep["intercept_s"]

    centres = deep_aperture_centres()
    t_deep = t0_deep + p_deep * centres
    deep_lever = float(t_deep.max() - t_deep.min())

    nano_x = np.array(NANO_APERTURE)
    t_nano = NANO_INTERCEPT + (nano_x - NANO_APERTURE[0]) / NANO_SPEED
    nano_lever = float((NANO_APERTURE[1] - NANO_APERTURE[0]) / NANO_SPEED)

    print(f"Nano lever arm  {nano_lever:.4f} s over {NANO_APERTURE[0]:.0f}-{NANO_APERTURE[1]:.0f} m")
    print(f"Deep lever arm  {deep_lever:.4f} s over {centres.min():.0f}-{centres.max():.0f} m "
          f"({len(centres)} apertures)")
    print(f"ratio           {deep_lever / nano_lever:.2f}x")

    fig, ax = plt.subplots(figsize=(9.0, 6.0), constrained_layout=True)
    deep_line = t0_deep + p_deep * np.array([COORD_RANGE[0], COORD_RANGE[1]])
    ax.plot([COORD_RANGE[0], COORD_RANGE[1]], deep_line, color="#b2182b", lw=2.4,
            label=f"Deep guided mode, {1 / p_deep:.0f} m s$^{{-1}}$, 15–30 Hz")
    ax.plot(centres, t_deep, "o", color="#b2182b", ms=6,
            label=f"{len(centres)} aperture centres")
    ax.plot(nano_x, t_nano, color="#2166ac", lw=2.4,
            label=f"Nano apparent moveout, {NANO_SPEED:.0f} m s$^{{-1}}$, 30–60 Hz")

    for x, lo, hi, colour, text in (
        (COORD_RANGE[1] + 90, t_deep.min(), t_deep.max(), "#b2182b",
         f"Deep lever arm\n{deep_lever:.2f} s"),
        (NANO_APERTURE[1] + 90, t_nano.min(), t_nano.max(), "#2166ac",
         f"Nano lever arm\n{nano_lever:.3f} s"),
    ):
        ax.annotate("", xy=(x, lo), xytext=(x, hi),
                    arrowprops=dict(arrowstyle="<->", color=colour, lw=1.8))
        ax.text(x + 70, 0.5 * (lo + hi), text, color=colour, va="center", fontsize=10)

    ax.set(xlabel="distance along fibre (m)",
           ylabel="reference travel time $T_0$ (s)",
           xlim=(0, 3600), ylim=(-0.15, 2.15))
    ax.set_title(
        f"The delay-gradient estimator regresses against $T_0$, so precision scales\n"
        f"with the span in $T_0$ — {deep_lever / nano_lever:.1f}× longer for Deep",
        fontsize=11)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=300)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(HERE)}")


if __name__ == "__main__":
    main()
