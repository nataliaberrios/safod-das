"""Poster-grade figures, built from the frozen analysis outputs.

The analysis scripts write six-panel diagnostic figures. Those are right for a
methods appendix and unreadable on a poster at one metre. These are the same
numbers in poster form: one message per figure, large type, direct labels.

Three figures, in the order a visitor reads them:

    P1  does the measurement work        recovery vs injected
    P2  what is the answer               sensitivity of the three observables
    P3  why                              the lever-arm/timing trade (make_lever_arm_figure.py)

Palette is the validated categorical set: blue #2a78d6, orange #eb6834,
aqua #1baf7a. Every series is direct-labelled, which is also the relief the
aqua slot's contrast warning requires.

    python make_poster_figures.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = HERE / "figures"

NANO = "#2a78d6"
OUTB = "#eb6834"
RETN = "#1baf7a"
INK = "#0b0b0b"
MUTED = "#52514e"

# Poster type: readable at ~1 m.
plt.rcParams.update({
    "font.size": 17, "axes.titlesize": 21, "axes.labelsize": 19,
    "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 16,
    "axes.edgecolor": "#b8b7b2", "axes.linewidth": 1.2,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.labelcolor": INK,
})


def summary_rows(leg: str) -> list[dict]:
    rows = [
        r for r in csv.DictReader(open(ANALYSIS / "deep_dvv_summary.csv"))
        if r["population"] == "heldout" and r["leg"] == leg
        and r["band"] == "15_30" and r["pass"] == "primary"
    ]
    return sorted(rows, key=lambda r: float(r["injected_dvv"]))


def figure_p1() -> None:
    """Recovery vs injected — the measurement works, and here is how small."""
    fig, ax = plt.subplots(figsize=(11.0, 8.5), constrained_layout=True)
    limit = 1.15

    ax.plot([-limit, limit], [-limit, limit], color="#9a9a95", lw=1.6, ls=(0, (6, 4)),
            zorder=1)
    ax.text(0.78, 0.92, "perfect recovery", color=MUTED, fontsize=16,
            rotation=41, ha="center", va="center")

    # One series only. Both legs recover almost identically, so plotting both
    # just hides one under the other; the leg comparison lives in P2.
    rows = summary_rows("outbound")
    x = np.array([float(r["injected_dvv"]) for r in rows]) * 100
    y = np.array([float(r["median_estimated_dvv"]) for r in rows]) * 100
    e = np.array([float(r["robust_scatter_1p4826_mad"]) for r in rows]) * 100
    ax.errorbar(x, y, yerr=e, fmt="o", ms=13, color=OUTB, capsize=0,
                elinewidth=2.6, mec="white", mew=2.0, zorder=3)

    ax.axhspan(-0.184, 0.184, color=OUTB, alpha=0.10, zorder=0)
    ax.text(-limit * 0.94, 0.235, "noise floor  ±0.18%",
            color=OUTB, fontsize=17, fontweight="bold", va="bottom")
    ax.text(-limit * 0.94, -limit * 0.78,
            "Deep outbound branch\n23 held-out bursts\nrecovery code never saw\nthe injected value",
            color=MUTED, fontsize=16, va="bottom", linespacing=1.5)

    ax.set(xlim=(-limit, limit), ylim=(-limit, limit),
           xlabel="injected velocity change (%)",
           ylabel="recovered velocity change (%)")
    ax.set_title("Known changes injected into real data are recovered blind\n"
                 "down to a few tenths of a percent", pad=14)
    ax.grid(alpha=0.18, lw=0.9)
    ax.set_aspect("equal")
    fig.savefig(OUT / "poster_p1_recovery.png", dpi=200)
    plt.close(fig)
    print("  poster_p1_recovery.png")


def figure_p2() -> None:
    """The answer: smallest reliably detected change, three observables."""
    comparison = list(csv.DictReader(open(ANALYSIS / "deep_dvv_nano_comparison.csv")))
    nano = next(r for r in comparison if r["observable"].startswith("Nano"))

    def row(prefix: str, leg: str) -> dict:
        return next(r for r in comparison
                    if r["observable"].startswith(prefix) and r["leg"] == leg
                    and r["band_hz"] == "15-30" and r["population"] == "heldout")

    def pair(r: dict) -> tuple[float, float]:
        return (float(r["null_threshold_dvv"]) * 100,
                float(r["reliable_level_detection"]) * 100)

    entries = [
        ("Nano\ncemented fiber", *pair(nano), NANO),
        ("Deep outbound\nwireline", *pair(row("Deep", "outbound")), OUTB),
        ("Deep return\nwireline", *pair(row("Deep", "return")), RETN),
    ]

    fig, ax = plt.subplots(figsize=(12.5, 5.0), constrained_layout=True)
    ypos = np.arange(len(entries))[::-1]

    for y, (name, thresh, reliable, colour) in zip(ypos, entries):
        ax.plot([thresh, reliable], [y, y], color=colour, lw=5, solid_capstyle="round",
                alpha=0.35, zorder=2)
        ax.plot(thresh, y, "o", ms=15, color=colour, mec="white", mew=2.2, zorder=3)
        ax.plot(reliable, y, "D", ms=15, color=colour, mec="white", mew=2.2, zorder=3)
        ax.text(reliable * 1.14, y, f"{reliable:g}%", color=colour, fontsize=21,
                fontweight="bold", va="center")
        ax.text(thresh * 0.86, y, f"{thresh:.2g}%", color=MUTED, fontsize=16,
                va="center", ha="right")

    ax.set_yticks(ypos)
    ax.set_yticklabels([e[0] for e in entries], fontsize=18)
    ax.set_xscale("log")
    ax.set_xlim(0.10, 2.4)
    ax.set_xticks([0.1, 0.2, 0.5, 1.0, 2.0])
    ax.set_xticklabels(["0.1%", "0.2%", "0.5%", "1%", "2%"])
    ax.set_xlabel("fractional change in apparent velocity")
    ax.set_title("The outbound wireline branch resolves a change twice as small\n"
                 "circle = noise floor      diamond = smallest reliably detected",
                 fontsize=20, pad=12)
    ax.set_ylim(-0.55, len(entries) - 0.45)
    ax.grid(axis="x", alpha=0.18, lw=0.9)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    fig.savefig(OUT / "poster_p2_sensitivity.png", dpi=200)
    plt.close(fig)
    print("  poster_p2_sensitivity.png")


def main() -> None:
    OUT.mkdir(exist_ok=True)
    print("poster figures:")
    figure_p1()
    figure_p2()
    print("  (P3 = fig01b_lever_arm.png, from make_lever_arm_figure.py)")


if __name__ == "__main__":
    main()
