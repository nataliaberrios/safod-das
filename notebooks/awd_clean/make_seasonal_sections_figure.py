#!/usr/bin/env python3
"""Create a useful seasonal figure showing the actual correlation sections."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE / "ambient_transfer"
DATES = [
    "2024-05-11", "2024-06-17", "2024-06-26", "2024-10-28",
    "2024-11-30", "2024-12-20", "2025-02-24", "2025-03-04",
]


def load(path):
    z = np.load(path)
    return z["top"], z["lags"], z["dist"]


def main():
    products = [(date, *load(OUT / f"fk_seasonal_{date}_negative.npz")) for date in DATES]
    aggregate = load(OUT / "fk_seasonal_aggregate_negative.npz")
    all_top = np.concatenate([top.ravel() for _, top, _, _ in products] + [aggregate[0].ravel()])
    limit = float(np.nanpercentile(np.abs(all_top), 98.5))
    limit = min(max(limit, 0.18), 0.35)

    fig, axes = plt.subplots(3, 3, figsize=(14, 11), constrained_layout=True, sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, (date, top, lags, dist) in zip(axes[:8], products):
        im = ax.imshow(
            top,
            extent=[lags[0], lags[-1], dist[-1], dist[0]],
            aspect="auto",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            interpolation="nearest",
        )
        # Black dashed is the Lellouch comparison; gold is the conditional ridge.
        ax.plot(dist / 3200.0, dist, "k--", lw=1.0, alpha=0.9)
        ax.plot(dist / 3075.0, dist, color="gold", lw=1.0, alpha=0.95)
        ax.set_title(date)
        ax.set_xlim(-0.25, 0.25)
        ax.set_ylim(700, 0)
        ax.grid(alpha=0.15, lw=0.4)
    top, lags, dist = aggregate
    ax = axes[8]
    im = ax.imshow(
        top,
        extent=[lags[0], lags[-1], dist[-1], dist[0]],
        aspect="auto",
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
        interpolation="nearest",
    )
    ax.plot(dist / 3200.0, dist, "k--", lw=1.2, alpha=0.9,
            label="3.2 km/s reference")
    ax.plot(dist / 3075.0, dist, color="gold", lw=1.2, alpha=0.95,
            label="3.075 km/s conditional ridge")
    ax.set_title("Across-day weighted aggregate")
    ax.set_xlim(-0.25, 0.25)
    ax.set_ylim(700, 0)
    ax.grid(alpha=0.15, lw=0.4)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.85)

    for ax in axes[6:]:
        ax.set_xlabel("Correlation lag (s)")
    for ax in axes[::3]:
        ax.set_ylabel("Receiver position from channel 0 (m)")
    cbar = fig.colorbar(im, ax=axes.tolist(), shrink=0.82, pad=0.015)
    cbar.set_label("Normalized correlation (display-clipped)")
    fig.suptitle(
        "Seasonal repeatability of the negative signed F–K correlation section",
        fontsize=15,
    )
    fig.text(
        0.5,
        0.008,
        "Black dashed: 3.2 km/s Lellouch comparison reference; gold: 3.075 km/s "
        "conditional ridge estimate. Panels use identical display limits; "
        "the F–K wedge is fixed and includes 2.5–4.5 km/s by construction.",
        ha="center",
        fontsize=9,
    )
    fig.savefig(OUT / "fk_seasonal_day_sections.png", dpi=350)
    fig.savefig(OUT / "fk_seasonal_day_sections.pdf")
    print("wrote", OUT / "fk_seasonal_day_sections.png")
    print("wrote", OUT / "fk_seasonal_day_sections.pdf")


if __name__ == "__main__":
    main()
