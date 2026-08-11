#!/usr/bin/env python3
"""Plot the completed seasonal ambient negative and positive signed F-K branches."""
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE / "ambient_transfer"


def load_branch(name):
    product = np.load(OUT / f"fk_seasonal_aggregate_{name}.npz")
    return {key: product[key] for key in product.files}


def main():
    negative = load_branch("negative")
    positive = load_branch("positive")
    report = json.loads((OUT / "fk_seasonal_aggregate.json").read_text())
    metrics = report["modes"]

    common = np.concatenate([negative["top"].ravel(), positive["top"].ravel()])
    limit = float(np.nanpercentile(np.abs(common), 98.5))

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), constrained_layout=True)
    for ax, name, data, moveout_sign in (
        (axes[0, 0], "Negative signed F-K branch", negative, 1.0),
        (axes[0, 1], "Positive signed F-K branch", positive, -1.0),
    ):
        image = ax.imshow(
            data["top"],
            extent=[data["lags"][0], data["lags"][-1], data["dist"][-1], data["dist"][0]],
            aspect="auto",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            interpolation="nearest",
        )
        # Signed 3.2-km/s physical moveout estimate, drawn from channel 0.
        moveout_distance = np.r_[0.0, data["dist"]]
        ax.plot(moveout_sign * moveout_distance / 3200.0, moveout_distance,
                "k--", lw=1.2, label="signed 3.2 km/s estimate")
        ax.scatter(0.0, 0.0, marker="*", s=110, color="#ffd84d",
                   edgecolor="black", linewidth=0.8, zorder=5, clip_on=False)
        ax.annotate("virtual source (channel 0)", xy=(0.0, 0.0),
                    xytext=(8, 10), textcoords="offset points", fontsize=8)
        mode = "negative" if "negative" in name.lower() else "positive"
        result = metrics[mode]
        if mode == "negative":
            ax.set_title(
                f"{name}: peak {result['peak_v_mps']/1000:.3f} km/s, "
                f"p={result['p_peak']:.4g}"
            )
        else:
            ax.set_title(f"{name}: opposite-slope estimate; signed-lag rerun pending")
        ax.set_xlim(-0.25, 0.25)
        ax.set_ylim(700, -15)
        ax.set_xlabel("Correlation lag (s)")
        ax.set_ylabel("Receiver position from channel 0 (m)")
        ax.grid(alpha=0.15, lw=0.4)

    cbar = fig.colorbar(image, ax=axes[0, :].tolist(), shrink=0.9, pad=0.02)
    cbar.set_label("Normalized correlation (common display scale)")

    ax = axes[1, 0]
    ax.plot(negative["vv"] / 1000.0, negative["scores"], color="tab:blue", lw=2, label="Negative branch (+lag target)")
    ax.plot(positive["vv"] / 1000.0, positive["scores"], color="tab:orange", lw=2, label="Positive-mask leakage (+lag target)")
    ax.axhline(metrics["negative"]["null95"], color="tab:blue", ls="--", lw=1, label="Negative null 95%")
    ax.axhline(metrics["positive"]["null95"], color="tab:orange", ls="--", lw=1, label="Positive null 95%")
    ax.axvline(3.2, color="0.25", ls=":", lw=1.2, label="3.2 km/s")
    ax.set_xlabel("Trial apparent velocity (km/s)")
    ax.set_ylabel("Median normalized correlation")
    ax.set_title("Velocity-score comparison")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = axes[1, 1]
    names = ["Negative target", "Positive-mask\nleakage"]
    values = [metrics["negative"]["peak_score"], metrics["positive"]["peak_score"]]
    nulls = [metrics["negative"]["null95"], metrics["positive"]["null95"]]
    x = np.arange(2)
    ax.bar(x, values, color=["tab:blue", "tab:orange"], alpha=0.8, label="Observed peak")
    ax.scatter(x, nulls, color="black", marker="_", s=450, linewidths=2, label="Permutation null 95%")
    for index, mode in enumerate(("negative", "positive")):
        ax.text(index, values[index] + 0.012, f"p={metrics[mode]['p_peak']:.4g}", ha="center", fontsize=9)
    ax.set_xticks(x, names)
    ax.set_ylabel("Peak velocity score")
    ax.set_title("Observed peak versus receiver-permutation null")
    ax.set_ylim(0, max(values) * 1.28)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle(
        "Eight-day ambient signed F-K branch comparison "
        f"({report['weighted_files']:,} weighted one-minute records)",
        fontsize=14,
    )
    fig.savefig(OUT / "fk_seasonal_signed_branch_comparison.png", dpi=350)
    fig.savefig(OUT / "fk_seasonal_signed_branch_comparison.pdf")


if __name__ == "__main__":
    main()
