import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


DEFAULT_EVENT_FILE = (
    "/oak/stanford/groups/ettore88/nberrios/event_cc/"
    "event20_base/EQ_20240525_045628_event20s.npz"
)
DEFAULT_PROCESSING_VERSION = "event20_base"
DEFAULT_FIG_DIR = "/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a single event cross-correlation gather from an event .npz file."
    )
    parser.add_argument("--event-file", default=DEFAULT_EVENT_FILE)
    parser.add_argument("--processing-version", default=DEFAULT_PROCESSING_VERSION)
    parser.add_argument("--fig-dir", default=DEFAULT_FIG_DIR)
    parser.add_argument("--vp", type=float, default=3200.0)
    parser.add_argument("--vs", type=float, default=1600.0)
    parser.add_argument("--dz-m", type=float, default=1.0)
    parser.add_argument("--xlim-show", type=float, default=1.0)
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def build_output_path(fig_dir: str, processing_version: str, event_file: str) -> str:
    stem = Path(event_file).stem
    return os.path.join(fig_dir, f"{processing_version}_{stem}.png")


def format_time_label(value: str) -> str:
    text = str(value).replace("T", " ")
    text = text.replace(".150000+00:00", "")
    text = text.replace("+00:00", " UTC")
    return text


def main() -> None:
    args = parse_args()

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
        }
    )

    d = np.load(args.event_file, allow_pickle=True)

    cc = d["cc_combo"]
    max_lag = float(d["max_lag"])
    ch_start = int(d["ch_start"])
    ch_end = int(d["ch_end"])
    src = int(d["src"])
    nstack = int(d["nstack_combo"])

    channels = np.arange(ch_start, ch_end)
    lags = np.linspace(-max_lag, max_lag, cc.shape[1])
    dist = np.abs(channels - src) * args.dz_m

    event_time = str(d["event_time"]) if "event_time" in d.files else "unknown"
    snippet_start = str(d["snippet_start"]) if "snippet_start" in d.files else "unknown"
    snippet_end = str(d["snippet_end"]) if "snippet_end" in d.files else "unknown"

    mask = np.abs(channels - src) > 10
    vsrc = np.abs(cc[mask]) if np.any(mask) else np.abs(cc)
    vlim = np.percentile(vsrc, 99) if vsrc.size else 1.0
    if vlim == 0:
        vlim = 1.0

    fig, ax = plt.subplots(figsize=(8.8, 6.8))

    im = ax.imshow(
        cc,
        extent=[lags[0], lags[-1], ch_start, ch_end],
        origin="lower",
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vlim,
        vmax=vlim,
        interpolation="nearest",
    )

    ax.plot(+dist / args.vp, channels, "--", color="k", lw=0.9, zorder=5, label=f"Vp = {args.vp:.0f} m/s")
    ax.plot(-dist / args.vp, channels, "--", color="k", lw=0.8, zorder=5)
    ax.plot(+dist / args.vs, channels, "--", color="gray", lw=0.9, zorder=5, label=f"Vs = {args.vs:.0f} m/s")
    ax.plot(-dist / args.vs, channels, "--", color="gray", lw=0.8, zorder=5)

    ax.axhline(src, color="lime", lw=0.8, ls=":")
    ax.set_xlim(-args.xlim_show, args.xlim_show)
    ax.set_xlabel("Lag (s)")
    ax.set_ylabel("Channel")
    ax.legend(loc="upper right", fontsize=8, frameon=True, borderpad=0.4, handlelength=2.4)

    fig.suptitle(
        f"{args.processing_version}   |   Nstack = {nstack}",
        y=0.975,
        fontsize=12,
        fontweight="semibold",
    )
    fig.text(
        0.125,
        0.94,
        f"Event:   {format_time_label(event_time)}\n"
        f"Snippet: {format_time_label(snippet_start)} to {format_time_label(snippet_end)}",
        ha="left",
        va="top",
        fontsize=9,
    )

    cbar = fig.colorbar(im, ax=ax, pad=0.03)
    cbar.set_label("CC amplitude")
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 2))
    cbar.formatter = formatter
    cbar.update_ticks()

    os.makedirs(args.fig_dir, exist_ok=True)
    out_png = build_output_path(args.fig_dir, args.processing_version, args.event_file)

    fig.tight_layout(rect=[0, 0, 1, 0.89])
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight")

    if args.show:
        plt.show()
    else:
        plt.close(fig)

    print(f"Saved {out_png}")


if __name__ == "__main__":
    main()
