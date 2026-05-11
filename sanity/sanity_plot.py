"""
Post-processing and visualization for the canonical Lellouch CC pipeline.

Loads a daily CC npz produced by sanity_cc.py, applies the Lellouch et al. 2019
adjacent-channel pre-shifted stacking step (R±10 with travel-time shift at the
average velocity), and produces three plots:

  1. cc_raw.png     : the raw (un-pair-stacked) CC vs channel/depth, full 5-20 Hz band.
                      This is what comes out of sanity_cc.py directly. The 3200 m/s
                      reference moveout is overlaid.
  2. cc_stacked.png : the Lellouch Fig 7c equivalent — adjacent-channel pre-shifted
                      stack at V_AVG, shown at the same band. If anything looks like
                      Fig 7c, it will be visible here.
  3. cc_fk.png      : F-K spectrum of the stacked CC. Tells us what apparent velocities
                      are present in the recovered correlations.

This script is intentionally a re-runnable post-processing pass — it does not need
the cluster, can run anywhere with the daily npz in hand, and you can iterate on
plotting choices without re-running CC.
"""
import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- USER PARAMETERS ----------------
DEFAULT_NPZ = os.environ.get("SANITY_NPZ", "")
DEFAULT_NPZ_DIR = os.environ.get(
    "SANITY_OUT",
    "/oak/stanford/groups/ettore88/nberrios/sanity_v1",
)

V_AVG       = float(os.environ.get("SANITY_VAVG", "3200.0"))   # Lellouch's pre-shift velocity
PAIR_HALF   = int(os.environ.get("SANITY_PAIR_HALF", "10"))    # R±10 -> 21-channel stack
DZ_M        = float(os.environ.get("SANITY_DZ", "1.0"))        # channel spacing in meters

OUT_DIR = os.environ.get(
    "SANITY_PLOT_OUT",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_out"),
)
os.makedirs(OUT_DIR, exist_ok=True)

XLIM = (-0.4, 0.4)   # Lellouch Fig 7c is roughly ±0.3 s
# -------------------------------------------------


def resolve_npz_path():
    if len(sys.argv) > 1:
        return sys.argv[1]
    if DEFAULT_NPZ:
        return DEFAULT_NPZ

    candidates = sorted(
        glob.glob(os.path.join(DEFAULT_NPZ_DIR, "sanity_cc_*.npz")),
        key=os.path.getmtime,
    )
    if not candidates:
        raise FileNotFoundError(
            f"No sanity_cc_*.npz files found in {DEFAULT_NPZ_DIR}. "
            "Set SANITY_NPZ or pass the npz path as an argument."
        )
    return candidates[-1]


def adjacent_pair_preshift_stack(cc, channels, isource_idx, dz_m, v_avg, pair_half, dt):
    """
    For each receiver R, sum C_{S,Z} over Z in [R-pair_half, R+pair_half] after
    shifting each by (Z-R)*dz_m / v_avg seconds.

    cc           : (nch, nlag) CC array
    channels     : (nch,) channel indices
    isource_idx  : index of the source channel in cc[ich, :] (relative to channels[0])
    dz_m         : channel spacing in meters
    v_avg        : average velocity for moveout pre-shift (m/s)
    pair_half    : half-width of the receiver stack (e.g., 10 -> R±10 = 21 channels)
    dt           : time step (s)

    Returns
    -------
    cc_stack : same shape as cc, but each row R is the pre-shifted sum over its window.
    """
    nch, nlag = cc.shape
    cc_stack = np.zeros_like(cc)
    nlag_half = nlag // 2
    for iR in range(nch):
        z0 = max(0, iR - pair_half)
        z1 = min(nch, iR + pair_half + 1)
        acc = np.zeros(nlag)
        for iZ in range(z0, z1):
            shift_s     = (iZ - iR) * dz_m / v_avg
            shift_samps = int(round(shift_s / dt))
            trace = cc[iZ, :]
            shifted = np.roll(trace, shift_samps)
            # Zero out the wraparound region
            if shift_samps > 0:
                shifted[:shift_samps] = 0
            elif shift_samps < 0:
                shifted[shift_samps:] = 0
            acc += shifted
        cc_stack[iR, :] = acc / max(1, (z1 - z0))
    return cc_stack


def plot_cc_panel(cc, lags, channels, source_channel, dz_m, title, out_png,
                  xlim=XLIM, mask_self=True):
    """Standard Lellouch-style: lags on x, channels (depth proxy) on y."""
    nch, _ = cc.shape
    # Mask zero-lag self-correlation if requested
    cc_disp = cc.copy()
    if mask_self:
        # exclude near the source channel from vmin/vmax computation
        mask = np.abs(channels - source_channel) > 5
    else:
        mask = np.ones(nch, dtype=bool)

    vlim = np.percentile(np.abs(cc_disp[mask]), 99)

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.imshow(
        cc_disp,
        extent=[lags[0], lags[-1], channels[-1], channels[0]],
        aspect="auto", origin="upper", cmap="RdBu_r",
        vmin=-vlim, vmax=vlim, interpolation="nearest",
    )
    # Reference moveout line — distance / V_AVG, where distance = |ch - source| * dz
    dist_m = np.abs(channels - source_channel) * dz_m
    ax.plot(+dist_m / V_AVG, channels, "--", color="k", lw=0.8, label=f"+{V_AVG:.0f} m/s")
    ax.plot(-dist_m / V_AVG, channels, "--", color="k", lw=0.8)
    ax.axhline(source_channel, color="lime", lw=0.8, ls=":", label=f"source ch {source_channel}")
    ax.set_xlim(*xlim)
    ax.set_xlabel("Lag (s)")
    ax.set_ylabel("Channel index (top → bottom of fiber)")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_fk(cc, channels, dt, dz_m, out_png, title):
    """F-K of the CC (channel vs lag)."""
    nch, nlag = cc.shape
    win_t = np.hanning(nlag)[None, :]
    win_x = np.hanning(nch)[:, None]
    Xw = cc * win_t * win_x

    F = np.fft.rfft(Xw, axis=1)
    FK = np.fft.fftshift(np.fft.fft(F, axis=0), axes=0)
    P = 10.0 * np.log10(np.maximum(np.abs(FK) ** 2, 1e-30))
    f = np.fft.rfftfreq(nlag, d=dt)
    k = np.fft.fftshift(np.fft.fftfreq(nch, d=dz_m))

    fmask = f <= 25.0
    f_d = f[fmask]
    P_d = P[:, fmask]

    fig, ax = plt.subplots(figsize=(8, 6))
    pcm = ax.pcolormesh(
        k, f_d, P_d.T, shading="auto", cmap="magma",
        vmin=np.percentile(P_d, 70), vmax=np.percentile(P_d, 99.5),
    )
    for label, v in [("3200 m/s", 3200.0), ("1600 m/s", 1600.0),
                     ("1500 m/s", 1500.0), ("500 m/s", 500.0)]:
        f_line = np.linspace(0.1, 25.0, 200)
        k_line = f_line / v
        ax.plot( k_line, f_line, "--", lw=0.8, label=label)
        ax.plot(-k_line, f_line, "--", lw=0.8, color=ax.lines[-1].get_color())
    ax.set_xlim(-0.05, 0.05)
    ax.set_ylim(0, 25)
    ax.axhspan(5, 20, color="cyan", alpha=0.15)
    ax.set_xlabel("k (cycles/m)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(pcm, ax=ax, label="Power (dB)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    npz_path = resolve_npz_path()
    print(f"Loading {npz_path}")
    d = np.load(npz_path, allow_pickle=True)
    cc       = d["cc"]
    lags     = d["lags"]
    channels = d["channels"]
    source_ch = int(d["source_channel"])
    isource_idx = int(d["isource_idx"])
    fs        = float(d["fs"])
    dt        = float(d["dt"])
    fmin      = float(d["fmin"])
    fmax      = float(d["fmax"])
    n_stack   = int(d["n_stack"])
    date_str  = str(d["date"])
    use_whiten = bool(d["use_whiten"]) if "use_whiten" in d.files else False
    ch_start = int(d["ch_start"]) if "ch_start" in d.files else int(channels[0])
    ch_end = int(d["ch_end"]) if "ch_end" in d.files else int(channels[-1]) + 1
    hour_mode = str(d["hour_mode"]) if "hour_mode" in d.files else "unknown"
    run_tag = f"{date_str}_hours{hour_mode}_ch{ch_start}-{ch_end}_src{source_ch}"
    whitening_label = "phase-only whitening" if use_whiten else "whitening off"

    print(f"Date={date_str}  cc shape={cc.shape}  fs={fs}  n_stack={n_stack}")
    print(f"Bandpass {fmin}-{fmax} Hz   source channel = {source_ch}")

    # 1. Raw CC (no pair-stacking)
    title_raw = (
        f"Sanity CC — {date_str}, raw (no pair-stack)\n"
        f"hours={hour_mode}, src ch {source_ch}, {fmin}-{fmax} Hz, n_stack={n_stack}, "
        f"running-AM 0.1 s, {whitening_label}"
    )
    plot_cc_panel(cc, lags, channels, source_ch, DZ_M, title_raw,
                  os.path.join(OUT_DIR, f"cc_raw_{run_tag}.png"))
    print(f"Wrote cc_raw_{run_tag}.png")

    # 2. Adjacent-channel pre-shift stack (Lellouch's R±10 @ V_AVG)
    print(f"Applying R±{PAIR_HALF} pre-shift stack at V_AVG={V_AVG:.0f} m/s...")
    cc_stack = adjacent_pair_preshift_stack(
        cc, channels, isource_idx, DZ_M, V_AVG, PAIR_HALF, dt,
    )
    title_stack = (
        f"Sanity CC — {date_str}, R±{PAIR_HALF} pre-shifted stack at {V_AVG:.0f} m/s\n"
        f"hours={hour_mode}, src ch {source_ch}, {fmin}-{fmax} Hz, n_stack={n_stack}"
    )
    plot_cc_panel(cc_stack, lags, channels, source_ch, DZ_M, title_stack,
                  os.path.join(OUT_DIR, f"cc_stacked_{run_tag}.png"))
    print(f"Wrote cc_stacked_{run_tag}.png")

    # 3. F-K of the stacked CC
    plot_fk(cc_stack, channels, dt, DZ_M,
            os.path.join(OUT_DIR, f"cc_fk_{run_tag}.png"),
            title=f"F-K of stacked CC — {date_str}")
    print(f"Wrote cc_fk_{run_tag}.png")
    print("Done.")


if __name__ == "__main__":
    main()
