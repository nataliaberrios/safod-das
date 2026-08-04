"""Clean six-figure advisor set for the SAFOD June 2026 AWD-DAS project.

This script only reformats accepted/reviewed reduced products. It does not
rerun the raw-data reconstruction or introduce a new estimator. Diagnostic
figures remain available under their original names; the ``core0*`` products
are the concise presentation set.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "figures" / "awd_2026"
OUT.mkdir(parents=True, exist_ok=True)

# Colorblind-safe, restrained journal palette.
NAVY = "#173F5F"
BLUE = "#20639B"
TEAL = "#3CAEA3"
ORANGE = "#F6A01A"
RED = "#C23B22"
GRAY = "#6B7280"
LIGHT = "#E8EEF3"
DARK = "#1F2933"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.linewidth": .8,
    "xtick.major.width": .7, "ytick.major.width": .7,
    "savefig.facecolor": "white", "figure.facecolor": "white",
})


def panel(ax, letter):
    ax.text(-0.12, 1.04, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", color=NAVY, va="bottom")
    ax.tick_params(direction="out", length=3, width=.7, colors=DARK)
    for spine in ax.spines.values():
        spine.set_color("#59636E")


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=.06)
    fig.savefig(OUT / f"{stem}.png", dpi=400, bbox_inches="tight", pad_inches=.06)
    plt.close(fig)


def fig1_experiment():
    m = pd.read_csv(HERE / "awd_manifest.csv")
    m["utc_time"] = pd.to_datetime(m["utc_time"], utc=True)
    burst = m.groupby("burst_id").agg(start=("utc_time", "min"), drops=("drop_id", "size"))
    hours = (burst.start - burst.start.iloc[0]).dt.total_seconds().to_numpy() / 3600
    fig = plt.figure(figsize=(7.2, 4.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=(1.25, 1), height_ratios=(1, 1.08))
    ax = fig.add_subplot(gs[:, 0]); panel(ax, "A")
    ax.vlines(hours, 0, burst.drops, color=BLUE, lw=.8, alpha=.6)
    ax.scatter(hours, burst.drops, s=22, color=NAVY, zorder=3)
    ax.set(xlabel="hours from first qualifying burst", ylabel="drops per burst",
           title="Repeated AWD sampling")
    ax.text(.03, .96, "49 bursts · 988 drops", transform=ax.transAxes,
            va="top", color=NAVY, fontweight="bold")
    ax.set_ylim(0, max(burst.drops) + 4); ax.grid(axis="y", color=LIGHT, lw=.7)
    ax2 = fig.add_subplot(gs[0, 1]); panel(ax2, "B")
    ax2.plot([0, 0], [0, 1], color=RED, lw=2)
    ax2.plot([0, 0], [1.1, 2.2], color=BLUE, lw=2)
    ax2.plot([0, 0], [2.3, 9.1], color=ORANGE, lw=2)
    ax2.scatter([-.08, -.08, -.08], [.5, 1.65, 5.7], color=[RED, BLUE, ORANGE], s=24)
    ax2.set_xlim(-.8, 1.0); ax2.set_ylim(-.15, 9.5)
    ax2.set_xticks([]); ax2.set_yticks([.5, 1.65, 5.7], ["AWD source", "Nano", "Deep"])
    ax2.set_title("Installation geometry", loc="left")
    ax2.text(.5, .5, "15 m lateral\noffset", transform=ax2.transAxes,
             ha="center", va="center", color=RED, fontsize=8)
    ax2.text(.82, .21, "cemented\n80–460 m", transform=ax2.transAxes,
             ha="center", color=BLUE, fontsize=8)
    ax2.text(.82, .66, "wireline\n≈3.4 km fiber", transform=ax2.transAxes,
             ha="center", color=ORANGE, fontsize=8)
    ax2.set_xlabel("schematic, not to scale")
    ax3 = fig.add_subplot(gs[1, 1]); panel(ax3, "C")
    ax3.axis("off")
    ax3.text(.02, .88, "What this data set supports", color=NAVY, fontweight="bold", fontsize=9)
    lines = [
        "• repeated-source coherence",
        "• installation-dependent mode visibility",
        "• drop- and burst-scale repeatability",
        "• an apparent-moveout sensitivity test",
    ]
    for j, line in enumerate(lines):
        ax3.text(.05, .70 - .18*j, line, color=DARK, fontsize=8.5)
    save(fig, "core01_experiment_geometry")


def fig2_nano_mode():
    d = np.load(HERE / "fk_dispersion.npz")
    f = d["frequency"]
    p = d["slowness"] * 1e3  # ms/m
    power = d["nano_80_440_m_power"]
    keep = (f >= 10) & (f <= 80)
    ff = f[keep]; pp = power[keep]
    csv = pd.read_csv(HERE / "nano_mode_identification.csv")
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.5), constrained_layout=True,
                           gridspec_kw={"width_ratios": (1.12, 1)})
    panel(ax[0], "A")
    im = ax[0].pcolormesh(p, ff, np.log10(np.maximum(pp, 1e-5)), shading="auto",
                          cmap="magma", vmin=-4, vmax=-.1)
    ridge = d["nano_80_440_m_ridge_slowness"] * 1e3
    mask = (f >= 10) & (f <= 80) & (ridge > 0) & (ridge < .9)
    ax[0].plot(ridge[mask], f[mask], color="#75D5C5", lw=1.8, label="ridge")
    ax[0].axvline(0, color="white", lw=.7, alpha=.8)
    ax[0].axvspan(1/2990*1e3, 1/2950*1e3, color="white", alpha=.14, lw=0)
    ax[0].set(xlim=(-.75, .75), ylim=(10, 80), xlabel="signed slowness (ms m⁻¹)",
              ylabel="frequency (Hz)", title="Nano phase coherence")
    ax[0].text(.03, .94, "positive sign = increasing fiber coordinate", transform=ax[0].transAxes,
               color="white", fontsize=7.5, va="top")
    cb = fig.colorbar(im, ax=ax[0], pad=.02, fraction=.046); cb.set_label("log₁₀ phase-beam power")
    ax[0].legend(frameon=False, loc="lower right")
    panel(ax[1], "B")
    x = csv.band_center_hz.to_numpy(); y = csv.bootstrap_speed_median_mps.to_numpy()
    lo = csv.bootstrap_speed_p2p5_mps.to_numpy(); hi = csv.bootstrap_speed_p97p5_mps.to_numpy()
    ax[1].axhspan(2950, 2990, color=TEAL, alpha=.16, label="2.95–2.99 km s⁻¹ interval")
    ax[1].errorbar(x, y, yerr=[y-lo, hi-y], fmt="o", color=NAVY, ms=4.5,
                   lw=1.1, capsize=2, label="burst-bootstrap median ±95%")
    ax[1].axhline(0, color="#999", lw=.6)
    ax[1].set(xlim=(12, 78), ylim=(2400, 3600), xlabel="frequency (Hz)",
              ylabel="apparent speed (m s⁻¹)", title="Dispersion estimate")
    ax[1].grid(axis="y", color=LIGHT, lw=.7); ax[1].legend(frameon=False, loc="upper right")
    save(fig, "core02_nano_mode_fk_dispersion")


def fig3_repeatability():
    d = pd.read_csv(HERE / "nano_drop_repeatability.csv")
    b = pd.read_csv(HERE / "nano_burst_repeatability_hierarchical.csv")
    c = pd.read_csv(HERE / "nano_stack_convergence.csv")
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.9), constrained_layout=True)
    panel(ax[0], "A")
    for vals, color, label in [(d.loo_noise_ncc, GRAY, "noise window"),
                               (d.loo_signal_ncc, BLUE, "moveout window")]:
        vals = np.sort(vals.to_numpy()); ecdf = np.arange(1, len(vals)+1)/len(vals)
        ax[0].plot(vals, ecdf, color=color, lw=1.5, label=label)
    ax[0].set(xlim=(-.05, 1.02), ylim=(0, 1.02), xlabel="NCC", ylabel="cumulative fraction",
              title="Individual drops")
    ax[0].grid(color=LIGHT, lw=.7); ax[0].legend(frameon=False, loc="lower right")
    ax[0].axvline(d.loo_signal_ncc.median(), color=BLUE, lw=.7, ls="--")
    ax[0].axvline(d.loo_noise_ncc.median(), color=GRAY, lw=.7, ls="--")
    panel(ax[1], "B")
    ax[1].plot(b.burst_id, b.loo_burst_noise_ncc, color=GRAY, lw=.7, alpha=.75, label="noise")
    ax[1].plot(b.burst_id, b.loo_burst_signal_ncc, "o-", color=BLUE, lw=1.0, ms=2.7, label="signal")
    ax[1].set(xlim=(0, 48), ylim=(-.02, 1.02), xlabel="burst index", ylabel="leave-one-burst-out NCC",
              title="Across-burst stability")
    ax[1].grid(color=LIGHT, lw=.7); ax[1].legend(frameon=False, loc="lower right")
    ax[1].text(.04, .08, "49/49 bursts: signal > noise", transform=ax[1].transAxes, fontsize=7.5, color=NAVY)
    panel(ax[2], "C")
    g = c.groupby("n_drops_per_substack").independent_substack_ncc
    ns = np.array(sorted(g.groups)); med = np.array([g.get_group(n).median() for n in ns])
    lo = np.array([g.get_group(n).quantile(.16) for n in ns]); hi = np.array([g.get_group(n).quantile(.84) for n in ns])
    ax[2].fill_between(ns, lo, hi, color=BLUE, alpha=.18, lw=0)
    ax[2].plot(ns, med, "o-", color=BLUE, lw=1.5, ms=4)
    ax[2].set_xscale("log", base=2); ax[2].set_xticks(ns); ax[2].set_xticklabels(ns)
    ax[2].set(xlim=(.8, 9), ylim=(-.02, 1.02), xlabel="drops per substack", ylabel="NCC",
              title="Within-burst convergence")
    ax[2].grid(color=LIGHT, lw=.7)
    save(fig, "core03_repeatability_hierarchy")


def fig4_installation_contrast():
    fk = np.load(HERE / "fk_dispersion.npz")
    deep = np.load(HERE / "deep_sliding_fk.npz")
    np_f = fk["frequency"]; np_p = fk["slowness"] * 1e3; np_power = fk["nano_80_440_m_power"]
    nmask = (np_f >= 25) & (np_f <= 60)
    ncurve = np.nanmedian(np_power[nmask], axis=0); ncurve /= np.nanmax(ncurve)
    dp = deep["slowness"] * 1e3; dcurve = np.nanmedian(deep["outbound_power"][0], axis=0); dcurve /= np.nanmax(dcurve)
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True, constrained_layout=True)
    panel(ax[0], "A"); ax[0].plot(np_p, ncurve, color=BLUE, lw=1.8)
    ax[0].axvspan(1/2990*1e3, 1/2950*1e3, color=TEAL, alpha=.18)
    ax[0].axvline(1/2975*1e3, color=TEAL, ls="--", lw=1)
    ax[0].set(xlim=(-.1, .9), ylim=(0, 1.05), xlabel="positive signed slowness (ms m⁻¹)",
              ylabel="normalized phase-beam power", title="Nano · cemented · 25–60 Hz")
    ax[0].text(.05, .90, "2.95–2.99 km s⁻¹", transform=ax[0].transAxes, color=TEAL, fontsize=8)
    ax[0].grid(color=LIGHT, lw=.7)
    panel(ax[1], "B"); ax[1].plot(dp, dcurve, color=ORANGE, lw=1.8)
    ax[1].axvspan(1/1560*1e3, 1/1400*1e3, color=ORANGE, alpha=.18)
    ax[1].axvline(1/1480*1e3, color=ORANGE, ls="--", lw=1)
    ax[1].set(xlim=(-.1, .9), xlabel="positive signed slowness (ms m⁻¹)",
              title="Deep · wireline · outbound · 3–15 Hz")
    ax[1].text(.05, .90, "1.4–1.56 km s⁻¹", transform=ax[1].transAxes, color=ORANGE, fontsize=8)
    ax[1].grid(color=LIGHT, lw=.7)
    save(fig, "core04_installation_mode_contrast")


def fig5_aperture_sensitivity():
    obs = np.load(HERE / "nano_frequency_observability.npz")
    inj = np.load(HERE / "nano_dvv_injection_recovery.npz")
    f = obs["frequency"]; x = obs["distance_centers"]; snr = obs["median_snr_db"]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.35), constrained_layout=True,
                           gridspec_kw={"width_ratios": (1.08, 1)})
    panel(ax[0], "A")
    im = ax[0].imshow(snr.T, origin="lower", aspect="auto",
                      extent=(x.min()-20, x.max()+20, f.min(), f.max()),
                      cmap="RdYlBu_r", vmin=-5, vmax=20)
    ax[0].contour(x, f, snr.T, levels=[0], colors=[DARK], linewidths=1.0)
    ax[0].axhspan(30, 60, color="white", alpha=.12, lw=0)
    ax[0].set(xlim=(40, 480), ylim=(5, 80), xlabel="Nano distance along fiber (m)", ylabel="frequency (Hz)",
              title="Frequency-dependent observability")
    ax[0].text(.04, .94, "0 dB contour = median detectability boundary", transform=ax[0].transAxes,
               color=DARK, fontsize=7.5, va="top")
    cb = fig.colorbar(im, ax=ax[0], pad=.02, fraction=.046); cb.set_label("median SNR (dB)")
    panel(ax[1], "B")
    levels = np.unique(inj["injected_dvv"])
    med=[]; lo=[]; hi=[]
    for q in levels:
        z=inj["estimated_dvv"][inj["injected_dvv"]==q]
        med.append(np.median(z)); lo.append(np.percentile(z,16)); hi.append(np.percentile(z,84))
    ax[1].fill_between([-0.012,0.012], [-.003246,-.003246], [.003246,.003246], color=GRAY, alpha=.12, label="95% null band")
    ax[1].plot([-0.012,.012],[-.012,.012], color="#888", ls="--", lw=.9, label="1:1")
    ax[1].errorbar(levels, med, yerr=[np.asarray(med)-lo, hi-np.asarray(med)], fmt="o", color=NAVY, ms=3.8, lw=1, capsize=2, label="median ±16–84%")
    ax[1].axvline(.01, color=RED, ls=":", lw=1); ax[1].axvline(-.01, color=RED, ls=":", lw=1)
    ax[1].axhline(0, color="#999", lw=.6)
    ax[1].set(xlim=(-.0115,.0115), ylim=(-.0115,.0115), xlabel="injected apparent Δv/v", ylabel="estimated apparent Δv/v",
              title="Blind injection–recovery")
    ax[1].grid(color=LIGHT, lw=.7); ax[1].legend(frameon=False, loc="upper left")
    ax[1].text(.97, .06, "1% is the smallest tested\n95%-correct-sign level", transform=ax[1].transAxes,
               ha="right", color=RED, fontsize=7.5)
    save(fig, "core05_aperture_detection_limit")


def fig6_pgsi():
    d = np.load(OUT / "pgsi_checkshot_allshots.npz", allow_pickle=True)
    keys = d["shot_keys"].astype(str); y = np.arange(len(keys))
    speed = d["velocity_mps"] / 1000; lo=d["velocity_lo_mps"]/1000; hi=d["velocity_hi_mps"]/1000
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.45), constrained_layout=True, gridspec_kw={"width_ratios":(1,1.12)})
    panel(ax[0], "A")
    ax[0].axvspan(2.95, 2.99, color=TEAL, alpha=.18, label="Nano 2.95–2.99 km s⁻¹")
    for j, k in enumerate(keys):
        color = TEAL if d["usable"][j] else GRAY
        alpha = .55 if (k == "1071") else 1
        ax[0].errorbar(speed[j], j, xerr=[[speed[j]-lo[j]], [hi[j]-speed[j]]], fmt="o", color=color, alpha=alpha, ms=5, capsize=2, lw=1.2)
    ax[0].axvline(2.975, color=TEAL, ls="--", lw=.9)
    ax[0].set(xlim=(2.0, 6.6), ylim=(-.6, len(keys)-.4), yticks=y, yticklabels=keys,
              xlabel="geometry-corrected apparent speed (km s⁻¹)", title="All-shot PGSI fits")
    ax[0].invert_yaxis(); ax[0].grid(axis="x", color=LIGHT, lw=.7)
    ax[0].text(.04, .06, "5/6 usable; 1040 is the only Nano match", transform=ax[0].transAxes, color=NAVY, fontsize=7.5)
    panel(ax[1], "B")
    depth=d["receiver_depth"]
    colors=[TEAL if k=="1040" else (RED if d["usable"][j] else GRAY) for j,k in enumerate(keys)]
    for j,k in enumerate(keys):
        good=d["good"][j]; obs=d["observed"][j]
        if good.sum()==0: continue
        t=(obs[good]-np.nanmin(obs[good]))*1000
        ax[1].plot(depth[good], t, ".", ms=2.7, color=colors[j], alpha=.55)
        # Reference slope in depth coordinates is a visualization of registration, not a ray-traced fit.
        ref=(depth[good]-depth[good].min())/(2975)*1000
        ax[1].plot(depth[good], ref, color=colors[j], lw=1.0, alpha=.7)
    ax[1].set(xlabel="PGSI measured depth (m)", ylabel="relative travel time (ms)",
              title="Conditional registration with Nano depth axis")
    ax[1].grid(color=LIGHT, lw=.7)
    ax[1].text(.04, .94, "relative slopes only; offsets removed", transform=ax[1].transAxes, fontsize=7.5, color=NAVY, va="top")
    save(fig, "core06_pgsi_calibration_registration")


def main():
    fig1_experiment(); fig2_nano_mode(); fig3_repeatability(); fig4_installation_contrast(); fig5_aperture_sensitivity(); fig6_pgsi()
    print("Wrote six core figures to", OUT)


if __name__ == "__main__":
    main()
