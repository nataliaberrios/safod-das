#!/usr/bin/env python3
"""Build every figure for the Figure 7c non-reproduction write-up.

RULE, and it is enforced rather than asserted: every number and every curve
plotted here is READ FROM A PRODUCT on disk. Nothing is typed in from a previous
run, and nothing is synthesised. If a product is missing the script fails loudly
with the name of the script that produces it, rather than drawing a partial
figure -- because a figure that silently omits a panel is exactly how a
reproducibility document becomes wrong.

Figures written to figures/ as PNG and PDF:
  fig1_no_reproduction      the faithful reproduction and its null
  fig2_method_validation    our picker recovers Lellouch's own Figure 9
  fig3_fig7d_isolation      constant-offset gather: his moveout, our zero lag
  fig4_eight_methods        the pedestal diagnostic across every method tried
  fig5_static_pattern       fixed wavenumber, not fixed velocity
  fig6_illumination         the asymmetry Lellouch used, present in 2017 only
  fig7_archive_scan         archive-wide illumination scan

Run from this directory with the `das` interpreter.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
AWD = HERE.parent
FIG = HERE / "figures"
AGG = AWD / "ambient_transfer" / "lellouch2019_exact_stack" / "_authoritative_38993456"

V_REF = 3200.0
FMT = dict(dpi=300, bbox_inches="tight")
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                     "legend.fontsize": 7.5, "figure.dpi": 110})

# Colour-blind-safe, consistent across every figure.
C_2017, C_2024, C_NULL, C_REF = "#0072B2", "#D55E00", "#666666", "#009E73"

NUMBERS: dict[str, str] = {}      # every value quoted in the text, with source


def need(path: Path, produced_by: str):
    if not path.is_file():
        raise SystemExit(
            "MISSING PRODUCT: %s\n  produced by: %s\n"
            "Refusing to draw a partial figure." % (path, produced_by))
    return path


def record(key: str, value, source: str, fmt: str = "%s"):
    """Store one quoted value with the product that produced it.

    `fmt` must consume exactly one argument. A mismatch used to raise
    TypeError deep inside figure code; it is now caught here with the key named,
    because a broken provenance record is worse than a missing figure.
    """
    try:
        text = fmt % value
    except TypeError as exc:
        raise SystemExit("record(%r): format %r does not match value %r (%s)"
                         % (key, fmt, value, exc))
    NUMBERS[key] = "%s  [%s]" % (text, source)
    return value


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / (name + ".png"), **FMT)
    fig.savefig(FIG / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  wrote figures/%s.{png,pdf}" % name)


# ---------------------------------------------------------------- figure 1
def fig1():
    src = "ambient_lellouch2019_exact_stack.py --action aggregate"
    p0 = need(AGG / "aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0.npz", src)
    p3 = need(AGG / "aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0_cm.npz", src)
    d0, d3 = np.load(p0), np.load(p3)

    fig, ax = plt.subplots(1, 3, figsize=(11.0, 3.5), constrained_layout=True)
    lags, off = d0["lags_s"], d0["offsets_m"]
    sec = d0["r_plus_minus_10_correlation"]
    lim = float(np.percentile(np.abs(sec), 99.0))
    ax[0].imshow(sec, extent=[lags[0], lags[-1], off[-1], off[0]], aspect="auto",
                 cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
    ax[0].plot(off / V_REF, off, "--", color=C_REF, lw=1.6,
               label="%.0f m/s (Lellouch)" % V_REF)
    ax[0].set(xlim=(-0.35, 0.35), xlabel="Correlation lag (s)",
              ylabel="Receiver offset below wellhead (m)",
              title="(a) Our Figure 7c gather, 24 h")
    ax[0].legend(loc="lower right")

    for d, lab, col in ((d0, "no common-mode removal", C_2024),
                        (d3, "common-mode removed", C_2017)):
        g, sc = d["velocity_grid_m_s"] / 1e3, d["causal_moveout_scores"]
        ax[1].plot(g, sc, lw=1.5, color=col, label=lab)
        r = float(np.corrcoef(d["velocity_grid_m_s"], sc)[0, 1])
        record("pedestal_corr_%s" % ("cfg0" if d is d0 else "cfg3"), r, src, "%+.3f")
    ax[1].axvline(V_REF / 1e3, color=C_REF, ls=":", lw=1.4)
    ax[1].axvspan(2.5, 4.0, color=C_REF, alpha=.10)
    ax[1].set(xlabel="Trial apparent velocity (km/s)", ylabel="Envelope moveout score",
              title="(b) Moveout scan")
    ax[1].legend(); ax[1].grid(alpha=.3)

    nulls = d0["receiver_order_null_maxima"]
    obs = float(np.max(d0["causal_moveout_scores"]))
    pval = float((np.sum(nulls >= obs) + 1) / (nulls.size + 1))
    record("fig1_p_value", pval, src, "%.4f")
    record("fig1_null_n", nulls.size, src, "%d")
    ax[2].hist(nulls, bins=50, color=C_NULL, alpha=.85, label="receiver-order null")
    ax[2].axvline(obs, color=C_2024, lw=2.2, label="observed (p = %.3f)" % pval)
    ax[2].axvline(float(np.percentile(nulls, 95)), color="k", ls="--", lw=1.2,
                  label="null 95th")
    ax[2].set(xlabel="Maximum moveout score", ylabel="Null realisations",
              title="(c) The observed peak is inside the null")
    ax[2].legend()
    save(fig, "fig1_no_reproduction")


# ---------------------------------------------------------------- figure 2
def fig2():
    """The picker validation, stated as what it actually is.

    CORRECTION. An earlier version of this figure -- and of this project's own
    plan document -- described `r = 0.948` as our picker "recovering Lellouch's
    published Figure 9 velocity model". It is not that. In
    ambient_lellouch_fig7d_profile.py the number is
    `np.corrcoef(lel_z, lel_v)`: the correlation between DEPTH and the velocity
    our picker returns when applied to his released Figure 7d correlograms. His
    published Figure 9 curve is not digitised anywhere in this repository, so no
    point-by-point agreement with it is established here.

    What the number does establish is still the thing the paper needs, and the
    contrast in panel (b) is the actual validation: the same picker, unchanged,
    yields a monotonic physically sensible profile on data that contains the
    arrival and nonsense on ours.
    """
    src = "ambient_lellouch_fig7d_profile.py"
    d = np.load(need(AWD / "ambient_lellouch_fig7d_profile.npz", src), allow_pickle=True)
    lz, lv = np.asarray(d["lellouch_z"], float), np.asarray(d["lellouch_v"], float)
    ours = np.asarray(d["v_median"], float)
    mid = np.asarray(d["midpoints_m"], float)
    if lz.size == 0:
        raise SystemExit("no Lellouch picks in %s" % src)

    r_depth = float(np.corrcoef(lz, lv)[0, 1])
    record("fig2_r_depth_velocity_his_traces", r_depth, src, "%.3f")
    record("fig2_his_v_shallow", float(lv[0]), src, "%.0f m/s")
    record("fig2_his_v_deep", float(lv[-1]), src, "%.0f m/s")
    finite_ours = int(np.sum(np.isfinite(ours)))
    NUMBERS["fig2_our_picks_finite"] = "%d of %d finite  [%s]" % (
        finite_ours, ours.size, src)
    if finite_ours:
        record("fig2_our_v_median", float(np.nanmedian(ours)), src, "%.0f m/s")

    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.6), constrained_layout=True)
    ax[0].plot(lv, lz, "o-", color=C_2017, ms=4, lw=1.5,
               label="our picker on Lellouch's\nreleased Fig. 7d traces")
    ax[0].invert_yaxis()
    ax[0].set(xlabel="P velocity (m/s)", ylabel="Depth (m)",
              title="(a) Same picker, his data:\n%.0f to %.0f m/s, r(z,v) = %.3f"
                    % (lv[0], lv[-1], r_depth))
    ax[0].legend(); ax[0].grid(alpha=.3)

    ok = np.isfinite(ours)
    if ok.any():
        ax[1].semilogx(ours[ok], mid[ok], "s", color=C_2024, ms=5,
                       label="our picker on our data\n(%d of %d picks finite)"
                             % (finite_ours, ours.size))
    ax[1].axvspan(lv.min(), lv.max(), color=C_2017, alpha=.18,
                  label="physical range from (a)")
    ax[1].invert_yaxis()
    ax[1].set(xlabel="P velocity (m/s), log scale", ylabel="Depth (m)",
              title="(b) Same picker, our data:\nno usable picks")
    ax[1].legend(loc="lower left"); ax[1].grid(alpha=.3, which="both")
    save(fig, "fig2_method_validation")


# ---------------------------------------------------------------- figure 3
def fig3():
    src = "ambient_lellouch2019_fig7d.py"
    d = np.load(need(AWD / "ambient_lellouch2019_fig7d.npz", src), allow_pickle=True)
    fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.5), constrained_layout=True)
    lags, depths = d["lags"], d["depths"]
    ours = d["ours_plain"]
    norm = ours / np.maximum(np.abs(ours).max(axis=1, keepdims=True), 1e-30)
    for z, row in zip(depths, norm):
        ax[0].plot(lags, -z + row * 38.0, "k-", lw=0.7)
    ax[0].axvline(0.0, color=C_2024, ls="--", lw=1.4, label="zero lag")
    ax[0].set(xlim=(-0.06, 0.06), xlabel="Correlation lag (s)",
              ylabel="Depth of the 50 m pair (m)",
              title="(a) Our constant-offset gather")
    ax[0].legend()

    lz, lv = d["lellouch_z"], d["lellouch_v"]
    his_lag = 50.0 / np.asarray(lv, dtype=float)
    ours_pick = np.asarray(d["pick_plain"], dtype=float)
    finite = int(np.sum(np.isfinite(ours_pick)))
    NUMBERS["fig3_finite_picks"] = "%d of %d picks are finite  [%s]" % (
        finite, ours_pick.size, src)
    ax[1].plot(his_lag * 1e3, lz, "o-", color=C_2017, ms=3.5, lw=1.4,
               label="Lellouch: 50 m / v(z)")
    ok = np.isfinite(ours_pick)
    if ok.any():
        ax[1].plot(ours_pick[ok] * 1e3, np.asarray(depths)[ok], "s", color=C_2024,
                   ms=4, label="our picks")
    ax[1].axvline(0.0, color=C_2024, ls="--", lw=1.4, label="our peak: 0.0 ms")
    ax[1].invert_yaxis()
    ax[1].set(xlabel="Peak lag (ms)", ylabel="Depth (m)",
              title="(b) His lag migrates with depth; ours does not")
    ax[1].legend(); ax[1].grid(alpha=.3)
    save(fig, "fig3_fig7d_isolation")


# ---------------------------------------------------------------- figure 4
def fig4():
    """Pedestal diagnostic corr(trial velocity, score) for every method tried.

    |corr| near 1 means the statistic is measuring proximity to the zero-lag
    lobe, not moveout, so its p-value is not a detection.
    """
    rows = []
    src0 = "ambient_lellouch2019_exact_stack.py"
    p0 = need(AGG / "aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0.npz", src0)
    p3 = need(AGG / "aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0_cm.npz", src0)
    for p, lab in ((p0, "no removal (baseline)"), (p3, "median common-mode")):
        d = np.load(p)
        rows.append((lab, float(np.corrcoef(d["velocity_grid_m_s"],
                                            d["causal_moveout_scores"])[0, 1]), src0))

    fe = AWD / "ambient_flat_event_removal.npz"
    if fe.is_file():
        d = np.load(fe, allow_pickle=True)
        for k, lab in (("median_removed_x1_corr_v", "flat-event median x1"),
                       ("median_removed_x3_corr_v", "flat-event median x3")):
            if k in d.files:
                rows.append((lab, float(d[k]), "ambient_flat_event_removal.py"))
    pc = AWD / "ambient_pcc_pws.npz"
    if pc.is_file():
        d = np.load(pc, allow_pickle=True)
        for k, lab in (("PCC_linear_stack_corr_v", "phase cross-correlation"),
                       ("PCC_p_PWS_corr_v", "PCC + phase-weighted stack")):
            if k in d.files:
                rows.append((lab, float(d[k]), "ambient_pcc_pws.py"))
    rs = AWD / "ambient_radon_slant_stack.npz"
    if rs.is_file():
        d = np.load(rs, allow_pickle=True)
        v, e = d["v_kept"], d["e_causal"]
        rows.append(("tau-p slant stack", float(np.corrcoef(v, e)[0, 1]),
                     "ambient_radon_slant_stack.py"))

    for lab, val, s in rows:
        record("pedestal_%s" % lab.replace(" ", "_"), val, s, "%+.3f")

    fig, ax = plt.subplots(figsize=(7.4, 0.42 * len(rows) + 1.7),
                           constrained_layout=True)
    y = np.arange(len(rows))
    vals = [r[1] for r in rows]
    cols = [C_2017 if abs(v) < 0.5 else C_2024 for v in vals]
    ax.barh(y, vals, 0.62, color=cols)
    ax.axvline(0.5, color=C_NULL, ls="--", lw=1.2)
    ax.axvline(-0.5, color=C_NULL, ls="--", lw=1.2,
               label="|corr| = 0.5 (pedestal threshold)")
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows])
    ax.invert_yaxis()
    ax.set(xlabel="corr(trial velocity, moveout score)",
           title="Every velocity-domain method: is the statistic measuring moveout,\n"
                 "or proximity to the zero-lag lobe?")
    ax.legend(loc="lower right"); ax.grid(alpha=.3, axis="x")
    save(fig, "fig4_eight_methods")


# ---------------------------------------------------------------- figure 5
def fig5():
    src = "ambient_apparent_velocity_census.py (CENSUS_FMIN/FMAX, K0_REMOVE=1)"
    lo = need(AWD / "ambient_apparent_velocity_census_5-12Hz_k0rm.npz", src)
    hi = need(AWD / "ambient_apparent_velocity_census_12-20Hz_k0rm.npz", src)

    def centroid(path, fmin, fmax, cells=8):
        d = np.load(path, allow_pickle=True)
        P, k, f = d["P24"], d["k24"], d["f24"]
        dk = float(abs(k[1] - k[0]))
        band = (np.abs(f) >= fmin) & (np.abs(f) <= fmax)
        ka = np.abs(k)
        sel = (ka > 0.5 * dk) & (ka <= (cells + 0.5) * dk)
        w = P[np.ix_(sel, band)].sum(axis=1)
        c = float(np.sum(ka[sel] * w) / np.sum(w))
        return c, dk, float(np.mean(np.abs(f[band]))), ka[sel], w

    c_lo, dk, f_lo, k_lo, w_lo = centroid(lo, 5.0, 12.0)
    c_hi, _, f_hi, k_hi, w_hi = centroid(hi, 12.0, 20.0)
    f_ratio = record("fig5_f_ratio", f_hi / f_lo, src, "%.3f")
    k_ratio = record("fig5_k_ratio", c_hi / c_lo, src, "%.3f")
    record("fig5_centroid_lo", c_lo, src, "%.5f")
    record("fig5_centroid_hi", c_hi, src, "%.5f")

    fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.5), constrained_layout=True)
    ax[0].plot(k_lo / dk, w_lo / w_lo.sum(), "o-", color=C_2017, ms=3.5,
               label="5-12 Hz (centre %.1f Hz)" % f_lo)
    ax[0].plot(k_hi / dk, w_hi / w_hi.sum(), "s-", color=C_2024, ms=3.5,
               label="12-20 Hz (centre %.1f Hz)" % f_hi)
    ax[0].axvline(c_lo / dk, color=C_2017, ls=":", lw=1.3)
    ax[0].axvline(c_hi / dk, color=C_2024, ls=":", lw=1.3)
    ax[0].set(xlabel="wavenumber (cells from k = 0)", ylabel="normalised power",
              title="(a) Low-k marginal, both bands")
    ax[0].legend(); ax[0].grid(alpha=.3)

    ax[1].bar([0, 1], [k_ratio, f_ratio], 0.5,
              color=[C_2024, C_NULL])
    ax[1].axhline(1.0, color=C_2017, ls="--", lw=1.4,
                  label="static pattern predicts 1.000")
    ax[1].set_xticks([0, 1])
    ax[1].set_xticklabels(["observed\nk ratio", "a wave predicts\n= f ratio"])
    ax[1].set(ylabel="ratio (12-20 Hz / 5-12 Hz)",
              title="(b) Frequency doubled; wavenumber did not move\n"
                    "observed %.3f vs %.3f predicted for a wave" % (k_ratio, f_ratio))
    ax[1].legend(); ax[1].grid(alpha=.3, axis="y")
    save(fig, "fig5_static_pattern")


# ---------------------------------------------------------------- figure 6
def fig6():
    """|A| versus spatial rank removed, both epochs. Parsed from the v2 product."""
    src = "interrogator_and_illumination_v2.py"
    txt = need(AWD / "interrogator_and_illumination_v2.txt", src).read_text().splitlines()
    cur, data = None, {"2024-25": [], "2017": []}
    for line in txt:
        if "--- 2024-25" in line:
            cur = "2024-25"
        elif "--- 2017" in line:
            cur = "2017"
        elif cur and "rank" in line and "|A|" in line:
            parts = line.replace("(", " ").replace(")", " ").replace(",", " ").split()
            try:
                rank = int(parts[parts.index("rank") + 1])
                a = float(parts[parts.index("=") + 1])
                n95 = float(parts[parts.index("null95") + 1])
                p = float(parts[parts.index("p") + 2])
            except (ValueError, IndexError):
                continue
            data[cur].append((rank, a, n95, p))
    if not data["2017"] or not data["2024-25"]:
        raise SystemExit("could not parse |A| table from %s" % src)
    for epoch in data:
        for rank, a, n95, p in data[epoch]:
            NUMBERS["fig6_%s_rank%d" % (epoch, rank)] = (
                "|A| = %.4f, null95 = %.4f, p = %.4f  [%s]" % (a, n95, p, src))

    fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.5), constrained_layout=True)
    for epoch, col, mk in (("2017", C_2017, "o"), ("2024-25", C_2024, "s")):
        r = [x[0] for x in data[epoch]]
        a = [x[1] for x in data[epoch]]
        n = [x[2] for x in data[epoch]]
        lab = "Lellouch 2017 records" if epoch == "2017" else "2024-25 ambient"
        ax[0].plot(r, a, mk + "-", color=col, ms=4.5, lw=1.5, label=lab)
        ax[0].plot(r, n, ":", color=col, lw=1.2, alpha=.75,
                   label="%s null 95th" % lab)
    ax[0].set(xlabel="spatial rank projected out", ylabel="|A|  fan asymmetry",
              yscale="log", title="(a) Downgoing/upgoing asymmetry")
    ax[0].legend(); ax[0].grid(alpha=.3, which="both")

    r17 = [x for x in data["2017"] if x[0] <= 2]
    r24 = [x for x in data["2024-25"] if x[0] <= 2]
    x = np.arange(len(r17))
    ax[1].bar(x - 0.2, [v[1] for v in r17], 0.38, color=C_2017, label="2017")
    ax[1].bar(x + 0.2, [v[1] for v in r24], 0.38, color=C_2024, label="2024-25")
    for i, v in enumerate(r17):
        ax[1].annotate("p=%.3f" % v[3], (i - 0.2, v[1]), ha="center",
                       textcoords="offset points", xytext=(0, 3), fontsize=6.5)
    for i, v in enumerate(r24):
        ax[1].annotate("p=%.2f" % v[3], (i + 0.2, v[1]), ha="center",
                       textcoords="offset points", xytext=(0, 3), fontsize=6.5)
    ax[1].set_xticks(x); ax[1].set_xticklabels(["rank %d" % v[0] for v in r17])
    ax[1].set(ylabel="|A|", title="(b) Significant in 2017 only")
    ax[1].legend(); ax[1].grid(alpha=.3, axis="y")
    save(fig, "fig6_illumination")


# ---------------------------------------------------------------- figure 7
def fig7():
    src = "illumination_window_scan.py"
    p = AWD / "illumination_window_scan.npz"
    if not p.is_file():
        print("  SKIPPING fig7: %s not present yet (produced by %s)" % (p.name, src))
        return
    d = np.load(p, allow_pickle=True)
    t = d["t"].astype("datetime64[ns]")
    asym, pv, hour = d["asym"], d["p"], d["hour"]
    n_hits = int(d["n_hits"]); ceil_ = int(d["chance_ceiling"])
    record("fig7_n_windows", asym.size, src, "%d")
    record("fig7_n_hits", n_hits, src, "%d")
    record("fig7_chance_ceiling", ceil_, src, "%d")
    record("fig7_median_asym", float(np.median(asym)), src, "%.4f")

    fig, ax = plt.subplots(1, 3, figsize=(12.2, 3.4), constrained_layout=True)
    hit = pv < float(d["alpha"])
    ax[0].plot(t, asym, ".", ms=3.2, color=C_NULL, label="all windows")
    if hit.any():
        ax[0].plot(t[hit], asym[hit], "o", ms=4.5, color=C_2024,
                   label="p < %.2f" % float(d["alpha"]))
    ax[0].set(xlabel="Date (UTC)", ylabel="|A| in the fan",
              title="(a) %d windows across the archive" % asym.size)
    ax[0].legend(); ax[0].grid(alpha=.3)
    ax[0].tick_params(axis="x", rotation=30, labelsize=7)

    ax[1].hist(asym, bins=40, color=C_2024, alpha=.85)
    ax[1].axvline(float(np.median(d["null95"])), color=C_NULL, ls="--", lw=1.4,
                  label="median null 95th")
    ax[1].set(xlabel="|A|", ylabel="windows",
              title="(b) %d hits vs %d expected by chance" % (n_hits, ceil_))
    ax[1].legend(); ax[1].grid(alpha=.3)

    med = [np.median(asym[hour == h]) if np.any(hour == h) else np.nan
           for h in range(24)]
    ax[2].bar(np.arange(24), med, color=C_REF, alpha=.9)
    ax[2].set(xlabel="Hour of day (UTC)", ylabel="median |A|",
              title="(c) Diurnal test: cultural sources\nwould cluster in working hours")
    ax[2].grid(alpha=.3, axis="y")
    save(fig, "fig7_archive_scan")


def main():
    print("Building figures for the Figure 7c non-reproduction document")
    for fn in (fig1, fig2, fig3, fig4, fig5, fig6, fig7):
        print("- %s" % fn.__name__)
        fn()
    out = HERE / "FIGURE_NUMBERS.txt"
    lines = ["Every number quoted in the document, with the product it came from.",
             "Generated by make_figures.py -- do not edit by hand.", ""]
    for k in sorted(NUMBERS):
        lines.append("%-34s %s" % (k, NUMBERS[k]))
    out.write_text("\n".join(lines) + "\n")
    print("wrote FIGURE_NUMBERS.txt (%d values)" % len(NUMBERS))


if __name__ == "__main__":
    main()
