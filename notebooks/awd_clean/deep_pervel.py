#!/usr/bin/env python3
"""Per-velocity test of the Deep-fibre low-velocity peak, and the curve shape.

WHY THIS EXISTS. The max-over-grid statistic is unusable here. The Deep ARM A
score curve carries a strong LOW-velocity trend, so the maximum of a permuted
curve is set by whatever the lowest velocity in the scan happens to be. Widening
the scan from 1500-6000 to 300-6000 m/s left the observed peak unchanged (3.3327
at 1675 m/s for source 400) while moving its p-value from 0.0002 to 0.1596 -- the
p-value was tracking the grid, not the data. Any conclusion from that statistic
is therefore void, in either direction.

The PER-VELOCITY null has no such dependence: at each trial velocity the observed
score is compared against permuted scores AT THAT SAME VELOCITY, so the grid
edges cannot influence it. It also answers the question that actually matters --
is there a coherent arrival at some specific velocity -- rather than "is the
biggest number on this curve big".

It also measures the curve SHAPE, because the pedestal correlation weakened from
-0.921 to -0.647 when the low-velocity range was added. A purely monotonic trend
would have strengthened toward -1; weakening implies the curve turns over, i.e. a
genuine local maximum rather than a slope. That is checked explicitly.

CORRECTION, 2026-08-20. The "1675 m/s" above is RETRACTED, and the "strong
low-velocity trend" it refers to was largely an artefact. The +-0.35 s lag window
cannot hold a 700 m offset below 1934 m/s, so on a 300-6000 m/s grid the far
offsets' gates fell outside the window and the median ran over the surviving NEAR
offsets -- the ones closest to the zero-lag lobe, which score high. The score
therefore rose as velocity fell for a geometric reason, which is the "pedestal".
Re-aggregated at +-2.5 s the peak is 1525 m/s (src 400) and 1550 m/s (src 800),
a 25 m/s spread where it had been 138, and matching an independent envelope-pick
regression of 1443-1537 m/s. `ex.moveout_scores` now refuses to score a velocity
whose gates do not fit, and this script refuses stale narrow-window products.
See AUDIT_2026-08-20.md.

Reads the stored aggregates only. No reprocessing.
Output: deep_pervel.{npz,png,txt}
"""
from __future__ import annotations

import glob, os, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ambient_lellouch2019_exact_stack as ex
import safod_geometry as geo

STEM = HERE / "deep_pervel"
NULLS = 2000
SEED = 20260820
V_REF = 3200.0


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    rng = np.random.default_rng(SEED)

    say("Per-velocity test of the Deep ARM A low-velocity peak")
    say("  the max-over-grid p moved 0.0002 -> 0.1596 for source 400 when the")
    say("  scan was widened, with the peak unchanged, so that statistic is void.")
    say("  Per-velocity nulls compare like with like and cannot depend on the")
    say("  grid edges. %d permutations per velocity." % NULLS)
    say("")

    # TWO ADMISSION TESTS, both added 2026-08-20 after the audits.
    #
    # 1  LAG WINDOW vs VELOCITY GRID. A product is usable only if its own lag
    #    window can hold its own largest offset at its own slowest trial
    #    velocity. Below that speed the far offsets' gates fall outside the
    #    window, the median runs over the surviving NEAR offsets, and those sit
    #    closest to the zero-lag lobe -- so the score climbs as velocity falls
    #    for a geometric reason. That is what put the Deep peak at 1675 m/s;
    #    re-run at +-2.5 s it is 1525-1550. Mixing a +-0.35 s product with a
    #    +-2.5 s one would compare a biased curve against a corrected one.
    #
    # 2  SOURCE GEOMETRY. Only near-vertical OUTBOUND sources give a moveout
    #    that is a vertical velocity. Channel 98 is surface lead-in; 1200 and
    #    1600 sit at 49.6 and 54.3 degrees; 1600's receivers additionally cross
    #    the fibre reversal at 1700, so its gather contains both limbs; 2000 and
    #    2400 are on the return limb entirely.
    files = sorted(glob.glob(str(HERE / "ambient_transfer" / "deep_exact_stack" /
                                 "aggregate_deepA_src*_ordered_r0.npz")))
    keep, rejected = [], []
    for f in files:
        src = os.path.basename(f).split("src")[1].split("_")[0]
        d = np.load(f, allow_pickle=True)
        g = d["velocity_grid_m_s"]
        if g.size < 200:
            rejected.append((src, "velocity grid is only %d points" % g.size))
            continue
        max_lag = float(np.abs(d["lags_s"]).max())
        need = float(np.max(np.abs(d["offsets_m"])) / g.min())
        if max_lag < need:
            rejected.append((src, "lag window +-%.2f s cannot hold %.0f m at %.0f m/s "
                                  "(needs %.2f s) -- STALE, re-aggregate with --max-lag"
                                  % (max_lag, np.max(np.abs(d["offsets_m"])), g.min(), need)))
            continue
        try:
            gd = geo.load()
            j = int(np.searchsorted(gd["channel"], int(src)))
            if not gd["near_vertical"][j]:
                rejected.append((src, geo.describe(int(src)).split(": ", 1)[1]))
                continue
            # THE RECEIVERS MUST BE NEAR-VERTICAL TOO, not just the source.
            # Checking only the source passed src 800, whose 700 m aperture runs
            # from channel 810 to 1153 -- and the near-vertical section ends at
            # 949, so 204 of its 354 channels are in the deviated hole (17-55
            # degrees). Along-fibre distance is not vertical depth there, so its
            # "velocity" is an apparent one and cannot join a depth comparison.
            cc = np.asarray(d["center_channels"]).ravel().astype(int)
            span = np.arange(cc.min() - 10, cc.max() + 11)   # R+-10 neighbours
            idx = np.clip(np.searchsorted(gd["channel"], span), 0,
                          gd["channel"].size - 1)
            nv = gd["near_vertical"][idx]
            if not nv.all():
                rejected.append((src, "receivers ch %d-%d leave the near-vertical "
                                      "section (%d of %d channels deviated); the "
                                      "moveout is apparent, not vertical"
                                      % (span.min(), span.max(),
                                         int((~nv).sum()), nv.size)))
                continue
        except Exception:
            pass                      # geometry unavailable; fall through
        keep.append((src, f))

    if rejected:
        say("=== products NOT used ===")
        for src, why in rejected:
            say("  src %-5s %s" % (src, why))
        say("")
    if not keep:
        raise SystemExit("no usable deep aggregates: every product was rejected above")
    say("=== products used ===")
    for src, _ in keep:
        say("  src %-5s %s" % (src, geo.describe(int(src)).split(": ", 1)[1]))
    say("")

    out = {}
    for src, f in keep:
        d = np.load(f, allow_pickle=True)
        g = d["velocity_grid_m_s"]
        gather = d["r_plus_minus_10_correlation"]
        lags = d["lags_s"]
        offs = d["offsets_m"]
        cs = d["causal_moveout_scores"]

        # per-velocity null: permute receiver order, rescore the WHOLE curve
        null = np.empty((NULLS, g.size))
        for i in range(NULLS):
            perm = rng.permutation(len(offs))
            null[i] = ex.moveout_scores(gather[perm], lags, offs, velocities=g)
        thresh = np.percentile(null, 95.0, axis=0)
        pv = (np.sum(null >= cs[None, :], axis=0) + 1.0) / (NULLS + 1.0)
        clears = np.flatnonzero(cs > thresh)

        k = int(np.argmax(cs))
        # shape test: is the maximum interior, with the curve falling on BOTH sides?
        interior = 0 < k < g.size - 1
        lo_side = float(cs[:k].min()) if k > 0 else np.nan
        hi_side = float(cs[k+1:].min()) if k < g.size-1 else np.nan
        turnover = interior and cs[k] > cs[0] and cs[k] > cs[-1]

        say("--- source channel %s ---" % src)
        say("  peak %.4f at %.0f m/s | per-velocity p there = %.4f"
            % (cs[k], g[k], pv[k]))
        say("  velocities clearing their OWN 95th percentile: %d of %d"
            % (clears.size, g.size))
        if clears.size:
            say("    range %.0f-%.0f m/s, best p = %.4f at %.0f m/s"
                % (g[clears].min(), g[clears].max(), pv.min(), g[int(np.argmin(pv))]))
        say("  curve shape: maximum is %s; higher than both ends: %s"
            % ("interior" if interior else "AT AN EDGE", turnover))
        say("  score at grid ends: %.4f (%.0f m/s) and %.4f (%.0f m/s)"
            % (cs[0], g[0], cs[-1], g[-1]))
        say("  per-velocity p at %.0f m/s (Lellouch): %.4f" % (V_REF, float(np.interp(V_REF, g, pv))))
        say("")
        out[src] = dict(g=g, cs=cs, thresh=thresh, pv=pv, peak_v=float(g[k]),
                        peak_p=float(pv[k]), n_clear=int(clears.size),
                        turnover=bool(turnover))

    say("=== reading ===")
    any_clear = [s for s, v in out.items() if v["n_clear"] > 0]
    turn = [s for s, v in out.items() if v["turnover"]]
    say("  source channels with >=1 velocity clearing its own null: %s"
        % (any_clear if any_clear else "NONE"))
    say("  source channels whose curve has an interior maximum: %s"
        % (turn if turn else "NONE"))
    say("")
    if not any_clear:
        say("  NO velocity, at any source channel, clears its own per-velocity")
        say("  null. The low-velocity feature is therefore not a detection: it is")
        say("  a trend in the statistic that receiver-order permutation reproduces.")
        say("  The earlier p = 0.0002 was an artefact of the 1500 m/s scan floor.")
    else:
        say("  At least one velocity clears its own per-velocity null, which the")
        say("  grid cannot explain. Report the velocity and source channel; this")
        say("  is a candidate arrival and needs the remaining gates (pedestal,")
        say("  causal dominance, and an input-level null) before it is a result.")

    fig, ax = plt.subplots(1, len(out), figsize=(5.2*len(out), 4.0), squeeze=False,
                           constrained_layout=True)
    for i, (src, v) in enumerate(sorted(out.items())):
        a = ax[0][i]
        a.plot(v["g"]/1e3, v["cs"], "k-", lw=1.6, label="observed")
        a.plot(v["g"]/1e3, v["thresh"], "--", color="crimson", lw=1.2,
               label="per-velocity null 95th")
        a.axvline(V_REF/1e3, color="#009E73", ls=":", lw=1.4, label="3200 m/s")
        a.axvline(v["peak_v"]/1e3, color="#0072B2", ls="-.", lw=1.2,
                  label="peak %.0f m/s" % v["peak_v"])
        a.set(xlabel="trial velocity (km/s)", ylabel="moveout score" if i==0 else "",
              title="src %s: %d/%d velocities clear" % (src, v["n_clear"], v["g"].size))
        a.legend(fontsize=7); a.grid(alpha=.3)
    fig.savefig(str(STEM)+".png", dpi=190)
    np.savez(str(STEM)+".npz", **{f"{s}_{k}": v[k] for s, v in out.items()
                                  for k in ("g","cs","thresh","pv")})
    Path(str(STEM)+".txt").write_text("\n".join(log)+"\n")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
