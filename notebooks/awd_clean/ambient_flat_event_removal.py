#!/usr/bin/env python3
"""Median flat-event removal on the Figure 7c gather, then rescore.

THE GAP THIS FILLS.  Every method tried so far attacks the pedestal BEFORE or
DURING correlation: F-K wedges (hard and tapered), median common-mode removal,
rank-k subspace projection, phase cross-correlation. The VSP wavefield-separation
literature also uses a fourth family that none of those cover -- MEDIAN FILTERING
ALONG THE OFFSET AXIS OF THE FINAL GATHER, applied to suppress a flat
(zero-moveout) event while preserving events that move out. It appears throughout
that literature (Rao & Wang 2016 combine f-k with median filtering; Radon-
constrained median filtering is a standard VSP separator).

WHY IT SHOULD WORK HERE WHEN THE OTHERS DID NOT.  The measured properties of the
contaminant are:
  - common to all channels     -> the per-sample median removes it pre-correlation
  - phase-coherent, not amplitude-driven -> PCC (amplitude-blind) does NOT remove it
  - BROAD in lag               -> no slowness mute separates it; muting faster than
                                  8000/6000/5000/4200 m/s leaves corr(v,score) at
                                  +0.957/+0.940/+0.934/+0.914
A broad lobe at the SAME LAG in every trace is, by definition, a flat event in the
gather. Slowness filters fail on it because it is broad; a median across offsets
does not care about its width, only that it is at a common lag. Conversely a real
3200 m/s arrival sits at t = x/3200, i.e. a DIFFERENT lag in every trace, so the
across-offset median at any given lag sees it in at most one or two traces and it
survives the subtraction.

The other reason to do this at the gather stage: the pre-correlation median leaves
a residual (a_i - median(a))*c(t) that is coherent and therefore ACCUMULATES under
stacking. Measured: corr(v,score) is -0.381 for one day but +0.951 once four days
are coherently stacked. A post-stack removal is applied after that accumulation
has already happened, so it is not defeated by it.

CONTROLS.  Same receiver-order permutation null and the same precondition gate as
the rest of the tree: no recovery is declared unless the pedestal is suppressed,
the peak lies inside 2500-4000 m/s and off the scan edge, the causal side
dominates at 3200 m/s, and familywise p < 0.05.

IMPORTANT CAVEAT, stated before the result.  Subtracting an across-offset median
is itself a filter, and it could in principle create apparent moveout from
nothing. That is exactly the failure that has recurred in this project, so the
null here is computed by applying THE IDENTICAL removal to each permuted gather.
The null therefore carries the operator, and a p-value from it cannot be inflated
by the removal itself.

Output: ambient_flat_event_removal.{npz,png,txt}
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert

HERE = Path(__file__).resolve().parent
STEM = HERE / "ambient_flat_event_removal"
AGG = HERE / "ambient_transfer" / "lellouch2019_exact_stack"

V_REF, V_LO, V_HI = 3200.0, 2500.0, 4000.0
VELOCITY_GRID = np.arange(1500.0, 6000.1, 25.0)
GATE_S = 0.012
NULL_COUNT = 4000
SEED = 20260814


def remove_flat(gather, passes=1):
    """Subtract the across-offset median trace. Repeat `passes` times."""
    g = gather.copy()
    for _ in range(max(1, passes)):
        g = g - np.median(g, axis=0, keepdims=True)
    return g


def score(gather, lags, offsets, v, sign=1.0):
    env = np.abs(hilbert(gather, axis=1))
    env = env / np.maximum(np.median(env, axis=1, keepdims=True), 1e-30)
    half = max(1, int(GATE_S / (lags[1] - lags[0])))
    vals = []
    for row, x in zip(env, offsets):
        k = int(np.argmin(np.abs(lags - sign * x / v)))
        lo, hi = max(0, k - half), min(len(lags), k + half + 1)
        vals.append(row[lo:hi].mean())
    return float(np.median(vals))


def curve(gather, lags, offsets, sign=1.0, grid=VELOCITY_GRID):
    return np.array([score(gather, lags, offsets, v, sign) for v in grid])


def evaluate(gather, lags, offsets, label, passes, log, rng):
    """Apply the removal, score it, and null it with the SAME removal applied."""
    g = remove_flat(gather, passes) if passes else gather
    cs = curve(g, lags, offsets)
    ac = curve(g, lags, offsets, sign=-1.0)
    k = int(np.argmax(cs))
    coarse = VELOCITY_GRID[::4]
    nulls = np.empty(NULL_COUNT)
    for i in range(NULL_COUNT):
        p = rng.permutation(len(offsets))
        gp = remove_flat(gather[p], passes) if passes else gather[p]
        nulls[i] = max(score(gp, lags, offsets, v) for v in coarse)
    n95 = float(np.percentile(nulls, 95))
    pval = float((np.sum(nulls >= cs[k]) + 1) / (NULL_COUNT + 1))
    corr_v = float(np.corrcoef(VELOCITY_GRID, cs)[0, 1])
    c32 = float(np.interp(V_REF, VELOCITY_GRID, cs))
    a32 = float(np.interp(V_REF, VELOCITY_GRID, ac))
    ratio = c32 / a32 if a32 else np.nan

    ped = abs(corr_v) < 0.5
    fan = V_LO <= VELOCITY_GRID[k] <= V_HI
    edge = k not in (0, len(cs) - 1)
    dom = ratio > 1.0
    ok = all((ped, fan, edge, dom, pval < 0.05))

    log.append("--- %s ---" % label)
    log.append("  corr(velocity, score) %+.3f   peak %.4f at %.0f m/s"
               % (corr_v, cs[k], VELOCITY_GRID[k]))
    log.append("  at 3200: causal %.4f acausal %.4f ratio %.2f" % (c32, a32, ratio))
    log.append("  null95 %.4f | observed %.4f | p = %.4f" % (n95, cs[k], pval))
    log.append("  preconditions: pedestal %s | in fan %s | off edge %s | causal %s | p<0.05 %s"
               % (ped, fan, edge, dom, pval < 0.05))
    log.append("  >>> RECOVERED <<<" if ok else "  not recovered")
    log.append("")
    return dict(g=g, cs=cs, ac=ac, k=k, n95=n95, p=pval, corr_v=corr_v,
                c32=c32, a32=a32, ratio=ratio, ok=ok, nulls=nulls)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--branch", choices=("base", "cm"), default="base")
    a = ap.parse_args()
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    suffix = "_cm" if a.branch == "cm" else ""
    cands = [str(AGG / "_authoritative_38993456" /
                 ("aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0%s.npz" % suffix))]
    cands += sorted(glob.glob(str(AGG / ("aggregate_2024-12-20_*_ordered_r0%s.npz" % suffix))))
    src = next((c for c in cands if Path(c).is_file()), None)
    if src is None:
        raise SystemExit("no aggregate for branch " + a.branch)
    d = np.load(src, allow_pickle=True)
    lags = d["lags_s"] if "lags_s" in d.files else d["lags"]
    gather = np.asarray(d["r_plus_minus_10_correlation"], dtype=float)
    offsets = np.asarray(d["offsets_m"], dtype=float)

    say("Median flat-event removal along the offset axis, then rescore")
    say("  source %s" % Path(src).name)
    say("  %d offsets %.0f-%.0f m, %d lags" % (len(offsets), offsets.min(), offsets.max(), len(lags)))
    say("  null: the IDENTICAL removal is applied to every permuted gather, so the")
    say("        operator is inside the null and cannot inflate the p-value")
    say("")

    rng = np.random.default_rng(SEED)
    res = {}
    for passes, lab in ((0, "no removal (reference)"), (1, "median removed x1"),
                        (2, "median removed x2"), (3, "median removed x3")):
        res[lab] = evaluate(gather, lags, offsets, lab, passes, log, rng)
        print("\n".join(log[-7:]), flush=True)

    best = min(res.items(), key=lambda kv: abs(kv[1]["corr_v"]))
    say("Best pedestal suppression: %s at corr = %+.3f" % (best[0], best[1]["corr_v"]))
    if any(v["ok"] for v in res.values()):
        say("At least one configuration satisfies every precondition.")
    else:
        say("No configuration satisfies every precondition. Flat-event removal does")
        say("not expose a 3200 m/s arrival in this gather.")

    fig, ax = plt.subplots(1, 3, figsize=(17, 5.5), constrained_layout=True)
    for a_, lab in ((ax[0], "no removal (reference)"), (ax[1], "median removed x1")):
        g = res[lab]["g"]; gn = g / np.maximum(np.abs(g).max(axis=1, keepdims=True), 1e-30)
        for o, row in zip(offsets, gn):
            a_.plot(lags, -o + row * 45, "k-", lw=0.8)
        a_.plot(offsets / V_REF, -offsets, "--", color="crimson", lw=1.5, label="3200 m/s")
        a_.set_xlim(-0.35, 0.35); a_.set_xlabel("lag (s)"); a_.legend(fontsize=7); a_.set_title(lab)
    ax[0].set_ylabel("depth below wellhead (m)")
    for lab, v in res.items():
        ax[2].plot(VELOCITY_GRID, v["cs"], lw=1.2, label="%s (p=%.3f)" % (lab, v["p"]))
    ax[2].axvline(V_REF, color="steelblue", ls=":")
    ax[2].axvspan(V_LO, V_HI, color="orange", alpha=.15)
    ax[2].set_xlabel("trial velocity (m/s)"); ax[2].set_ylabel("moveout score")
    ax[2].legend(fontsize=6); ax[2].grid(alpha=.3); ax[2].set_title("scan after flat-event removal")
    fig.savefig(str(STEM) + ".png", dpi=190)
    np.savez(str(STEM) + ".npz", lags=lags, offsets=offsets, grid=VELOCITY_GRID,
             **{("%s_%s" % (k.replace(" ", "_").replace("(", "").replace(")", ""), f)): v[f]
                for k, v in res.items() for f in ("cs", "ac", "p", "corr_v")})
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
