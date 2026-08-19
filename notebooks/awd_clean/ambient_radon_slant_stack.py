#!/usr/bin/env python3
"""Linear Radon (tau-p) analysis of the Figure 7c correlation gather.

WHY tau-p RATHER THAN F-K.  The VSP wavefield-separation literature is explicit
that f-k separation fails where up- and down-going energy overlap near ZERO
apparent velocity, producing leakage and ringing (e.g. adaptive masking-filter
work in J. Geophys. Eng.; Radon-constrained median filtering). That is precisely
this dataset's problem:

  - The correlation gather is dominated by a broad zero-lag lobe, i.e. energy at
    slowness p ~ 0 (infinite apparent velocity).
  - The array aperture is 900 x 1.0209 m = 918.9 m, so dk = 0.00109 cycles/m and
    the 2.5-4.5 km/s fan spans k = 0.00111-0.00200 at 5 Hz -- UNDER ONE resolution
    cell wide, and only 1.4 cells from k = 0. F-K cannot separate the target from
    the p ~ 0 pedestal at the low end of the band.
  - Measured consequence: the moveout score is ~97 % pedestal,
    corr(trial velocity, score) = +0.976, and its maximum is pinned to whatever
    velocity the scan grid stops at.

The tau-p transform scans slowness DIRECTLY, so p ~ 0 can be muted explicitly
rather than being the axis where the transform degenerates. Radon "exploits
linear coherence of events with different apparent velocities and is
theoretically a better separator than either f-k or median filtering."

WHAT THIS DOES

  1. Forward linear Radon of the R+-10 correlation gather:
         R(p, tau) = sum_x d(x, tau + p*x)
     over the paper's 50-700 m offsets. Unlike the existing moveout score -- which
     samples ONE lag per offset for each trial velocity and takes a median -- the
     slant stack retains the tau axis, so a p ~ 0 lobe at tau ~ 0 and a genuine
     3.2 km/s arrival at tau = x/3200 occupy DIFFERENT cells and can be told apart.
  2. Mutes a declared band around p = 0 and reports the residual energy along the
     3.2 km/s trajectory.
  3. Tests it with the same receiver-order permutation null used everywhere else
     in this tree, so the result is comparable to the config 0/3 p-values.

PRECONDITIONS, checked and reported. This script refuses to declare a recovery
unless ALL of:
    (a) the p-mute actually suppresses the pedestal
        -- corr(trial velocity, energy) must fall below 0.5 in magnitude;
    (b) the peak slowness lies inside the 2500-4000 m/s fan, not at a scan edge;
    (c) the causal side dominates the acausal side at 3200 m/s;
    (d) familywise p < 0.05 against the receiver-order null.
Three automated verdicts were printed from broken inputs on 2026-08-14; the
precondition gate exists so that cannot happen again.

Reads the already-computed aggregates -- no raw data, no reprocessing.

Output: ambient_radon_slant_stack.{npz,png,txt}
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
STEM = HERE / "ambient_radon_slant_stack"
AGG = HERE / "ambient_transfer" / "lellouch2019_exact_stack"

V_REF = 3200.0
V_LO, V_HI = 2500.0, 4000.0          # the physically admissible fan
NULL_COUNT = 10000
SEED = 20260814
# slowness grid: 1/6000 .. 1/1200 s/m, plus the p~0 region so the pedestal is
# visible rather than hidden by the grid choice
P_GRID = np.linspace(0.0, 1.0 / 1200.0, 241)


def slant_stack(gather, lags, offsets, p_grid):
    """R(p, tau) = sum_x d(x, tau + p*x), linear interpolation in lag."""
    n_tau = len(lags)
    out = np.zeros((len(p_grid), n_tau))
    for i, p in enumerate(p_grid):
        acc = np.zeros(n_tau)
        for row, x in zip(gather, offsets):
            shifted = np.interp(lags + p * x, lags, row, left=0.0, right=0.0)
            acc += shifted
        out[i] = acc / len(offsets)
    return out


def envelope(a):
    return np.abs(hilbert(a, axis=-1))


def energy_along_p(radon, lags, causal=True):
    """Envelope energy per slowness, restricted to one lag side."""
    m = (lags > 0) if causal else (lags < 0)
    return envelope(radon[:, m]).mean(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--branch", choices=("base", "cm"), default="cm",
                    help="cm = median-removed (the pedestal-suppressed branch)")
    ap.add_argument("--p-mute", type=float, default=1.0 / 8000.0,
                    help="mute |p| below this slowness (s/m); default 1/8000, i.e. "
                         "everything faster than 8 km/s is treated as pedestal")
    a = ap.parse_args()

    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    suffix = "_cm" if a.branch == "cm" else ""
    cands = ([str(AGG / "_authoritative_38993456" /
                  ("aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0%s.npz" % suffix))]
             + sorted(glob.glob(str(AGG / ("aggregate_2024-12-20_*_ordered_r0%s.npz" % suffix)))))
    src = next((c for c in cands if Path(c).is_file()), None)
    if src is None:
        raise SystemExit("no aggregate for branch %s" % a.branch)

    d = np.load(src, allow_pickle=True)
    lags = d["lags_s"] if "lags_s" in d.files else d["lags"]
    gather = d["r_plus_minus_10_correlation"]
    offsets = np.asarray(d["offsets_m"], dtype=float)

    say("Linear Radon (tau-p) analysis of the Figure 7c gather")
    say("  source   : %s" % Path(src).name)
    say("  branch   : %s" % ("median-removed (config 3)" if a.branch == "cm" else "baseline (config 0)"))
    say("  gather   : %d offsets %.0f-%.0f m, %d lags %.3f..%.3f s"
        % (len(offsets), offsets.min(), offsets.max(), len(lags), lags[0], lags[-1]))
    say("  p mute   : |p| < %.3e s/m  (faster than %.0f m/s)" % (a.p_mute, 1.0 / a.p_mute))
    say("  target   : p = %.3e s/m (%.0f m/s)" % (1.0 / V_REF, V_REF))
    say("")

    radon = slant_stack(gather, lags, offsets, P_GRID)
    keep = P_GRID >= a.p_mute
    p_kept = P_GRID[keep]
    v_kept = 1.0 / np.maximum(p_kept, 1e-12)

    e_all = energy_along_p(radon, lags, causal=True)
    e_causal = e_all[keep]
    e_acausal = energy_along_p(radon, lags, causal=False)[keep]

    # (a) has the pedestal actually gone?
    corr_before = float(np.corrcoef(1.0 / np.maximum(P_GRID[1:], 1e-12), e_all[1:])[0, 1])
    corr_after = float(np.corrcoef(v_kept, e_causal)[0, 1])
    say("--- precondition (a): pedestal suppression ---")
    say("  corr(apparent velocity, energy) before mute : %+.3f" % corr_before)
    say("  corr(apparent velocity, energy) after  mute : %+.3f" % corr_after)
    pedestal_gone = abs(corr_after) < 0.5
    say("  pedestal suppressed (|corr| < 0.5): %s" % pedestal_gone)

    k = int(np.argmax(e_causal))
    v_peak = float(v_kept[k])
    e_ref = float(np.interp(V_REF, v_kept[::-1], e_causal[::-1]))
    a_ref = float(np.interp(V_REF, v_kept[::-1], e_acausal[::-1]))
    say("")
    say("--- observation ---")
    say("  peak causal energy %.5f at %.0f m/s" % (e_causal[k], v_peak))
    say("  at %.0f m/s: causal %.5f, acausal %.5f, ratio %.2f"
        % (V_REF, e_ref, a_ref, e_ref / a_ref if a_ref else np.nan))

    rng = np.random.default_rng(SEED)
    nulls = np.empty(NULL_COUNT)
    for i in range(NULL_COUNT):
        perm = rng.permutation(len(offsets))
        r = slant_stack(gather[perm], lags, offsets, p_kept)
        nulls[i] = energy_along_p(r, lags, causal=True).max()
    null95 = float(np.percentile(nulls, 95))
    p_value = float((np.sum(nulls >= e_causal[k]) + 1) / (NULL_COUNT + 1))
    say("")
    say("--- receiver-order permutation null (%d) ---" % NULL_COUNT)
    say("  null 95th %.5f | observed %.5f | p = %.4f" % (null95, e_causal[k], p_value))

    in_fan = V_LO <= v_peak <= V_HI
    not_edge = (k not in (0, len(e_causal) - 1))
    causal_dom = (e_ref / a_ref) > 1.0 if a_ref else False
    say("")
    say("--- preconditions ---")
    say("  (a) pedestal suppressed        : %s" % pedestal_gone)
    say("  (b) peak in %.0f-%.0f m/s fan   : %s (peak %.0f m/s)" % (V_LO, V_HI, in_fan, v_peak))
    say("      peak not at a scan edge    : %s" % not_edge)
    say("  (c) causal dominates at %.0f    : %s (%.2f)" % (V_REF, causal_dom, e_ref / a_ref if a_ref else np.nan))
    say("  (d) familywise p < 0.05        : %s (%.4f)" % (p_value < 0.05, p_value))
    say("")
    if all((pedestal_gone, in_fan, not_edge, causal_dom, p_value < 0.05)):
        say("  RECOVERED: a causal arrival inside the physical fan clears the")
        say("  receiver-order null with the pedestal demonstrably suppressed.")
    else:
        say("  NOT RECOVERED. Failing preconditions are listed above; a p-value")
        say("  alone is not sufficient and is not reported as a result when the")
        say("  pedestal or fan conditions fail.")

    fig, ax = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    im = ax[0].imshow(envelope(radon), aspect="auto", origin="lower", cmap="magma",
                      extent=[lags[0], lags[-1], P_GRID[0] * 1e3, P_GRID[-1] * 1e3])
    ax[0].axhline(a.p_mute * 1e3, color="cyan", ls="--", lw=1, label="p mute")
    ax[0].axhline(1e3 / V_REF, color="w", ls=":", lw=1.2, label="3200 m/s")
    ax[0].set_xlabel("tau (s)"); ax[0].set_ylabel("slowness (ms/m)")
    ax[0].legend(fontsize=7); ax[0].set_title("tau-p panel (envelope)")
    fig.colorbar(im, ax=ax[0])
    ax[1].plot(v_kept, e_causal, "k-", label="causal")
    ax[1].plot(v_kept, e_acausal, color="gray", ls="--", lw=0.9, label="acausal")
    ax[1].axhline(null95, color="crimson", ls="--", lw=1, label="null 95th")
    ax[1].axvline(V_REF, color="steelblue", ls=":", label="3200 m/s")
    ax[1].axvspan(V_LO, V_HI, color="orange", alpha=.15)
    ax[1].set_xlabel("apparent velocity (m/s)"); ax[1].set_ylabel("Radon envelope energy")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3); ax[1].set_title("energy vs slowness, p~0 muted")
    ax[2].hist(nulls, bins=60, color="steelblue", alpha=.8)
    ax[2].axvline(e_causal[k], color="k", lw=2, label="observed")
    ax[2].axvline(null95, color="crimson", ls="--", label="null 95th")
    ax[2].set_xlabel("null max energy"); ax[2].legend(fontsize=8)
    ax[2].set_title("receiver-order null, p = %.4f" % p_value)
    fig.savefig(str(STEM) + ".png", dpi=190)

    np.savez(str(STEM) + ".npz", radon=radon, p_grid=P_GRID, lags=lags,
             offsets=offsets, v_kept=v_kept, e_causal=e_causal,
             e_acausal=e_acausal, nulls=nulls, p_value=p_value,
             branch=a.branch, p_mute=a.p_mute)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
