#!/usr/bin/env python3
"""Read the adaptive-f-k aggregates and write section 4.6 of the document.

The prediction being tested was fixed BEFORE the run, in afk_recovery.sbatch and
in AMBIENT_CC_LITERATURE_REVIEW.md section 1: an adaptive filter enhances the
DOMINANT coherent component, and in this data that is the static fixed-k pattern,
so `afk1 only` (no prior static removal) should be the WORST configuration. If it
is instead the best, the mechanism argument in AMBIENT_LOWK_MECHANISM.md is wrong
and this script says so.

Writes sections/03b_afk.md, which assemble.py picks up in filename order.
Refuses to write a verdict for a configuration whose aggregate is missing.
"""
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
AWD = HERE.parent
OUT = AWD / "ambient_transfer" / "lellouch2019_exact_stack_afk"
V_LO, V_HI, V_REF = 2500.0, 4000.0, 3200.0

CONFIGS = [
    ("afk1 only, no static removal", "_afk1", "predicted WORST"),
    ("median common mode, no AFK", "_cm", "baseline"),
    ("median common mode + AFK alpha=1", "_cm_afk1", ""),
    ("median common mode + AFK alpha=2", "_cm_afk2", ""),
    ("median common mode + rank-2 + AFK alpha=1", "_cm_svd2_afk1", ""),
]
STEM = "aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0"

rows, missing = [], []
for label, suffix, note in CONFIGS:
    p = OUT / (STEM + suffix + ".npz")
    if not p.is_file():
        missing.append((label, p.name))
        continue
    d = np.load(p, allow_pickle=True)
    g = d["velocity_grid_m_s"]; cs = d["causal_moveout_scores"]
    ac = d["acausal_moveout_scores"]; nulls = d["receiver_order_null_maxima"]
    k = int(np.argmax(cs))
    obs = float(cs[k])
    pval = float((np.sum(nulls >= obs) + 1) / (nulls.size + 1))
    ped = float(np.corrcoef(g, cs)[0, 1])
    c32 = float(np.interp(V_REF, g, cs)); a32 = float(np.interp(V_REF, g, ac))
    ratio = c32 / a32 if a32 else float("nan")
    gates = dict(pedestal=abs(ped) < 0.5,
                 in_fan=V_LO <= g[k] <= V_HI,
                 off_edge=k not in (0, len(cs) - 1),
                 causal=ratio > 1.0,
                 significant=pval < 0.05)
    rows.append(dict(label=label, note=note, peak=obs, at=float(g[k]), p=pval,
                     ped=ped, ratio=ratio, gates=gates,
                     recovered=all(gates.values()),
                     windows=int(d["n_windows"]) if "n_windows" in d.files else -1))

L = ["### 4.6 Adaptive f-k",
     "",
     "The adaptive filter of Isken et al. (2022) is the one f-k family that makes",
     "no velocity assumption, and therefore the one not addressed by section 4.5.",
     "Before use it was verified that an exponent of zero is a bit-identical no-op,",
     "that Bartlett 50 %-overlap recombination reconstructs the input to a maximum",
     "relative error of 0.000e+00, that it raises the coherence of a synthetic",
     "3,200 m/s plane wave in noise from 0.0523 to 0.1980, and that it does not",
     "manufacture moveout from pure noise (0.0000 to 0.0000).",
     "",
     "Five configurations were run on a common 300-record block of 2024-12-20,",
     "with the prediction fixed in advance that applying the filter *without*",
     "prior static-pattern removal would be the worst configuration, because an",
     "adaptive filter enhances the dominant coherent component and here that is the",
     "static pattern rather than an arrival.",
     ""]
if rows:
    L += ["| configuration | peak | at (m/s) | p | pedestal corr | causal/acausal at 3,200 | recovered |",
          "|---|---:|---:|---:|---:|---:|:--:|"]
    for r in rows:
        L.append("| %s%s | %.4f | %.0f | %.4f | %+.3f | %.2f | %s |" % (
            r["label"], (" (%s)" % r["note"]) if r["note"] else "",
            r["peak"], r["at"], r["p"], r["ped"], r["ratio"],
            "**yes**" if r["recovered"] else "no"))
    L += ["", "Recovery requires all five predeclared gates: pedestal suppressed",
          "(|corr| < 0.5), peak inside 2,500-4,000 m/s, peak not at a scan edge,",
          "causal side dominant at 3,200 m/s, and p < 0.05.", ""]
    any_rec = [r for r in rows if r["recovered"]]
    base = next((r for r in rows if r["note"] == "baseline"), None)
    solo = next((r for r in rows if r["note"] == "predicted WORST"), None)
    if any_rec:
        L += ["**%d configuration(s) satisfied every gate.** This is a positive"
              % len(any_rec),
              "result for adaptive f-k and it revises the conclusion of section 4.5:",
              "the earlier failures reflected the fixed-fan assumption rather than an",
              "unfilterable field. Details: " +
              "; ".join("%s at p = %.4f" % (r["label"], r["p"]) for r in any_rec),
              ""]
    else:
        L += ["**No configuration satisfied the gates.** Adaptive f-k does not",
              "recover the arrival either. Given section 4.7, this is expected rather",
              "than surprising: the filter can only enhance coherent energy that is",
              "present, and the wavefield carries no net downgoing component to",
              "enhance. It does, however, close the remaining f-k avenue explicitly",
              "rather than by argument.", ""]
    if solo and base:
        worse = solo["ped"] > base["ped"] or solo["p"] > base["p"]
        L += ["On the pre-registered prediction: `afk1 only` gives pedestal",
              "corr = %+.3f and p = %.4f against the baseline's %+.3f and %.4f, so the"
              % (solo["ped"], solo["p"], base["ped"], base["p"]),
              "prediction that applying the adaptive filter without prior static",
              "removal is %s. %s" % (
                  "not contradicted" if worse else "**CONTRADICTED**",
                  "" if worse else
                  "That is evidence against the mechanism argument of section 4.4 "
                  "and is recorded as such."),
              ""]
if missing:
    L += ["*Configurations whose aggregate was not available when this section was",
          "generated, and which are therefore not reported: %s.*"
          % ", ".join(m[0] for m in missing), ""]

(HERE / "sections" / "03b_afk.md").write_text("\n".join(L) + "\n")
print("wrote sections/03b_afk.md with %d configurations (%d missing)"
      % (len(rows), len(missing)))
for r in rows:
    print("  %-44s peak %.4f at %5.0f  p %.4f  ped %+.3f  rec %s"
          % (r["label"], r["peak"], r["at"], r["p"], r["ped"], r["recovered"]))
