#!/usr/bin/env python3
"""Is the ambient field illuminated from the surface?  Lellouch's own diagnostic.

WHY THIS IS THE TEST THAT DECIDES WHETHER RECOVERY IS POSSIBLE.

Lellouch et al. (2019, JGR 10.1029/2019JB017533) did not assume a downgoing P
existed -- they inferred it.  Their stated evidence that the dominant ambient
sources are at the surface is the AMPLITUDE ASYMMETRY between downgoing and
upgoing energy: "the strong signals of these downgoing waves compared to the
upgoing ones suggest that the dominant ambient field sources originate at the
surface" (citing Zhou & Paulssen 2017).  Behm (2016, Geophys. Prospect.
10.1111/1365-2478.12424) reaches the same conclusion independently for borehole
arrays in a producing oil field: the interferograms "clearly indicate an origin
of the ambient seismic energy from above the arrays, thus suggesting surface
activities as sources", and there body-wave retrieval worked from as little as
30 s of noise.

That last number is the important one for us.  Under good illumination, 30 s is
enough.  Our stacks go the other way: 6 days give min p = 0.1345 and a coherent
96 h stack gives p = 0.9184 -- MORE data, WORSE result.  A signal that never
converges under stacking is not a weak signal, it is an absent one.  The
mechanism for absence would be missing illumination.

THE MEASUREMENT, and why it needs no sign convention.  A one-sided illumination
produces energy travelling predominantly ONE way along the fibre, i.e. an
imbalance between +k and -k at positive frequency.  Define

    A = (E(+k) - E(-k)) / (E(+k) + E(-k))        over the body-wave fan

|A| near 1 means strongly one-directional; |A| near 0 means balanced.  Which sign
corresponds to "downgoing" depends on the channel-ordering convention, and this
project has an explicit constant for that (BRANCH_LAG_SIGN in
ambient_signed_fk_v2.py) precisely because it is easy to get wrong -- so this
script deliberately reports only |A|, which is convention-free.  The question
"is there a preferred direction at all" does not require knowing which one.

CONTROL.  |A| is not zero for a finite sample even when the field is balanced,
because of estimation noise.  The null is built by the same statistic on a
phase-randomised version of the same spectrum, which preserves the power
spectrum and destroys any directional preference, repeated to give a
distribution.  A measured |A| inside that distribution is not evidence of
illumination.

Both epochs are measured, from the SAME census products, so the 2017 pre-event
noise (verified genuine pre-event by lellouch2017_window_audit.py: arrivals at
4.100 s and 4.916 s of 5.00 s, both after the 2.5 s cut) acts as the positive
control -- Lellouch reported asymmetry there, so this arm should show it.

Reads census .npz only. Output: ambient_directional_asymmetry.{png,txt}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = HERE / "ambient_directional_asymmetry"
V_LO, V_HI = 2500.0, 4000.0
NULL_COUNT = 400
SEED = 20260819


def fan_asymmetry(P, k, f, fmin, fmax, v_lo=V_LO, v_hi=V_HI):
    """|A| = |E(+k) - E(-k)| / (E(+k) + E(-k)) over the body-wave fan.

    POSITIVE FREQUENCIES ONLY.  For real input, P(-k,-f) = P(k,f), so a budget
    taken over all f would make +k and -k identical by construction and force
    A = 0.  Restricting to f > 0 makes the two signs physically distinct
    directions.  Getting this wrong would produce a guaranteed null result that
    looks like a finding.
    """
    K, F = np.meshgrid(k, f, indexing="ij")
    pos_f = (F >= fmin) & (F <= fmax)          # strictly positive frequencies
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(np.abs(K) > 0, F / np.abs(K), np.inf)
    fan = pos_f & (v >= v_lo) & (v <= v_hi)
    up = fan & (K > 0)
    dn = fan & (K < 0)
    e_up, e_dn = float(P[up].sum()), float(P[dn].sum())
    total = e_up + e_dn
    if total <= 0 or up.sum() < 4 or dn.sum() < 4:
        return None
    return dict(asym=abs(e_up - e_dn) / total, e_plus=e_up, e_minus=e_dn,
                cells=int(up.sum() + dn.sum()),
                fan_share=100.0 * total / float(P[pos_f].sum()))


def null_asymmetry(P, k, f, fmin, fmax, rng, n=NULL_COUNT):
    """Same statistic on spectra whose +k/-k assignment is randomised.

    The power values inside the fan are shuffled between the +k and -k halves.
    This preserves the total fan energy and the marginal distribution of cell
    powers while destroying any directional preference, which is exactly the
    null hypothesis "no preferred propagation direction".
    """
    K, F = np.meshgrid(k, f, indexing="ij")
    pos_f = (F >= fmin) & (F <= fmax)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(np.abs(K) > 0, F / np.abs(K), np.inf)
    fan = pos_f & (v >= V_LO) & (v <= V_HI)
    up, dn = fan & (K > 0), fan & (K < 0)
    pool = np.concatenate([P[up].ravel(), P[dn].ravel()])
    n_up = int(up.sum())
    out = np.empty(n)
    for i in range(n):
        s = rng.permutation(pool)
        a, b = float(s[:n_up].sum()), float(s[n_up:].sum())
        out[i] = abs(a - b) / (a + b) if (a + b) > 0 else 0.0
    return out


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    say("Directional asymmetry of the ambient field in the body-wave fan")
    say("  This is Lellouch et al. (2019)'s own diagnostic for surface")
    say("  illumination: downgoing energy strongly exceeding upgoing.")
    say("  Reported as |A| only -- convention-free, see module docstring.")
    say("  Fan %.0f-%.0f m/s, positive frequencies only." % (V_LO, V_HI))
    say("")

    sources = [("5-12 Hz", "ambient_apparent_velocity_census_5-12Hz_k0rm.npz", 5.0, 12.0),
               ("12-20 Hz", "ambient_apparent_velocity_census_12-20Hz_k0rm.npz", 12.0, 20.0)]
    rng = np.random.default_rng(SEED)
    rows = []
    for label, name, lo, hi in sources:
        p = HERE / name
        if not p.is_file():
            say("  (missing, skipped: %s)" % name)
            continue
        d = np.load(p, allow_pickle=True)
        say("=== %s ===" % label)
        for epoch, Pk, kk, ff in (("2017 pre-event", "P17", "k17", "f17"),
                                  ("2024-25 ambient", "P24", "k24", "f24")):
            if Pk not in d.files:
                say("  %-16s not in product" % epoch)
                continue
            res = fan_asymmetry(d[Pk], d[kk], d[ff], lo, hi)
            if res is None:
                say("  %-16s too few fan cells to measure" % epoch)
                continue
            nulls = null_asymmetry(d[Pk], d[kk], d[ff], lo, hi, rng)
            n95 = float(np.percentile(nulls, 95))
            pval = float((np.sum(nulls >= res["asym"]) + 1) / (len(nulls) + 1))
            say("  %-16s |A| = %.4f  (null 95th %.4f, p = %.4f)  fan holds %.2f%% "
                "of in-band energy over %d cells"
                % (epoch, res["asym"], n95, pval, res["fan_share"], res["cells"]))
            rows.append((label, epoch, res["asym"], n95, pval, res["fan_share"]))
        say("")

    if not rows:
        raise SystemExit("no census products found; run the census first")

    say("=== interpretation ===")
    got = {}
    for label, epoch, a, n95, pv, share in rows:
        got.setdefault(epoch, []).append((label, a, pv, share))
    for epoch, items in got.items():
        sig = [i for i in items if i[2] < 0.05]
        say("  %-16s significant asymmetry in %d of %d bands"
            % (epoch, len(sig), len(items)))
    e24 = got.get("2024-25 ambient", [])
    e17 = got.get("2017 pre-event", [])
    if e17 and e24:
        s17 = sum(1 for i in e17 if i[2] < 0.05)
        s24 = sum(1 for i in e24 if i[2] < 0.05)
        say("")
        if s17 > 0 and s24 == 0:
            say("  THE MISSING INGREDIENT. The 2017 pre-event noise carries a")
            say("  significant directional preference in the body-wave fan and the")
            say("  2024-25 ambient field does not. Lellouch's inference of surface")
            say("  sources rests on exactly that asymmetry, so this says the")
            say("  2024-25 field is NOT illuminated from the surface in this band.")
            say("  A downgoing P cannot be recovered from a field that has no net")
            say("  downgoing energy: no filter can create a propagation direction")
            say("  that is not in the data. This is consistent with Behm (2016),")
            say("  where 30 s sufficed under good illumination, and with our stacks")
            say("  getting WORSE with more data (96 h coherent stack p = 0.9184).")
        elif s24 > 0:
            say("  The 2024-25 field DOES carry a significant directional")
            say("  preference in the body-wave fan. Illumination is therefore not")
            say("  the blocker, and recovery should be pursued: the obstacle is the")
            say("  static fixed-k pattern (AMBIENT_LOWK_MECHANISM.md), which is a")
            say("  spatial-domain problem with a spatial-domain remedy.")
        elif s17 == 0:
            say("  NEITHER epoch shows significant asymmetry. The 2017 arm is the")
            say("  positive control -- Lellouch reported asymmetry in this data --")
            say("  so a null there means THIS MEASUREMENT lacks the sensitivity to")
            say("  decide, most likely because only ~5 s of 2017 noise exists.")
            say("  No conclusion about 2024-25 illumination may be drawn from it.")
    say("")
    say("  LIMITS: 2017 is ~5 s of pre-event noise from two records, so its")
    say("  estimate is noisy; the fan is the frozen 2500-4000 m/s selection; and")
    say("  |A| measures a preferred direction ALONG THE FIBRE, which in a near-")
    say("  vertical borehole is close to but not exactly vertical.")

    if rows:
        fig, ax = plt.subplots(figsize=(8.5, 5), constrained_layout=True)
        labels = ["%s\n%s" % (r[0], r[1]) for r in rows]
        x = np.arange(len(rows))
        ax.bar(x, [r[2] for r in rows], 0.55,
               color=["tab:blue" if "2017" in r[1] else "tab:red" for r in rows],
               label="|A| observed")
        ax.plot(x, [r[3] for r in rows], "k_", ms=26, mew=2, label="null 95th")
        for i, r in enumerate(rows):
            ax.annotate("p=%.3f" % r[4], (i, r[2]), textcoords="offset points",
                        xytext=(0, 4), ha="center", fontsize=8)
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
        ax.set(ylabel="|A|  directional asymmetry in the fan",
               title="Surface illumination test (Lellouch 2019's diagnostic)\n"
                     "high = one-directional field; at/below null = balanced")
        ax.legend(fontsize=8); ax.grid(alpha=.3, axis="y")
        fig.savefig(str(STEM) + ".png", dpi=190)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
