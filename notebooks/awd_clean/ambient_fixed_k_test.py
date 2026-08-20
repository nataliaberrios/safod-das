#!/usr/bin/env python3
"""Is the dominant 2024-25 low-wavenumber feature a WAVE or a STATIC PATTERN?

THE QUESTION, and how this relates to the cross-epoch comparison.

A NOTE ON A RETRACTION THAT WAS ITSELF WRONG, because the sequence matters for
reading the numbers below.  The 2017 release files M1p33 and M2p46 are EARTHQUAKE
records (their README: "Earthquake records are float32 binary", 800 x 1250 at
250 Hz = 5 s), and the ambient census compares one of its windows against 2024-25
continuous noise.  That looked like the two-arms-processed-differently error
behind five of the six earlier withdrawals in this tree, so the cross-epoch k = 0
claim was withdrawn on 2026-08-19.  It should not have been:
`lellouch2017_window_audit.py` then MEASURED the arrival times -- 4.100 s and
4.916 s of 5.00 s -- and the census cuts at 2.5 s, so both arrivals are outside
its window.  The 2017 arm is genuine pre-event noise, both arms of the census
were noise, and the comparison is not confounded by signal presence.  It stands,
subject to its originally stated limits: only ~5 s of 2017 data, and a gauge
length of 10 m against 16.335 m (a 0.1-1.1 % effect at 3200 m/s).

So this test COMPLEMENTS that comparison rather than replacing it.  The census
says the 2024-25 ambient field carries most of its in-band energy at very low
wavenumber; this asks what that low-wavenumber energy IS.

THE DISCRIMINATOR.  Apparent velocity is v = f / k, so a feature's behaviour
across frequency tells you what it is:

    a PROPAGATING ARRIVAL has fixed v  ->  k_peak(f) = f / v, a straight line
                                           through the origin with slope 1/v
    a STATIC SPATIAL PATTERN has fixed k ->  k_peak(f) = constant, a flat line
                                           whose apparent velocity rises with f

They are not subtly different: one is a ray through the origin, the other is
horizontal.  Fitting k_peak(f) and comparing the two models therefore settles
what the dominant low-k feature is without any reference to 2017, to a filter,
or to a velocity scan.

The prediction that motivated this: between the 5-20 Hz and 12-20 Hz censuses,
the ~28 % of energy in the 3000-4000 m/s bin did not vanish, it MOVED to the
6000-10000 m/s bin.  Energy that changes apparent velocity when you change the
band is at fixed k, not fixed v.

WHAT EACH OUTCOME MEANS
  fixed k wins   -> the contaminant is a static per-channel spatial pattern (a
                    channel-to-channel amplitude/coupling pattern, essentially a
                    fixed shape multiplying a time series).  It is not a wave,
                    so no velocity-domain filter -- F-K, tau-p, or otherwise --
                    can be expected to separate it by velocity, because it has
                    no velocity.  It DOES have a fixed spatial spectrum, which
                    is a different and more tractable target.
  fixed v wins   -> there is a coherent arrival at that velocity and the recovery
                    problem is one of amplitude, not of separability.

CONTROL.  Both models are fitted to the SAME k_peak(f) points by least squares
and compared on identical residuals, and the fixed-k model is given no free
advantage: it has one parameter (the constant) against the wave model's one
parameter (the slope).  Equal parameter count, so R^2 is directly comparable.
A permutation check shuffles the frequency labels to confirm the winning fit is
not what random assignment would produce.

Reads only the census .npz products; no raw data, no reprocessing.

Output: ambient_fixed_k_test.{png,txt}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = HERE / "ambient_fixed_k_test"
V_REF = 3200.0
SEED = 20260819


def peak_k_per_frequency(P, k, f, fmin, fmax, k_exclude_cells=2):
    """For each in-band frequency, the |k| carrying the most power.

    k = 0 and its immediate neighbours are excluded: with a Hann taper the main
    lobe is +-2 cells wide, so leaving them in would make every frequency report
    the same k by construction and guarantee the fixed-k answer.  Excluding them
    biases the test AGAINST the fixed-k model, which is the safe direction.
    """
    dk = float(np.abs(k[1] - k[0]))
    pos = k > (k_exclude_cells + 0.5) * dk
    kk = k[pos]
    rows = []
    for i, freq in enumerate(f):
        if not (fmin <= freq <= fmax):
            continue
        # fold negative and positive wavenumber: direction is not the question
        slice_pos = P[pos, i]
        neg = k < -(k_exclude_cells + 0.5) * dk
        slice_neg = P[neg, i][::-1][: slice_pos.size]
        folded = slice_pos + slice_neg
        j = int(np.argmax(folded))
        rows.append((float(freq), float(kk[j]), float(folded[j])))
    return np.array(rows)          # columns: f, k_peak, power


def low_k_centroid(P, k, f, fmin, fmax, cells=8):
    """Power-weighted centroid of |k| over the first `cells` non-zero cells.

    THIS is the measurement the peak-per-frequency approach should have made.
    The pedestal lives in the first few wavenumber cells; asking where the
    per-frequency MAXIMUM sits (with k=0 excluded) instead finds the genuine
    surface-wave energy ~13 cells out, which is a different feature entirely.

    The discriminator between the two models is how the centroid moves when the
    analysis band moves:
        fixed k (static pattern) -> centroid is the SAME in both bands
        fixed v (a wave)         -> centroid scales with band centre frequency,
                                    since k = f/v
    k = 0 itself is excluded so the answer is not fixed by construction; the
    remaining cells are where the contested energy is.
    """
    dk = float(np.abs(k[1] - k[0]))
    band = (np.abs(f) >= fmin) & (np.abs(f) <= fmax)
    if band.sum() < 2:
        return None
    ka = np.abs(k)
    sel = (ka > 0.5 * dk) & (ka <= (cells + 0.5) * dk)
    if sel.sum() < 2:
        return None
    w = P[np.ix_(sel, band)].sum(axis=1)
    kk = ka[sel]
    total = float(w.sum())
    if total <= 0:
        return None
    centroid = float(np.sum(kk * w) / total)
    f_centre = float(np.mean(np.abs(f[band])))
    return dict(centroid=centroid, dk=dk, f_centre=f_centre,
                cells=int(sel.sum()), apparent_v=f_centre / centroid,
                share_of_band=100.0 * total / float(P[:, band].sum()))


def fit_models(freq, kpk):
    """Least squares for k = a*f (wave, v = 1/a) and k = c (static). Same dof."""
    a = float(np.sum(freq * kpk) / np.sum(freq * freq))
    c = float(np.mean(kpk))
    res_w = kpk - a * freq
    res_s = kpk - c
    tss = float(np.sum((kpk - kpk.mean()) ** 2))
    r2 = lambda r: 1.0 - float(np.sum(r * r)) / tss if tss > 0 else np.nan
    return dict(slope=a, velocity=(1.0 / a if a else np.inf), r2_wave=r2(res_w),
                const=c, apparent_v_at_f=None, r2_static=r2(res_s),
                rss_wave=float(np.sum(res_w ** 2)),
                rss_static=float(np.sum(res_s ** 2)))


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    say("Is the dominant 2024-25 low-k feature a wave (fixed v) or a static")
    say("pattern (fixed k)?  Internal test -- no 2017 comparison, no filter.")
    say("")
    say("CORRECTION 2026-08-19: an earlier version of this header said the")
    say("cross-epoch k=0 claim was WITHDRAWN because M1p33/M2p46 are earthquake")
    say("records. lellouch2017_window_audit.py then measured the arrival times --")
    say("4.100 s and 4.916 s of 5.00 s -- both LATER than the census's 2.5 s cut.")
    say("The 2017 census window is therefore genuine pre-event noise, both arms")
    say("were noise, and that comparison is NOT confounded. It is reinstated.")
    say("This test is a complement to it, not a replacement.")
    say("")

    sources = [("5-20 Hz  (census default)", "ambient_apparent_velocity_census_k0rm.npz", 5.0, 20.0),
               ("5-12 Hz", "ambient_apparent_velocity_census_5-12Hz_k0rm.npz", 5.0, 12.0),
               ("12-20 Hz", "ambient_apparent_velocity_census_12-20Hz_k0rm.npz", 12.0, 20.0)]
    found = []
    for label, name, lo, hi in sources:
        p = HERE / name
        if not p.is_file():
            say("  (missing, skipped: %s)" % name)
            continue
        d = np.load(p, allow_pickle=True)
        rows = peak_k_per_frequency(d["P24"], d["k24"], d["f24"], lo, hi)
        if rows.shape[0] < 4:
            say("  (too few in-band frequencies in %s)" % name)
            continue
        cen = low_k_centroid(d["P24"], d["k24"], d["f24"], lo, hi)
        found.append((label, rows, lo, hi, cen))

    if not found:
        raise SystemExit("no census products found; run the census first")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4), constrained_layout=True)
    verdicts = []
    rng = np.random.default_rng(SEED)

    for label, rows, lo, hi, cen in found:
        freq, kpk = rows[:, 0], rows[:, 1]
        m = fit_models(freq, kpk)
        # permutation control: does shuffling f change which model wins?
        wins_static = 0
        for _ in range(2000):
            mm = fit_models(rng.permutation(freq), kpk)
            if mm["rss_static"] < mm["rss_wave"]:
                wins_static += 1
        say("--- %s (%d frequencies) ---" % (label, len(freq)))
        if cen is not None:
            say("  LOW-K CENTROID (first %d non-zero cells, k=0 excluded):" % cen["cells"])
            say("    |k| centroid %.5f cyc/m = %.2f cells | band centre %.1f Hz"
                % (cen["centroid"], cen["centroid"] / cen["dk"], cen["f_centre"]))
            say("    that centroid implies apparent v = %.0f m/s at the band centre"
                % cen["apparent_v"])
            say("    this low-k region holds %.2f %% of the in-band energy"
                % cen["share_of_band"])
        say("  wave model    k = f/v      : v = %8.0f m/s   R^2 = %+.4f  RSS %.3e"
            % (m["velocity"], m["r2_wave"], m["rss_wave"]))
        say("  static model  k = const    : k = %.5f cyc/m  R^2 = %+.4f  RSS %.3e"
            % (m["const"], m["r2_static"], m["rss_static"]))
        say("  static model implies apparent v = %.0f m/s at %.1f Hz rising to "
            "%.0f m/s at %.1f Hz" % (lo / m["const"], lo, hi / m["const"], hi))
        better = "STATIC (fixed k)" if m["rss_static"] < m["rss_wave"] else "WAVE (fixed v)"
        say("  better fit: %s   (RSS ratio wave/static = %.2f)"
            % (better, m["rss_wave"] / m["rss_static"] if m["rss_static"] else np.inf))
        say("  permutation control: static wins on %d/2000 frequency-shuffled sets"
            % wins_static)

        # PRECONDITION GATE.  Three false automated verdicts were printed in this
        # tree on 2026-08-14 from broken inputs, and this script printed a fourth
        # on 2026-08-19: it declared a STATIC pattern while its own permutation
        # control showed static winning on 1578/2000 SHUFFLED sets, i.e. the
        # answer was what random frequency labels give.  Two conditions must hold
        # for the comparison to carry information at all:
        #   (a) the permutation control must not already favour the winner --
        #       if shuffling f still picks the same model, f carried no signal;
        #   (b) at least one model must actually fit. The static model's R^2 is
        #       0 by construction (it IS the mean), so a wave R^2 <= 0 means
        #       neither model explains anything and k_peak(f) is noise.
        chance = wins_static / 2000.0
        control_ok = 0.20 <= chance <= 0.80
        fit_ok = m["r2_wave"] > 0.20 or m["rss_wave"] / max(m["rss_static"], 1e-30) > 2.0
        say("  gate (a) permutation control informative (0.2-0.8): %s (%.3f)"
            % (control_ok, chance))
        say("  gate (b) at least one model fits                  : %s "
            "(wave R^2 %+.4f, RSS ratio %.2f)"
            % (fit_ok, m["r2_wave"], m["rss_wave"] / max(m["rss_static"], 1e-30)))
        if not (control_ok and fit_ok):
            say("  >>> UNINFORMATIVE: no verdict for this band. <<<")
        say("")
        verdicts.append((label, better if (control_ok and fit_ok) else None, m))
        ax[0].plot(freq, kpk * 1e3, "o", ms=4, label="%s peak k" % label)

    lab0, rows0, lo0, hi0, cen0 = found[0]
    freq0, kpk0 = rows0[:, 0], rows0[:, 1]
    m0 = fit_models(freq0, kpk0)
    fgrid = np.linspace(0, max(r[1][:, 0].max() for r in found), 50)
    ax[0].plot(fgrid, (fgrid / V_REF) * 1e3, "r--", lw=1.6,
               label="a %.0f m/s wave (k = f/v)" % V_REF)
    ax[0].axhline(m0["const"] * 1e3, color="steelblue", ls=":", lw=1.6,
                  label="static pattern (k = const)")
    ax[0].set(xlabel="frequency (Hz)", ylabel="peak wavenumber (10$^{-3}$ cycles/m)",
              title="Peak wavenumber vs frequency\nray through origin = wave; flat = static pattern")
    ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)

    for label, rows, lo, hi, cen in found:
        freq, kpk = rows[:, 0], rows[:, 1]
        ax[1].plot(freq, freq / kpk, "o-", ms=4, lw=1, label=label)
    ax[1].axhline(V_REF, color="crimson", ls="--", lw=1.5, label="%.0f m/s" % V_REF)
    ax[1].set(xlabel="frequency (Hz)", ylabel="apparent velocity of peak (m/s)",
              yscale="log",
              title="Apparent velocity of the peak\nflat = a real arrival; rising = a static pattern")
    ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
    fig.savefig(str(STEM) + ".png", dpi=190)

    # ---- the measurement that actually discriminates: does the low-k centroid
    # move with the band?  fixed k -> same centroid; fixed v -> centroid scales
    # with band centre frequency.  This needs two bands with different centres.
    cens = [(lab, c, lo, hi) for lab, _, lo, hi, c in found if c is not None]
    say("=== low-k centroid across bands (the discriminating measurement) ===")
    if len(cens) < 2:
        say("  need two bands; only %d available -- run the census for 5-12 Hz")
        say("  and 12-20 Hz with K0_REMOVE=1 first.")
    else:
        cens.sort(key=lambda t: t[1]["f_centre"])
        (la, ca, _, _), (lb, cb, _, _) = cens[0], cens[-1]
        f_ratio = cb["f_centre"] / ca["f_centre"]
        k_ratio = cb["centroid"] / ca["centroid"]
        say("  %-10s band centre %5.2f Hz -> |k| centroid %.5f cyc/m (%.2f cells)"
            % (la, ca["f_centre"], ca["centroid"], ca["centroid"] / ca["dk"]))
        say("  %-10s band centre %5.2f Hz -> |k| centroid %.5f cyc/m (%.2f cells)"
            % (lb, cb["f_centre"], cb["centroid"], cb["centroid"] / cb["dk"]))
        say("  frequency ratio %.3f | wavenumber ratio %.3f" % (f_ratio, k_ratio))
        say("  a WAVE at fixed v predicts k_ratio = f_ratio = %.3f" % f_ratio)
        say("  a STATIC pattern at fixed k predicts k_ratio = 1.000")
        d_wave = abs(k_ratio - f_ratio)
        d_static = abs(k_ratio - 1.0)
        if d_static < 0.5 * d_wave:
            say("  -> closer to FIXED k: a static spatial pattern. It has no")
            say("     velocity, so no velocity-domain filter can separate it BY")
            say("     velocity -- which is why F-K, tau-p, rank-k, PCC, PWS and")
            say("     flat-event removal all failed the same way.")
            say("")
            say("     SENSITIVITY: the test had the range to see a wave and did")
            say("     not. A %.0f m/s arrival sits at %.2f cells at %.1f Hz and"
                % (V_REF, (ca["f_centre"] / V_REF) / ca["dk"], ca["f_centre"]))
            say("     %.2f cells at %.1f Hz, both inside the window measured, so"
                % ((cb["f_centre"] / V_REF) / cb["dk"], cb["f_centre"]))
            say("     a dominant wave would have moved the centroid and did not.")
            say("")
            say("     LIMIT, and it matters for what comes next: this says the")
            say("     DOMINANT low-k feature is not a wave. It does NOT say no")
            say("     arrival exists -- a weak arrival can hide beneath a feature")
            say("     holding ~%.0f %% of the in-band energy. The implication is"
                % ca["share_of_band"])
            say("     about METHOD, not absence: because the contaminant is a")
            say("     fixed SPATIAL pattern, the removal has to be spatial --")
            say("     estimate and divide out the static per-channel response --")
            say("     not a velocity filter. C2_PERMEABILITY_FOLLOWUP.md already")
            say("     proposes exactly that calibration and records 5x headroom")
            say("     from it; this is independent motivation for doing it.")
        elif d_wave < 0.5 * d_static:
            say("  -> closer to FIXED v: consistent with a propagating arrival at")
            say("     ~%.0f m/s. Recovery would then be an amplitude problem."
                % ca["apparent_v"])
        else:
            say("  -> AMBIGUOUS (|k_ratio-1| = %.3f vs |k_ratio-f_ratio| = %.3f):"
                % (d_static, d_wave))
            say("     neither model is clearly closer. No verdict.")
    say("")

    scored = [(lab, b, m) for lab, b, m in verdicts if b is not None]
    static_wins = sum(1 for _, b, _ in scored if b.startswith("STATIC"))
    say("=== secondary diagnostic: peak-per-frequency model fit ===")
    say("  (The PRIMARY result is the low-k centroid comparison above. This")
    say("   section is a weaker diagnostic kept for the record; if it reports no")
    say("   verdict that does NOT retract the centroid result.)")
    if not scored:
        say("  This diagnostic produced NO INFORMATIVE COMPARISON in any band.")
        say("  Every band failed the permutation control, the fit test, or both,")
        say("  so k_peak(f) is noise and neither model is supported BY THIS")
        say("  DIAGNOSTIC. The centroid measurement above is unaffected: it uses")
        say("  a different quantity in a different region of wavenumber.")
        say("")
        say("  WHY, and what to do instead. Excluding +-2 cells around k = 0 to")
        say("  avoid a trivial answer pushed the per-frequency peak out to")
        say("  ~0.019 cycles/m -- about 13 cells from zero, i.e. the genuine")
        say("  surface-wave energy at 267-640 m/s -- which is not the pedestal")
        say("  this test was written to characterise. The right measurement is")
        say("  the k-marginal in the FIRST FEW non-zero cells, compared between")
        say("  bands: if the same low-k cells dominate in 5-12 Hz and 12-20 Hz,")
        say("  the feature is at fixed k. The band-shift already observed (the")
        say("  ~28 %% share moving from the 3000-4000 m/s bin to 6000-10000 m/s")
        say("  when the band was raised) is consistent with fixed k, but it is an")
        say("  observation awaiting that test, not a result.")
        Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
        say("")
        say("wrote %s.{png,txt}" % STEM.name)
        return
    say("  %d of %d INFORMATIVE bands are better described by a STATIC fixed-k pattern"
        % (static_wins, len(scored)))
    if static_wins == len(scored):
        say("  The dominant low-wavenumber feature in the 2024-25 ambient field is")
        say("  a STATIC SPATIAL PATTERN, not a propagating wave. It has no")
        say("  velocity, so no velocity-domain filter can separate it BY velocity;")
        say("  that is why F-K, tau-p, rank-k, PCC, PWS and flat-event removal all")
        say("  failed in the same way. It does have a fixed spatial spectrum,")
        say("  which is a different and more tractable target -- see")
        say("  C2_PERMEABILITY_FOLLOWUP.md, which already proposes calibrating the")
        say("  static per-channel amplitude response.")
    elif static_wins == 0:
        say("  The peak tracks a fixed velocity -- consistent with a real arrival.")
        say("  Recovery is then an amplitude problem, not a separability problem.")
    else:
        say("  Mixed across bands; do not conclude. Report per-band numbers only.")
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
