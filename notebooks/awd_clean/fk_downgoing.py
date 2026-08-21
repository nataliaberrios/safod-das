#!/usr/bin/env python3
"""Directional f-k separation: keep only DOWNGOING energy, with the honest control.

WHAT THIS IS, AND WHY IT IS NOT THE F-K FILTERING THIS PROJECT REJECTED.

`Ambient_FK_QC_workflow.ipynb` rejected a fixed 2.5-4.5 km/s VELOCITY FAN on the
2024-25 record: it produced an apparent ridge that survived a pre-filter
channel-scramble, i.e. the wedge was manufacturing the result. That verdict
stands and is not reopened here.

This is a different and much weaker operation. A velocity fan keeps a narrow
wedge and therefore imposes a velocity; a DIRECTIONAL filter keeps a half-plane
in wavenumber and imposes only a propagation SENSE. It is standard VSP
wavefield separation, it assumes no velocity, and on these gathers most of the
ringing is upgoing/acausal energy that it removes without touching the arrival.

Because a smaller assumption is still an assumption, the same control applies:

    THE INPUT-LEVEL NULL. The receiver order is scrambled BEFORE the filter runs,
    then the identical filter is applied. A directional filter cannot know which
    ordering is real, so if it produces comparable apparent moveout on scrambled
    input, whatever it produces on real input is the operator and not the data.
    This is precisely the gate the velocity fan failed, and it is the reason a
    post-filter null is not sufficient: the receiver-order null used elsewhere in
    this tree permutes the FINISHED gather, so an operator applied before the
    gather is formed sits outside its own null.

SIGN CONVENTION, now fixed ANALYTICALLY and cross-checked against the data.
`ambient_signed_fk_v2` carries a BRANCH_LAG_SIGN constant precisely because this
is easy to invert. numpy's `fft2` uses exp(-2i*pi*(k z + f t)), so a wave
u(z,t) = w(t - z/v) travelling towards increasing `positions` transforms to
G(f)*delta(k + f/v): k and f carry OPPOSITE signs. Downgoing is therefore the
k*f < 0 quadrant (`DOWNGOING_KF_SIGN` below; verified to 45:1 on a synthetic pure
downgoing gather by `audit_test_fk.py`). The empirical pick is still computed and
the two are compared out loud, because a disagreement means the gather is not
what the caller thinks it is.

TAPER, AND WHAT IT COSTS. A brick wall in k rings, so the quadrant edge is
tapered with a raised cosine over `--taper-cells` cells in BOTH k and f;
`--taper-cells 0` reproduces the brick wall. The cost is not cosmetic and is
reported per run: the taper is a fixed number of CELLS, and cell size is set by
the aperture, so at the production geometry (367 traces x 2.0419 m, +-0.35 s)
six cells span 0 to ~13 Hz in k and 0 to ~8.6 Hz in f. Across most of the 5-20 Hz
band the mask is therefore still inside its own transition, the two panels are
far from disjoint, and the down/up ratio is compressed towards 1. Read the
"separation contrast" table the script prints before quoting a ratio.

Operates on gathers already saved by the section scripts; no reprocessing.

Output: fk_downgoing[_<tag>].{npz,png,txt}
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import safod_geometry as geo

STEM = HERE / "fk_downgoing"
INK, MUTED = "#444444", "#6b6b6b"
C1, C2, C3 = "#0072B2", "#D55E00", "#009E73"
NULLS = 200
SEED = 20260820

# Sign of k*f occupied by energy propagating towards INCREASING `positions`,
# under numpy's exp(-2i*pi*(k z + f t)) convention. See the module docstring.
DOWNGOING_KF_SIGN = -1
KF_NEG, KF_POS = "kf<0", "kf>0"
# accepted for backwards compatibility with products written before 2026-08-20,
# when the two quadrants were mislabelled as half-planes in k
_LEGACY = {"negative-k": KF_NEG, "positive-k": KF_POS}


def _hermitian(mask):
    """Force m(-k,-f) == m(k,f), including the Nyquist row and column.

    An fftshift'd axis of EVEN length carries -k_Nyquist but not +k_Nyquist, so
    the conjugate partner of the Nyquist bin is ITSELF. A quadrant mask built
    pointwise therefore hands a bin and its own partner opposite values, the
    masked spectrum is not Hermitian, and `np.real()` silently throws the
    imaginary residue away -- measured at 1.9e-3 of the panel amplitude by
    `audit_test_fk.py` (A4/A5), spread over the k-Nyquist row and f-Nyquist
    column. Averaging each bin with its partner removes it, splits those bins
    50/50 (propagation direction is undefined at Nyquist), and preserves
    m + (1 - m) == 1 so the round-trip identity still holds exactly.
    """
    m = np.fft.ifftshift(mask)
    nk, nf = m.shape
    partner = m[(-np.arange(nk)) % nk][:, (-np.arange(nf)) % nf]
    return np.fft.fftshift(0.5 * (m + partner))


def separation_contrast(shape, dz, dt, taper_cells, v_ref, band=(5.0, 20.0)):
    """How disjoint are the two panels, as a function of frequency?

    The taper is a fixed number of CELLS and cell size is 1/aperture, so on a
    short array the transition zone can cover the whole signal band. For a wave
    at `v_ref` the mask value on its own dispersion line is the fraction of that
    wave kept by the downgoing panel; 1.0 is perfect separation and 0.5 is none.
    Returns (frequencies, kept_fraction).
    """
    nk, nf = shape
    dk, df = 1.0 / (nk * abs(dz)), 1.0 / (nf * abs(dt))
    fr = np.linspace(band[0], band[1], 16)
    kk = fr / v_ref                                  # |k| on the moveout line
    if taper_cells > 0:
        sk = np.sin(0.5 * np.pi * np.clip(kk / (taper_cells * dk), 0.0, 1.0))
        sf = np.sin(0.5 * np.pi * np.clip(fr / (taper_cells * df), 0.0, 1.0))
        return fr, 0.5 * (1.0 + sk * sf)
    return fr, np.ones_like(fr)


def separate(gather, lags, positions, taper_cells=6.0, keep="auto", v_ref=1675.0,
             quiet=False):
    """Return (downgoing, upgoing, which_quadrant_was_kept).

    `positions` is depth or along-fibre distance, uniformly sampled and
    INCREASING. A decreasing axis gives a negative dz, which negates the k axis
    and swaps the two quadrants; callers must sort by depth first (the section
    scripts do).
    """
    positions = np.asarray(positions, dtype=float)
    dz = float(np.median(np.diff(positions)))
    if dz <= 0:
        raise SystemExit(
            "positions must INCREASE: median spacing is %.4f. A descending axis "
            "negates k and swaps the downgoing and upgoing quadrants. Sort the "
            "gather by depth before calling separate()." % dz)
    dt = float(lags[1] - lags[0])
    F = np.fft.fftshift(np.fft.fft2(gather))
    k = np.fft.fftshift(np.fft.fftfreq(gather.shape[0], dz))
    f = np.fft.fftshift(np.fft.fftfreq(gather.shape[1], dt))
    K, Fq = np.meshgrid(k, f, indexing="ij")

    # QUADRANTS, not half-planes. For real-valued input the 2-D FFT obeys
    # F(-k,-f) = conj(F(k,f)), so masking to the k > 0 HALF-PLANE and taking the
    # real part returns exactly half the original signal -- all the information is
    # still present and nothing is separated. (An earlier version of this script
    # did that; its "raw", "downgoing" and "upgoing" panels were identical, which
    # is the signature.)
    # Direction lives in the SIGN PRODUCT: energy propagating one way has k and f
    # of the same sign, the other way has them opposite. Both conjugate quadrants
    # must be kept together so the inverse transform is real.
    if taper_cells > 0:
        dk = float(np.abs(k[1] - k[0]))
        df = float(np.abs(f[1] - f[0]))
        wk = np.clip(K / (taper_cells * dk), -1.0, 1.0)
        wf = np.clip(Fq / (taper_cells * df), -1.0, 1.0)
        sk = np.sin(0.5 * np.pi * wk)      # smooth sign of k
        sf = np.sin(0.5 * np.pi * wf)      # smooth sign of f
        ramp = 0.5 * (1.0 + sk * sf)       # -> 1 where signs agree, 0 where they oppose
    else:
        ramp = ((K * Fq) > 0).astype(float)
    # Hermitian-symmetrise before use: without this the k-Nyquist row and
    # f-Nyquist column are assigned inconsistently and np.real() below drops the
    # difference on the floor.
    ramp = _hermitian(ramp)
    mask_pos, mask_neg = ramp, 1.0 - ramp

    def back(m):
        z = np.fft.ifft2(np.fft.ifftshift(F * m))
        return np.real(z), float(np.max(np.abs(z.imag)))

    (a, ia), (b, ib) = back(mask_pos), back(mask_neg)
    scale = max(1e-30, float(np.max(np.abs(gather))))
    # sanity: the two quadrants must sum back to the input, must DIFFER, and the
    # inverse transform must be real to machine precision
    recon = float(np.max(np.abs(a + b - gather))) / scale
    diff = float(np.max(np.abs(a - b))) / scale
    imag = max(ia, ib) / scale
    if not quiet:
        print("  separation check: |a+b-raw|/|raw| = %.2e   |a-b|/|raw| = %.3f"
              "   max|imag|/|raw| = %.2e" % (recon, diff, imag), flush=True)
    if imag > 1e-9:
        raise SystemExit("masked spectrum is not Hermitian (imag/raw = %.2e); "
                         "the inverse transform is discarding energy." % imag)
    if diff < 1e-3:
        raise SystemExit("the two wavefields are identical -- the mask is not "
                         "separating. Check the quadrant logic.")
    analytic = KF_NEG if DOWNGOING_KF_SIGN < 0 else KF_POS
    # `a` is the k*f > 0 quadrant, `b` is k*f < 0
    sep = positions - positions[0]

    def along(x):
        env = np.abs(hilbert(x, axis=1))
        tot = 0.0
        for row, s in zip(env, sep):
            j = int(np.argmin(np.abs(lags - s / v_ref)))
            lo, hi = max(0, j - 3), min(len(lags), j + 4)
            tot += row[lo:hi].max()
        return tot

    empirical = KF_POS if along(a) >= along(b) else KF_NEG
    if keep == "auto":
        if empirical != analytic and not quiet:
            print("  WARNING: the quadrant carrying most energy along t = z/%.0f "
                  "(%s) is NOT the analytic downgoing quadrant (%s). Keeping the "
                  "ANALYTIC one. Either the arrival is not downgoing at this "
                  "velocity, or `positions` is not increasing with depth."
                  % (v_ref, empirical, analytic), flush=True)
        keep = analytic
    keep = _LEGACY.get(keep, keep)
    if keep not in (KF_NEG, KF_POS):
        raise SystemExit("keep must be 'auto', %r or %r, not %r"
                         % (KF_NEG, KF_POS, keep))
    down, up = (a, b) if keep == KF_POS else (b, a)
    return down, up, keep


def moveout_energy(gather, lags, positions, v, sep=None):
    """Envelope energy along t = separation / v. One number per gather.

    `sep` is the separation from the VIRTUAL SOURCE. If it is not supplied it
    falls back to `positions - positions[0]`, which is short by one channel
    spacing because the first receiver is not at the source: the section scripts
    save the true `separation` and main() now passes it.
    """
    env = np.abs(hilbert(gather, axis=1))
    env = env / np.maximum(np.median(env, axis=1, keepdims=True), 1e-30)
    sep = (positions - positions[0]) if sep is None else np.asarray(sep, float)
    vals = []
    for row, s in zip(env, sep):
        j = int(np.argmin(np.abs(lags - s / v)))
        lo, hi = max(0, j - 3), min(len(lags), j + 4)
        vals.append(row[lo:hi].mean())
    return float(np.median(vals))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help=".npz from a section script")
    ap.add_argument("--v-ref", type=float, default=1675.0)
    ap.add_argument("--taper-cells", type=float, default=6.0)
    ap.add_argument("--keep", default="auto", choices=("auto", KF_NEG, KF_POS),
                    help="which quadrant to call downgoing; 'auto' uses the "
                         "analytic convention and warns if the data disagree")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    rng = np.random.default_rng(SEED)

    d = np.load(a.input, allow_pickle=True)
    gather = d["gather"]; lags = d["lags"]
    pos = d["tvd"] if "tvd" in d.files else d["offsets"]
    axis = "true vertical depth (m)" if "tvd" in d.files else "distance along fibre (m)"
    if "hours" in d.files:
        hours = float(d["hours"])
    elif "n_windows" in d.files:
        # older section products predate the 'hours' key. 30 s windows at 15 s
        # step: total span = (N-1)*15 + 30 seconds. Reported rather than left NaN,
        # because every figure must state the duration stacked.
        nw_ = int(d["n_windows"]); hours = ((nw_ - 1) * 15.0 + 30.0) / 3600.0
    else:
        hours = float("nan")
    t0 = str(d["t_first"]) if "t_first" in d.files else ""
    t1 = str(d["t_last"]) if "t_last" in d.files else ""
    fibre = "Nano" if "nano" in Path(a.input).name.lower() else "Deep"
    # true separation from the virtual source when the section script saved it;
    # positions - positions[0] is short by one channel because the first
    # receiver sits one channel below the source
    sep_true = np.asarray(d["separation"], float) if "separation" in d.files else None

    say("Directional f-k separation (downgoing quadrant, k*f < 0)")
    say("  input   : %s" % Path(a.input).name)
    say("  gather  : %d traces x %d lags, %.2f h stacked" % (*gather.shape, hours))
    say("  taper   : %.0f cells in k AND f (0 = brick wall)" % a.taper_cells)
    say("")

    # WHAT THE TAPER COSTS, in band, before any result is quoted.
    dz_ = float(np.median(np.diff(np.asarray(pos, float))))
    dt_ = float(lags[1] - lags[0])
    fr, kept = separation_contrast(gather.shape, dz_, dt_, a.taper_cells, a.v_ref)
    say("  separation contrast on the %.0f m/s line (1.0 = disjoint, 0.5 = none):"
        % a.v_ref)
    say("    " + "  ".join("%.0fHz %.2f" % (x, y) for x, y in
                           zip(fr[::3], kept[::3])))
    if kept.min() < 0.75:
        say("    NOTE: below ~%.0f Hz the mask is still inside its own transition,"
            % float(fr[np.argmax(kept >= 0.75)] if (kept >= 0.75).any() else fr[-1]))
        say("    so the two panels are NOT disjoint there and the down/up ratio")
        say("    below is compressed towards 1. Calibrate it against")
        say("    audit_test_fk.py's known-input table before interpreting it.")
    say("")

    down, up, which = separate(gather, lags, pos, a.taper_cells, keep=a.keep,
                               v_ref=a.v_ref)
    say("  downgoing quadrant kept: %s" % which)
    e_raw = moveout_energy(gather, lags, pos, a.v_ref, sep_true)
    e_dn = moveout_energy(down, lags, pos, a.v_ref, sep_true)
    e_up = moveout_energy(up, lags, pos, a.v_ref, sep_true)
    # ALSO report absolute amplitude, because the moveout metric normalises each
    # trace by its own median envelope and is therefore blind to uniform scaling.
    amp = lambda g_: float(np.median(np.abs(g_).max(axis=1)))
    say("  moveout energy at %.0f m/s : raw %.4f | downgoing %.4f | upgoing %.4f"
        % (a.v_ref, e_raw, e_dn, e_up))
    say("  absolute amplitude        : raw %.4e | downgoing %.4e | upgoing %.4e"
        % (amp(gather), amp(down), amp(up)))
    say("  moveout-score ratio  down/up = %.3f  <- NOT an energy ratio"
        % (e_dn / e_up if e_up else np.nan))
    # THE INTERPRETABLE NUMBER. `moveout_energy` normalises every trace by its
    # own median envelope, so its ratio is scale-free and, calibrated against
    # known input by audit_test_fk.py, is NON-MONOTONIC in the true down:up
    # ratio -- a 1:1 synthetic scores 6.06 there, not 1.0. The plain quadrant
    # energy ratio has none of that ambiguity, so quote this one.
    e_ratio = float(np.sum(down ** 2) / max(1e-300, np.sum(up ** 2)))
    say("  QUADRANT ENERGY RATIO  down/up = %.3f  <- quote this one" % e_ratio)
    say("")

    # ---- THE CONTROL: scramble receivers BEFORE filtering ----
    say("=== input-level null: scramble the receiver order, THEN filter ===")
    say("  a directional filter cannot know the true ordering, so comparable")
    say("  moveout on scrambled input means the operator, not the data.")
    null = np.empty(NULLS)
    for i in range(NULLS):
        perm = rng.permutation(gather.shape[0])
        dn_s, _, _ = separate(gather[perm], lags, pos, a.taper_cells,
                              keep=which, v_ref=a.v_ref, quiet=True)
        null[i] = moveout_energy(dn_s, lags, pos, a.v_ref, sep_true)
    n95 = float(np.percentile(null, 95.0))
    pval = float((np.sum(null >= e_dn) + 1) / (NULLS + 1))
    say("  scrambled-then-filtered: median %.4f, 95th %.4f" % (np.median(null), n95))
    say("  observed downgoing      : %.4f" % e_dn)
    say("  p = %.4f  (%d realisations)" % (pval, NULLS))
    say("")
    if pval < 0.05:
        say("  PASSES the input-level null. The downgoing filter is not")
        say("  manufacturing the moveout: scrambling the input first destroys it.")
    else:
        say("  FAILS the input-level null. The filter produces comparable moveout")
        say("  from scrambled input, so its output on real data is not evidence.")
        say("  This is the same failure mode as the 2.5-4.5 km/s velocity fan.")

    stem = Path(str(STEM) + (("_" + a.tag) if a.tag else ""))
    np.savez_compressed(str(stem) + ".npz", down=down, up=up, lags=lags,
                        positions=pos, kept=which, e_raw=e_raw, e_down=e_dn,
                        e_up=e_up, null=null, p_value=pval, hours=hours,
                        contrast_f=fr, contrast_kept=kept,
                        energy_ratio=e_ratio,
                        taper_cells=float(a.taper_cells))

    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "font.size": 9.5,
        "axes.titlesize": 10, "axes.labelsize": 9.5, "axes.edgecolor": "#b0b0b0",
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "axes.labelcolor": INK, "legend.frameon": False, "legend.fontsize": 8})
    lab, utc_note = (geo.figure_label(t0, t1, hours, fibre,
                            extra="directional f-k separation, downgoing kept "
                                  "(%s quadrant); input-level null p = %.4f"
                                  % (which, pval))
           if t0 else ("%s fibre, %.1f h stacked" % (fibre, hours), ""))
    fig, ax = plt.subplots(1, 4, figsize=(19.5, 6.3), constrained_layout=True,
                           sharey=True)
    fig.suptitle(lab, fontsize=11)
    fig.text(0.995, 0.002, utc_note, ha="right", va="bottom",
             fontsize=6.5, color="#9a9a9a")
    ext = [lags[0], lags[-1], pos[-1], pos[0]]
    sep = (pos - pos[0]) if sep_true is None else sep_true
    for i, (dat, ttl) in enumerate(((gather, "(a) Raw gather"),
                                    (down, "(b) DOWNGOING only"),
                                    (up, "(c) Upgoing only"))):
        tn = dat / np.maximum(np.abs(dat).max(axis=1, keepdims=True), 1e-30)
        lim = float(np.percentile(np.abs(tn), 97.0))
        ax[i].imshow(tn, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                     interpolation="nearest", extent=ext)
        ax[i].plot(sep / a.v_ref, pos, "-", color=C3, lw=1.4)
        ax[i].set(xlim=(-0.35, 0.35), xlabel="correlation lag (s)", title=ttl)
        ax[i].grid(False)
    ax[0].set_ylabel(axis)
    ax[3].hist(null, bins=30, color="#8a8a8a", alpha=.85,
               label="scrambled then filtered")
    ax[3].axvline(e_dn, color=C2, lw=2.2, label="observed (p = %.3f)" % pval)
    ax[3].axvline(n95, color="k", ls="--", lw=1.2, label="null 95th")
    ax[3].set(xlabel="moveout energy at %.0f m/s" % a.v_ref, ylabel="realisations",
              title="(d) Input-level null")
    ax[3].legend(); ax[3].grid(alpha=.3)
    fig.savefig(str(stem) + ".png", bbox_inches="tight")
    Path(str(stem) + ".txt").write_text("\n".join(log) + "\n")
    say("wrote %s.{npz,png,txt}" % stem.name)


if __name__ == "__main__":
    main()
