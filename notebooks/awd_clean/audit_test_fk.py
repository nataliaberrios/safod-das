#!/usr/bin/env python3
"""Synthetic audit of fk_downgoing.separate() -- directional f-k separation.

Pure numpy/scipy, no $OAK data, runs in a few seconds on a login node.

Checks (checklist letters A and B of the 2026-08-20 numerical audit):

  A1  known down + known up, equal amplitude -> recovered panels match the truth
  A2  up/down are not SWAPPED (the correlation matrix must be diagonal-dominant)
  A3  a PURE upgoing gather -> the "downgoing" panel is near zero
  A4  the inverse transform is REAL: imaginary residual at machine precision
  A5  Nyquist / DC rows: is the masked spectrum Hermitian, i.e. is any energy
      doubled or dropped at k = k_Nyq or f = f_Nyq?
  B   round trip: down + up == input to machine precision

Run:  python audit_test_fk.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fk_downgoing as fk

RESULTS = []


def record(name, verdict, detail):
    RESULTS.append((name, verdict, detail))
    print("%-6s %-46s %s" % (verdict, name, detail), flush=True)


def ricker(t, f0=12.0):
    a = (np.pi * f0 * t) ** 2
    return (1.0 - 2.0 * a) * np.exp(-a)


# The geometry the production runs actually used (fk_downgoing_deep_dense.txt):
# 367 traces x 701 lags, 2.0419 m channel spacing, lags +-0.35 s.
PROD = dict(n_ch=367, n_t=701, dz=2.0419, dt=0.7 / 700.0, v=1675.0)


def make_gather(n_ch=367, n_t=701, dz=2.0419, dt=0.001, v=1675.0,
                amp_down=1.0, amp_up=1.0, f0=12.0):
    """A causal downgoing arrival at t = +z/v and its acausal mirror at -z/v.

    This is the shape of a virtual-source gather, so `along()`'s trajectory
    t = z/v_ref coincides exactly with the downgoing wave -- the fair test of
    the automatic quadrant pick.
    """
    z = np.arange(n_ch) * dz
    lags = (np.arange(n_t) - n_t // 2) * dt
    T = lags[None, :]
    Z = z[:, None]
    down = amp_down * ricker(T - Z / v, f0)
    up = amp_up * ricker(T + Z / v, f0)
    return down, up, z, lags


def norm_corr(a, b):
    a = a.ravel() - a.mean()
    b = b.ravel() - b.mean()
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d else np.nan


# --------------------------------------------------------------- A1/A2/B
def test_separation():
    down_t, up_t, z, lags = make_gather()
    gather = down_t + up_t
    d_out, u_out, which = fk.separate(gather, lags, z, taper_cells=6.0,
                                      keep="auto", v_ref=1675.0)

    # B: round trip
    resid = float(np.max(np.abs(d_out + u_out - gather)))
    scale = float(np.max(np.abs(gather)))
    record("B  round-trip down+up == input",
           "PASS" if resid / scale < 1e-12 else "FAIL",
           "max|d+u-raw|/|raw| = %.3e" % (resid / scale))

    cdd = norm_corr(d_out, down_t)
    cdu = norm_corr(d_out, up_t)
    cud = norm_corr(u_out, down_t)
    cuu = norm_corr(u_out, up_t)
    print("    corr matrix   down_out.down_true=%+.4f  down_out.up_true=%+.4f" % (cdd, cdu))
    print("                  up_out  .down_true=%+.4f  up_out  .up_true=%+.4f" % (cud, cuu))
    record("A1 recovered downgoing matches truth",
           "PASS" if cdd > 0.9 else "FAIL",
           "corr(down_out, down_true) = %+.4f (kept=%s)" % (cdd, which))
    record("A2 panels not swapped",
           "PASS" if (cdd > cdu and cuu > cud) else "FAIL",
           "diag %+.4f/%+.4f vs off-diag %+.4f/%+.4f" % (cdd, cuu, cdu, cud))
    return which


# --------------------------------------------------------------- A3
def _bandpass(x, lo, hi, fs):
    from scipy.signal import butter, sosfiltfilt
    return sosfiltfilt(butter(4, [lo, hi], btype="bandpass", fs=fs,
                              output="sos"), x, axis=1)


def test_pure_upgoing(which_from_mixed):
    """A directional f-k filter cannot separate energy at k ~ 0 or f ~ 0, so
    leakage must be quoted with the band it was measured in. 13-20 Hz is where
    `separation_contrast` says the mask is disjoint at this geometry; 5-20 Hz is
    the band the production runs actually used.

    Leakage is quoted as an ENERGY fraction. The peak-amplitude fraction is 20x
    larger (0.42 vs 0.021 at 13-20 Hz) because the leaked energy is concentrated
    in edge ringing: `separate()` applies NO window in z or t before fft2, so the
    rectangular aperture and the arrival running off the +-0.35 s lag edge both
    leak. Both numbers are reported; the energy one is the physical statement.
    """
    fs = 1.0 / PROD["dt"]
    for lo, hi, strict in ((13.0, 20.0, True), (5.0, 20.0, False)):
        _, up_t, z, lags = make_gather(amp_down=0.0, amp_up=1.0, **PROD)
        up_b = _bandpass(up_t, lo, hi, fs)
        d_out, _, _ = fk.separate(up_b, lags, z, taper_cells=6.0,
                                  keep=which_from_mixed, v_ref=PROD["v"],
                                  quiet=True)
        e1 = float(np.sum(d_out ** 2) / np.sum(up_b ** 2))
        m1 = float(np.max(np.abs(d_out)) / np.max(np.abs(up_b)))
        down_t, _, z2, lags2 = make_gather(amp_down=1.0, amp_up=0.0, **PROD)
        dn_b = _bandpass(down_t, lo, hi, fs)
        _, u2, _ = fk.separate(dn_b, lags2, z2, taper_cells=6.0,
                               keep=which_from_mixed, v_ref=PROD["v"], quiet=True)
        e2 = float(np.sum(u2 ** 2) / np.sum(dn_b ** 2))
        m2 = float(np.max(np.abs(u2)) / np.max(np.abs(dn_b)))
        worst = max(e1, e2)
        record("A3 %s Hz: wrong-panel ENERGY leakage" % ("%.0f-%.0f" % (lo, hi)),
               ("PASS" if worst < 0.10 else "FAIL") if strict
               else ("PASS" if worst < 0.20 else "FAIL"),
               "pure up -> down %.4f, pure down -> up %.4f (peak-amplitude "
               "fraction %.2f / %.2f -- edge ringing)%s"
               % (e1, e2, m1, m2, "" if strict else "  [production band]"))


# --------------------------------------------------------------- A4/A5
def _mask(gather, lags, positions, taper_cells):
    """Re-derive fk_downgoing's mask so the imaginary part can be inspected."""
    dz = float(np.median(np.diff(positions)))
    dt = float(lags[1] - lags[0])
    F = np.fft.fftshift(np.fft.fft2(gather))
    k = np.fft.fftshift(np.fft.fftfreq(gather.shape[0], dz))
    f = np.fft.fftshift(np.fft.fftfreq(gather.shape[1], dt))
    K, Fq = np.meshgrid(k, f, indexing="ij")
    if taper_cells > 0:
        dk = float(np.abs(k[1] - k[0])); df = float(np.abs(f[1] - f[0]))
        sk = np.sin(0.5 * np.pi * np.clip(K / (taper_cells * dk), -1.0, 1.0))
        sf = np.sin(0.5 * np.pi * np.clip(Fq / (taper_cells * df), -1.0, 1.0))
        ramp = 0.5 * (1.0 + sk * sf)
    else:
        ramp = ((K * Fq) > 0).astype(float)
    return F, ramp, k, f


def test_real_and_nyquist():
    down_t, up_t, z, lags = make_gather()
    gather = down_t + up_t
    for taper in (6.0, 0.0):
        F, ramp, k, f = _mask(gather, lags, z, taper)
        cplx = np.fft.ifft2(np.fft.ifftshift(F * ramp))
        imag = float(np.max(np.abs(cplx.imag))) / float(np.max(np.abs(cplx.real)))
        record("A4 masked ifft2 is real (taper=%g)" % taper,
               "PASS" if imag < 1e-10 else "FAIL",
               "max|imag|/max|real| = %.3e" % imag)

        # A5: Hermitian test of the MASK itself. For real output we need
        # m(-k,-f) == m(k,f) at every bin INCLUDING the Nyquist row/col, where
        # -k_nyq aliases back to +k_nyq and the partner index is not the mirror.
        nk, nf = ramp.shape
        m_un = np.fft.ifftshift(ramp)              # back to fft ordering
        herm = m_un[(-np.arange(nk)) % nk][:, (-np.arange(nf)) % nf]
        worst = float(np.max(np.abs(m_un - herm)))
        # where does it break?
        bad = np.argwhere(np.abs(m_un - herm) > 1e-12)
        rows = sorted(set(bad[:, 0].tolist())) if bad.size else []
        cols = sorted(set(bad[:, 1].tolist())) if bad.size else []
        record("A5 mask Hermitian-symmetric (taper=%g)" % taper,
               "PASS" if worst < 1e-12 else "FAIL",
               "max|m(k,f)-m(-k,-f)| = %.3f; %d bad bins on fft rows %s cols %s"
               % (worst, bad.shape[0], rows[:4], cols[:4]))

        # how much ENERGY sits in the mis-assigned bins?
        e_tot = float(np.sum(np.abs(F) ** 2))
        if bad.size:
            F_un = np.fft.ifftshift(F)
            e_bad = float(np.sum(np.abs(F_un[bad[:, 0], bad[:, 1]]) ** 2))
        else:
            e_bad = 0.0
        record("A5b energy in non-Hermitian bins (taper=%g)" % taper,
               "PASS" if e_bad / e_tot < 1e-6 else "WARN",
               "%.3e of total spectral energy" % (e_bad / e_tot))


# --------------------------------------------------------------- sign convention
def test_sign_convention():
    """Which quadrant does a DOWNGOING (+z) wave actually occupy?

    numpy's fft2 uses exp(-2i*pi*(k z + f t)).  For u(z,t)=w(t - z/v) the
    stationary condition gives k = -f/v, i.e. k and f have OPPOSITE signs.
    Verify numerically, and report which sign product fk_downgoing's `keep`
    label corresponds to.
    """
    down_t, _, z, lags = make_gather(amp_down=1.0, amp_up=0.0)
    dz = float(np.median(np.diff(z))); dt = float(lags[1] - lags[0])
    F = np.fft.fftshift(np.fft.fft2(down_t))
    k = np.fft.fftshift(np.fft.fftfreq(down_t.shape[0], dz))
    f = np.fft.fftshift(np.fft.fftfreq(down_t.shape[1], dt))
    K, Fq = np.meshgrid(k, f, indexing="ij")
    p = float(np.sum(np.abs(F[(K * Fq) > 0]) ** 2))
    n = float(np.sum(np.abs(F[(K * Fq) < 0]) ** 2))
    record("sign: downgoing occupies sign(k)*sign(f) < 0",
           "PASS" if n > p else "FAIL",
           "energy k*f>0 = %.3e, k*f<0 = %.3e  (ratio %.1f)"
           % (p, n, n / max(p, 1e-300)))
    print("    -> k*f>0 is physically UPGOING; fk_downgoing.DOWNGOING_KF_SIGN "
          "= %d encodes this." % fk.DOWNGOING_KF_SIGN)
    record("sign: DOWNGOING_KF_SIGN agrees with the measurement",
           "PASS" if (fk.DOWNGOING_KF_SIGN < 0) == (n > p) else "FAIL",
           "constant says %s, measurement says %s"
           % (fk.KF_NEG if fk.DOWNGOING_KF_SIGN < 0 else fk.KF_POS,
              fk.KF_NEG if n > p else fk.KF_POS))


# --------------------------------------------------------------- geometry sanity
def test_descending_positions():
    """separate() takes dz = median(diff(positions)). On the return limb TVD
    DECREASES with channel index; the section scripts argsort so it increases,
    but a caller passing raw descending TVD would get a NEGATIVE dz, which
    negates k and swaps the quadrants. That must be refused, not absorbed."""
    down_t, up_t, z, lags = make_gather()
    g = down_t + up_t
    try:
        fk.separate(g, lags, z[::-1].copy(), taper_cells=6.0, keep="auto")
        record("descending positions are refused", "FAIL",
               "separate() accepted a negative dz and silently swapped quadrants")
    except SystemExit as e:
        record("descending positions are refused", "PASS", str(e).split(".")[0])


def test_ratio_calibration():
    """CALIBRATION. What down/up ratio does the pipeline report for a gather
    whose TRUE ratio is known? The production runs quote 1.56-1.97 and the Nano
    run 0.985; those numbers mean nothing until this table exists.

    Run at the exact production geometry in the 5-20 Hz analysis band, with an
    incoherent noise background so `moveout_energy`'s per-trace median-envelope
    normalisation behaves as it does on real data (with a noiseless synthetic
    the median envelope is set by the arrival itself and the metric saturates).
    """
    fs = 1.0 / PROD["dt"]
    rng = np.random.default_rng(4)
    print("    true down:up   ->  reported e_down/e_up   (taper 6 | brick wall)")
    rows = []
    for ad, au, lab in ((1.0, 0.0, "pure down"), (1.0, 0.25, "4 : 1"),
                        (1.0, 1.0, "1 : 1"), (0.25, 1.0, "1 : 4"),
                        (0.0, 1.0, "pure up")):
        d_t, u_t, z, lags = make_gather(amp_down=ad, amp_up=au, **PROD)
        noise = rng.standard_normal(d_t.shape) * 0.5
        g = _bandpass(d_t + u_t + noise, 5.0, 20.0, fs)
        out = []
        for taper in (6.0, 0.0):
            dn, up, _ = fk.separate(g, lags, z, taper_cells=taper,
                                    keep=fk.KF_NEG, v_ref=PROD["v"], quiet=True)
            e_d = fk.moveout_energy(dn, lags, z, PROD["v"])
            e_u = fk.moveout_energy(up, lags, z, PROD["v"])
            out.append(e_d / e_u if e_u else np.inf)
        rows.append((lab, out[0], out[1]))
        print("      %-10s   ->  %8.3f              |  %8.3f"
              % (lab, out[0], out[1]), flush=True)
    ratios = [r[1] for r in rows]
    monotone = all(ratios[i] >= ratios[i + 1] for i in range(len(ratios) - 1))
    record("calibration: reported ratio decreases with true ratio",
           "PASS" if monotone else "FAIL",
           "pure-down %.2f, 4:1 %.2f, 1:1 %.2f, 1:4 %.2f, pure-up %.2f"
           % tuple(ratios))
    record("calibration: 1:1 input reports ~1",
           "PASS" if 0.7 < ratios[2] < 1.4 else "FAIL",
           "balanced input reports down/up = %.3f" % ratios[2])
    record("calibration: production 1.56-1.97 is inside the usable range",
           "PASS" if ratios[0] > 1.97 else "FAIL",
           "pure downgoing tops out at %.2f, so 1.97 corresponds to roughly "
           "%s of the dynamic range" % (ratios[0],
           "%.0f%%" % (100.0 * (1.97 - ratios[-1]) /
                       max(1e-9, ratios[0] - ratios[-1]))))
    print("      -> production reported 1.56-1.97 (Deep) and 0.985 (Nano).")


def test_contrast_table():
    """The taper is a fixed number of CELLS and a cell is 1/aperture, so on a
    short array the transition can swallow the signal band."""
    fr, kept = fk.separation_contrast((PROD["n_ch"], PROD["n_t"]), PROD["dz"],
                                      PROD["dt"], 6.0, PROD["v"])
    print("    kept fraction on the %.0f m/s line, production geometry:" % PROD["v"])
    print("      " + "  ".join("%.0fHz %.2f" % (x, y)
                               for x, y in zip(fr[::2], kept[::2])))
    record("taper leaves the mask disjoint across the 5-20 Hz band",
           "PASS" if kept.min() > 0.9 else "FAIL",
           "worst in-band kept fraction %.3f at %.1f Hz (1.0 = disjoint, "
           "0.5 = no separation at all)" % (kept.min(), fr[int(np.argmin(kept))]))


def main():
    print("=" * 78)
    print("audit_test_fk.py -- synthetic audit of fk_downgoing.separate()")
    print("=" * 78)
    which = test_separation()
    print()
    test_pure_upgoing(which)
    print()
    test_real_and_nyquist()
    print()
    test_sign_convention()
    print()
    test_descending_positions()
    print()
    test_contrast_table()
    print()
    test_ratio_calibration()
    print()
    n_fail = sum(1 for _, v, _ in RESULTS if v == "FAIL")
    print("%d checks, %d FAIL, %d WARN"
          % (len(RESULTS), n_fail, sum(1 for _, v, _ in RESULTS if v == "WARN")))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
