"""
Shared dv/v estimators. Imported by dvv_hrsn.py (G3) and dvv_measure.py (G4/G5).

Two estimators, deliberately different in what they assume:

  STRETCHING (coda) -- high sensitivity, path-integrated, NOT depth-local. Coda
  recorded at 200 m and 800 m share most of their scattering volume, so a
  per-channel coda dv/v is not a measurement of the medium at that channel. Use it
  for the depth-integrated value and for contrast between bins, never as a profile.

  DIFFERENTIAL TIMING (direct arrival) -- lower sensitivity but genuinely local, and
  immune to origin-time error because differencing along the array cancels the
  unknown source instant. This is interval velocity, the same quantity a check shot
  measures.

SIGN CONVENTION, stated once because it is easy to get backwards. We resample b
onto t*(1+eps) and correlate against a. If eps > 0 maximises the correlation, then
compressing b's time axis matches a, so b's arrivals sit at LATER times than a's,
so b has larger travel times, so b is SLOWER. Therefore

    dv/v (from a to b) = -eps

Errors are bootstrapped over the independent units (stations for HRSN, channels for
DAS) rather than taken from the correlation-peak curvature, which assumes white
noise and is optimistic in practice.
"""
import numpy as np
from scipy.signal import butter, sosfiltfilt


def bandpass(x, fs, band, order=4):
    sos = butter(order, list(band), btype='band', fs=fs, output='sos')
    return sosfiltfilt(sos, np.asarray(x, float), axis=-1)


def _norm(x):
    x = x - x.mean()
    n = np.sqrt(np.sum(x * x))
    return x / n if (np.isfinite(n) and n > 0) else None


def bulk_align(a, b, fs, max_lag_s=2.0):
    """Shift b to align with a. Returns (b_shifted, lag_s).

    NECESSARY, NOT OPTIONAL. Catalog origin times carry 0.1-0.5 s of error, and at
    5-20 Hz a 0.2 s misalignment drives the correlation of two identical waveforms
    to zero. Measuring coda CC or stretch without removing this first returns ~0 for
    genuine repeaters -- which is exactly what happened in the first pass of
    dvv_hrsn.py and coda_window_survey.py, on pairs that correlate at 0.99 when a
    lag search is used.

    The bulk shift is bookkeeping (origin time), not signal. The medium change lives
    in the STRETCH measured after alignment.
    """
    n = min(a.size, b.size)
    a, b = a[:n], b[:n]
    aa, bb = a - a.mean(), b - b.mean()
    if not np.any(aa) or not np.any(bb):
        return b, 0.0
    cc = np.correlate(bb, aa, 'full')
    mid = n - 1
    pad = int(max_lag_s * fs)
    lo = max(mid - pad, 0)
    seg = cc[lo:mid + pad + 1]
    if seg.size == 0:
        return b, 0.0
    k = int(np.argmax(np.abs(seg))) + lo - mid
    return np.roll(b, -k), k / fs


def stretch_cc(a, b, fs, eps, lag_search_s=0.0):
    """Correlation between a and b resampled onto t*(1+eps). Scalar per eps.

    lag_search_s > 0 takes the maximum over residual lags at each trial eps, which
    makes the estimate robust to imperfect bulk alignment. Keep it small (a few
    tenths of a second): a wide lag search can trade shift against stretch and
    flatten the eps peak.
    """
    n = min(a.size, b.size)
    t = np.arange(n) / fs
    an = _norm(a[:n])
    if an is None:
        return np.full(np.size(eps), np.nan)
    pad = int(lag_search_s * fs)
    out = np.full(np.size(eps), np.nan)
    for i, e in enumerate(np.atleast_1d(eps)):
        bi = np.interp(t * (1.0 + e), t, b[:n], left=0.0, right=0.0)
        bn = _norm(bi)
        if bn is None:
            continue
        if pad <= 0:
            out[i] = float(np.dot(an, bn))
        else:
            cc = np.correlate(bn, an, 'full')
            mid = n - 1
            seg = cc[max(mid - pad, 0):mid + pad + 1]
            out[i] = float(seg[int(np.argmax(np.abs(seg)))]) if seg.size else np.nan
    return out


def stretch_dvv(a, b, fs, eps_max=0.04, n_eps=161, refine=True,
                lag_search_s=0.2):
    """dv/v from coda stretching. Returns (dvv, cc_at_peak, eps_grid, cc_curve).

    Grid search then parabolic refinement on the three points around the peak. A
    peak at the edge of the grid is returned as NaN rather than clipped, so an
    unconverged fit cannot masquerade as a large measurement.
    """
    eps = np.linspace(-eps_max, eps_max, n_eps)
    cc = stretch_cc(a, b, fs, eps, lag_search_s=lag_search_s)
    if np.all(np.isnan(cc)):
        return np.nan, np.nan, eps, cc
    k = int(np.nanargmax(cc))
    if k == 0 or k == len(eps) - 1:
        return np.nan, float(cc[k]), eps, cc
    if refine:
        y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]
        den = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        e_hat = eps[k] + shift * (eps[1] - eps[0])
    else:
        e_hat = eps[k]
    return -float(e_hat), float(cc[k]), eps, cc


def stretch_dvv_bootstrap(traces_a, traces_b, fs, n_boot=200, seed=0, **kw):
    """dv/v with a bootstrap error over independent units (stations or channels).

    Each unit is measured separately; the reported value is the median and the
    error is the bootstrap std of the median. Units that fail to converge are
    dropped and counted.
    """
    vals, ccs = [], []
    for a, b in zip(traces_a, traces_b):
        d, c, _, _ = stretch_dvv(a, b, fs, **kw)
        if np.isfinite(d):
            vals.append(d); ccs.append(c)
    vals = np.array(vals)
    if vals.size == 0:
        return np.nan, np.nan, 0, np.nan
    rng = np.random.default_rng(seed)
    boot = [np.median(rng.choice(vals, vals.size, replace=True))
            for _ in range(n_boot)]
    return (float(np.median(vals)), float(np.std(boot)), int(vals.size),
            float(np.median(ccs)))


def sub_sample_delay(a, b, fs, f_band, max_lag_s=0.5):
    """Delay of b relative to a by cross-spectral phase slope.

    Poupinet, Ellsworth & Frechet 1984 (doi:10.1029/JB089iB07p05719). The phase of
    the cross spectrum is linear in frequency with slope 2*pi*delay; fitting the
    slope over a band where coherence is high beats picking a correlation peak,
    which is quantised at the sample interval (10 ms at 100 Hz is 72 degrees of
    phase at 20 Hz -- far too coarse for a 1 ms target).

    Returns (delay_s, weighted_coherence). Positive delay means b arrives later.
    """
    n = min(a.size, b.size)
    a, b = a[:n], b[:n]
    if not np.any(a) or not np.any(b):
        return np.nan, 0.0
    # integer alignment first, so the phase never wraps within the band
    cc = np.correlate(b - b.mean(), a - a.mean(), 'full')
    mid = n - 1
    pad = int(max_lag_s * fs)
    seg = cc[max(mid - pad, 0):mid + pad + 1]
    if seg.size == 0:
        return np.nan, 0.0
    k0 = int(np.argmax(np.abs(seg))) + max(mid - pad, 0) - mid
    b_sh = np.roll(b, -k0)
    A = np.fft.rfft(a * np.hanning(n))
    B = np.fft.rfft(b_sh * np.hanning(n))
    f = np.fft.rfftfreq(n, 1.0 / fs)
    m = (f >= f_band[0]) & (f <= f_band[1])
    if m.sum() < 4:
        return np.nan, 0.0
    X = B[m] * np.conj(A[m])
    ph = np.unwrap(np.angle(X))
    w = np.abs(X)
    w = w / w.sum() if w.sum() > 0 else np.ones(m.sum()) / m.sum()
    # weighted least squares slope through the origin-free linear phase
    fm = f[m]
    fbar = np.sum(w * fm)
    pbar = np.sum(w * ph)
    slope = (np.sum(w * (fm - fbar) * (ph - pbar)) /
             max(np.sum(w * (fm - fbar) ** 2), 1e-30))
    resid = ph - (pbar + slope * (fm - fbar))
    coh = float(np.exp(-np.var(resid) / 2))          # phase-stability proxy
    # SIGN. numpy's FFT uses exp(-2i.pi.f.t), so for b(t) = a(t - tau) the cross
    # spectrum B.conj(A) has phase -2.pi.f.tau and therefore tau = -slope/(2.pi).
    # This returned +slope/(2.pi), which SUBTRACTS the sub-sample residual instead
    # of adding it, so the error was exactly -2x the fractional part of the true
    # delay -- up to one whole sample, 10 ms at 100 Hz, against a stated 1 ms
    # target. Verified: a pair delayed by exactly +3.00 ms returned -2.73 ms, and
    # the -5.73 ms error was constant across true delays of +3, +13 and -7 ms
    # because all three share that fractional part.
    return float(k0 / fs - slope / (2 * np.pi)), coh


def interval_dvv(delays, depths, tau_of_z, z_lo, z_hi, half=15.0):
    """Interval dv/v over [z_lo, z_hi] from per-depth delays.

    dv/v = -( dt(z_lo) - dt(z_hi) ) / tau,  tau = interval travel time.

    The difference cancels any depth-independent term, which is why this is immune
    to origin-time error: an unknown constant added to every delay drops out.

    `half` is the half-width (m) of the depth band averaged at each edge. Shrink it
    below the 15 m default when an edge sits near the wellhead, so the band cannot
    reach into the uncemented lead-in.
    """
    m_lo = (depths >= z_lo - half) & (depths <= z_lo + half)
    m_hi = (depths >= z_hi - half) & (depths <= z_hi + half)
    if m_lo.sum() < 2 or m_hi.sum() < 2:
        return np.nan, np.nan
    d_lo = np.nanmedian(delays[m_lo])
    d_hi = np.nanmedian(delays[m_hi])
    tau = abs(tau_of_z(z_hi) - tau_of_z(z_lo))
    if not np.isfinite(tau) or tau <= 0:
        return np.nan, np.nan
    sd = np.sqrt(np.nanvar(delays[m_lo]) / m_lo.sum() +
                 np.nanvar(delays[m_hi]) / m_hi.sum())
    # SIGN. From dt(z) = dt0 - (dv/v).tau(z) with z_hi deeper than z_lo, so
    # tau_hi > tau_lo:  dt_lo - dt_hi = +(dv/v).|tau_hi - tau_lo|, hence
    # dv/v = +(d_lo - d_hi)/tau. This returned the negative, i.e. the exact
    # opposite of slope_dvv on identical input. Verified against a synthetic with
    # dv/v = -0.500%: slope_dvv gave -0.5000% (correct), interval_dvv gave
    # +0.4985%. slope_dvv is the reference; this function was wrong.
    return float((d_lo - d_hi) / tau), float(sd / tau)


def slope_dvv(delays, taus, weights=None, n_boot=400, seed=0,
              robust_iters=2, clip=3.5, block=1):
    """Interval dv/v as the slope of delay against travel time. Preferred estimator.

    interval_dvv uses only two edge bands and throws away everything between them.
    But dt accumulates along the ray, so within an interval of uniform fractional
    velocity change dv/v the delay is exactly linear in the one-way travel time:

        dt(z) = dt_0 - (dv/v) * tau(z)      =>   dv/v = -d(dt)/d(tau)

    Fitting the slope over every channel in the interval uses N points instead of
    two bands, and the intercept absorbs the origin-time difference for free -- so
    this inherits interval_dvv's immunity to the unknown source instant.

    The residual scatter about the line is itself diagnostic: a genuine uniform
    change is linear, whereas a cycle-skipped channel or a coupling defect shows up
    as an outlier rather than as a change in slope.

    Errors are bootstrapped over channels rather than taken from the fit
    covariance, which would assume the per-channel delays are independent and
    identically distributed. They are neither -- neighbouring channels share
    wavefield.

    `block` -- BOOTSTRAP BLOCK LENGTH IN CHANNELS. Read this before trusting an
    error bar. A plain bootstrap that resamples channels one at a time makes the
    very independence assumption the paragraph above rejects, so it understates the
    error by roughly sqrt(block). On this fiber the 16.335 m gauge length spans
    16.335 / 1.021 = 16 channels, which all average the same strain, so the correct
    call for DAS delays is `block=16`.

    The default stays 1 only so that results already published from this function
    (g5_shallow_dvv) remain reproducible. It is NOT the right value for DAS. An
    earlier measurement, ray_parameter_test.py, divided by sqrt(N) with N ~ 700 and
    produced error bars about 4x too small; the corroborating symptom was a median
    |z| of 1.75 where 0.674 was expected. Pass block=16 for any new work.

    Returns (dvv, err, n_used, rms_residual_seconds).
    """
    d = np.asarray(delays, float)
    t = np.asarray(taus, float)
    ok = np.isfinite(d) & np.isfinite(t)
    if weights is not None:
        w_all = np.asarray(weights, float)
        ok &= np.isfinite(w_all)
    d, t = d[ok], t[ok]
    w = (np.asarray(weights, float)[ok] if weights is not None
         else np.ones(d.size))
    n = d.size
    if n < 6 or np.ptp(t) <= 0:
        return np.nan, np.nan, int(n), np.nan

    def fit(dd, tt, ww):
        sw = ww.sum()
        if sw <= 0:
            return np.nan
        tb, db = np.sum(ww * tt) / sw, np.sum(ww * dd) / sw
        den = np.sum(ww * (tt - tb) ** 2)
        return np.sum(ww * (tt - tb) * (dd - db)) / den if den > 0 else np.nan

    # Iteratively reject channels whose residual is more than `clip` robust sigma
    # off the line. A cycle-skipped channel is displaced by a whole period -- tens
    # of ms against a target of well under one -- so a single skip inside an
    # interval can dominate a least-squares slope entirely. Rejecting on the
    # residual rather than on |dt| is what makes this work: the delays themselves
    # legitimately trend with depth, so a fixed cut on |dt| either keeps the skips
    # or throws away real signal at the ends of the interval.
    s = fit(d, t, w)
    if not np.isfinite(s):
        return np.nan, np.nan, int(n), np.nan
    for _ in range(max(int(robust_iters), 0)):
        sw = w.sum()
        tb, db = np.sum(w * t) / sw, np.sum(w * d) / sw
        r = d - (db + s * (t - tb))
        sig = 1.4826 * np.median(np.abs(r - np.median(r)))
        if not np.isfinite(sig) or sig <= 0:
            break
        keep = np.abs(r - np.median(r)) <= clip * sig
        if keep.all() or keep.sum() < 6:
            break
        d, t, w = d[keep], t[keep], w[keep]
        s2 = fit(d, t, w)
        if not np.isfinite(s2):
            break
        s = s2
    n = d.size

    sw = w.sum()
    tb, db = np.sum(w * t) / sw, np.sum(w * d) / sw
    resid = d - (db + s * (t - tb))
    rms = float(np.sqrt(np.sum(w * resid ** 2) / sw))

    rng = np.random.default_rng(seed)
    boot = []
    blk = max(int(block), 1)
    if blk > 1:
        # Moving-block bootstrap. Delays are ordered along the fiber (taus increase
        # monotonically with channel), so contiguous index blocks are contiguous in
        # depth and one block spans one gauge length of shared strain.
        nblk = max(int(np.ceil(n / blk)), 1)
        starts_max = max(n - blk, 0)
    for _ in range(n_boot):
        if blk > 1:
            st = rng.integers(0, starts_max + 1, nblk)
            k = np.concatenate([np.arange(x, min(x + blk, n)) for x in st])[:n]
        else:
            k = rng.integers(0, n, n)
        if k.size < 6:
            continue
        b = fit(d[k], t[k], w[k])
        if np.isfinite(b):
            boot.append(b)
    err = float(np.std(boot)) if len(boot) > 10 else np.nan
    return -float(s), err, int(n), rms


def seasonal_phase_days(t1, t2):
    """Separation in seasonal phase (days), 0-182.6, independent of elapsed time.

    A seasonal signal scales with this; secular or instrumental drift scales with
    elapsed time. Across the confirmed pair set the two are nearly uncorrelated,
    which makes them separable.
    """
    d1 = t1.dayofyear if hasattr(t1, 'dayofyear') else t1.timetuple().tm_yday
    d2 = t2.dayofyear if hasattr(t2, 'dayofyear') else t2.timetuple().tm_yday
    d = abs(d1 - d2) % 365.25
    return float(min(d, 365.25 - d))
