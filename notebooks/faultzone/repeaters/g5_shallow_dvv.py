"""
G5: how shallow is the velocity change? Differential direct-arrival timing on DAS.

WHY THIS IS THE DELIVERABLE, and why it only became well-posed after G3.

G3 measured coda dv/v on the HRSN borehole seismometers -- the better instrument,
sitting at 251-284 m -- across the seven confirmed repeater pairs and found NULL:
0/10 pairs above a control-derived floor of 1.94%. Li & Ben-Zion 2023 put the
seasonal signal in the top few tens of metres with peak sensitivity near 17 m. A
null at 250 m is therefore not a failure; it localises the change ABOVE the depth
where the best instrument sits, and the only sensor in the ground that samples
that interval is the shallow end of this fibre.

G0 established that the shallow end is usable: the wellhead is channel ~23, the
uncemented lead-in is channels 0-23 (35x baseline noise), and from channel 23 down
the noise is flat within +/-10%. Channel 40 is the ~17 m peak-sensitivity depth.
The CH_LO = 100 cut used everywhere else in this directory is depth 79 m -- it was
discarding the entire target interval. This script deliberately does not use it.

WHY DIFFERENTIAL TIMING AND NOT CODA STRETCHING. Coda recorded at 30 m and at
800 m shares nearly all of its scattering volume, so a per-channel coda dv/v is not
a measurement of the medium at that channel; it is the same path-integrated number
copied 900 times. It cannot produce a depth profile no matter how it is binned.
The direct arrival can: the wave travels UP the borehole, so the arrival at depth z
has sampled only the medium BELOW z. Differencing delays between two depths
isolates the interval between them. That is interval velocity -- the same quantity
a check shot measures -- and it is what dvv_core.slope_dvv estimates.

    dt(z) = dt_0 - (dv/v) * tau(z)        =>       dv/v = -d(dt)/d(tau)

The intercept dt_0 absorbs the origin-time difference between the two occurrences,
so no absolute timing is needed. This is the one property that makes the
measurement possible at all with catalog origin times that carry 0.1-0.5 s of
error.

THE CONFOUNDER THAT MATTERS, AND THE CONTROL THAT SETTLES IT. A few metres of
hypocentre separation between the two occurrences also produces a differential
delay across the array. It is not distinguishable from a velocity change at any
single depth -- but it has a different SHAPE. A source-position difference
perturbs the ray's apparent slowness over the WHOLE path, so dt is linear in tau
from the bottom of the array to the top. A change confined above 100 m accumulates
delay only in the shallow interval and leaves dt FLAT below it. So the deep
intervals are the control, and G3 independently predicts they must read zero.

  shallow slope, flat below   -> shallow-confined velocity change     (the result)
  uniform slope at all depths -> hypocentre separation                (artifact)
  no slope anywhere           -> no resolvable change                 (upper bound)

Three further controls, all run here rather than deferred:
  (a) random non-repeating pairs -- the method's noise floor on the same channels
  (b) acausal: time-reverse B before timing. Cannot produce a coherent delay
      profile, so it measures the false-positive rate of the profile shape itself
  (c) SNR equalisation: shallow and deep channels do not have equal SNR, and the
      slope estimator's scatter depends on SNR. Additive noise is used to degrade
      whichever band is cleaner until the two match, and the contrast must survive

Inputs are all cached: cache_all/*.npz (206 events, 900 ch, 25 s at 100 Hz),
g0_refine.npz (registration + the 2005 check-shot travel-time curve),
hrsn_control.csv (the confirmed pairs and the random controls).
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dvv_core import bandpass, slope_dvv, interval_dvv, sub_sample_delay  # noqa: E402
from dvv_core import seasonal_phase_days                                  # noqa: E402

CACHE = os.path.join(HERE, 'cache_all')
PRE_S = 5.0                       # cache starts 5 s before catalog origin

BANDS = [(5.0, 30.0), (5.0, 20.0)]   # primary, then the band used elsewhere
PICK_BAND = (5.0, 30.0)

WIN = (-0.20, 0.80)               # timing window about the per-channel P arrival
NOISE_WIN = (-4.5, -1.0)          # pre-event, relative to catalog origin
MAX_DELAY_S = 0.020               # cycle-skip guard: half period at 25 Hz
MIN_COH = 0.60                    # phase-stability floor per channel
SMOOTH_CH = 21                    # median-filter width for the moveout curve

# Depth intervals. The first is the target; the last two are the artifact control.
INTERVALS = [(2.0, 75.0), (75.0, 250.0), (250.0, 550.0), (550.0, 850.0)]
# Finer profile, reported for shape. Individually less significant, by design.
PROFILE = [(2.0, 25.0), (25.0, 50.0), (50.0, 100.0), (100.0, 200.0),
           (200.0, 400.0), (400.0, 600.0), (600.0, 850.0)]
SHALLOW, DEEP = (2.0, 75.0), (250.0, 850.0)

FIBRE_IN_WELL = 864.0


# --------------------------------------------------------------------------
# registration and the travel-time curve
# --------------------------------------------------------------------------
def registration():
    """Channel->depth and depth->one-way vertical travel time.

    tau below the shallowest check-shot level (20.3 m) is extrapolated at that
    level's own average velocity, which is the constant-velocity continuation
    consistent with the measurement rather than a new assumption. That top layer
    is slow -- 20.3 m in 26.9 ms, i.e. 753 m/s -- and it is slowness, not
    thickness, that gives the shallow interval its leverage: the top 20 m alone
    carries more one-way travel time than the 250 m below it.
    """
    d = np.load(os.path.join(HERE, 'g0_refine.npz'), allow_pickle=True)
    wh, dx = float(d['wellhead']), float(d['dx'])
    top, bot = int(d['top']), int(d['bot'])
    zcs, tcs = np.asarray(d['zcs'], float), np.asarray(d['tcs'], float)
    v0 = zcs[0] / tcs[0]

    def tau(z):
        z = np.asarray(z, float)
        deep = np.interp(np.clip(z, zcs[0], zcs[-1]), zcs, tcs)
        return np.where(z < zcs[0], np.maximum(z, 0.0) / v0, deep)

    n_ch = 900
    depth = (np.arange(n_ch) - wh) * dx
    ch_hi = min(bot, int(wh + FIBRE_IN_WELL / dx))
    ok = np.zeros(n_ch, bool)
    ok[top + 2:ch_hi] = True                 # +2: clear of the coupling transition
    return dict(wh=wh, dx=dx, depth=depth, tau=tau, ok=ok, v0=v0,
                ch_lo=top + 2, ch_hi=ch_hi, zcs=zcs, tcs=tcs)


# --------------------------------------------------------------------------
# arrival picking and moveout
# --------------------------------------------------------------------------
def stalta(x, fs, sta=0.30, lta=2.50):
    """Leading-STA over trailing-LTA. Peaks at an onset, not at the maximum."""
    e = np.asarray(x, float) ** 2
    n = e.size
    c = np.concatenate([[0.0], np.cumsum(e)])
    i = np.arange(n)
    ns, nl = max(int(sta * fs), 1), max(int(lta * fs), 1)
    j = np.minimum(i + ns, n)
    k = np.maximum(i - nl, 0)
    S = (c[j] - c[i]) / np.maximum(j - i, 1)
    L = (c[i] - c[k]) / np.maximum(i - k, 1)
    r = S / np.maximum(L, 1e-30)
    r[:nl] = 0.0
    return r


def pick_p(stack, fs, t0=PRE_S - 0.5, t1=PRE_S + 8.0):
    """First arrival on the aligned stack.

    Take the STA/LTA maximum inside the search range, then walk BACK to the
    earliest sample within 3 s of it that still clears a detection level. The
    maximum alone frequently lands on S, which is larger; the walk-back targets P,
    which is what section 2.7 of METHODS_STATUS says to weight toward -- a vertical
    fibre is nearly blind to vertically-incident S anyway.
    """
    r = stalta(stack, fs)
    i0, i1 = int(t0 * fs), min(int(t1 * fs), r.size)
    if i1 <= i0 + 2:
        return np.nan, np.nan
    seg = r[i0:i1]
    kmax = i0 + int(np.argmax(seg))
    lvl = max(3.0, 0.40 * r[kmax])
    back = max(kmax - int(3.0 * fs), i0)
    cand = np.where(r[back:kmax + 1] > lvl)[0]
    k = back + int(cand[0]) if cand.size else kmax
    return k / fs, float(r[kmax])


def _seg(x, fs, t_c, w=WIN):
    i0 = int(round((t_c + w[0]) * fs))
    i1 = i0 + int(round((w[1] - w[0]) * fs))
    if i0 < 0 or i1 > x.size:
        return None
    return x[i0:i1]


def _xcorr_lag(a, b, fs, pad_s):
    """Lag of b relative to a, parabolically refined. Positive = b later."""
    n = min(a.size, b.size)
    a, b = a[:n] - a[:n].mean(), b[:n] - b[:n].mean()
    na, nb = np.sqrt(np.sum(a * a)), np.sqrt(np.sum(b * b))
    if na <= 0 or nb <= 0:
        return np.nan, 0.0
    c = np.correlate(b, a, 'full') / (na * nb)
    mid = n - 1
    pad = max(int(pad_s * fs), 1)
    lo = max(mid - pad, 0)
    seg = c[lo:mid + pad + 1]
    if seg.size < 3:
        return np.nan, 0.0
    k = int(np.argmax(seg))
    peak = float(seg[k])
    sub = 0.0
    if 0 < k < seg.size - 1:
        den = seg[k - 1] - 2 * seg[k] + seg[k + 1]
        if den != 0:
            sub = 0.5 * (seg[k - 1] - seg[k + 1]) / den
    return (k + lo - mid + sub) / fs, peak


def _median_smooth(y, w):
    """Running median that ignores NaN and leaves all-NaN neighbourhoods NaN."""
    n = y.size
    out = np.full(n, np.nan)
    h = w // 2
    for i in range(n):
        s = y[max(i - h, 0):min(i + h + 1, n)]
        s = s[np.isfinite(s)]
        if s.size:
            out[i] = np.median(s)
    return out


def moveout(Xb, fs, reg, n_iter=2):
    """Per-channel direct-arrival time (s from cache start), and the aligned stack.

    Starts from the moveout PREDICTED by the 2005 check shot -- for an upgoing ray
    the arrival time is a constant minus tau(z), because the deepest channel is
    reached first -- then refines each channel against the aligned stack. Starting
    from the prediction rather than from flat avoids the cycle-skipping that
    channel-to-channel tracking suffers over a 900-channel aperture.
    """
    ok = reg['ok']
    trel = -reg['tau'](reg['depth'])
    ref = int(np.median(np.where(ok)[0]))
    tarr = trel - trel[ref]                     # relative for now

    nrm = np.sqrt(np.sum(Xb ** 2, axis=1))
    use = ok & (nrm > 0) & np.isfinite(nrm)

    def stack_on(tarr_abs, t_ref):
        acc = np.zeros(Xb.shape[1])
        cnt = 0
        for k in np.where(use)[0]:
            sh = int(round((tarr_abs[k] - t_ref) * fs))
            acc += np.roll(Xb[k] / nrm[k], -sh)
            cnt += 1
        return acc / max(cnt, 1)

    t_ref, snr_pk = pick_p(stack_on(tarr + PRE_S, PRE_S), fs)
    if not np.isfinite(t_ref):
        return None, None, np.nan
    tarr = tarr + t_ref                          # absolute

    for _ in range(n_iter):
        s = stack_on(tarr, t_ref)
        ss = _seg(s, fs, t_ref, (WIN[0] - 0.25, WIN[1] + 0.25))
        if ss is None:
            break
        r = np.full(Xb.shape[0], np.nan)
        for k in np.where(use)[0]:
            xk = _seg(Xb[k], fs, tarr[k], (WIN[0] - 0.25, WIN[1] + 0.25))
            if xk is None:
                continue
            lag, pk = _xcorr_lag(ss, xk, fs, 0.15)
            if np.isfinite(lag) and pk > 0.3:
                r[k] = lag
        r = _median_smooth(r, SMOOTH_CH)
        tarr = tarr + np.where(np.isfinite(r), r, 0.0)

    return tarr, stack_on(tarr, t_ref), t_ref


# --------------------------------------------------------------------------
# per-pair delay profile
# --------------------------------------------------------------------------
def load_event(tag, band):
    f = os.path.join(CACHE, f'{tag}.npz')
    if not os.path.exists(f):
        return None, None
    d = np.load(f)
    fs = float(d['fs'])
    return bandpass(d['X'].astype(np.float64), fs, band), fs


def channel_snr(Xb, fs, tarr):
    """Arrival-window RMS over pre-event RMS, per channel."""
    i0, i1 = int((PRE_S + NOISE_WIN[0]) * fs), int((PRE_S + NOISE_WIN[1]) * fs)
    noise = np.sqrt(np.mean(Xb[:, max(i0, 0):i1] ** 2, axis=1))
    sig = np.full(Xb.shape[0], np.nan)
    for k in range(Xb.shape[0]):
        s = _seg(Xb[k], fs, tarr[k]) if np.isfinite(tarr[k]) else None
        if s is not None:
            sig[k] = np.sqrt(np.mean(s ** 2))
    return sig / np.where(noise > 0, noise, np.nan), noise, sig


def equalise_snr(Xb, fs, snr, noise, reg, target, seed):
    """Add band-limited noise so channels above `target` SNR are degraded to it.

    Only degradation is possible, so the cleaner band is brought down to the
    noisier one. The noise is filtered into the timing band, otherwise it would sit
    largely outside the band the phase-slope estimator actually uses and the
    equalisation would be cosmetic.
    """
    rng = np.random.default_rng(seed)
    Y = Xb.copy()
    hit = np.where(np.isfinite(snr) & (snr > target) & reg['ok'])[0]
    for k in hit:
        want = noise[k] * snr[k] / target        # required total noise RMS
        add = np.sqrt(max(want ** 2 - noise[k] ** 2, 0.0))
        if add <= 0:
            continue
        z = bandpass(rng.normal(0.0, 1.0, Xb.shape[1]), fs, PICK_BAND)
        z_rms = np.sqrt(np.mean(z ** 2))
        if z_rms > 0:
            Y[k] += z * (add / z_rms)
    return Y, hit.size


def delay_profile(Xa, Xb, fs, reg, band, acausal=False):
    """Per-channel delay of B relative to A on the direct arrival.

    The arrival curve is measured on A and reused for B after a bulk shift. Two
    repeating events have the same moveout by construction, so re-picking on B
    would only inject pick noise into the differential -- the quantity being
    measured is the difference between them, and it must not be re-estimated
    independently on each side.
    """
    tarr, stack_a, t_ref = moveout(Xa, fs, reg)
    if tarr is None:
        return None

    # bulk shift from the array stack: origin-time bookkeeping, not signal
    sb = np.zeros(Xb.shape[1])
    cnt = 0
    for k in np.where(reg['ok'])[0]:
        n = np.sqrt(np.sum(Xb[k] ** 2))
        if n > 0:
            sb += np.roll(Xb[k] / n, -int(round((tarr[k] - t_ref) * fs)))
            cnt += 1
    sb /= max(cnt, 1)
    bulk, _ = _xcorr_lag(stack_a, sb, fs, 2.0)
    if not np.isfinite(bulk):
        bulk = 0.0

    dt = np.full(Xa.shape[0], np.nan)
    coh = np.zeros(Xa.shape[0])
    for k in np.where(reg['ok'])[0]:
        a = _seg(Xa[k], fs, tarr[k])
        b = _seg(Xb[k], fs, tarr[k] + bulk)
        if a is None or b is None:
            continue
        if acausal:
            b = b[::-1]
        d, c = sub_sample_delay(a, b, fs, band, max_lag_s=0.15)
        if np.isfinite(d):
            dt[k], coh[k] = d, c

    good = np.isfinite(dt) & (coh > MIN_COH) & reg['ok']
    if good.sum() >= 20:
        # Cycle-skip rejection about the robust centre, not about zero. The fixed
        # guard alone is not enough: at 5-30 Hz a full period is ~59 ms and a
        # half-period skip ~29 ms, but the per-channel precision actually achieved
        # is nearer 2 ms, so a surviving 15 ms outlier is a 7-sigma point and a
        # least-squares slope over a 65 ms interval cannot absorb it. The MAD cut
        # adapts to the precision the data really has instead of to a period.
        dt = dt - np.median(dt[good])
        good &= np.abs(dt) < MAX_DELAY_S
        if good.sum() >= 20:
            c0 = np.median(dt[good])
            s0 = 1.4826 * np.median(np.abs(dt[good] - c0))
            if np.isfinite(s0) and s0 > 0:
                good &= np.abs(dt - c0) < 4.0 * s0
    return dict(dt=dt, coh=coh, good=good, tarr=tarr, bulk=bulk, t_ref=t_ref)


def intervals_from(prof, reg, spans, weights=None):
    """slope_dvv on every span, plus the two-band interval_dvv as a cross-check."""
    out = []
    dt, good = prof['dt'], prof['good']
    z, tau = reg['depth'], reg['tau'](reg['depth'])
    w = prof['coh'] ** 2 if weights is None else weights
    for zl, zh in spans:
        m = good & (z >= zl) & (z <= zh)
        if m.sum() < 6:
            out.append(dict(z_lo=zl, z_hi=zh, dvv=np.nan, err=np.nan,
                            n=int(m.sum()), rms=np.nan, dvv_edge=np.nan,
                            tau_s=np.nan))
            continue
        d, e, n, rms = slope_dvv(dt[m], tau[m], w[m])
        half = min(15.0, max(0.30 * (zh - zl), 5.0), max(zl - 1.0, 5.0))
        de, _ = interval_dvv(np.where(good, dt, np.nan), z, reg['tau'],
                             zl, zh, half=half)
        out.append(dict(z_lo=zl, z_hi=zh, dvv=d, err=e, n=n, rms=rms,
                        dvv_edge=de,
                        tau_s=float(reg['tau'](zh) - reg['tau'](zl))))
    return out


# --------------------------------------------------------------------------
def main():
    reg = registration()
    print('G5 -- shallow differential-timing dv/v on DAS\n')
    print(f'registration: wellhead ch {reg["wh"]:.0f}, usable ch '
          f'{reg["ch_lo"]}-{reg["ch_hi"]}, depth '
          f'{reg["depth"][reg["ch_lo"]]:.0f}-{reg["depth"][reg["ch_hi"]-1]:.0f} m')
    print(f'top-layer velocity from check shot: {reg["v0"]:.0f} m/s '
          f'(0-{reg["zcs"][0]:.1f} m)')
    print('\ninterval one-way travel time, and the delay a 1% change would make:')
    for zl, zh in INTERVALS:
        t = float(reg['tau'](zh) - reg['tau'](zl))
        print(f'  {zl:6.0f}-{zh:6.0f} m : tau {1000*t:7.2f} ms  '
              f'-> 1% = {1000*0.01*t:6.3f} ms')
    print()

    R = pd.read_csv(os.path.join(HERE, 'hrsn_control.csv'))
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    ev['time'] = pd.to_datetime(ev['time'], utc=True, format='mixed')
    for c in ('t_i', 't_j'):
        R[c] = pd.to_datetime(R[c], utc=True, format='mixed')
    cand = R[R.is_cand & (R.hrsn > 0.78)].copy()
    ctrl = R[~R.is_cand].copy()
    print(f'{len(cand)} confirmed repeater pairs, {len(ctrl)} random controls\n',
          flush=True)

    rows, prof_store = [], {}
    for label, D, is_rep in [('repeater', cand, True), ('control', ctrl, False)]:
        for _, r in D.iterrows():
            ti, tj = ev.tag[int(r.i)], ev.tag[int(r.j)]
            tag = f'{r.t_i:%Y-%m-%d}/{r.t_j:%Y-%m-%d}'
            for band in BANDS:
                Xa, fs = load_event(ti, band)
                Xb, _ = load_event(tj, band)
                if Xa is None or Xb is None:
                    break
                prof = delay_profile(Xa, Xb, fs, reg, band)
                if prof is None:
                    break
                variants = {'main': prof}

                if is_rep:
                    variants['acausal'] = delay_profile(Xa, Xb, fs, reg, band,
                                                        acausal=True)
                    snr, noise, _ = channel_snr(Xa, fs, prof['tarr'])
                    zz = reg['depth']
                    ms = reg['ok'] & (zz >= SHALLOW[0]) & (zz <= SHALLOW[1])
                    md = reg['ok'] & (zz >= DEEP[0]) & (zz <= DEEP[1])
                    t_snr = float(np.nanmin([np.nanmedian(snr[ms]),
                                             np.nanmedian(snr[md])]))
                    if np.isfinite(t_snr) and t_snr > 0:
                        Ya, na_ = equalise_snr(Xa, fs, snr, noise, reg, t_snr, 1)
                        snb, nob, _ = channel_snr(Xb, fs, prof['tarr'])
                        Yb, nb_ = equalise_snr(Xb, fs, snb, nob, reg, t_snr, 2)
                        variants['snr_eq'] = delay_profile(Ya, Yb, fs, reg, band)
                    else:
                        na_ = nb_ = 0

                for vname, p in variants.items():
                    if p is None:
                        continue
                    for iv in intervals_from(p, reg, INTERVALS + PROFILE):
                        rows.append(dict(
                            kind=label, is_rep=is_rep, variant=vname,
                            band=f'{band[0]:.0f}-{band[1]:.0f}',
                            i=int(r.i), j=int(r.j), pair=tag,
                            dt_days=r.dt_days,
                            phase_days=seasonal_phase_days(r.t_i, r.t_j),
                            hrsn_cc=r.hrsn, das_cc=r.das, **iv))
                if is_rep and band == BANDS[0]:
                    prof_store[tag] = dict(dt=prof['dt'], good=prof['good'],
                                           coh=prof['coh'])
            if is_rep:
                sh = [x for x in rows if x['pair'] == tag and x['variant'] == 'main'
                      and x['band'] == f'{BANDS[0][0]:.0f}-{BANDS[0][1]:.0f}'
                      and (x['z_lo'], x['z_hi']) == SHALLOW]
                v = sh[0] if sh else dict(dvv=np.nan, err=np.nan, n=0)
                print(f'  {tag}  dt={r.dt_days:5.0f}d   shallow dv/v '
                      f'{100*v["dvv"]:+7.3f} +/- {100*v["err"]:.3f} %  '
                      f'({v["n"]} ch)', flush=True)

    D = pd.DataFrame(rows)
    if D.empty:
        print('no measurements'); return
    D.to_csv(os.path.join(HERE, 'g5_shallow_dvv.csv'), index=False)

    b0 = f'{BANDS[0][0]:.0f}-{BANDS[0][1]:.0f}'
    M = D[(D.band == b0) & (D.variant == 'main')]
    rep, con = M[M.is_rep], M[~M.is_rep]

    def agg(sub):
        """Combine pairs. Returns (mean, err, ivw_err, scatter_err, n, sigma_ms).

        Two errors are carried because they fail in opposite directions and the
        disagreement between them is diagnostic. The inverse-variance error trusts
        the per-pair bootstrap; if the per-channel delays are correlated (they
        share a wavefield) it is optimistic. The scatter of the pair values around
        their mean makes no such assumption but is itself noisy at n=10. The
        reported error is the LARGER of the two, so a claim never rests on the
        more flattering of two defensible numbers.
        """
        s = sub.dropna(subset=['dvv', 'err'])
        s = s[s.err > 0]
        if len(s) < 2:
            return (np.nan,) * 4 + (len(s), np.nan)
        w = 1.0 / s.err ** 2
        m = float(np.sum(w * s.dvv) / np.sum(w))
        e_ivw = float(np.sqrt(1.0 / np.sum(w)))
        e_sc = float(s.dvv.std(ddof=1) / np.sqrt(len(s)))
        sig = float(1000 * np.nanmedian(s.rms)) if 'rms' in s else np.nan
        return m, max(e_ivw, e_sc), e_ivw, e_sc, len(s), sig

    def A(sub):
        """(mean, err) only -- the two numbers most call sites want."""
        r = agg(sub)
        return r[0], r[1]

    print(f'\n{"="*86}\nDEPTH INTERVALS -- stacked over repeater pairs '
          f'({b0} Hz)\n{"="*86}')
    print(f'{"interval (m)":>16}{"tau ms":>9}{"dv/v %":>10}{"err %":>9}'
          f'{"[1/var":>8}{"scat]":>8}{"npair":>7}{"sigma ms":>10}'
          f'{"floor %":>10}')
    summary = {}
    for zl, zh in INTERVALS:
        s = rep[(rep.z_lo == zl) & (rep.z_hi == zh)]
        c = con[(con.z_lo == zl) & (con.z_hi == zh)].dvv.dropna()
        m, e, ei, es, n, sig = agg(s)
        cf = float(np.percentile(np.abs(c), 90)) if len(c) > 3 else np.nan
        tau_s = float(reg['tau'](zh) - reg['tau'](zl))
        summary[(zl, zh)] = dict(m=m, e=e, ei=ei, es=es, n=n, floor=cf,
                                 sig=sig, tau=tau_s)
        print(f'{zl:7.0f}-{zh:<8.0f}{1000*tau_s:9.2f}{100*m:10.3f}{100*e:9.3f}'
              f'{100*ei:8.3f}{100*es:8.3f}{n:7d}{sig:10.2f}{100*cf:10.3f}')
    print('  sigma ms = median per-channel delay residual about the fitted line;')
    print('  it is the precision that actually limits every number in this table.')

    print(f'\nFINE PROFILE ({b0} Hz) -- and the common-mode test')
    print(f'{"interval (m)":>16}{"tau ms":>9}{"dv/v %":>10}{"err %":>9}'
          f'{"scatter %":>11}{"npair":>7}   verdict')
    fine = {}
    for zl, zh in PROFILE:
        m, e, ei, es, n, sig = agg(rep[(rep.z_lo == zl) & (rep.z_hi == zh)])
        fine[(zl, zh)] = (m, e, es, n)
        # A seasonal change must vary from pair to pair: the pairs span 39-440 day
        # baselines at different seasonal phases, so they sample different amounts
        # of the cycle. A value that is large but IDENTICAL across all ten pairs
        # cannot be seasonal -- it is something the processing does to every pair
        # alike. Pair-to-pair scatter far below the mean is that signature.
        tag = ''
        if np.isfinite(m) and np.isfinite(es) and n >= 5:
            if abs(m) > 3 * e and abs(es * np.sqrt(n)) < 0.35 * abs(m):
                tag = '<- COMMON MODE, not seasonal'
            elif abs(m) > 3 * e:
                tag = '<- varies across pairs'
        print(f'{zl:7.0f}-{zh:<8.0f}'
              f'{1000*float(reg["tau"](zh)-reg["tau"](zl)):9.2f}'
              f'{100*m:10.3f}{100*e:9.3f}{100*es*np.sqrt(max(n,1)):11.3f}'
              f'{n:7d}   {tag}')

    # Does the shallow value track seasonal phase, as a seasonal signal must?
    sh = rep[(rep.z_lo == SHALLOW[0]) & (rep.z_hi == SHALLOW[1])].dropna(
        subset=['dvv'])
    if len(sh) >= 5:
        print('\nSEASONAL vs SYSTEMATIC, shallow interval')
        for k, nm in [('phase_days', 'seasonal phase'), ('dt_days', 'elapsed')]:
            rr = float(np.corrcoef(sh[k], sh.dvv)[0, 1])
            print(f'  r(dv/v, {nm:>14}) = {rr:+.2f}')
        print('  a seasonal change tracks phase; a processing systematic '
              'tracks neither')

    MIN_CTRL_N = 5           # a control on 2 pairs is not a control
    print('\nCONTROLS')
    ctrl_res = {}
    for v in ('acausal', 'snr_eq'):
        s = D[(D.band == b0) & (D.variant == v) & D.is_rep &
              (D.z_lo == SHALLOW[0]) & (D.z_hi == SHALLOW[1])]
        m, e = A(s)
        n = int(agg(s)[4])
        ctrl_res[v] = (m, e, n)
        note = '' if n >= MIN_CTRL_N else '   <- too few pairs to conclude'
        print(f'  ({v:8s}) shallow {SHALLOW[0]:.0f}-{SHALLOW[1]:.0f} m: '
              f'dv/v {100*m:+.3f} +/- {100*e:.3f} %  (n={n}){note}')
    sd = summary[SHALLOW]
    n_con = int(len(con[(con.z_lo == SHALLOW[0]) &
                        (con.z_hi == SHALLOW[1])].dvv.dropna()))
    print(f'  (random  ) shallow control floor, 90th pct |dv/v|: '
          f'{100*sd["floor"]:.3f} %  (n={n_con} of {len(ctrl)} pairs converged)')
    if n_con < 0.3 * len(ctrl):
        print('             most random pairs yield NO shallow measurement at '
              'all --\n             the coherence gate refuses to time '
              'dissimilar waveforms, which is\n             the control working, '
              'but it leaves the floor poorly sampled.')
    for bb in BANDS[1:]:
        bl = f'{bb[0]:.0f}-{bb[1]:.0f}'
        m, e = A(D[(D.band == bl) & (D.variant == 'main') & D.is_rep &
                   (D.z_lo == SHALLOW[0]) & (D.z_hi == SHALLOW[1])])
        print(f'  (band {bl:>6}) shallow: dv/v {100*m:+.3f} +/- {100*e:.3f} %'
              f'   <- must agree with {b0} Hz if the signal is in the medium')

    # ---------------------------------------------------------------- verdict
    deep = [summary[k] for k in INTERVALS if k[0] >= 250.0]
    dm = np.array([d['m'] for d in deep]); de_ = np.array([d['e'] for d in deep])
    shal = summary[SHALLOW]
    ac_m, ac_e, ac_n = ctrl_res['acausal']
    sq_m, sq_e, sq_n = ctrl_res['snr_eq']

    thr_s = max(2 * shal['e'], shal['floor'] if np.isfinite(shal['floor']) else 0)
    shallow_sig = np.isfinite(shal['m']) and abs(shal['m']) > thr_s
    deep_flat = bool(np.all(np.abs(dm) < np.maximum(2 * de_, 1e-12)))
    # An undersampled control cannot clear a result, but it must not condemn one
    # either -- the first run called the acausal test CONTAMINATED off two pairs.
    ac_clean = (ac_n < MIN_CTRL_N) or (not np.isfinite(ac_m)) or abs(ac_m) < thr_s
    ac_known = ac_n >= MIN_CTRL_N
    sq_holds = (sq_n >= MIN_CTRL_N and np.isfinite(sq_m)
                and abs(sq_m) > 0.5 * abs(shal['m']))
    band_agree = True
    for bb in BANDS[1:]:
        bl = f'{bb[0]:.0f}-{bb[1]:.0f}'
        m2, e2 = A(D[(D.band == bl) & (D.variant == 'main') & D.is_rep &
                     (D.z_lo == SHALLOW[0]) & (D.z_hi == SHALLOW[1])])
        if np.isfinite(m2) and abs(m2 - shal['m']) > 2 * np.hypot(e2, shal['e']):
            band_agree = False

    print(f'\n{"="*74}\nG5 VERDICT\n{"="*74}')
    print(f'  shallow {SHALLOW[0]:.0f}-{SHALLOW[1]:.0f} m : '
          f'{100*shal["m"]:+.3f} +/- {100*shal["e"]:.3f} %   '
          f'threshold {100*thr_s:.3f} %  -> '
          f'{"RESOLVED" if shallow_sig else "not resolved"}')
    for k in INTERVALS:
        if k[0] >= 250.0:
            print(f'  deep {k[0]:.0f}-{k[1]:.0f} m (artifact control): '
                  f'{100*summary[k]["m"]:+.3f} +/- {100*summary[k]["e"]:.3f} %')
    print(f'  acausal control  : {100*ac_m:+.3f} % (n={ac_n}, '
          f'{"clean" if ac_clean and ac_known else "undersampled" if not ac_known else "CONTAMINATED"})')
    print(f'  SNR-equalised    : {100*sq_m:+.3f} % (n={sq_n}, '
          f'{"contrast survives" if sq_holds else "does NOT survive"})')
    ba_msg = ('yes' if band_agree else
              'NO -- the bands disagree beyond their errors, which a\n'
              '                     medium change cannot do and a timing '
              'artifact can')
    print(f'  band agreement   : {ba_msg}')

    if shallow_sig and not band_agree:
        print('\n  -> REJECTED BY THE BAND CHECK. A real change in the medium '
              'shifts the\n     arrival by the same TIME in every passband. Two '
              'bands disagreeing\n     beyond their errors means the estimator is '
              'reading something\n     frequency-dependent -- cycle skipping, or '
              'the window catching a\n     different part of the waveform in each '
              'band.')
    elif shallow_sig and deep_flat and ac_clean and sq_holds:
        print('\n  -> SHALLOW-CONFINED VELOCITY CHANGE.\n'
              '     Delay accumulates above 75 m and is flat below 250 m, which\n'
              '     hypocentre separation cannot produce -- a source shift perturbs\n'
              '     the whole ray. Consistent with the G3 null on HRSN at 251-284 m\n'
              '     and with Li & Ben-Zion 2023 peak sensitivity near 17 m.')
    elif shallow_sig and not deep_flat:
        print('\n  -> SLOPE AT ALL DEPTHS. The deep intervals are supposed to be the\n'
              '     control and they are not flat, so this is most likely hypocentre\n'
              '     separation rather than a shallow medium change. Relocate the\n'
              '     pairs, or restrict to the pairs with the smallest separation,\n'
              '     before claiming a depth profile.')
    elif shallow_sig and not ac_clean:
        print('\n  -> REJECTED BY THE ACAUSAL CONTROL. A time-reversed template\n'
              '     cannot carry a real delay profile, so a comparable shallow slope\n'
              '     on it means the estimator is reading window shape, not medium.')
    elif shallow_sig and not sq_holds:
        print('\n  -> REJECTED BY SNR EQUALISATION. The shallow/deep contrast\n'
              '     disappears once the two bands have the same SNR, so it was the\n'
              '     SNR gradient along the fibre, not a velocity gradient.')
    else:
        lim = 100 * max(2 * shal['e'],
                        shal['floor'] if np.isfinite(shal['floor']) else 0)
        print(f'\n  -> NO RESOLVED CHANGE, even in the shallow interval.\n'
              f'     This is still a publishable number: with the G3 null at\n'
              f'     251-284 m it becomes a two-depth UPPER BOUND on any seasonal\n'
              f'     velocity change at SAFOD -- < {lim:.2f}% above 75 m and\n'
              f'     < 1.94% at HRSN depth, over baselines of 39-440 days.')

    np.savez(os.path.join(HERE, 'g5_shallow_dvv.npz'),
             depth=reg['depth'], tau=reg['tau'](reg['depth']), ok=reg['ok'],
             intervals=np.array(INTERVALS), profile=np.array(PROFILE),
             pairs=np.array(list(prof_store)),
             dt=np.array([prof_store[k]['dt'] for k in prof_store]),
             good=np.array([prof_store[k]['good'] for k in prof_store]),
             coh=np.array([prof_store[k]['coh'] for k in prof_store]))

    # ----------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 6.0))
    z = reg['depth']

    a = ax[0]
    for tag, p in prof_store.items():
        g = p['good']
        if g.sum() > 20:
            a.plot(1000 * p['dt'][g], z[g], lw=0.8, alpha=0.75, label=tag)
    a.axvline(0, color='k', lw=0.8)
    a.axhspan(SHALLOW[0], SHALLOW[1], color='C1', alpha=0.12)
    a.text(0.02, 0.02, 'target interval', transform=a.transAxes, color='C1',
           fontsize=8)
    a.invert_yaxis()
    a.set(xlabel='differential delay B-A (ms)', ylabel='depth below surface (m)',
          xlim=(-1000 * MAX_DELAY_S, 1000 * MAX_DELAY_S),
          title='A  Delay profile per repeater pair\n'
                'shallow-only change = slope above, flat below')
    a.legend(fontsize=6); a.grid(alpha=0.3)

    a = ax[1]
    zc = [0.5 * (zl + zh) for zl, zh in PROFILE]
    for lab, col, sub in [('repeaters', 'C3', rep), ('random pairs', '0.6', con)]:
        m = [agg(sub[(sub.z_lo == zl) & (sub.z_hi == zh)]) for zl, zh in PROFILE]
        a.errorbar([100 * x[0] for x in m], zc,
                   xerr=[100 * x[1] for x in m], fmt='o-', color=col, ms=5,
                   capsize=3, lw=1.4, label=lab)
    a.axvline(0, color='k', lw=0.8)
    a.axhline(250, color='C0', ls='--', lw=1.2)
    a.text(a.get_xlim()[0], 245, 'HRSN sensor depth (G3 null)', color='C0',
           fontsize=7, va='bottom')
    a.invert_yaxis()
    a.set(xlabel='interval dv/v (%)', ylabel='depth (m)',
          title='B  Depth profile of dv/v\n(stacked over pairs)')
    a.legend(fontsize=8); a.grid(alpha=0.3)

    a = ax[2]
    labs, vals, errs = [], [], []
    for lab, s in [('repeaters', rep[(rep.z_lo == SHALLOW[0]) &
                                     (rep.z_hi == SHALLOW[1])]),
                   ('acausal', D[(D.variant == 'acausal') & (D.band == b0) &
                                 (D.z_lo == SHALLOW[0]) & (D.z_hi == SHALLOW[1])]),
                   ('SNR-equalised', D[(D.variant == 'snr_eq') & (D.band == b0) &
                                       (D.z_lo == SHALLOW[0]) &
                                       (D.z_hi == SHALLOW[1])]),
                   ('random pairs', con[(con.z_lo == SHALLOW[0]) &
                                        (con.z_hi == SHALLOW[1])])]:
        r = agg(s)
        m, e, n = r[0], r[1], r[4]
        labs.append(f'{lab}\n(n={n})'); vals.append(100 * m); errs.append(100 * e)
    a.bar(range(len(labs)), vals, yerr=errs, capsize=4,
          color=['C3', '0.5', 'C0', '0.75'])
    a.axhline(0, color='k', lw=0.8)
    a.set_xticks(range(len(labs))); a.set_xticklabels(labs, fontsize=8)
    a.set(ylabel='dv/v (%)',
          title=f'C  Shallow interval {SHALLOW[0]:.0f}-{SHALLOW[1]:.0f} m\n'
                'against every control')
    a.grid(alpha=0.3, axis='y')

    fig.suptitle('G5: differential direct-arrival timing on shallow DAS channels '
                 '— how deep does the change reach?', fontsize=12)
    fig.tight_layout()
    p = os.path.join(HERE, 'g5_shallow_dvv.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
