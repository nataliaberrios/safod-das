"""
D0: can this array measure a depth-localized dv/v at all? Gate before any science.

THE MEASUREMENT BEING GATED. For a repeater pair, the per-channel delay dt(z)
accumulates along the ray, so inside an interval of uniform fractional velocity
change dv/v the delay is linear in one-way travel time tau:

    dt(z) = dt_0 - (dv/v) * tau(z)        =>   dv/v = -d(dt)/d(tau)

The intercept dt_0 absorbs the origin-time difference, so this is immune to the
unknown source instant. tau(z) = p * z is taken from the event's OWN measured
apparent slowness (moveout_test.slant_scan), so no velocity model is assumed.

WHY THIS IS WORTH MEASURING. Every published dv/v at Parkfield -- including the
2002-2022 record in doi:10.1029/2023jb028084 -- is volume-averaged over an unknown
scattering volume, and depth is inferred only from frequency dependence. A vertical
aperture localises the change to a known interval. Takagi & Okada 2012
(doi:10.1029/2012gl051342) did this on KiK-net vertical arrays, but with one or two
sensors per borehole, so they had to compare five stations of differing borehole
depth to get any depth information at all. This fiber has 900 receiver depths in one
hole.

TWO THINGS THAT MUST NOT BE GOT WRONG, both learned the hard way here.

1. BOTH EVENTS GET THE SAME ALIGNMENT SHIFT. Aligning each event with its own
   fitted slowness would remove exactly the differential moveout being measured.
   Applying one common shift to both preserves dt(z) identically while still
   flattening the wavefront for any per-bin stacking. The mean of the two fitted
   slownesses is used.

2. THE BOOTSTRAP IS BLOCKED AT ONE GAUGE LENGTH. 16.335 m / 1.021 m = 16 channels
   share the same averaged strain, so resampling channels singly understates the
   error by ~4x. dvv_core.slope_dvv(..., block=16). ray_parameter_test.py made
   precisely this mistake and its error bars were 4x too small.

--------------------------------------------------------------------------------
PASS CRITERIA, FIXED BEFORE RUNNING.

  G1  sigma(dv/v) from the random-pair null must be <= 0.1 %.
      Requirement comes from tau ~ 864 m / 2500 m/s = 0.35 s: resolving 0.1 %
      needs 0.35 ms of timing, and the predicted budget is 0.23 ms (Cramer-Rao
      calibrated against the measured 5.31 ms at the old 100 Hz configuration).
      The margin is only ~1.5x, which is why this is measured rather than assumed.

  G2  SYNTHETIC RECOVERY UNBIASED. Inject a known uniform dv/v of 0.1, 0.3 and
      1.0 % into one event of a pair and recover it within error. This is the
      single most important test in the plan: it validates the estimator against
      ground truth rather than against my expectations. Injection reuses
      moveout_test.align with slowness -eps*p, which applies exactly the
      shift -eps*p*z that a uniform fractional velocity change produces.

  G3  ACAUSAL FLOOR. Time-reverse one event of each repeater pair. A reversed
      record cannot contain a real interval arrival, so any apparent dv/v is pure
      method output and must be consistent with the null.

FAIL IS AN ACCEPTABLE OUTCOME. If sigma exceeds 0.1 % the deliverable becomes a
depth-localized UPPER BOUND on dv/v in 0-864 m, which no one has published for
Parkfield either. Report the achieved precision and stop; do not tune thresholds
until something passes.
--------------------------------------------------------------------------------
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.signal import decimate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import moveout_test as MT                                       # noqa: E402
from dvv_core import sub_sample_delay, slope_dvv, bulk_align     # noqa: E402

CACHE_HF = os.path.join(HERE, 'cache_hf')
FS_TARGET = 500.0
BAND = (5.0, 80.0)          # frozen: 5-80 Hz, the usable band per hf_snr_test
GAUGE_CH = 16               # 16.335 m gauge / 1.021 m spacing
WIN_P = (-0.2, 1.5)         # about the picked P, seconds
WIN_CODA = (2.0, 12.0)      # coda, from coda_window_survey's usable range
DELAY_MAX_LAG = 0.05        # per-channel RESIDUAL search only -- valid because
                            # prep_pair() removes the bulk inter-event lag first
MIN_COH = 0.3               # phase-stability floor from sub_sample_delay;
                            # 0.6 kept only 0-20 of ~700 channels
INJECT = [0.001, 0.003, 0.010]
WINDOWS = {'P': WIN_P, 'coda': WIN_CODA}

# ------------------------------------------------------------------------------
# G4 -- THE CONFOUND THAT DECIDES WHETHER THIS MEASUREMENT MEANS ANYTHING.
#
# Source offset and velocity change are EXACTLY degenerate in the direct arrival:
# both appear as a fractional change in apparent slowness, delta_p / p. For a
# source displaced by delta at hypocentral distance R and incidence i,
#
#       apparent dv/v  =  tan(i) * delta / R
#
# Computed for these 10 pairs with delta = one source radius -- i.e. for events
# that are genuine repeaters BY DEFINITION, not failed ones:
#
#       incidence 13.8-65.3 deg, contamination 0.088-1.408 %, median 0.958 %
#       only 1 of 10 pairs falls below the 0.1 % requirement
#
# So the direct-P estimator measures source position, not velocity, and cannot
# separate them. Averaging helps only as 1/sqrt(N) with random offset direction,
# which needs ~100 pairs; there are 10.
#
# THE ESCAPE, AND WHY IT IS THE PUBLISHED ONE. Poupinet et al. 1984 and Takagi &
# Okada 2012 both use CODA, not the direct arrival. Source separations here are
# 12-48 m against a 150 m wavelength at 20 Hz / 3 km/s, so the coda waveform is
# nearly unchanged -- which is exactly why these pairs correlate at 0.99 in the
# coda. The coda's angular spectrum is set by the scatterer distribution rather
# than by the precise source point, so tan(i) sensitivity is suppressed.
#
# THAT IS AN ARGUMENT, AND ARGUMENTS ARE WHAT KILLED THE PREVIOUS FOUR DIRECTIONS
# HERE. So it is measured instead. The predicted contamination is known per pair,
# so the test is a regression:
#
#   G4  r(measured dv/v, predicted contamination) must be STRONG in the P window
#       and WEAK in the coda window.
#
#       strong in P  -> the confound is real and the estimator is reading geometry
#       weak in coda -> the coda estimator is not, and is the one to use
#       strong in BOTH -> the measurement is geometry everywhere. Report the
#                         contamination as an upper bound on dv/v and stop.
#       weak in BOTH -> either the confound model is wrong or nothing is resolved;
#                         check against the null before believing anything
#
# There is no threshold to tune here: the sign and relative size across the two
# windows carry the information, and all four outcomes are stated in advance.
# ------------------------------------------------------------------------------


def load(tag):
    """Load one event at the common rate. Some files are 5000 Hz, not 500."""
    f = os.path.join(CACHE_HF, f'{tag}.npz')
    if not os.path.exists(f):
        return None
    d = np.load(f)
    X, fs = d['X'].astype(np.float64), float(d['fs'])
    if fs > FS_TARGET * 1.5:
        q = int(round(fs / FS_TARGET))
        if abs(fs / q - FS_TARGET) > 1:
            return None
        X = decimate(X, q, axis=-1, ftype='fir', zero_phase=True)
        fs = fs / q
    if abs(fs - FS_TARGET) > 1:
        return None
    return X, fs


def prep_pair(tag_a, tag_b):
    """Both events, bandpassed, on common channels, with ONE shared alignment."""
    ra, rb = load(tag_a), load(tag_b)
    if ra is None or rb is None:
        return None
    MT.BAND = BAND
    MT.CACHE = CACHE_HF
    ga, _ = MT.prep(tag_a)
    gb, _ = MT.prep(tag_b)
    if ga is None or gb is None:
        return None
    Aa, za, fsa = ga
    Ab, zb, fsb = gb
    # prep() drops dead channels independently per event, so intersect on depth
    zc, ia, ib = np.intersect1d(za, zb, return_indices=True)
    if zc.size < 100:
        return None
    Aa, Ab = Aa[ia], Ab[ib]
    if abs(fsa - fsb) > 1:
        return None

    pa, tpa, _sa, _S, _ps = MT.slant_scan(Aa, zc, fsa)
    pb, tpb, _sb, _S2, _p2 = MT.slant_scan(Ab, zc, fsb)
    pbar = 0.5 * (pa + pb)
    # SAME shift to both -- see docstring point 1
    Ba = MT.align(Aa, zc, fsa, pbar)
    Bb = MT.align(Ab, zc, fsb, pbar)

    # BULK-ALIGN THE TWO EVENTS TO EACH OTHER. Aligning channels for moveout does
    # NOT align event B to event A: catalog origin times carry 0.1-0.5 s of error,
    # and the per-channel residual search below is only +/-50 ms. Omitting this
    # returned 0-20 usable channels out of ~700 with 46 ms rms. dvv_core.bulk_align
    # says "NECESSARY, NOT OPTIONAL" and names dvv_hrsn and coda_window_survey as
    # the two scripts that forgot it; this was the third.
    #
    # The shift is rigid -- identical for every channel -- so it lands entirely in
    # the intercept of dt(z) = a + b*z and cannot bias the slope b, which is the
    # measurement. It is bookkeeping, not signal.
    #
    # Stored as a sample offset applied to B's WINDOW INDICES rather than by
    # np.roll: rolling a cut window folds its tail back in at the front, which is a
    # bug already made once in beta_similarity.py.
    sa_ = Ba.mean(axis=0)
    sb_ = Bb.mean(axis=0)
    _, lag_s = bulk_align(sa_, sb_, fsa, max_lag_s=2.0)
    shift = int(round(lag_s * fsa))
    return dict(A=Ba, B=Bb, z=zc, fs=fsa, p=pbar, tp=0.5 * (tpa + tpb),
                pa=pa, pb=pb, shift=shift, lag_s=float(lag_s))


def per_channel_delays(P, invert_eps=None, reverse_b=False, win=WIN_P):
    """dt(z) between the two events, channel by channel. Returns (dt, tau, coh)."""
    A, B, z, fs, p = P['A'], P['B'], P['z'], P['fs'], P['p']
    if invert_eps is not None:
        # inject a uniform fractional velocity change: shift -eps*p*z, which is
        # exactly align() with slowness -eps*p
        B = MT.align(B, z, fs, -invert_eps * p)
    if reverse_b:
        B = B[:, ::-1]
    i0 = int((P['tp'] + win[0]) * fs)
    i1 = int((P['tp'] + win[1]) * fs)
    # B's window is offset by the rigid inter-event lag; slicing rather than
    # rolling so nothing wraps around
    sh = int(P.get('shift', 0))
    j0, j1 = i0 + sh, i1 + sh
    if i0 < 0 or i1 > A.shape[1] or j0 < 0 or j1 > B.shape[1]:
        return None, None, None
    dt = np.full(A.shape[0], np.nan)
    coh = np.zeros(A.shape[0])
    for k in range(A.shape[0]):
        d_, c_ = sub_sample_delay(A[k, i0:i1], B[k, j0:j1], fs, BAND,
                                  max_lag_s=DELAY_MAX_LAG)
        dt[k], coh[k] = d_, c_
    tau = p * z                      # one-way travel time from measured slowness
    return dt, tau, coh


def measure(P, **kw):
    dt, tau, coh = per_channel_delays(P, **kw)
    if dt is None:
        return dict(dvv=np.nan, err=np.nan, n=0, rms=np.nan)
    m = np.isfinite(dt) & (coh > MIN_COH)
    if m.sum() < 20:
        return dict(dvv=np.nan, err=np.nan, n=int(m.sum()), rms=np.nan)
    dvv, err, n, rms = slope_dvv(dt[m], tau[m], weights=coh[m],
                                 block=GAUGE_CH, n_boot=400)
    return dict(dvv=dvv, err=err, n=n, rms=rms, tau_span=float(np.ptp(tau[m])))


def predicted_contamination(ev, pa, i_idx, j_idx, z_arr=0.43):
    """Apparent dv/v a one-source-radius offset would produce: tan(i)*delta/R."""
    out = []
    for k in (i_idx, j_idx):
        e = ev.iloc[k]
        X = pa['km_safod'].get(e.tag, np.nan)
        inc = np.degrees(np.arctan2(X, e.depth - z_arr))
        M0 = 10 ** (1.5 * e.mag + 9.1)
        r = (7 * M0 / (16 * 3e6)) ** (1 / 3)
        R = np.sqrt(X ** 2 + (e.depth - z_arr) ** 2) * 1000.0
        out.append((inc, r, R))
    inc = 0.5 * (out[0][0] + out[1][0])
    r = max(out[0][1], out[1][1])
    R = 0.5 * (out[0][2] + out[1][2])
    return float(np.tan(np.radians(inc)) * r / R), float(inc), float(r), float(R)


def offset_from_slope(dpp, inc_deg, R):
    """Invert delta_p/p = tan(i) * delta / R for the relative source offset.

    THE SAME-PATCH TEST. This is the identical measurement as the dv/v above, read
    the other way round: what was 'contamination' for a medium change IS the signal
    for source geometry. With 500 Hz + aligned stacking the block-bootstrapped
    sigma(delta_p/p) is ~2.7e-4, which over these incidences (15-65 deg) and
    distances (2.9-8.5 km) resolves a relative offset of 1.2-13.7 m, median 1.9 m
    -- against source radii of 14-48 m. So the array can say whether two ruptures
    overlap, which waveform similarity alone cannot (Gao, Kao & Wang 2021,
    doi:10.1029/2021gl092815).

    A velocity change contributes to delta_p/p as well and is NOT separable per
    pair. It is separable at the population level, because it is common to family
    and control pairs while the source offset is not -- which is why the
    family-vs-control comparison, not the single number, is the test.
    """
    t = np.tan(np.radians(inc_deg))
    if not np.isfinite(t) or abs(t) < 1e-6:
        return np.nan
    return float(abs(dpp) * R / t)


def main():
    mt = pd.read_csv(os.path.join(HERE, 'moveout_test.csv'))
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    pa = pd.read_csv(os.path.join(HERE, 'phaseA_events.csv')
                     ).drop_duplicates('tag').set_index('tag')
    have = {f[:-4] for f in os.listdir(CACHE_HF) if f.endswith('.npz')}
    mt = mt[mt.i.map(lambda k: ev.tag[k] in have) &
            mt.j.map(lambda k: ev.tag[k] in have)].copy()
    rep = mt[mt.is_cand]
    ctl = mt[~mt.is_cand]
    print(f'{len(rep)} repeater pairs, {len(ctl)} control pairs usable at 500 Hz')
    print(f'band {BAND[0]:.0f}-{BAND[1]:.0f} Hz, bootstrap block {GAUGE_CH} ch\n',
          flush=True)

    rows = []
    cache = {}

    def get(i, j):
        key = (i, j)
        if key not in cache:
            cache[key] = prep_pair(ev.tag[i], ev.tag[j])
        return cache[key]

    for wname, win in WINDOWS.items():
        print(f'--- repeater pairs, {wname} window '
              f'{win[0]:+.1f} to {win[1]:+.1f} s ---', flush=True)
        print(f'{"events":>22}{"dv/v %":>10}{"err %":>9}{"nch":>6}'
              f'{"rms ms":>9}{"predContam %":>14}{"i deg":>8}')
        for _, r in rep.iterrows():
            P = get(int(r.i), int(r.j))
            if P is None:
                continue
            m = measure(P, win=win)
            contam, inc, rr_, RR_ = predicted_contamination(ev, pa,
                                                            int(r.i), int(r.j))
            lab = f'{ev.tag[int(r.i)][3:11]}/{ev.tag[int(r.j)][3:11]}'
            print(f'{lab:>22}{100*m["dvv"]:10.4f}{100*m["err"]:9.4f}{m["n"]:6d}'
                  f'{1e3*m.get("rms",np.nan):9.3f}{100*contam:14.3f}{inc:8.1f}',
                  flush=True)
            rows.append(dict(kind='repeater', win=wname, i=int(r.i), j=int(r.j),
                             lab=lab, contam=contam, inc=inc, rad=rr_, R=RR_,
                             offset_m=offset_from_slope(m['dvv'], inc, RR_),
                             **m))

        for _, r in rep.iterrows():
            P = get(int(r.i), int(r.j))
            if P is None:
                continue
            m = measure(P, reverse_b=True, win=win)
            rows.append(dict(kind='acausal', win=wname, i=int(r.i), j=int(r.j),
                             lab='rev', contam=np.nan, inc=np.nan,
                             rad=np.nan, R=np.nan, offset_m=np.nan, **m))
        for _, r in ctl.iterrows():
            P = get(int(r.i), int(r.j))
            if P is None:
                continue
            m = measure(P, win=win)
            c_, inc_, rr_, RR_ = predicted_contamination(ev, pa,
                                                         int(r.i), int(r.j))
            rows.append(dict(kind='null', win=wname, i=int(r.i), j=int(r.j),
                             lab='ctl', contam=c_, inc=inc_, rad=rr_, R=RR_,
                             offset_m=offset_from_slope(m['dvv'], inc_, RR_),
                             **m))
        print(flush=True)

    D0 = pd.DataFrame(rows)
    print('--- G1 null / G3 acausal, by window ---')
    print(f'{"window":>8}{"null n":>8}{"null med %":>12}{"null sd %":>11}'
          f'{"acausal med %":>15}')
    sigma_by_win = {}
    for wname in WINDOWS:
        nu = D0[(D0.kind == 'null') & (D0.win == wname)].dvv.dropna()
        ac = D0[(D0.kind == 'acausal') & (D0.win == wname)].dvv.dropna()
        sd = float(nu.std()) if len(nu) > 3 else np.nan
        sigma_by_win[wname] = sd
        print(f'{wname:>8}{len(nu):8d}{100*nu.median():12.4f}{100*sd:11.4f}'
              f'{100*ac.median() if len(ac) else np.nan:15.4f}')
    # the P window is the historical/default configuration; keep its sigma as the
    # headline so the pass/fail below is not silently chosen after the fact
    sigma_null = sigma_by_win.get('P', np.nan)
    nu = list(D0[(D0.kind == 'null') & (D0.win == 'P')].dvv.dropna())
    ac = list(D0[(D0.kind == 'acausal') & (D0.win == 'P')].dvv.dropna())

    print('\n--- G4 does the measurement track the geometric confound? ---')
    print(f'{"window":>8}{"n":>5}{"r(dvv,contam)":>16}{"slope":>10}'
          f'{"|dvv| med %":>13}{"contam med %":>14}')
    g4 = {}
    for wname in WINDOWS:
        s = D0[(D0.kind == 'repeater') & (D0.win == wname)].dropna(
            subset=['dvv', 'contam'])
        if len(s) < 4:
            continue
        rr = float(np.corrcoef(np.abs(s.dvv), s.contam)[0, 1])
        sl = float(np.polyfit(s.contam, np.abs(s.dvv), 1)[0])
        g4[wname] = rr
        print(f'{wname:>8}{len(s):5d}{rr:16.3f}{sl:10.2f}'
              f'{100*np.median(np.abs(s.dvv)):13.4f}'
              f'{100*s.contam.median():14.3f}')
    print('  interpretation fixed in advance (see header): strong in P and weak in')
    print('  coda -> use coda. Strong in both -> the number is geometry; report a')
    print('  bound. Weak in both -> check against the null before believing it.')

    # ---------------------------------------------------------------- same-patch
    # The identical slope, read as source geometry instead of as a medium change.
    # This is the test Ellsworth's "repeaters or neighbours" question actually
    # asks, and the one waveform similarity cannot answer (Gao et al. 2021).
    #
    # A velocity change also enters the slope and is NOT separable per pair. It IS
    # separable here, because it is common to family and control pairs while the
    # source offset is not. So the statistic is the family-vs-control CONTRAST,
    # never a single pair's number.
    print('\n--- SAME-PATCH TEST: inferred relative source offset ---')
    print(f'{"window":>8}{"group":>10}{"n":>5}{"offset med m":>14}'
          f'{"offset IQR m":>14}{"rupture rad m":>15}{"overlap?":>10}')
    same_patch = {}
    for wname in WINDOWS:
        for grp, kind in (('family', 'repeater'), ('control', 'null')):
            s = D0[(D0.kind == kind) & (D0.win == wname)].dropna(
                subset=['offset_m', 'rad'])
            if len(s) < 3:
                continue
            med = float(np.median(s.offset_m))
            q1, q3 = np.percentile(s.offset_m, [25, 75])
            rad = float(np.median(s.rad))
            print(f'{wname:>8}{grp:>10}{len(s):5d}{med:14.1f}'
                  f'{q3-q1:14.1f}{rad:15.1f}'
                  f'{"YES" if med < rad else "no":>10}')
            same_patch[(wname, grp)] = (med, len(s), s.offset_m.values)
    print('\n  contrast (this is the statistic, not the individual numbers):')
    for wname in WINDOWS:
        f_ = same_patch.get((wname, 'family'))
        c_ = same_patch.get((wname, 'control'))
        if not (f_ and c_):
            continue
        try:
            from scipy.stats import mannwhitneyu
            u, pv = mannwhitneyu(f_[2], c_[2], alternative='less')
        except Exception:
            pv = np.nan
        print(f'    {wname:>6}: family {f_[0]:7.1f} m vs control {c_[0]:7.1f} m'
              f'   ratio {f_[0]/max(c_[0],1e-9):5.2f}   Mann-Whitney p={pv:.4f}')
    print('  family offsets must be BOTH smaller than controls AND below the')
    print('  rupture radius for the same-patch assumption to hold. If family and')
    print('  control are indistinguishable, waveform similarity is not selecting')
    print('  co-located ruptures and recurrence-based creep inherits that.')

    print('\n--- G2 synthetic recovery (inject into event B) ---', flush=True)
    print(f'{"injected %":>12}{"recovered %":>14}{"err %":>9}{"bias %":>9}'
          f'{"npairs":>8}')
    recov = {}
    for eps in INJECT:
        got, ers = [], []
        for _, r in rep.iterrows():
            P = get(int(r.i), int(r.j))
            if P is None:
                continue
            m = measure(P, invert_eps=eps)          # default win is WIN_P
            base = next((x['dvv'] for x in rows if x['kind'] == 'repeater'
                         and x['win'] == 'P'
                         and x['i'] == int(r.i) and x['j'] == int(r.j)), np.nan)
            if np.isfinite(m['dvv']) and np.isfinite(base):
                got.append(m['dvv'] - base)      # differential removes any real dv/v
                ers.append(m['err'])
        if got:
            recov[eps] = (float(np.median(got)), float(np.median(ers)), len(got))
            print(f'{100*eps:12.3f}{100*recov[eps][0]:14.4f}'
                  f'{100*recov[eps][1]:9.4f}{100*(recov[eps][0]-eps):9.4f}'
                  f'{recov[eps][2]:8d}')

    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(HERE, 'interval_dvv_gate.csv'), index=False)

    print('\n' + '=' * 74)
    print('VERDICT against criteria fixed before running')
    g1 = np.isfinite(sigma_null) and sigma_null <= 0.001
    print(f'  G1 null sigma <= 0.100 % : {100*sigma_null:.4f} %   '
          f'{"PASS" if g1 else "FAIL"}')
    g2 = True
    for eps in INJECT:
        if eps not in recov:
            g2 = False; continue
        bias, er, _ = recov[eps]
        ok = abs(bias - eps) <= max(2 * er, 0.2 * eps)
        g2 &= ok
        print(f'  G2 recover {100*eps:.1f} %      : bias '
              f'{100*(bias-eps):+.4f} %  {"PASS" if ok else "FAIL"}')
    g3 = bool(ac) and abs(np.median(ac)) <= 2 * (sigma_null if np.isfinite(sigma_null) else np.inf)
    print(f'  G3 acausal ~ null        : {"PASS" if g3 else "FAIL"}')
    rP, rC = g4.get('P', np.nan), g4.get('coda', np.nan)
    usable = None
    if np.isfinite(rP) and np.isfinite(rC):
        if rP > 0.5 and rC < 0.3:
            usable = 'coda'
            verdict = 'confound real in P, suppressed in coda -> USE CODA'
        elif rP > 0.5 and rC >= 0.3:
            verdict = 'confound present in BOTH -> the number is geometry, report a bound'
        elif rP <= 0.5 and rC < 0.3:
            verdict = ('confound not detected in either -- check |dvv| against the '
                       'null before believing it')
        else:
            verdict = 'unexpected pattern (weak in P, strong in coda) -- investigate'
        print(f'  G4 r(P)={rP:+.3f}  r(coda)={rC:+.3f} : {verdict}')
    print()
    if g1 and g2 and g3 and usable:
        print(f'  -> GATE PASSED using the {usable} window. Freeze band, channel')
        print('     range and window, then proceed to D1/D2.')
    else:
        print('  -> GATE NOT PASSED. The deliverable is a depth-localized UPPER')
        print('     BOUND on dv/v over 0-864 m. Report the achieved precision and')
        print('     the geometric contamination alongside it.')
        print('     Do not retune thresholds to force a pass.')

    fig, ax = plt.subplots(1, 4, figsize=(20, 4.6))
    for kind, c in [('repeater', 'C3'), ('null', '0.6'), ('acausal', 'C0')]:
        v = D[(D.kind == kind) & (D.win == 'P') & np.isfinite(D.dvv)].dvv * 100
        if len(v):
            ax[0].hist(v, bins=18, alpha=.6, color=c, label=f'{kind} (n={len(v)})')
    ax[0].set(xlabel='dv/v (%)', ylabel='pairs', title='A  P window: signal vs null')
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    if recov:
        xs = [100 * e for e in recov]
        ys = [100 * recov[e][0] for e in recov]
        es = [100 * recov[e][1] for e in recov]
        ax[1].errorbar(xs, ys, yerr=es, fmt='o-', color='C2')
        lim = [0, max(xs) * 1.1]
        ax[1].plot(lim, lim, 'k--', lw=1, label='y = x')
        ax[1].set(xlabel='injected dv/v (%)', ylabel='recovered (%)',
                  title='B  G2 synthetic recovery')
        ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    rp = D[(D.kind == 'repeater') & (D.win == 'P') & np.isfinite(D.dvv)]
    if len(rp):
        ax[2].errorbar(range(len(rp)), 100 * rp.dvv, yerr=100 * rp.err,
                       fmt='o', color='C3')
        ax[2].axhline(0, color='k', lw=1)
        if np.isfinite(sigma_null):
            ax[2].axhspan(-100 * sigma_null, 100 * sigma_null, color='0.85',
                          label='null 1 sigma')
            ax[2].legend(fontsize=8)
        ax[2].set(xlabel='pair', ylabel='dv/v (%)', title='C  per-pair, P window')
        ax[2].grid(alpha=.3)
    # D: the decisive panel. If |dv/v| rides the predicted contamination line, the
    # measurement is reading source geometry rather than the medium.
    for wname, c in [('P', 'C3'), ('coda', 'C0')]:
        s = D[(D.kind == 'repeater') & (D.win == wname)].dropna(
            subset=['dvv', 'contam'])
        if len(s):
            ax[3].scatter(100 * s.contam, 100 * np.abs(s.dvv), color=c,
                          label=f'{wname}  r={g4.get(wname, np.nan):+.2f}')
    lim = ax[3].get_xlim()
    ax[3].plot(lim, lim, 'k--', lw=1, label='|dv/v| = contamination')
    ax[3].set(xlabel='predicted geometric contamination (%)',
              ylabel='|measured dv/v| (%)',
              title='D  G4: is this geometry or medium?')
    ax[3].legend(fontsize=8); ax[3].grid(alpha=.3)
    fig.tight_layout()
    p = os.path.join(HERE, 'interval_dvv_gate.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
