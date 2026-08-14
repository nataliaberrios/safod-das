"""
C2 PHASE 0: are the tube-wave amplitude candidates significant at all?

`tube_wave_gate.py` reported C2 PASS on the criterion

    c2 = any(F['drops'].size >= 3 for F in fibers.values())

-- "at least three channels below -2 sigma of a linear log-amplitude trend on
either fiber."  That is not a statistical test.  For Gaussian residuals the
expected count below -2 sigma is 2.3% of channels, and the wireline fit spans
1603 channels, so ~37 are expected by chance and three is guaranteed.  The gate
reported seven.

This script asks the questions the gate did not, in the order that could retire
the claim soonest.  It changes nothing about the measurement -- same band, same
window, same trajectory search, same linear detrend -- so that any difference in
conclusion comes from the statistics and not from reprocessing.

  A  REPRODUCE          re-derive the gate's amplitude profile and candidate list
  B  COUNT vs CHANCE    exact binomial, and effective N from the measured
                        residual autocorrelation rather than an assumed
                        gauge-length correlation span
  C  IS SIGMA HONEST    sigma is estimated from residuals that contain the
                        features being detected; recompute robustly, and
                        decompose the residual variance by wavelength to see
                        whether sigma is set by scatter or by unmodelled trend
  D  SURROGATE NULL     phase-randomised surrogates preserve the residual
                        autocorrelation and destroy localisation, so they give
                        the null distribution of "most extreme deficit" that a
                        Gaussian assumption does not
  E  STEP OR DIP        the physics predicts amplitude STAYS lower below a
                        permeable fracture -- a step.  The gate detected
                        excursions.  Test each candidate for a sustained offset.
  F  DISTINCT FEATURES  collapse candidates separated by less than the measured
                        correlation length
  G  SPLIT-SAMPLE       odd vs even bursts.  A property of the rock appears in
                        both halves; a processing artefact need not.

Reads the stored epoch stacks only.  No raw data access.  Writes a text report,
an npz of the profiles, and a figure.  Modifies nothing that exists.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt
from scipy import stats

FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
NPZ = os.path.join(FIG_DIR, 'epoch_stacks_paired_deep_all.npz')
OUT_NPZ = os.path.join(FIG_DIR, 'c2_phase0_significance.npz')
OUT_PNG = os.path.join(FIG_DIR, 'c2_phase0_significance.png')

# Everything below is copied from tube_wave_gate.py so the measurement is the
# same one.  Do not tune these -- the point is to re-judge the gate's numbers.
DX_CEM, DX_WIRE = 1.26606202, 2.0419
FS, PRE_S = 1000.0, 0.5
BAND = (5.0, 20.0)
V_TUBE_GRID = np.arange(900.0, 1900.0, 10.0)
Z_LO, Z_HI = 130.0, 3400.0
WIN_S = 0.060

N_SURROGATE = 2000
SEED = 20260813
SIGMA_K = 2.0            # the gate's threshold, in sigma
STEP_HALFWIN_M = 150.0   # half-window each side for the step test
LONG_WAVELENGTH_M = 200.0  # boundary for the variance decomposition


def bandpass(x, lo=BAND[0], hi=BAND[1], fs=FS):
    return sosfiltfilt(butter(4, [lo, hi], btype='band', fs=fs, output='sos'),
                       x, axis=-1)


def semblance(sec, zrel, i_ref, vgrid, dtgrid, nw):
    """Same statistic as the gate, vectorised over the window gather."""
    out = np.zeros((vgrid.size, dtgrid.size))
    nt = sec.shape[1]
    n = zrel.size
    rows = np.arange(n)[:, None]
    off = np.arange(nw)[None, :]
    for iv, v in enumerate(vgrid):
        sh = (zrel / v * FS).astype(int)
        for it, dt in enumerate(dtgrid):
            idx = i_ref + int(dt * FS) + sh
            if idx.min() < 0 or idx.max() + nw >= nt:
                continue
            g = sec[rows, idx[:, None] + off]
            den = n * np.sum(g ** 2)
            out[iv, it] = np.sum(np.sum(g, axis=0) ** 2) / den if den > 0 else 0.0
    return out


def amplitude_profile(mean, z, ch0, a, b, v, t0, i0, nw):
    """RMS amplitude in the window following the tube trajectory, per channel."""
    amp = np.full(z.size, np.nan)
    for c in range(a, b + 1):
        k = i0 + int((t0 + (z[c] - Z_LO) / v) * FS)
        if 0 <= k and k + nw < mean.shape[1]:
            amp[c] = np.sqrt(np.mean(mean[c - ch0, k:k + nw] ** 2))
    return amp


def detrend_log(amp, z):
    m = np.isfinite(amp) & (amp > 0)
    la = np.log(amp[m])
    zz = z[m]
    trend = np.polyval(np.polyfit(zz, la, 1), zz)
    return zz, la - trend, la, trend


def correlation_length(resid, dz):
    """Lag at which the residual autocorrelation first falls below 1/e."""
    r = resid - resid.mean()
    ac = np.correlate(r, r, mode='full')[r.size - 1:]
    ac /= ac[0]
    below = np.where(ac < 1.0 / np.e)[0]
    n_lag = int(below[0]) if below.size else r.size
    return max(n_lag, 1), n_lag * dz, ac


def phase_randomised(resid, rng, n):
    """Surrogates with the same power spectrum, hence the same autocorrelation,
    but no localised features.  The correct null for "is this dip special?"."""
    F = np.fft.rfft(resid)
    mag = np.abs(F)
    out = np.empty((n, resid.size))
    for i in range(n):
        ph = rng.uniform(-np.pi, np.pi, mag.size)
        ph[0] = 0.0
        if resid.size % 2 == 0:
            ph[-1] = 0.0
        s = np.fft.irfft(mag * np.exp(1j * ph), n=resid.size)
        out[i] = s - s.mean()
    return out


def deficit_stats(resid, k=SIGMA_K):
    """Statistics a real localised amplitude loss would make extreme."""
    sig = 1.4826 * stats.median_abs_deviation(resid, scale=1.0)
    return dict(n_below=int(np.sum(resid < -k * sig)),
                min_resid=float(resid.min()),
                sigma_robust=float(sig))


def best_step(resid, min_side=25):
    """Largest sustained downward offset from any single breakpoint.

    A permeable fracture removes energy from the tube wave, so amplitude below it
    stays lower.  This is the statistic that shape predicts; the gate's threshold
    crossing is not.
    """
    n = resid.size
    cs = np.concatenate([[0.0], np.cumsum(resid)])
    tot = cs[-1]
    ks = np.arange(min_side, n - min_side)
    above = cs[ks] / ks
    below = (tot - cs[ks]) / (n - ks)
    step = above - below            # positive = drops going down the hole
    j = int(np.argmax(step))
    return float(step[j]), int(ks[j])


def variance_by_wavelength(resid, dz, cut_m=LONG_WAVELENGTH_M):
    """Fraction of residual variance at wavelengths longer than cut_m.

    If most of sigma comes from long wavelengths, the linear detrend is
    inadequate and a threshold in units of sigma is not a scatter threshold.
    """
    F = np.fft.rfft(resid - resid.mean())
    freq = np.fft.rfftfreq(resid.size, d=dz)      # cycles per metre
    p = np.abs(F) ** 2
    long_m = (freq > 0) & (freq < 1.0 / cut_m)
    tot = p[freq > 0].sum()
    return float(p[long_m].sum() / tot) if tot > 0 else np.nan


def cluster(depths, min_gap_m):
    """Collapse candidates closer together than the correlation length."""
    if depths.size == 0:
        return []
    groups, cur = [], [depths[0]]
    for x in depths[1:]:
        if x - cur[-1] <= min_gap_m:
            cur.append(x)
        else:
            groups.append(cur)
            cur = [x]
    groups.append(cur)
    return groups


def main():
    rng = np.random.default_rng(SEED)
    i0 = int(PRE_S * FS)
    nw = int(WIN_S * FS)

    print(__doc__.split('Reads the stored')[0].rstrip())
    print('\n' + '=' * 74)
    print(f'band {BAND} Hz | window {WIN_S*1e3:.0f} ms | '
          f'range {Z_LO:.0f}-{Z_HI:.0f} m | {N_SURROGATE} surrogates | seed {SEED}')
    print('=' * 74 + '\n')

    d = np.load(NPZ)
    n_common = d['n_common']
    good = n_common > 0
    w = n_common[good].astype(float)
    n_ep = int(good.sum())
    print(f'{n_ep} epochs, {int(w.sum())} paired drops\n')

    # odd/even epoch split for stage G, indexed within the good epochs
    halves = {'odd': np.arange(n_ep) % 2 == 1, 'even': np.arange(n_ep) % 2 == 0}

    results = {}
    for name, key, dx in [('cemented', 'nano_stacks', DX_CEM),
                          ('wireline', 'deep_stacks', DX_WIRE)]:
        n_ch_full = d[key].shape[1]
        z = np.arange(n_ch_full) * dx
        zhi = min(Z_HI, z[-1])
        a, b = int(Z_LO / dx), int(zhi / dx)
        # slice channels before widening to float64 -- the file is 2.7 GB
        ep = d[key][good][:, a:b + 1, :].astype(np.float64)
        if key == 'deep_stacks':                # OptaSense strain -> strain rate
            ep = np.gradient(ep, 1.0 / FS, axis=-1)
        ep = bandpass(ep)
        mean = np.tensordot(w, ep, axes=(0, 0)) / w.sum()

        print(f'--- {name}: channels {a}-{b} ({b-a+1}), '
              f'{z[a]:.0f}-{z[b]:.0f} m along fiber ---')

        # ---- A  reproduce the gate ------------------------------------------
        dtg = np.arange(0.0, 0.60, 0.004)
        sm = semblance(mean, z[a:b + 1] - Z_LO, i0, V_TUBE_GRID, dtg, nw)
        iv, it = np.unravel_index(np.argmax(sm), sm.shape)
        v, t0, s_tube = V_TUBE_GRID[iv], dtg[it], sm[iv, it]
        print(f'A  trajectory: V = {v:.0f} m/s, t0 = {t0*1e3:.1f} ms, '
              f'semblance = {s_tube:.3f}')

        amp = amplitude_profile(mean, z, a, a, b, v, t0, i0, nw)
        zz, resid, la, trend = detrend_log(amp, z)
        N = resid.size
        sigma_naive = float(np.std(resid))
        thr_naive = -SIGMA_K * sigma_naive
        drops_naive = zz[resid < thr_naive]
        print(f'A  gate reproduction: sigma = {sigma_naive:.2f}, '
              f'{drops_naive.size} channels below -{SIGMA_K:.0f} sigma')
        if drops_naive.size:
            print('A  candidate depths: '
                  + ', '.join(f'{x:.0f}' for x in drops_naive) + ' m')

        # ---- B  count against chance ----------------------------------------
        p_gauss = stats.norm.cdf(-SIGMA_K)
        exp_gauss = N * p_gauss
        n_lag, corr_len_m, ac = correlation_length(resid, dx)
        n_eff = N / n_lag
        # exact binomial, two-sided in the sense of "is the count anomalous"
        p_lo = stats.binom.cdf(drops_naive.size, N, p_gauss)
        p_hi = stats.binom.sf(drops_naive.size - 1, N, p_gauss)
        print(f'B  N = {N}, expected below -{SIGMA_K:.0f} sigma if Gaussian = '
              f'{exp_gauss:.0f}, observed = {drops_naive.size}')
        print(f'B  binomial P(<= observed) = {p_lo:.3g}, '
              f'P(>= observed) = {p_hi:.3g}')
        print(f'B  residual correlation length = {n_lag} channels '
              f'= {corr_len_m:.0f} m  ->  effective N = {n_eff:.0f}, '
              f'expected = {n_eff*p_gauss:.1f}')

        # ---- C  is sigma honest --------------------------------------------
        st = deficit_stats(resid)
        drops_rob = zz[resid < -SIGMA_K * st['sigma_robust']]
        frac_long = variance_by_wavelength(resid, dx)
        sk, ku = float(stats.skew(resid)), float(stats.kurtosis(resid))
        print(f'C  sigma naive = {sigma_naive:.2f} '
              f'(amplitude factor {np.exp(sigma_naive):.2f}x), '
              f'sigma robust = {st["sigma_robust"]:.2f}')
        print(f'C  robust threshold gives {drops_rob.size} candidates '
              f'(naive gave {drops_naive.size})')
        print(f'C  skewness = {sk:+.2f}, excess kurtosis = {ku:+.2f}  '
              f'(0, 0 = Gaussian)')
        print(f'C  fraction of residual variance at wavelengths '
              f'> {LONG_WAVELENGTH_M:.0f} m = {frac_long*100:.0f}%')

        # ---- D  surrogate null ---------------------------------------------
        obs = deficit_stats(resid)
        obs_step, obs_k = best_step(resid)
        sur = phase_randomised(resid, rng, N_SURROGATE)
        s_nbelow = np.empty(N_SURROGATE, int)
        s_min = np.empty(N_SURROGATE)
        s_step = np.empty(N_SURROGATE)
        for i in range(N_SURROGATE):
            ss = deficit_stats(sur[i])
            s_nbelow[i] = ss['n_below']
            s_min[i] = ss['min_resid']
            s_step[i] = best_step(sur[i])[0]
        p_nbelow = float((s_nbelow >= obs['n_below']).mean())
        p_min = float((s_min <= obs['min_resid']).mean())
        p_step = float((s_step >= obs_step).mean())
        print(f'D  count below -{SIGMA_K:.0f} sigma: observed {obs["n_below"]}, '
              f'surrogate median {np.median(s_nbelow):.0f} '
              f'[{np.percentile(s_nbelow,5):.0f}, '
              f'{np.percentile(s_nbelow,95):.0f}]  p = {p_nbelow:.3f}')
        print(f'D  deepest residual:  observed {obs["min_resid"]:.2f}, '
              f'surrogate 5th pct {np.percentile(s_min,5):.2f}  p = {p_min:.3f}')
        print(f'D  largest sustained step: observed {obs_step:.3f} at '
              f'{zz[obs_k]:.0f} m, surrogate 95th pct '
              f'{np.percentile(s_step,95):.3f}  p = {p_step:.3f}')

        # ---- E  step or dip -------------------------------------------------
        print(f'E  sustained offset across each candidate '
              f'(+/-{STEP_HALFWIN_M:.0f} m, positive = stays lower below):')
        half = int(STEP_HALFWIN_M / dx)
        guard = max(n_lag, 1)
        step_rows = []
        for zc in drops_naive:
            j = int(np.argmin(np.abs(zz - zc)))
            lo_a, lo_b = max(0, j - half), max(0, j - guard)
            hi_a, hi_b = min(N, j + guard), min(N, j + half)
            if lo_b - lo_a < 5 or hi_b - hi_a < 5:
                print(f'E    {zc:6.0f} m  too close to the range edge')
                continue
            up, dn = resid[lo_a:lo_b], resid[hi_a:hi_b]
            offset = float(np.median(up) - np.median(dn))
            # standard error with the correlation length folded in
            se = np.sqrt(np.var(up) / (up.size / guard)
                         + np.var(dn) / (dn.size / guard))
            step_rows.append((zc, offset, se, offset / se if se > 0 else np.nan))
            print(f'E    {zc:6.0f} m  offset = {offset:+.3f} '
                  f'+/- {se:.3f}  ({offset/se if se>0 else np.nan:+.1f} sigma)')

        # ---- F  distinct features ------------------------------------------
        groups = cluster(drops_naive, corr_len_m)
        print(f'F  {drops_naive.size} channels collapse to {len(groups)} '
              f'features at the {corr_len_m:.0f} m correlation length: '
              + ', '.join(f'{np.mean(g):.0f} m x{len(g)}' for g in groups))

        # ---- G  split-sample ------------------------------------------------
        print('G  odd/even burst reproducibility:')
        half_profiles = {}
        for hname, hmask in halves.items():
            hw = w[hmask]
            hmean = np.tensordot(hw, ep[hmask], axes=(0, 0)) / hw.sum()
            hamp = amplitude_profile(hmean, z, a, a, b, v, t0, i0, nw)
            hzz, hres, _, _ = detrend_log(hamp, z)
            hsig = float(np.std(hres))
            hdr = hzz[hres < -SIGMA_K * hsig]
            half_profiles[hname] = dict(z=hzz, resid=hres, sigma=hsig, drops=hdr)
            print(f'G    {hname:4s}: {int(hmask.sum())} epochs, '
                  f'sigma = {hsig:.2f}, {hdr.size} candidates'
                  + ('  at ' + ', '.join(f'{x:.0f}' for x in hdr[:8])
                     if hdr.size else ''))
        # do the two halves agree channel by channel?
        ro, re = half_profiles['odd']['resid'], half_profiles['even']['resid']
        nmin = min(ro.size, re.size)
        rho = float(np.corrcoef(ro[:nmin], re[:nmin])[0, 1])
        # a feature counts as reproduced only if BOTH halves put a candidate
        # within one correlation length of it
        d_odd, d_even = half_profiles['odd']['drops'], half_profiles['even']['drops']
        shared = []
        if d_odd.size and d_even.size:
            allz = np.sort(np.concatenate([d_odd, d_even]))
            for g in cluster(allz, corr_len_m):
                zc = float(np.mean(g))
                if (np.abs(d_odd - zc).min() <= corr_len_m
                        and np.abs(d_even - zc).min() <= corr_len_m):
                    shared.append(g)
        print(f'G    residual profiles correlate rho = {rho:.3f} between halves')
        print(f'G    features present in BOTH halves: {len(shared)}'
              + ('  at ' + ', '.join(f'{np.mean(g):.0f} m' for g in shared)
                 if shared else ''))

        results[name] = dict(
            z=zz, resid=resid, log_amp=la, trend=trend, amp=amp[np.isfinite(amp)],
            v_tube=v, t0_tube=t0, semblance=s_tube, N=N,
            sigma_naive=sigma_naive, sigma_robust=st['sigma_robust'],
            drops_naive=drops_naive, drops_robust=drops_rob,
            exp_gauss=exp_gauss, n_lag=n_lag, corr_len_m=corr_len_m, n_eff=n_eff,
            frac_long=frac_long, skew=sk, kurtosis=ku,
            p_nbelow=p_nbelow, p_min=p_min, p_step=p_step,
            obs_step=obs_step, obs_step_z=float(zz[obs_k]),
            n_features=len(groups), n_shared=len(shared), rho_halves=rho,
            step_rows=np.array(step_rows) if step_rows else np.zeros((0, 4)),
            surr_nbelow=s_nbelow, surr_min=s_min, surr_step=s_step,
            half_odd_resid=ro, half_even_resid=re,
            half_odd_z=half_profiles['odd']['z'],
            half_even_z=half_profiles['even']['z'])
        print()
        del ep, mean

    # ---- verdict -----------------------------------------------------------
    print('=' * 74)
    print('PHASE 0 VERDICT')
    print('=' * 74)
    for name, R in results.items():
        signif = (R['p_nbelow'] < 0.05 or R['p_min'] < 0.05 or R['p_step'] < 0.05)
        print(f'\n{name}:')
        print(f'  candidates                {R["drops_naive"].size} channels '
              f'-> {R["n_features"]} distinct features')
        print(f'  expected by chance        {R["exp_gauss"]:.0f} channels '
              f'(Gaussian, N = {R["N"]})')
        print(f'  surrogate p-values        count {R["p_nbelow"]:.3f} | '
              f'depth {R["p_min"]:.3f} | step {R["p_step"]:.3f}')
        print(f'  reproduced in both halves {R["n_shared"]} features '
              f'(profile rho = {R["rho_halves"]:.3f})')
        print(f'  sigma from long lambda    {R["frac_long"]*100:.0f}% of variance '
              f'> {LONG_WAVELENGTH_M:.0f} m')
        print(f'  -> {"SURVIVES Phase 0" if signif else "does NOT survive Phase 0"}')

    any_survive = any(
        (R['p_nbelow'] < 0.05 or R['p_min'] < 0.05 or R['p_step'] < 0.05)
        for R in results.values())
    print(f'\nC2 {"SURVIVES -- proceed to Phase 1" if any_survive else "RETIRED by Phase 0"}')
    if not any_survive:
        print('No candidate is more extreme than a surrogate with the same')
        print('autocorrelation and no localised feature.  The gate criterion')
        print('("at least three channels below -2 sigma") was satisfied by')
        print('chance, and in fact the observed count is BELOW the Gaussian')
        print('expectation, because sigma is set by unmodelled long-wavelength')
        print('structure rather than by channel-to-channel scatter.')

    np.savez_compressed(OUT_NPZ, **{
        f'{n}__{k}': v for n, R in results.items() for k, v in R.items()})
    print('\nSaved', OUT_NPZ)

    # ---- figure ------------------------------------------------------------
    fig, ax = plt.subplots(1, 4, figsize=(19, 6))
    for name, col in [('cemented', 'C0'), ('wireline', 'C1')]:
        R = results[name]
        ax[0].plot(R['log_amp'], R['z'], col, lw=0.7, alpha=0.6)
        ax[0].plot(R['trend'], R['z'], col, lw=2.0, ls='--')
        ax[1].plot(R['resid'], R['z'], col, lw=1.0, label=name)
        ax[1].axvline(-SIGMA_K * R['sigma_naive'], color=col, ls=':', lw=1.2)
        for zc in R['drops_naive']:
            ax[1].plot(-SIGMA_K * R['sigma_naive'], zc, 'v', color=col, ms=6)
        ax[2].hist(R['surr_nbelow'], bins=np.arange(-0.5, 60.5, 2), alpha=0.45,
                   color=col, label=f'{name} surrogate')
        ax[2].axvline(R['drops_naive'].size, color=col, lw=2.5,
                      label=f'{name} observed = {R["drops_naive"].size}')
        ax[2].axvline(R['exp_gauss'], color=col, ls='--', lw=1.5)
        ax[3].plot(R['half_odd_resid'][:R['half_even_resid'].size],
                   R['half_even_resid'][:R['half_odd_resid'].size],
                   '.', color=col, ms=2, alpha=0.4,
                   label=f'{name}  rho = {R["rho_halves"]:.2f}')
    ax[0].invert_yaxis()
    ax[0].set(xlabel='log RMS amplitude', ylabel='distance along fiber (m)',
              title='Amplitude and the linear detrend')
    ax[1].axvline(0, c='0.5', lw=1)
    ax[1].invert_yaxis()
    ax[1].set(xlabel='log-amplitude residual', ylabel='distance along fiber (m)',
              title=f'Residual, threshold, candidates')
    ax[1].legend(fontsize=8)
    ax[2].set(xlabel=f'channels below -{SIGMA_K:.0f} sigma',
              ylabel='surrogates',
              title='Observed vs null (dashed = Gaussian expectation)')
    ax[2].legend(fontsize=7)
    lim = 3.0
    ax[3].plot([-lim, lim], [-lim, lim], '0.5', lw=1)
    ax[3].set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel='odd-burst residual',
              ylabel='even-burst residual', title='Split-sample reproducibility')
    ax[3].legend(fontsize=8)
    fig.suptitle('C2 Phase 0 — are the tube-wave amplitude candidates significant?',
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=140)
    print('Saved', OUT_PNG)


if __name__ == '__main__':
    main()
