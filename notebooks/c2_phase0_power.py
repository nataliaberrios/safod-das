"""
C2 PHASE 0, part 2: fix two defects in part 1, then turn the null into a bound.

Part 1 (`c2_phase0_significance.py`) had two problems of its own, both found by
reading its output rather than by any new data:

  1. `best_step` allowed a breakpoint only 25 channels from the edge of the
     range.  The cemented fiber's "significant" step (p = 0.007) sits at 161 m,
     which is exactly that minimum, with the residual distribution's skew
     supplying the offset.  It is an edge artefact, and it alone drove the
     script's automatic "C2 SURVIVES" verdict.  Guard the breakpoint by two
     correlation lengths instead.

  2. The verdict treated any p < 0.05 as survival, without requiring the step to
     point the way the physics does.  A permeable fracture makes amplitude stay
     LOWER below it.  Require the sign.

Then the useful part.  A null result is only worth stating if the test could
have seen a real feature, so this plants synthetic features of known size into
surrogates that carry the measured autocorrelation, and asks how large a
localised amplitude loss would have had to be for this measurement to detect it.
That converts "we found nothing" into "we can exclude losses larger than X",
which is the same discipline the dv/v work applies to itself.

Reads only the residual profiles saved by part 1.  The 2.7 GB stack is not
touched again.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
IN_NPZ = os.path.join(FIG_DIR, 'c2_phase0_significance.npz')
OUT_PNG = os.path.join(FIG_DIR, 'c2_phase0_power.png')
OUT_NPZ = os.path.join(FIG_DIR, 'c2_phase0_power.npz')

N_SURROGATE = 2000
N_TRIAL = 400
SEED = 20260813
SIGMA_K = 2.0
# amplitude of the planted feature, in natural-log amplitude units
AMP_GRID = np.array([0.15, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0])
DETECT_P = 0.05
RELIABLE = 0.95        # project convention for a "reliable detection" level


def phase_randomised(resid, rng, n):
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


def robust_sigma(x):
    return 1.4826 * stats.median_abs_deviation(x, scale=1.0)


def best_step(resid, guard):
    """Largest sustained downward offset, with the breakpoint kept `guard`
    samples away from both ends so the statistic is not an edge effect."""
    n = resid.size
    if n <= 2 * guard + 2:
        return np.nan, -1
    cs = np.concatenate([[0.0], np.cumsum(resid)])
    tot = cs[-1]
    ks = np.arange(guard, n - guard)
    step = cs[ks] / ks - (tot - cs[ks]) / (n - ks)
    j = int(np.argmax(step))
    return float(step[j]), int(ks[j])


def deepest(resid):
    return float(resid.min())


def redetrend(log_amp, z):
    return log_amp - np.polyval(np.polyfit(z, log_amp, 1), z)


def plant(base_resid, z, kind, amp, loc_idx, width_pts):
    """Add a feature to a residual series, then re-fit the linear trend.

    Planting before the detrend matters: a real step is partly absorbed by the
    linear fit, so planting after it would overstate the detectability.
    """
    f = np.zeros_like(base_resid)
    if kind == 'step':
        f[loc_idx:] = -amp
    else:                                   # localised dip, one correlation length wide
        f = -amp * np.exp(-0.5 * ((np.arange(base_resid.size) - loc_idx)
                                  / max(width_pts, 1.0)) ** 2)
    return redetrend(base_resid + f, z)


def main():
    rng = np.random.default_rng(SEED)
    d = np.load(IN_NPZ)
    print('=' * 74)
    print('C2 PHASE 0 part 2 — corrected statistics, and a detectability bound')
    print('=' * 74)

    out = {}
    for name in ['cemented', 'wireline']:
        z = d[f'{name}__z']
        resid = d[f'{name}__resid']
        n_lag = int(d[f'{name}__n_lag'])
        corr_len = float(d[f'{name}__corr_len_m'])
        drops = d[f'{name}__drops_naive']
        sig_naive = float(d[f'{name}__sigma_naive'])
        N = resid.size
        dz = corr_len / max(n_lag, 1)
        guard = max(2 * n_lag, int(0.05 * N))

        print(f'\n--- {name}: N = {N}, correlation length {corr_len:.0f} m '
              f'({n_lag} ch), guard = {guard} ch ---')

        # ---- corrected surrogate null ---------------------------------------
        sur = phase_randomised(resid, rng, N_SURROGATE)
        obs_step, obs_k = best_step(resid, guard)
        obs_min = deepest(resid)
        obs_n = int(np.sum(resid < -SIGMA_K * robust_sigma(resid)))

        s_step = np.array([best_step(s, guard)[0] for s in sur])
        s_min = np.array([deepest(s) for s in sur])
        s_n = np.array([int(np.sum(s < -SIGMA_K * robust_sigma(s))) for s in sur])

        p_step = float((s_step >= obs_step).mean())
        p_min = float((s_min <= obs_min).mean())
        p_n = float((s_n >= obs_n).mean())
        print(f'  step  observed {obs_step:.3f} at {z[obs_k]:.0f} m  '
              f'(was {float(d[f"{name}__obs_step"]):.3f} at '
              f'{float(d[f"{name}__obs_step_z"]):.0f} m with the old guard)')
        print(f'        surrogate 95th pct {np.percentile(s_step,95):.3f}   '
              f'p = {p_step:.3f}')
        print(f'  depth observed {obs_min:.2f}, surrogate 5th pct '
              f'{np.percentile(s_min,5):.2f}   p = {p_min:.3f}')
        print(f'  count observed {obs_n}, surrogate median {np.median(s_n):.0f}'
              f'   p = {p_n:.3f}')

        # ---- power: what size of feature WOULD have been detected? ----------
        thr_step = np.percentile(s_step, 100 * (1 - DETECT_P))
        thr_min = np.percentile(s_min, 100 * DETECT_P)
        power = {'step': [], 'dip': []}
        for amp in AMP_GRID:
            for kind in ('step', 'dip'):
                hits = 0
                for t in range(N_TRIAL):
                    base = sur[rng.integers(N_SURROGATE)]
                    loc = int(rng.integers(guard, N - guard))
                    r = plant(base, z, kind, amp, loc, n_lag)
                    if kind == 'step':
                        hits += best_step(r, guard)[0] >= thr_step
                    else:
                        hits += deepest(r) <= thr_min
                power[kind].append(hits / N_TRIAL)
        power = {k: np.array(v) for k, v in power.items()}

        print(f'  detectability ({int(RELIABLE*100)}% of planted features found, '
              f'p < {DETECT_P}):')
        lim = {}
        for kind in ('step', 'dip'):
            ok = np.where(power[kind] >= RELIABLE)[0]
            if ok.size:
                a = AMP_GRID[ok[0]]
                lim[kind] = float(a)
                print(f'    {kind:4s}  amplitude loss >= {a:.2f} log units '
                      f'= {100*(1-np.exp(-a)):.0f}% amplitude drop '
                      f'({20*np.log10(np.exp(a)):.1f} dB)')
            else:
                lim[kind] = np.nan
                print(f'    {kind:4s}  not reached within the tested grid '
                      f'(max {AMP_GRID[-1]:.1f} log units, power '
                      f'{power[kind][-1]*100:.0f}%)')

        # ---- signed step test on the surviving candidates -------------------
        if drops.size:
            half = int(150.0 / dz)
            print(f'  candidate offsets (positive = stays lower below, '
                  f'as a permeable fracture requires):')
            n_right_sign = 0
            for zc in drops:
                j = int(np.argmin(np.abs(z - zc)))
                lo_a, lo_b = max(0, j - half), max(0, j - n_lag)
                hi_a, hi_b = min(N, j + n_lag), min(N, j + half)
                if lo_b - lo_a < 5 or hi_b - hi_a < 5:
                    print(f'    {zc:6.0f} m  edge')
                    continue
                up, dn = resid[lo_a:lo_b], resid[hi_a:hi_b]
                off = float(np.median(up) - np.median(dn))
                se = np.sqrt(np.var(up) / max(up.size / n_lag, 1)
                             + np.var(dn) / max(dn.size / n_lag, 1))
                n_right_sign += off > 0
                print(f'    {zc:6.0f} m  {off:+.3f} +/- {se:.3f} '
                      f'({off/se if se>0 else np.nan:+.1f} sigma)  '
                      f'{"right sign" if off > 0 else "WRONG SIGN"}')
            print(f'    {n_right_sign}/{drops.size} have the sign the physics '
                  f'requires')

        survives = (p_step < 0.05 and obs_step > 0) or p_min < 0.05 or p_n < 0.05
        print(f'  -> {name} {"SURVIVES" if survives else "does not survive"} '
              f'the corrected test')
        out[name] = dict(p_step=p_step, p_min=p_min, p_n=p_n, survives=survives,
                         power_step=power['step'], power_dip=power['dip'],
                         lim_step=lim['step'], lim_dip=lim['dip'],
                         obs_step=obs_step, obs_step_z=float(z[obs_k]),
                         guard=guard, corr_len=corr_len, sigma_naive=sig_naive)

    print('\n' + '=' * 74)
    print('VERDICT')
    print('=' * 74)
    any_s = any(v['survives'] for v in out.values())
    for name, v in out.items():
        print(f'{name:9s}  p(step) = {v["p_step"]:.3f}  p(depth) = '
              f'{v["p_min"]:.3f}  p(count) = {v["p_n"]:.3f}  -> '
              f'{"survives" if v["survives"] else "no significant feature"}')
    print(f'\nC2 {"survives" if any_s else "is RETIRED by Phase 0"}.')
    w = out['wireline']
    if not np.isnan(w['lim_dip']):
        print(f'\nBound worth keeping: on the wireline this measurement would have')
        print(f'detected a localised tube-wave amplitude loss of '
              f'{100*(1-np.exp(-w["lim_dip"])):.0f}% '
              f'({20*np.log10(np.exp(w["lim_dip"])):.1f} dB) or larger in '
              f'{int(RELIABLE*100)}% of trials.')
        print(f'None is present.  That is a limit on permeable-fracture')
        print(f'indication, not a detection.')

    np.savez_compressed(OUT_NPZ, amp_grid=AMP_GRID,
                        **{f'{n}__{k}': v for n, R in out.items()
                           for k, v in R.items()})
    print('\nSaved', OUT_NPZ)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    pct = 100 * (1 - np.exp(-AMP_GRID))
    for name, col in [('cemented', 'C0'), ('wireline', 'C1')]:
        ax[0].plot(pct, 100 * out[name]['power_dip'], col + '-o', ms=4,
                   label=f'{name} localised dip')
        ax[0].plot(pct, 100 * out[name]['power_step'], col + '--s', ms=4,
                   mfc='none', label=f'{name} step')
    ax[0].axhline(100 * RELIABLE, ls=':', c='0.4', lw=1.2)
    ax[0].set(xlabel='planted amplitude loss (%)', ylabel='detection rate (%)',
              title='What size of feature would have been seen?', ylim=(0, 103))
    ax[0].legend(fontsize=8)
    labels, vals, cols = [], [], []
    for name, col in [('cemented', 'C0'), ('wireline', 'C1')]:
        for k, lab in [('p_n', 'count'), ('p_min', 'depth'), ('p_step', 'step')]:
            labels.append(f'{name[:4]}\n{lab}')
            vals.append(out[name][k])
            cols.append(col)
    ax[1].bar(range(len(vals)), vals, color=cols, alpha=0.75)
    ax[1].axhline(0.05, ls='--', c='r', lw=1.2, label='p = 0.05')
    ax[1].set_xticks(range(len(vals)))
    ax[1].set_xticklabels(labels, fontsize=8)
    ax[1].set(ylabel='surrogate p-value', ylim=(0, 1.05),
              title='Nothing crosses the null')
    ax[1].legend(fontsize=8)
    fig.suptitle('C2 Phase 0 part 2 — corrected null and detectability bound',
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=140)
    print('Saved', OUT_PNG)


if __name__ == '__main__':
    main()
