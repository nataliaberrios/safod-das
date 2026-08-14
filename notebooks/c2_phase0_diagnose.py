"""
C2 PHASE 0, part 3: why the test has no power, and what would give it some.

Part 2 produced two results that need explaining rather than reporting:

  1. The cemented "significant step" (p = 0.007) moved from 161 m to 168 m when
     the edge guard was widened -- i.e. it stayed pinned to whatever the first
     allowed breakpoint was.  A statistic that always lands on the boundary is
     measuring the boundary.  The likely cause is that log-amplitude is not
     linear in depth near the source, where geometric spreading dominates.
     Test: drop the shallowest 300 m and see whether the step survives.

  2. Neither fiber reached 95% detection of a planted feature anywhere on the
     grid -- not even a 95% amplitude loss (3.0 log units, 26 dB).  So the
     measurement cannot detect a permeable fracture of any plausible size, and
     "no significant candidates" was never going to mean anything.

The reason for (2) is that the channel-to-channel log-amplitude scatter is ~1.2,
a factor of 3.3 in amplitude, which is far larger than the effect being sought.
But part 1 also found the odd- and even-burst residual profiles correlate at
rho = 0.922 on the wireline.  Reproducible scatter is not noise -- it is static
per-channel response (coupling, sensitivity, gauge-length position).  That can be
calibrated out; random noise cannot.

So this script splits sigma into its static and random parts using the two burst
halves, and re-runs the power analysis at the random-only level to show what a
calibrated amplitude measurement could achieve.  That is the number that decides
whether Phase 1 is worth doing.

  o = s + n_o,  e = s + n_e   (s static, n independent between halves)
  d = (o - e)/2               isolates noise; var(d) = sigma_n^2 / 2
  full 46-epoch stack has noise variance sigma_n^2 / 2 = var(d)

so the difference series gives both the noise level of the full stack and its
autocorrelation, with no model assumed.

Reads only the profiles saved by part 1.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
IN_NPZ = os.path.join(FIG_DIR, 'c2_phase0_significance.npz')
OUT_PNG = os.path.join(FIG_DIR, 'c2_phase0_diagnose.png')
OUT_NPZ = os.path.join(FIG_DIR, 'c2_phase0_diagnose.npz')

N_SURROGATE = 1500
N_TRIAL = 400
SEED = 20260813
SIGMA_K = 2.0
AMP_GRID = np.array([0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
DETECT_P = 0.05
RELIABLE = 0.95
SHALLOW_CUT_M = 300.0


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


def corr_len(resid):
    r = resid - resid.mean()
    ac = np.correlate(r, r, mode='full')[r.size - 1:]
    ac = ac / ac[0]
    below = np.where(ac < 1.0 / np.e)[0]
    return max(int(below[0]) if below.size else r.size, 1)


def best_step(resid, guard):
    n = resid.size
    if n <= 2 * guard + 2:
        return np.nan, -1
    cs = np.concatenate([[0.0], np.cumsum(resid)])
    tot = cs[-1]
    ks = np.arange(guard, n - guard)
    step = cs[ks] / ks - (tot - cs[ks]) / (n - ks)
    j = int(np.argmax(step))
    return float(step[j]), int(ks[j])


def redetrend(x, z):
    return x - np.polyval(np.polyfit(z, x, 1), z)


def plant(base, z, kind, amp, loc, width):
    f = np.zeros_like(base)
    if kind == 'step':
        f[loc:] = -amp
    else:
        f = -amp * np.exp(-0.5 * ((np.arange(base.size) - loc)
                                  / max(width, 1.0)) ** 2)
    return redetrend(base + f, z)


def power_curve(series, z, rng, guard, width, label):
    """Detection rate vs planted amplitude, for surrogates carrying `series`'
    autocorrelation and level."""
    sur = phase_randomised(series, rng, N_SURROGATE)
    s_step = np.array([best_step(s, guard)[0] for s in sur])
    s_min = np.array([s.min() for s in sur])
    thr_step = np.percentile(s_step, 100 * (1 - DETECT_P))
    thr_min = np.percentile(s_min, 100 * DETECT_P)
    pw = {'step': [], 'dip': []}
    for amp in AMP_GRID:
        for kind in ('step', 'dip'):
            hits = 0
            for _ in range(N_TRIAL):
                base = sur[rng.integers(N_SURROGATE)]
                loc = int(rng.integers(guard, series.size - guard))
                r = plant(base, z, kind, amp, loc, width)
                hits += (best_step(r, guard)[0] >= thr_step if kind == 'step'
                         else r.min() <= thr_min)
            pw[kind].append(hits / N_TRIAL)
    pw = {k: np.array(v) for k, v in pw.items()}
    print(f'  {label}: sigma = {np.std(series):.3f} '
          f'(amplitude factor {np.exp(np.std(series)):.2f}x)')
    lim = {}
    for kind in ('step', 'dip'):
        ok = np.where(pw[kind] >= RELIABLE)[0]
        a = AMP_GRID[ok[0]] if ok.size else np.nan
        lim[kind] = float(a) if ok.size else np.nan
        if ok.size:
            print(f'    {kind:4s} {int(RELIABLE*100)}% detection at '
                  f'{a:.2f} log units = {100*(1-np.exp(-a)):.0f}% amplitude '
                  f'loss ({20*np.log10(np.exp(a)):.1f} dB)')
        else:
            print(f'    {kind:4s} never reaches {int(RELIABLE*100)}% '
                  f'(max power {pw[kind][-1]*100:.0f}% at '
                  f'{100*(1-np.exp(-AMP_GRID[-1])):.0f}% loss)')
    return pw, lim


def main():
    rng = np.random.default_rng(SEED)
    d = np.load(IN_NPZ)
    print('=' * 74)
    print('C2 PHASE 0 part 3 — diagnosis: what limits this measurement')
    print('=' * 74)

    out = {}
    for name in ['cemented', 'wireline']:
        z = d[f'{name}__z']
        resid = d[f'{name}__resid']
        n_lag = int(d[f'{name}__n_lag'])
        N = resid.size
        dz = float(d[f'{name}__corr_len_m']) / max(n_lag, 1)
        print(f'\n{"="*70}\n{name}\n{"="*70}')

        # ---- 1  is the cemented step just near-source curvature? ------------
        print(f'1  step statistic vs where the range is allowed to start')
        for cut in (0.0, 150.0, SHALLOW_CUT_M, 600.0):
            m = z >= z[0] + cut
            if m.sum() < 200:
                continue
            zz, rr = z[m], redetrend(d[f'{name}__log_amp'][m], z[m])
            g = max(2 * corr_len(rr), int(0.10 * rr.size))
            st, k = best_step(rr, g)
            sur = phase_randomised(rr, rng, 600)
            p = float((np.array([best_step(s, g)[0] for s in sur]) >= st).mean())
            at_edge = (k <= g + 2)
            print(f'   drop shallowest {cut:5.0f} m -> step {st:.3f} at '
                  f'{zz[k]:6.0f} m, p = {p:.3f}'
                  f'{"   <-- pinned to the first allowed breakpoint" if at_edge else ""}')

        # ---- 2  split sigma into static and random ---------------------------
        ro, re = d[f'{name}__half_odd_resid'], d[f'{name}__half_even_resid']
        n = min(ro.size, re.size)
        ro, re, zc = ro[:n], re[:n], z[:n]
        rho = float(np.corrcoef(ro, re)[0, 1])
        diff = 0.5 * (ro - re)                 # noise only; static cancels
        sig_total = float(np.std(resid))
        sig_noise = float(np.std(diff))         # = noise sigma of the full stack
        sig_static = float(np.sqrt(max(sig_total ** 2 - sig_noise ** 2, 0.0)))
        print(f'\n2  variance decomposition from the two burst halves')
        print(f'   half-to-half correlation rho          {rho:+.3f}')
        print(f'   total residual sigma                  {sig_total:.3f} '
              f'(factor {np.exp(sig_total):.2f}x in amplitude)')
        print(f'   STATIC per-channel response sigma     {sig_static:.3f} '
              f'(factor {np.exp(sig_static):.2f}x)  <- calibratable')
        print(f'   RANDOM measurement noise sigma        {sig_noise:.3f} '
              f'(factor {np.exp(sig_noise):.2f}x)  <- irreducible here')
        print(f'   static fraction of variance           '
              f'{100*sig_static**2/sig_total**2:.0f}%')
        print(f'   potential improvement in threshold    '
              f'{sig_total/sig_noise:.1f}x')

        # ---- 3  power now, and power if the static part were removed --------
        guard = max(2 * n_lag, int(0.10 * N))
        print(f'\n3  detectability (guard {guard} ch, dip width {n_lag} ch)')
        pw_now, lim_now = power_curve(resid, z, rng, guard, n_lag,
                                     'AS MEASURED   ')
        gn = max(2 * corr_len(diff), int(0.10 * diff.size))
        pw_cal, lim_cal = power_curve(diff, zc, rng, gn, max(corr_len(diff), 1),
                                     'NOISE-LIMITED ')

        out[name] = dict(rho=rho, sig_total=sig_total, sig_static=sig_static,
                         sig_noise=sig_noise, pw_now_dip=pw_now['dip'],
                         pw_now_step=pw_now['step'], pw_cal_dip=pw_cal['dip'],
                         pw_cal_step=pw_cal['step'],
                         lim_now_dip=lim_now['dip'], lim_cal_dip=lim_cal['dip'],
                         lim_now_step=lim_now['step'], lim_cal_step=lim_cal['step'])

    print('\n' + '=' * 74)
    print('WHAT THIS MEANS FOR C2')
    print('=' * 74)
    for name, v in out.items():
        print(f'\n{name}:')
        print(f'  as measured, {int(RELIABLE*100)}% detection needs '
              + (f'{100*(1-np.exp(-v["lim_now_dip"])):.0f}% amplitude loss'
                 if not np.isnan(v['lim_now_dip'])
                 else 'more than a 95% amplitude loss — i.e. it cannot detect one'))
        print(f'  with static channel response removed, it would need '
              + (f'{100*(1-np.exp(-v["lim_cal_dip"])):.0f}%'
                 if not np.isnan(v['lim_cal_dip']) else '>95%'))
        print(f'  {100*v["sig_static"]**2/v["sig_total"]**2:.0f}% of the scatter '
              f'is static and therefore removable '
              f'({v["sig_total"]/v["sig_noise"]:.1f}x headroom)')

    np.savez_compressed(OUT_NPZ, amp_grid=AMP_GRID,
                        **{f'{n}__{k}': v for n, R in out.items()
                           for k, v in R.items()})
    print('\nSaved', OUT_NPZ)

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    pct = 100 * (1 - np.exp(-AMP_GRID))
    for name, col in [('cemented', 'C0'), ('wireline', 'C1')]:
        v = out[name]
        ax[0].semilogx(pct, 100 * v['pw_now_dip'], col + '-o', ms=4,
                       label=f'{name} as measured')
        ax[0].semilogx(pct, 100 * v['pw_cal_dip'], col + '--s', ms=4, mfc='none',
                       label=f'{name} noise-limited')
    ax[0].axhline(100 * RELIABLE, ls=':', c='0.4', lw=1.2)
    ax[0].set(xlabel='localised amplitude loss (%)', ylabel='detection rate (%)',
              ylim=(0, 103),
              title='Detectability now vs with channel response removed')
    ax[0].legend(fontsize=8)
    xs = np.arange(2)
    wd = 0.35
    for i, (name, col) in enumerate([('cemented', 'C0'), ('wireline', 'C1')]):
        v = out[name]
        ax[1].bar(xs + i * wd, [v['sig_static'], v['sig_noise']], wd,
                  color=col, alpha=0.8, label=name)
    ax[1].set_xticks(xs + wd / 2)
    ax[1].set_xticklabels(['static per-channel\nresponse (removable)',
                           'random noise\n(irreducible)'], fontsize=9)
    ax[1].set(ylabel='log-amplitude sigma',
              title='What the scatter is made of')
    ax[1].legend(fontsize=8)
    fig.suptitle('C2 Phase 0 part 3 — why the amplitude test has no power',
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=140)
    print('Saved', OUT_PNG)


if __name__ == '__main__':
    main()
