"""Redraw the Phase 0 part-3 figure from the saved npz, showing BOTH statistics.

The first version plotted only the localised-dip curves, while the headline
number -- an 18% detectable amplitude loss on the wireline -- comes from the step
statistic.  A figure that omits the curve its caption depends on is a trap for
whoever reads this next.  No recomputation; this only redraws.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIG = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
d = np.load(os.path.join(FIG, 'c2_phase0_diagnose.npz'))
amp = d['amp_grid']
pct = 100 * (1 - np.exp(-amp))

fig, ax = plt.subplots(1, 3, figsize=(17, 5))
for name, col in [('cemented', 'C0'), ('wireline', 'C1')]:
    for i, stat in enumerate(['step', 'dip']):
        ax[i].semilogx(pct, 100 * d[f'{name}__pw_now_{stat}'], col + '-o', ms=4,
                       label=f'{name} as measured')
        ax[i].semilogx(pct, 100 * d[f'{name}__pw_cal_{stat}'], col + '--s', ms=4,
                       mfc='none', label=f'{name} static response removed')
for i, stat in enumerate(['step  (what a permeable fracture makes)',
                          'localised dip  (what the gate looked for)']):
    ax[i].axhline(95, ls=':', c='0.4', lw=1.2)
    ax[i].set(xlabel='localised amplitude loss (%)',
              ylabel='detection rate (%)', ylim=(0, 103), title=stat)
    ax[i].legend(fontsize=8, loc='upper left')
w = float(d['wireline__lim_cal_step'])
ax[0].axvline(100 * (1 - np.exp(-w)), c='C1', lw=1, alpha=0.5)
ax[0].annotate(f'{100*(1-np.exp(-w)):.0f}%', (100 * (1 - np.exp(-w)), 50),
               color='C1', fontsize=9, ha='right')

xs, wd = np.arange(2), 0.35
for i, (name, col) in enumerate([('cemented', 'C0'), ('wireline', 'C1')]):
    ax[2].bar(xs + i * wd, [float(d[f'{name}__sig_static']),
                            float(d[f'{name}__sig_noise'])], wd,
              color=col, alpha=0.85, label=f'{name}  (rho = {float(d[f"{name}__rho"]):.2f})')
ax[2].set_xticks(xs + wd / 2)
ax[2].set_xticklabels(['static per-channel\nresponse (removable)',
                       'random noise\n(irreducible)'], fontsize=9)
ax[2].set(ylabel='log-amplitude sigma', title='What the scatter is made of')
ax[2].legend(fontsize=8)
fig.suptitle('C2 Phase 0 — the amplitude test had no power, and 96% of why is removable',
             fontsize=13)
fig.tight_layout()
out = os.path.join(FIG, 'c2_phase0_diagnose.png')
fig.savefig(out, dpi=140)
print('Saved', out)
