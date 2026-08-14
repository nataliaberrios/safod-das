"""
GATE TEST: is the tube wave worth building a thesis on?

Three claims were made for it. None is verified. This script tests all three and
is meant to be decisive -- if the answers are no, we stay on the current
two-fiber thesis and this costs one job.

  C1  DETECTABLE ON BOTH FIBERS?
      The tube wave travels in the borehole fluid. The wireline hangs in that
      fluid; the cemented fiber is bonded outside the casing. So the wireline
      should see it more strongly. If it does, that is a physical explanation
      for the wireline's poor shallow direct-P correlation (0.11) -- its record
      would be dominated by a different arrival, not badly coupled.
      PASS: semblance >= 0.15 on at least one fiber, and the two differ.

  C2  DOES IT IMAGE PERMEABLE STRUCTURE?
      Tube waves lose energy at permeable fractures; amplitude-vs-depth steps
      down across them. That is the classic hydrophone-VSP fracture indicator.
      PASS: amplitude decay is not a smooth exponential -- there are localised
      drops well above the channel-to-channel scatter.

  C3  DOES IT IMPROVE THE dv/v FLOOR?
      dv/v ~ dt/t. The tube wave is ~2.3x slower than the direct P, so the same
      timing scatter over the same depth gives a proportionally better floor.
      The current floor is 0.287%, 5.7x short of the 0.05% target.
      PASS: tube-wave floor beats the direct-P floor on the same channels.

Everything reads the stored epoch stacks. No raw data access.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt, hilbert

FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
NPZ = os.path.join(FIG_DIR, 'epoch_stacks_paired_deep_all.npz')

DX_CEM, DX_WIRE = 1.26606202, 2.0419
FS, PRE_S = 1000.0, 0.5
# Corrected after wireline_tube_look.py: the tube wave is unmistakable at
# 5-20 Hz across the full 3400 m down-leg and is already fading by 20-50 Hz.
# The original 10-40 Hz band sat on its high-frequency edge.
BAND = (5.0, 20.0)
V_TUBE_GRID = np.arange(900.0, 1900.0, 10.0)
V_P = 2975.0                 # measured direct-P velocity
# And it lives DEEP. The first pass used epoch_stacks_paired.npz, which holds
# only the top 460 wireline channels (~939 m) -- the tube wave's clearest
# expression is beyond that. Use the full-depth file and the whole down-leg.
Z_LO, Z_HI = 130.0, 3400.0
WIN_S = 0.060                # analysis window following each arrival


def bandpass(x, lo=BAND[0], hi=BAND[1], fs=FS):
    return sosfiltfilt(butter(4, [lo, hi], btype='band', fs=fs, output='sos'), x, axis=-1)


def semblance(sec, zrel, i_ref, vgrid, dtgrid, nw):
    out = np.zeros((vgrid.size, dtgrid.size))
    nt = sec.shape[1]
    for iv, v in enumerate(vgrid):
        sh = (zrel / v * FS).astype(int)
        for it, dt in enumerate(dtgrid):
            idx = i_ref + int(dt * FS) + sh
            if idx.min() < 0 or idx.max() + nw >= nt:
                continue
            g = np.stack([sec[k, i:i + nw] for k, i in enumerate(idx)])
            den = zrel.size * np.sum(g ** 2)
            out[iv, it] = np.sum(np.sum(g, axis=0) ** 2) / den if den > 0 else 0.0
    return out


def epoch_lag_scatter(ep, z, v, t0, i0, nw):
    """Burst-to-burst arrival-time scatter along a given moveout, per channel.

    Cross-correlates each epoch's windowed arrival against the all-epoch mean and
    takes the spread of the peak lags -- the same quantity the direct-P
    repeatability products report, so the two are directly comparable.
    """
    n_ep, n_ch = ep.shape[0], ep.shape[1]
    scat = np.full(n_ch, np.nan)
    for c in range(n_ch):
        a = i0 + int((t0 + z[c] / v) * FS)
        if a < 0 or a + nw >= ep.shape[2]:
            continue
        w = ep[:, c, a:a + nw]
        w = w - w.mean(axis=1, keepdims=True)
        ref = w.mean(axis=0)
        if np.all(ref == 0):
            continue
        lags = []
        for e in range(n_ep):
            if np.all(w[e] == 0):
                continue
            cc = np.correlate(w[e], ref, mode='full')
            k = int(np.argmax(cc))
            if 0 < k < cc.size - 1:
                y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]
                den = y0 - 2 * y1 + y2
                k = k + (0.5 * (y0 - y2) / den if den != 0 else 0.0)
            lags.append((k - (nw - 1)) / FS)
        if len(lags) > 4:
            scat[c] = np.std(lags)
    return scat


def main():
    d = np.load(NPZ)
    n_common = d['n_common']
    good = n_common > 0
    w = n_common[good].astype(float)
    i0 = int(PRE_S * FS)
    nw = int(WIN_S * FS)
    print(f'{good.sum()} epochs, {int(w.sum())} paired drops, band {BAND} Hz\n')

    fibers = {}
    for name, key, dx, diff in [('cemented', 'nano_stacks', DX_CEM, False),
                                ('wireline', 'deep_stacks', DX_WIRE, True)]:
        ep = d[key][good].astype(np.float64)
        if diff:                       # OptaSense strain -> strain rate
            ep = np.gradient(ep, 1.0 / FS, axis=-1)
        ep = bandpass(ep)
        mean = np.tensordot(w, ep, axes=(0, 0)) / w.sum()
        fibers[name] = dict(ep=ep, mean=mean, z=np.arange(ep.shape[1]) * dx, dx=dx)

    dtg = np.arange(0.0, 0.60, 0.004)
    print('=== C1  detectable on both fibers? ===')
    for name, F in fibers.items():
        zhi = min(Z_HI, F['z'][-1])
        a, b = int(Z_LO / F['dx']), int(zhi / F['dx'])
        sec, zrel = F['mean'][a:b + 1], F['z'][a:b + 1] - Z_LO
        sm = semblance(sec, zrel, i0, V_TUBE_GRID, dtg, nw)
        iv, it = np.unravel_index(np.argmax(sm), sm.shape)
        F.update(v_tube=V_TUBE_GRID[iv], t0_tube=dtg[it], s_tube=sm[iv, it],
                 sm=sm, a=a, b=b)
        print(f'  {name:9s} V_tube = {F["v_tube"]:6.0f} m/s   '
              f't0 = {F["t0_tube"]*1e3:6.1f} ms   semblance = {F["s_tube"]:.3f}')
    c1 = max(F['s_tube'] for F in fibers.values()) >= 0.15
    print(f'  -> C1 {"PASS" if c1 else "FAIL"}\n')

    print('=== C2  does tube-wave amplitude image structure? ===')
    for name, F in fibers.items():
        z, v, t0 = F['z'], F['v_tube'], F['t0_tube']
        amp = np.full(z.size, np.nan)
        for c in range(F['a'], F['b'] + 1):
            a = i0 + int((t0 + (z[c] - Z_LO) / v) * FS)
            if 0 <= a and a + nw < F['mean'].shape[1]:
                amp[c] = np.sqrt(np.mean(F['mean'][c, a:a + nw] ** 2))
        F['amp'] = amp
        m = np.isfinite(amp)
        la = np.log(amp[m]); zz = z[m]
        trend = np.polyval(np.polyfit(zz, la, 1), zz)
        resid = la - trend
        F['resid'], F['zres'] = resid, zz
        # localised drops = residual excursions well below the scatter
        thr = -2.0 * np.std(resid)
        drops = zz[resid < thr]
        F['drops'] = drops
        print(f'  {name:9s} log-amplitude residual sigma = {np.std(resid):.2f}, '
              f'{drops.size} channels below -2 sigma')
        if drops.size:
            print(f'            candidate depths: '
                  f'{", ".join(f"{x:.0f}" for x in drops[:8])}'
                  f'{" ..." if drops.size > 8 else ""} m')
    c2 = any(F['drops'].size >= 3 for F in fibers.values())
    print(f'  -> C2 {"PASS (candidates present)" if c2 else "FAIL"}\n')

    print('=== C3  does it improve the dv/v floor? ===')
    for name, F in fibers.items():
        z = F['z']
        s_tube = epoch_lag_scatter(F['ep'], z - Z_LO, F['v_tube'], F['t0_tube'], i0, nw)
        s_p = epoch_lag_scatter(F['ep'], z - Z_LO, V_P, 0.004, i0, nw)
        m = (z >= Z_LO) & (z <= min(Z_HI, z[-1]))
        fl_tube = s_tube / (z / F['v_tube'])
        fl_p = s_p / (z / V_P)
        F.update(fl_tube=fl_tube, fl_p=fl_p, mask=m)
        with np.errstate(invalid='ignore'):
            bt, bp = np.nanmedian(fl_tube[m]), np.nanmedian(fl_p[m])
        F['bt'], F['bp'] = bt, bp
        print(f'  {name:9s} median floor {Z_LO:.0f}-{Z_HI:.0f} m:  '
              f'tube {bt*100:6.3f} %   direct-P {bp*100:6.3f} %   '
              f'ratio {bp/bt if bt else np.nan:.2f}x')
    c3 = any(F['bt'] < F['bp'] for F in fibers.values())
    print(f'  -> C3 {"PASS" if c3 else "FAIL"}\n')

    print(f'VERDICT  C1 {"PASS" if c1 else "FAIL"} | C2 {"PASS" if c2 else "FAIL"} '
          f'| C3 {"PASS" if c3 else "FAIL"}')
    print('Pivot the thesis only if at least C1 and C3 pass.')

    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    for name, F, col in [('cemented', fibers['cemented'], 'C0'),
                         ('wireline', fibers['wireline'], 'C1')]:
        ax[0].plot(V_TUBE_GRID, F['sm'].max(axis=1), col, lw=1.8,
                   label=f'{name} — {F["v_tube"]:.0f} m/s ({F["s_tube"]:.2f})')
        ax[1].plot(F['resid'], F['zres'], col, lw=1.2, label=name)
        m = F['mask']
        ax[2].semilogx(F['fl_tube'][m] * 100, F['z'][m], col, lw=1.8, label=f'{name} tube')
        ax[2].semilogx(F['fl_p'][m] * 100, F['z'][m], col, lw=1.0, ls=':',
                       label=f'{name} direct-P')
    ax[0].axhline(0.15, ls='--', c='0.5', lw=1)
    ax[0].set(xlabel='trial velocity (m/s)', ylabel='peak semblance',
              title='C1  Tube wave present?')
    ax[0].legend(fontsize=8)
    ax[1].axvline(0, c='0.5', lw=1)
    ax[1].invert_yaxis()
    ax[1].set(xlabel='log-amplitude residual', ylabel='distance along fiber (m)',
              title='C2  Localised amplitude loss = permeable zones?')
    ax[1].legend(fontsize=8)
    ax[2].invert_yaxis()
    ax[2].set(xlabel='dv/v floor (%)', ylabel='distance along fiber (m)',
              title='C3  Tube wave vs direct P')
    ax[2].legend(fontsize=8)
    fig.suptitle('Tube-wave gate test — does it support a thesis?', fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, 'tube_wave_gate.png')
    fig.savefig(out, dpi=140)
    print('\nSaved', out)


if __name__ == '__main__':
    main()
