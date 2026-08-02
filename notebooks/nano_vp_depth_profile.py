"""
Interval P-wave velocity vs depth on the SAFOD cemented fiber, from the 2026 AWD
active-source survey.

Builds on nano_velocity_scan.py, which established that the AWD record contains
one dominant coherent arrival at ~2975 m/s over 130-530 m (semblance 0.58 at
20-50 Hz), plus slower guided/tube-wave energy near 1000-1500 m/s. That single
whole-window velocity is an average; this script resolves its depth dependence.

Method: run the same semblance in overlapping depth windows, parameterised so
the interval velocity is what is being solved for rather than trading off
against a global intercept:

    t(z) = t_top + (z - z_top) / V_int

with t_top searched about the prediction from the whole-window velocity. The
peak of semblance over (V_int, t_top) gives the interval velocity at that window.

Scientific point: Lellouch et al. 2019 (10.1029/2019JB017533) published a Vp
profile for this same cemented fiber from passive earthquake arrivals in 2017.
This is an independent active-source profile on the identical sensor, so the two
are directly comparable. Note the depth axis here is distance along fiber with
an uncalibrated absolute offset -- the profile SHAPE is what is comparable.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt

FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
NPZ = os.path.join(FIG_DIR, 'epoch_stacks_paired.npz')

DX = 1.26606202
PRE_S = 0.5
BAND = (20.0, 50.0)         # best direct-P SNR on this fiber

Z_MIN_M, Z_MAX_M = 130.0, 530.0
WIN_M, STEP_M = 100.0, 25.0
V_REF = 2975.0              # whole-window result, used only to centre searches
V_GRID = np.arange(1500.0, 5500.0, 25.0)
TTOP_HALF_S = 0.040
SEMB_WIN_S = 0.040

# For context on the poster: the dv/v product assumes this value.
V_ASSUMED_IN_DVV = 3000.0


def bandpass(x, fs, lo, hi):
    sos = butter(4, [lo, hi], btype='band', fs=fs, output='sos')
    return sosfiltfilt(sos, x, axis=-1)


def window_semblance(sec, zrel, fs, i_top, vgrid, dt_grid, nw):
    """Semblance over (V_int, t_top) for one depth window."""
    nt = sec.shape[1]
    out = np.zeros((vgrid.size, dt_grid.size))
    for iv, v in enumerate(vgrid):
        shifts = (zrel / v * fs).astype(int)
        for it, dt in enumerate(dt_grid):
            idx = i_top + int(dt * fs) + shifts
            if idx.min() < 0 or idx.max() + nw >= nt:
                continue
            g = np.empty((zrel.size, nw), dtype=np.float32)
            for k, i in enumerate(idx):
                g[k] = sec[k, i:i + nw]
            num = np.sum(np.sum(g, axis=0) ** 2)
            den = zrel.size * np.sum(g ** 2)
            out[iv, it] = num / den if den > 0 else 0.0
    return out


def main():
    d = np.load(NPZ)
    nano = d['nano_stacks']
    fs = float(d['fs'])
    n_common = d['n_common']
    good = n_common > 0
    w = n_common[good].astype(float)
    raw = np.tensordot(w, nano[good], axes=(0, 0)) / w.sum()
    print(f'stacked {good.sum()} epochs / {int(w.sum())} drops, band {BAND}')

    sec_all = bandpass(raw, fs, *BAND)
    i0 = int(PRE_S * fs)
    nw = int(SEMB_WIN_S * fs)
    dt_grid = np.arange(-TTOP_HALF_S, TTOP_HALF_S, 0.002)

    centers = np.arange(Z_MIN_M + WIN_M / 2, Z_MAX_M - WIN_M / 2 + 1, STEP_M)
    vp, vz, qual = [], [], []
    print(f'\n{"z_center":>9} {"V_int":>7} {"semblance":>10}')
    for zc in centers:
        z_top, z_bot = zc - WIN_M / 2, zc + WIN_M / 2
        c0, c1 = int(z_top / DX), int(z_bot / DX)
        sec = sec_all[c0:c1 + 1]
        zrel = np.arange(c0, c1 + 1) * DX - z_top
        i_top = i0 + int(z_top / V_REF * fs)

        sm = window_semblance(sec, zrel, fs, i_top, V_GRID, dt_grid, nw)
        iv, it = np.unravel_index(np.argmax(sm), sm.shape)
        v, s = V_GRID[iv], sm[iv, it]
        # a peak pinned to the edge of the search grid is not a measurement
        edge = iv in (0, V_GRID.size - 1)
        print(f'{zc:9.1f} {v:7.0f} {s:10.3f}'
              + ('   <- at grid edge, rejected' if edge else ''))
        if not edge:
            vp.append(v)
            vz.append(zc)
            qual.append(s)

    vp, vz, qual = np.array(vp), np.array(vz), np.array(qual)
    if vp.size == 0:
        raise RuntimeError('no interval velocities resolved')

    print(f'\ninterval Vp over {vz[0]:.0f}-{vz[-1]:.0f} m: '
          f'median {np.median(vp):.0f} m/s, range {vp.min():.0f}-{vp.max():.0f}')
    print(f'whole-window value was {V_REF:.0f} m/s; '
          f'dv/v analysis assumes {V_ASSUMED_IN_DVV:.0f} m/s '
          f'({100*(np.median(vp)-V_ASSUMED_IN_DVV)/V_ASSUMED_IN_DVV:+.1f}% vs measured)')

    np.savez(os.path.join(FIG_DIR, 'nano_vp_depth_profile.npz'),
             vp=vp, depth=vz, semblance_peak=qual, band=np.array(BAND),
             win_m=WIN_M, step_m=STEP_M, v_ref=V_REF, dx=DX,
             n_drops=int(w.sum()), z_min=Z_MIN_M, z_max=Z_MAX_M)

    fig, ax = plt.subplots(1, 2, figsize=(10, 7),
                           gridspec_kw={'width_ratios': [2, 1]})
    ax[0].plot(vp, vz, 'C0-o', lw=1.8, ms=4)
    ax[0].axvline(V_ASSUMED_IN_DVV, ls='--', c='C3', lw=1.2,
                  label=f'assumed in dv/v ({V_ASSUMED_IN_DVV:.0f} m/s)')
    ax[0].axvline(V_REF, ls=':', c='k', lw=1.2,
                  label=f'whole-window ({V_REF:.0f} m/s)')
    ax[0].invert_yaxis()
    ax[0].set_xlabel('interval $V_P$ (m/s)')
    ax[0].set_ylabel('distance along fiber (m), uncalibrated offset')
    ax[0].set_title(f'Active-source $V_P$, SAFOD cemented fiber\n'
                    f'{int(w.sum())} weight drops, {BAND[0]:.0f}-{BAND[1]:.0f} Hz, '
                    f'{WIN_M:.0f} m windows')
    ax[0].grid(alpha=0.3)
    ax[0].legend(fontsize=8)

    ax[1].plot(qual, vz, 'C2-o', lw=1.4, ms=4)
    ax[1].invert_yaxis()
    ax[1].set_xlabel('peak semblance')
    ax[1].set_title('Fit quality')
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, 'nano_vp_depth_profile.png')
    fig.savefig(out, dpi=140)
    print('\nSaved', out)


if __name__ == '__main__':
    main()
