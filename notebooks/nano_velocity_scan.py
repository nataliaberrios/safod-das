"""
Slant-stack (semblance) velocity scan of the stacked AWD record on the SAFOD
cemented fiber.

Why this instead of picking: a first-break picker has to be told what to look
for, and my first two attempts locked onto the wrong arrival -- once onto a late
coda 1-3 s after the drop, then onto something moving at ~1100 m/s, which for a
fluid-filled borehole is far more likely a tube wave than direct P. A slant
stack makes no such choice. It sums the section along every trial moveout and
reports which velocities carry coherent energy, so multiple arrivals show up as
separate peaks and we can see what is actually in the record before committing
to an interpretation.

For each trial velocity V and intercept t0, semblance is computed over
    t(z) = t0 + z / V
Peaks in the (V, t0) plane are the arrivals. Direct P in Santa Margarita
sediments should land somewhere around 1800-3500 m/s; a tube wave sits near
1300-1500 m/s and is typically much stronger on borehole DAS.

Outputs a semblance panel plus the record section with the top peaks overlaid,
so the interpretation is visible rather than asserted.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt, hilbert

FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
NPZ = os.path.join(FIG_DIR, 'epoch_stacks_paired.npz')

DX = 1.26606202
PRE_S = 0.5
BANDS = [(10.0, 40.0), (20.0, 50.0)]

Z_MIN_M, Z_MAX_M = 130.0, 530.0
V_GRID = np.arange(400.0, 6000.0, 25.0)
T0_MAX_S = 0.25             # intercept search range
SEMB_WIN_S = 0.040          # semblance window length


def bandpass(x, fs, lo, hi):
    sos = butter(4, [lo, hi], btype='band', fs=fs, output='sos')
    return sosfiltfilt(sos, x, axis=-1)


def semblance(sec, z, fs, i0, vgrid, t0grid, win_s):
    """Classic semblance over (V, t0): coherent energy / total energy."""
    nw = max(1, int(win_s * fs))
    nt = sec.shape[1]
    out = np.zeros((vgrid.size, t0grid.size))
    for iv, v in enumerate(vgrid):
        shifts = (z / v * fs).astype(int)
        for it, t0 in enumerate(t0grid):
            idx = i0 + int(t0 * fs) + shifts
            if idx.min() < 0 or idx.max() + nw >= nt:
                continue
            # gather the window following each trace's predicted arrival
            g = np.empty((z.size, nw), dtype=np.float32)
            for k, i in enumerate(idx):
                g[k] = sec[k, i:i + nw]
            num = np.sum(np.sum(g, axis=0) ** 2)
            den = z.size * np.sum(g ** 2)
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
    print(f'stacked {good.sum()} epochs / {int(w.sum())} drops')

    i0 = int(PRE_S * fs)
    c_lo, c_hi = int(Z_MIN_M / DX), int(Z_MAX_M / DX)
    z = np.arange(c_lo, c_hi + 1) * DX - Z_MIN_M      # relative to window top
    zabs = np.arange(c_lo, c_hi + 1) * DX
    t0grid = np.arange(0.0, T0_MAX_S, 0.002)

    fig, axes = plt.subplots(2, len(BANDS), figsize=(7 * len(BANDS), 10))
    for ib, band in enumerate(BANDS):
        sec = bandpass(raw, fs, *band)[c_lo:c_hi + 1]
        sm = semblance(sec, z, fs, i0, V_GRID, t0grid, SEMB_WIN_S)

        prof = sm.max(axis=1)
        order = np.argsort(prof)[::-1]
        # report distinct peaks, not neighbours of one peak
        peaks, seen = [], []
        for i in order:
            if all(abs(V_GRID[i] - V_GRID[j]) > 300 for j in seen):
                peaks.append(i)
                seen.append(i)
            if len(peaks) == 3:
                break
        print(f'\nband {band[0]:.0f}-{band[1]:.0f} Hz -- top coherent moveouts:')
        for i in peaks:
            it = int(np.argmax(sm[i]))
            print(f'   V = {V_GRID[i]:6.0f} m/s   t0 = {t0grid[it]*1e3:6.1f} ms'
                  f'   semblance {prof[i]:.3f}')

        ax = axes[0, ib] if len(BANDS) > 1 else axes[0]
        im = ax.imshow(sm.T, aspect='auto', origin='lower', cmap='magma',
                       extent=[V_GRID[0], V_GRID[-1], 0, t0grid[-1] * 1e3])
        ax.set_xlabel('trial velocity (m/s)')
        ax.set_ylabel('intercept $t_0$ (ms)')
        ax.set_title(f'Semblance, {band[0]:.0f}-{band[1]:.0f} Hz')
        plt.colorbar(im, ax=ax, label='semblance')
        for i in peaks:
            ax.plot(V_GRID[i], t0grid[int(np.argmax(sm[i]))] * 1e3,
                    'co', ms=7, mfc='none', mew=1.6)

        ax = axes[1, ib] if len(BANDS) > 1 else axes[1]
        clip = np.percentile(np.abs(sec), 99)
        ax.imshow(sec.T, aspect='auto', cmap='gray_r',
                  extent=[zabs[0], zabs[-1], (sec.shape[1] - i0) / fs, -PRE_S],
                  vmin=-clip, vmax=clip)
        for i, col in zip(peaks, ['C0', 'C1', 'C2']):
            it = int(np.argmax(sm[i]))
            ax.plot(zabs, t0grid[it] + z / V_GRID[i], col, lw=1.3,
                    label=f'{V_GRID[i]:.0f} m/s')
        ax.set_ylim(0.6, -0.05)
        ax.set_xlabel('distance along fiber (m)')
        ax.set_ylabel('time after drop (s)')
        ax.set_title(f'Record section, {band[0]:.0f}-{band[1]:.0f} Hz')
        ax.legend(loc='lower right', fontsize=8)

        np.savez(os.path.join(FIG_DIR,
                              f'nano_velocity_scan_{band[0]:.0f}_{band[1]:.0f}.npz'),
                 semblance=sm, v_grid=V_GRID, t0_grid=t0grid,
                 z_abs=zabs, band=np.array(band), n_drops=int(w.sum()))

    fig.suptitle('SAFOD cemented fiber: what moveouts are actually in the '
                 'AWD record?', fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, 'nano_velocity_scan.png')
    fig.savefig(out, dpi=140)
    print('\nSaved', out)


if __name__ == '__main__':
    main()
