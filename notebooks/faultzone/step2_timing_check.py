"""
Step 2, part 1: is the DAS timing reference what we think it is?

Before any registration can be attempted, one thing has to be checked. The
slant-stack put the cemented-fiber arrival at t0 = 4 ms at 130 m of fiber. The
digitised check shot says a P wave from a surface source reaches 130 m depth at
about 100 ms. Nothing launched at the surface arrives at 130 m in 4 ms -- that
would be 32 km/s. So one of these is true:

  (A) the stacks are not aligned to the drop instant, and t=0 in the stored
      window is offset by an unknown constant;
  (B) the arrival being tracked did not start at the surface at the drop time;
  (C) the channel at "130 m of fiber" is nowhere near 130 m depth.

This matters more than the 7.7% velocity discrepancy, and it is cheap to test:
just look at when energy actually appears, per channel, relative to the stored
drop sample, and compare it against the check-shot prediction.

No claim is made here about which of A/B/C holds. The point is to measure the
offset and see whether it is constant with depth -- a constant offset is a clock
problem, a depth-varying one is not.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt, hilbert

FZ = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone'
AWD = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'

DX_CEM = 1.26606202
FS, PRE_S = 1000.0, 0.5
BAND = (20.0, 50.0)
AWD_OFFSET = 15.0        # AWD source ~10-20 m from the wellhead


def bandpass(x, lo, hi):
    return sosfiltfilt(butter(4, [lo, hi], btype='band', fs=FS, output='sos'),
                       x, axis=-1)


def main():
    ck = np.load(os.path.join(FZ, 'checkshot_traveltime.npz'))
    z_ck, t_ck, off_ck = ck['depth'], ck['traveltime'], float(ck['hoffset'])
    print(f'check shot: {z_ck.min():.0f}-{z_ck.max():.0f} m, HOFFSET {off_ck} m')

    d = np.load(os.path.join(AWD, 'epoch_stacks_paired.npz'))
    n = d['n_common']; good = n > 0; w = n[good].astype(float)
    cem = np.tensordot(w, d['nano_stacks'][good], axes=(0, 0)) / w.sum()
    cem = bandpass(cem, *BAND)
    z = np.arange(cem.shape[0]) * DX_CEM
    i0 = int(PRE_S * FS)
    print(f'cemented stack {cem.shape}, {int(w.sum())} drops, drop sample {i0}')

    # Envelope onset per channel: first time the envelope exceeds a multiple of
    # the pre-drop noise level. Deliberately crude and assumption-free -- this is
    # a timing check, not a pick for velocity work.
    env = np.abs(hilbert(cem, axis=-1))
    noise = np.median(env[:, :i0 - 20], axis=1)
    t_on = np.full(z.size, np.nan)
    for c in range(z.size):
        if not np.isfinite(noise[c]) or noise[c] <= 0:
            continue
        hit = np.where(env[c, i0:] > 5.0 * noise[c])[0]
        if hit.size:
            t_on[c] = hit[0] / FS

    # Check-shot prediction, corrected from its 45.72 m offset to the AWD's ~15 m
    # (both are surface sources; only the slant path differs).
    def predict(zq):
        tq = np.interp(zq, z_ck, t_ck)
        s_ck = np.sqrt(zq ** 2 + off_ck ** 2)
        s_awd = np.sqrt(zq ** 2 + AWD_OFFSET ** 2)
        return tq * s_awd / np.maximum(s_ck, 1e-9)

    sel = (z >= 130) & (z <= 530) & np.isfinite(t_on)
    zz, tt = z[sel], t_on[sel]
    tp = predict(zz)
    resid = tt - tp
    print(f'\n{sel.sum()} channels with an onset in 130-530 m')
    print(f'  observed onset  {tt.min()*1e3:7.1f} - {tt.max()*1e3:7.1f} ms')
    print(f'  check predicts  {tp.min()*1e3:7.1f} - {tp.max()*1e3:7.1f} ms')
    print(f'  residual (obs-pred): median {np.median(resid)*1e3:+7.1f} ms, '
          f'slope {np.polyfit(zz, resid, 1)[0]*1e3:+.4f} ms/m')

    print('\n  depth   observed   predicted   residual')
    for zq in [150, 200, 250, 300, 350, 400, 450, 500]:
        j = np.argmin(np.abs(zz - zq))
        if abs(zz[j] - zq) < 5:
            print(f'  {zq:5d} m {tt[j]*1e3:9.1f} {tp[j]*1e3:11.1f} '
                  f'{resid[j]*1e3:+11.1f} ms')

    # Is the offset constant with depth? Constant -> clock. Sloped -> not.
    sl, ic = np.polyfit(zz, resid, 1)
    spread = np.std(resid - (sl * zz + ic))
    print(f'\nresidual fit: {ic*1e3:+.1f} ms {sl*1e3:+.4f} ms/m, '
          f'scatter about fit {spread*1e3:.1f} ms')
    if abs(sl) * (zz.max() - zz.min()) < 0.010:
        print('=> residual is FLAT with depth: consistent with a constant timing '
              'offset (case A)')
    else:
        print('=> residual VARIES with depth: not a simple clock offset '
              '(case B or C)')

    np.savez(os.path.join(FZ, 'step2_timing_check.npz'),
             depth=zz, t_onset=tt, t_predicted=tp, residual=resid,
             awd_offset=AWD_OFFSET)

    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    ax[0].plot(tt * 1e3, zz, 'C0-', lw=2, label='DAS envelope onset')
    ax[0].plot(tp * 1e3, zz, 'k--', lw=2, label='check shot, offset-corrected')
    ax[0].invert_yaxis()
    ax[0].set(xlabel='time after stored drop sample (ms)',
              ylabel='distance along fiber (m)',
              title='A  Do they agree in absolute time?')
    ax[0].legend(fontsize=9)
    ax[0].grid(alpha=0.3)

    ax[1].plot(resid * 1e3, zz, 'C3-', lw=2)
    ax[1].axvline(0, c='0.4', lw=1)
    ax[1].plot((sl * zz + ic) * 1e3, zz, 'k:', lw=1.5,
               label=f'{ic*1e3:+.0f} ms {sl*1e3:+.3f} ms/m')
    ax[1].invert_yaxis()
    ax[1].set(xlabel='observed - predicted (ms)',
              title='B  Constant offset, or depth-varying?')
    ax[1].legend(fontsize=9)
    ax[1].grid(alpha=0.3)
    fig.suptitle('Is the DAS timing reference the drop instant?', fontsize=12)
    fig.tight_layout()
    out = os.path.join(FZ, 'figures', 'step2_timing_check.png')
    fig.savefig(out, dpi=140)
    print('\nSaved', out)


if __name__ == '__main__':
    main()
