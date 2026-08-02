"""Was there a real tube wave, and is it contaminating the wireline result?

The velocity scans were saved with their full semblance planes, so the slow
arrival can be quantified without recomputing anything. Two questions:
  1. How strong is the low-velocity (fluid-column) energy vs the direct P?
  2. Does the window used for the repeatability metrics contain it?
"""
import os, glob
import numpy as np

FIG = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
TUBE_LO, TUBE_HI = 800.0, 1800.0      # plausible fluid-column / Stoneley range
P_LO = 2200.0                          # direct-P and faster

for f in sorted(glob.glob(os.path.join(FIG, 'nano_velocity_scan_*.npz'))):
    d = np.load(f)
    sm, v, t0 = d['semblance'], d['v_grid'], d['t0_grid']
    band = d['band']
    prof = sm.max(axis=1)

    tm = (v >= TUBE_LO) & (v <= TUBE_HI)
    pm = v >= P_LO
    it_t = int(np.argmax(prof[tm])); it_p = int(np.argmax(prof[pm]))
    v_t, s_t = v[tm][it_t], prof[tm][it_t]
    v_p, s_p = v[pm][it_p], prof[pm][it_p]
    # intercept at each peak
    t0_t = t0[int(np.argmax(sm[np.where(tm)[0][it_t]]))]
    t0_p = t0[int(np.argmax(sm[np.where(pm)[0][it_p]]))]

    print(f'\n=== {os.path.basename(f)}  ({band[0]:.0f}-{band[1]:.0f} Hz) ===')
    print(f'  direct P : V = {v_p:6.0f} m/s  semblance {s_p:.3f}  t0 {t0_p*1e3:5.1f} ms')
    print(f'  slow     : V = {v_t:6.0f} m/s  semblance {s_t:.3f}  t0 {t0_t*1e3:5.1f} ms')
    print(f'  slow/P semblance ratio: {s_t/s_p:.2f}')

    # Does the 0.1 s window that follows the direct P also catch the slow arrival?
    # At depth z the two are separated by z*(1/v_t - 1/v_p).
    print('   separation of the two arrivals vs the 0.10 s metric window:')
    for z in (200, 400, 600, 800):
        sep = z * (1 / v_t - 1 / v_p)
        flag = 'INSIDE window' if sep < 0.10 else 'clear of window'
        print(f'      z={z:4d} m : {sep*1e3:6.1f} ms   {flag}')
