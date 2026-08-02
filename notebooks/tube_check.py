"""What is actually at low velocity in the saved slant-stack scans?"""
import glob, os
import numpy as np
FIG='/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
for f in sorted(glob.glob(os.path.join(FIG,'nano_velocity_scan_*.npz'))):
    d=np.load(f); sm=d['semblance']; v=d['v_grid']; t0=d['t0_grid']
    print(f"\n=== {os.path.basename(f)}  band {d['band']} ===")
    prof=sm.max(axis=1)
    for lo,hi,lab in [(800,1800,'tube-wave window 0.8-1.8 km/s'),
                      (2500,3500,'direct-P window 2.5-3.5 km/s')]:
        m=(v>=lo)&(v<=hi)
        i=np.argmax(prof[m])
        vv=v[m][i]; s=prof[m][i]
        it=int(np.argmax(sm[np.where(v==vv)[0][0]]))
        print(f"  {lab:32s} peak V={vv:5.0f} m/s  semblance={s:.3f}  t0={t0[it]*1e3:6.1f} ms")
    # how much of the record's coherent energy sits slow?
    slow=prof[(v>=800)&(v<=1800)].max(); fast=prof[(v>=2500)&(v<=3500)].max()
    print(f"  slow/fast semblance ratio = {slow/fast:.2f}")
