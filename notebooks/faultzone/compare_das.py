"""Check-shot velocity over the exact interval the DAS measured, and what a
mismatch would imply for depth registration."""
import numpy as np
d = np.load('checkshot_traveltime.npz')
z, t, s = d['depth'], d['traveltime'], d['slant']

def vfit(lo, hi):
    m = (z >= lo) & (z <= hi)
    if m.sum() < 6: return np.nan, 0
    A = np.vstack([t[m], np.ones(m.sum())]).T
    return np.linalg.lstsq(A, s[m], rcond=None)[0][0], m.sum()

DAS_V, DAS_LO, DAS_HI = 2975.0, 130.0, 530.0
v, n = vfit(DAS_LO, DAS_HI)
print(f'check shot over {DAS_LO:.0f}-{DAS_HI:.0f} m : {v:.0f} m/s  (n={n})')
print(f'AWD DAS      over {DAS_LO:.0f}-{DAS_HI:.0f} m(fiber): {DAS_V:.0f} m/s')
print(f'ratio check/DAS = {v/DAS_V:.4f}   ({100*(v/DAS_V-1):+.1f} %)')
print()
print('interval Vp profile (200 m windows):')
for lo in range(0, 2400, 200):
    vv, nn = vfit(lo, lo+200)
    if nn: print(f'  {lo:5d}-{lo+200:5d} m : {vv:6.0f} m/s  (n={nn})')
print()
# What depth interval on the check-shot curve has slope == the DAS velocity?
print('where does the check shot run at the DAS velocity?')
for lo in range(0, 2400, 100):
    vv, nn = vfit(lo, lo+400)
    if nn and abs(vv-DAS_V) < 250:
        print(f'  {lo:5d}-{lo+400:5d} m : {vv:6.0f} m/s')
