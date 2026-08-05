"""Where exactly does the repeater candidate search fail?

ddrt_candidates.py returned zero sequences at a 150 m diameter. That could mean
the threshold is too tight, the DAS coverage misses one member of every pair, the
magnitude filter is too strict, or there genuinely are no repeaters in the
window. Those have different remedies, so decompose the count rather than tune.
"""
import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

HERE = os.path.dirname(os.path.abspath(__file__))
DDRT = ('/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/'
        'repeating/DDRT_DAS.txt')
SAFOD_LAT, SAFOD_LON = 35.982, -120.544

cols = ['yr','mo','day','hr','mi','sec','lat','lon','depth','mag']
df = pd.read_csv(DDRT, sep=r'\s+', header=None, usecols=range(10), names=cols,
                 skip_blank_lines=True).dropna().reset_index(drop=True)
df['time'] = pd.to_datetime(dict(
    year=df.yr.astype(int), month=df.mo.astype(int), day=df.day.astype(int),
    hour=df.hr.astype(int), minute=df.mi.astype(int),
    second=df.sec.astype(float).round().astype(int).clip(0,59)), utc=True)
klat = 111.19; klon = 111.19*np.cos(np.radians(SAFOD_LAT))
X = np.column_stack([(df.lon-SAFOD_LON)*klon, (df.lat-SAFOD_LAT)*klat, df.depth])*1000.
D = squareform(pdist(X))
days = set(open(os.path.join(HERE,'das_days_all.txt')).read().split())
cov = df.time.dt.strftime('%Y-%m-%d').isin(days).values
dt = np.abs(df.time.values[:,None] - df.time.values[None,:]) / np.timedelta64(1,'D')
dm = np.abs(df.mag.values[:,None] - df.mag.values[None,:])
iu = np.triu_indices(len(df), 1)

print(f'{len(df)} events, {cov.sum()} on covered days\n')
print(f'{"sep_m":>7}{"pairs":>8}{"both_cov":>10}{"+dmag<.3":>10}{"+>30d":>8}'
      f'{"+both_cov,dmag,30d":>20}')
for sep in [50, 100, 150, 250, 500, 1000]:
    m = D[iu] <= sep
    both = m & cov[iu[0]] & cov[iu[1]]
    dmg  = m & (dm[iu] < 0.3)
    t30  = m & (dt[iu] > 30)
    allc = both & (dm[iu] < 0.3) & (dt[iu] > 30)
    print(f'{sep:>7}{m.sum():>8}{both.sum():>10}{dmg.sum():>10}{t30.sum():>8}'
          f'{allc.sum():>20}')

print('\nclosest pairs overall (any filter):')
order = np.argsort(D[iu])[:12]
for k in order:
    i, j = iu[0][k], iu[1][k]
    print(f'  {D[iu][k]:7.0f} m  {df.time[i]:%Y-%m-%d} M{df.mag[i]:.2f} / '
          f'{df.time[j]:%Y-%m-%d} M{df.mag[j]:.2f}  '
          f'dt={dt[i,j]:6.1f} d  dmag={dm[i,j]:.2f}  '
          f'cov={int(cov[i])}{int(cov[j])}')
