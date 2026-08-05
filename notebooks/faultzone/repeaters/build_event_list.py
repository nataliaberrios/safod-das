"""
Stage 1: every DDRT event that falls on a day the cemented fiber was recording.

This is the input to a similarity-first repeater search. The previous approach
was backwards: it used catalog LOCATION as the prior for which events might
repeat, then checked waveforms. Three of four co-located, magnitude-matched pairs
came back at CC ~ 0, because location at +/-150-250 m is simply a poor predictor
of whether two events rupture the same patch.

The field does it the other way round (REDPy, Waldhauser-Schaff): correlate
everything and let waveform similarity define the sequences. Location then
becomes a check on the result rather than a filter on the input.

No location or magnitude filtering here on purpose -- the whole point is to let
similarity decide.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DDRT = ('/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/'
        'repeating/DDRT_DAS.txt')
DAYS = os.path.join(HERE, 'das_days_all.txt')
SAFOD_LAT, SAFOD_LON = 35.982, -120.544

cols = ['yr', 'mo', 'day', 'hr', 'mi', 'sec', 'lat', 'lon', 'depth', 'mag']
df = pd.read_csv(DDRT, sep=r'\s+', header=None, usecols=range(10), names=cols,
                 skip_blank_lines=True).dropna().reset_index(drop=True)
df['time'] = pd.to_datetime(dict(
    year=df.yr.astype(int), month=df.mo.astype(int), day=df.day.astype(int),
    hour=df.hr.astype(int), minute=df.mi.astype(int),
    second=df.sec.astype(float).round().astype(int).clip(0, 59)), utc=True)

klat = 111.19
klon = 111.19 * np.cos(np.radians(SAFOD_LAT))
n = (df.lat - SAFOD_LAT) * klat
e = (df.lon - SAFOD_LON) * klon
s = np.radians(45.0)
df['r'] = n * np.cos(s) + e * np.sin(s)
df['t'] = -n * np.sin(s) + e * np.cos(s)
df['km_safod'] = np.sqrt(n ** 2 + e ** 2)

days = set(open(DAYS).read().split())
cov = df[df.time.dt.strftime('%Y-%m-%d').isin(days)].copy()
cov = cov.sort_values('time').reset_index(drop=True)
cov['idx'] = np.arange(len(cov))
cov['tag'] = cov.time.dt.strftime('ev_%Y%m%dT%H%M%S')

out = os.path.join(HERE, 'all_events.csv')
cov[['idx', 'tag', 'time', 'mag', 'depth', 'lat', 'lon', 't', 'r',
     'km_safod']].to_csv(out, index=False)

print(f'{len(df)} DDRT events; {len(cov)} on covered days')
print(f'  span {cov.time.min():%Y-%m-%d} to {cov.time.max():%Y-%m-%d}')
print(f'  mag {cov.mag.min():.2f} to {cov.mag.max():.2f}')
print(f'  depth {cov.depth.min():.2f} to {cov.depth.max():.2f} km')
print(f'  pairs to correlate: {len(cov)*(len(cov)-1)//2}')
print(f'wrote {out}')
