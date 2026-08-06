"""Did ANY catalog earthquake occur at our sequence locations during the 2017 DAS run?

The 2017 acquisition spans 2017-06-21 to 2017-06-30 -- about ten days (SEG-Y file
names on scratch disks 2 and 3). If no catalog event occurred at a confirmed
sequence location in that window, the 8.4 TB on scratch contains nothing about
those sequences and does not need to be rescued for this project.

Free to run, and decides an 8.4 TB question.
"""
import numpy as np, pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

LAT, LON = 35.982, -120.544
T0, T1 = '2017-06-18', '2017-07-05'          # window + margin

c = Client('NCEDC', timeout=120)
cat = c.get_events(starttime=UTCDateTime(T0), endtime=UTCDateTime(T1),
                   latitude=LAT, longitude=LON, maxradius=0.15,
                   minmagnitude=-1.0)
rows = []
for ev in cat:
    o = ev.preferred_origin() or ev.origins[0]
    m = ev.preferred_magnitude() or (ev.magnitudes[0] if ev.magnitudes else None)
    rows.append(dict(t=pd.Timestamp(o.time.datetime, tz='UTC'),
                     lat=o.latitude, lon=o.longitude,
                     depth=o.depth/1000 if o.depth else np.nan,
                     mag=m.mag if m else np.nan))
E = pd.DataFrame(rows)
print(f'events within ~15 km of SAFOD, {T0} to {T1}: {len(E)}')
if len(E):
    print(E.sort_values('t').to_string(index=False))

# our sequence locations: the 8 HRSN-confirmed pairs + the student's 10
locs = []
h = pd.read_csv('hrsn_extended.csv').dropna(subset=['hrsn'])
h = h[h.hrsn > 0.90]
ev = pd.read_csv('correlate_all_events.csv')
ev['time'] = pd.to_datetime(ev.time, utc=True, format='mixed')
for _, r in h.iterrows():
    for t in (r.t_i, r.t_j):
        m = ev[(ev.time - pd.Timestamp(t)).abs() < pd.Timedelta('60s')]
        if len(m):
            locs.append(('mine', m.iloc[0].lat, m.iloc[0].lon, m.iloc[0].depth))
s = pd.read_csv('/home/groups/ettore88/ettore/research/projects/SAFOD/Catalogs/'
                'RepeatingEq.csv', sep=r'\s+')
for _, r in s.iterrows():
    locs.append(('student', r.latitude, r.longitude, r.depth))

print(f'\nchecking {len(locs)} sequence locations against the 2017 window')
kl, kn = 111.19, 111.19*np.cos(np.radians(36.0))
hits = 0
for src, la, lo, dz in locs:
    if not len(E):
        break
    d = np.sqrt(((E.lat-la)*kl)**2 + ((E.lon-lo)*kn)**2 + (E.depth-dz)**2)
    if (d < 0.5).any():
        hits += 1
        k = d.idxmin()
        print(f'  HIT ({src}): {E.t[k]:%Y-%m-%d %H:%M} M{E.mag[k]:.1f} '
              f'{1000*d[k]:.0f} m from a sequence location')
print(f'\nsequence locations with a 2017 event within 500 m: {hits}/{len(locs)}')
print('\nVERDICT')
if hits == 0:
    print('  No catalog event occurred at any of our sequence locations during the')
    print('  2017 DAS window. The 8.4 TB adds nothing to the repeater analysis.')
    print('  It may still matter for instrument comparison or the tap test -- but')
    print('  the tap test was on disk 1, which has already purged.')
else:
    print(f'  {hits} location(s) had a 2017 event. Worth rescuing that subset.')
