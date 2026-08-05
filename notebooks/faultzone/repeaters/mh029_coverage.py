"""How many of the confirmed-repeater events does SF.MH029 actually have?

MH029 is not continuously archived, so its value is decided by coverage of the
specific events, not by uptime. Every event it recorded is an independent
same-hole confirmation at 1000 Hz and three components.
"""
import pandas as pd, numpy as np
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

c = Client('NCEDC', timeout=180)
d = pd.read_csv('hrsn_extended.csv')
h = d.dropna(subset=['hrsn'])
g = h[h.hrsn > 0.90].copy()
ev = sorted(set(pd.to_datetime(g.t_i, utc=True, format='mixed')) |
            set(pd.to_datetime(g.t_j, utc=True, format='mixed')))
print(f'{len(ev)} events in pairs above HRSN CC 0.90\n')
got = {}
for t in ev:
    t0 = UTCDateTime(t.strftime('%Y-%m-%dT%H:%M:%S')) - 5
    try:
        st = c.get_waveforms('SF','MH029','*','*', t0, t0+25)
        n = len(st)
        s = f'{n} ch, {st[0].stats.npts} samp' if n else 'EMPTY'
    except Exception as e:
        n, s = 0, 'no data'
    got[str(t)] = n > 0
    print(f'  {t:%Y-%m-%d %H:%M:%S}  {"YES" if n else "no ":3s}  {s}')
print(f'\nMH029 coverage: {sum(got.values())}/{len(ev)} confirmed-repeater events')
pairs = 0
for _, r in g.iterrows():
    ti = pd.to_datetime(r.t_i, utc=True, format='mixed')
    tj = pd.to_datetime(r.t_j, utc=True, format='mixed')
    if got.get(str(ti), False) and got.get(str(tj), False):
        pairs += 1
print(f'PAIRS with BOTH events on MH029: {pairs}/{len(g)}')
print('\n(a pair needs both events to be usable as an independent confirmation)')
