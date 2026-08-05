"""Does SF.MH029 (SAFOD Main Hole, 2555 m) record continuously during the DAS period?

Matters because a same-hole geophone would be an independent confirmation channel
for the repeaters and a calibration reference for DAS amplitude -- but only if it
records continuously rather than on trigger, and only if data actually exists
across the 413 covered days.

Test: request fixed windows on quiet random days (no catalog event nearby). A
continuously recording instrument returns full-length traces; a triggered one
returns nothing.
"""
import numpy as np, pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

for dc in ['NCEDC', 'IRIS']:
    try:
        c = Client(dc, timeout=180)
        inv = c.get_stations(network='SF', station='MH029', level='response',
                             starttime='2024-05-01', endtime='2025-10-01')
    except Exception as e:
        print(f'[{dc}] metadata failed: {str(e)[:70]}'); continue
    print(f'\n[{dc}] SF.MH029 channels')
    for net in inv:
        for sta in net:
            print(f'  lat {sta.latitude:.5f} lon {sta.longitude:.5f} '
                  f'elev {sta.elevation}')
            for ch in sta:
                print(f'   {ch.code:5s} depth {ch.depth:8.1f} m  '
                      f'{ch.sample_rate:7.1f} Hz  az {ch.azimuth} dip {ch.dip}  '
                      f'{str(ch.start_date)[:10]}..{str(ch.end_date)[:10]}')
    break

c = Client('NCEDC', timeout=180)
print('\nCONTINUITY TEST -- 10 min on quiet days, no event targeted')
days = ['2024-06-15','2024-09-20','2024-12-05','2025-02-11','2025-05-01','2025-08-14']
ok = 0
for d in days:
    t0 = UTCDateTime(d + 'T03:00:00')
    try:
        st = c.get_waveforms('SF', 'MH029', '*', '*', t0, t0 + 600)
        if len(st):
            tot = sum(tr.stats.npts / tr.stats.sampling_rate for tr in st)
            print(f'  {d}: {len(st)} trace(s), {tot:7.1f} s of 600 s, '
                  f'{st[0].stats.sampling_rate:g} Hz, ch '
                  f'{sorted({tr.stats.channel for tr in st})}')
            ok += 1
        else:
            print(f'  {d}: empty')
    except Exception as e:
        print(f'  {d}: no data ({str(e)[:55]})')
print(f'\n{ok}/{len(days)} quiet windows returned data '
      f'-> {"CONTINUOUS" if ok >= len(days)-1 else "NOT continuous / gappy"}')

print('\nEVENT TEST -- windows around three confirmed repeaters')
for t in ['2024-07-08T08:30:36','2025-04-06T14:52:18','2024-05-13T09:30:29']:
    t0 = UTCDateTime(t) - 5
    try:
        st = c.get_waveforms('SF','MH029','*','*', t0, t0 + 25)
        print(f'  {t}: {len(st)} trace(s)'
              + (f', {st[0].stats.npts} samples @ {st[0].stats.sampling_rate:g} Hz'
                 if len(st) else ''))
    except Exception as e:
        print(f'  {t}: no data ({str(e)[:55]})')
