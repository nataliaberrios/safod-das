"""
Is any HRSN geophone co-located with the SAFOD fiber?

This decides whether 'cross-correlate the geophone against each DAS channel' is a
true INSTRUMENT comparison or just a measurement of path differences. HRSN
stations are borehole instruments, but if they sit kilometres away then a DAS
channel and a geophone see the same earthquake through different structure, and
their CC is dominated by path, not by the sensors. Co-located and depth-matched is
a completely different, and much more useful, experiment.

Also reports sensor emplacement depth, since the fiber spans 0-919 m: a geophone
inside that depth range can be compared channel-to-channel.
"""
import numpy as np
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

SAFOD_LAT, SAFOD_LON = 35.982, -120.544

c = Client('NCEDC', timeout=120)
inv = c.get_stations(network='BP,NC,PB', latitude=SAFOD_LAT, longitude=SAFOD_LON,
                     maxradius=0.25, level='channel',
                     starttime='2024-05-01', endtime='2025-10-01')
rows = []
for net in inv:
    for sta in net:
        d, _, _ = gps2dist_azimuth(SAFOD_LAT, SAFOD_LON, sta.latitude, sta.longitude)
        dep = sorted({ch.depth for ch in sta})
        cods = sorted({ch.code for ch in sta})
        rows.append((d, net.code, sta.code, sta.latitude, sta.longitude,
                     sta.elevation, dep, cods, sta.site.name))
rows.sort()
print(f'{"net":>4}{"sta":>7}{"dist_km":>9}{"elev":>7}{"sensor_depth_m":>16}   '
      f'channels / site')
for d, n, s, la, lo, el, dep, cods, name in rows:
    ds = ','.join(f'{x:.0f}' for x in dep)
    print(f'{n:>4}{s:>7}{d/1000:9.2f}{el:7.0f}{ds:>16}   '
          f'{",".join(cods[:6])} | {name[:34]}')

print('\nclosest station:')
d, n, s, la, lo, el, dep, cods, name = rows[0]
print(f'  {n}.{s} at {d/1000:.2f} km, sensor depth {dep}, "{name}"')
print(f'\nfiber spans 0-919 m depth at SAFOD. A geophone is directly comparable to')
print(f'a DAS channel only if it is BOTH within a few hundred metres laterally AND')
print(f'inside that depth range.')
