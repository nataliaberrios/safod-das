"""Is there ANY seismometer in the SAFOD hole during the DAS period?

The earlier check (station_geometry.py) restricted to networks BP, NC, PB. This
repeats it with no network restriction and against two FDSN data centres, and also
asks what existed historically -- SAFOD did host downhole instruments (the 2005
PGSI 80-level array is in awd_clean/pgsi_reference/), so "no geophone in the hole"
is only true for a stated time window.
"""
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
LAT, LON = 35.982, -120.544

for dc in ['NCEDC', 'IRIS']:
    for label, t0, t1 in [('DAS period 2024-05..2025-10', '2024-05-01', '2025-10-01'),
                          ('all time 1990..2026',          '1990-01-01', '2026-08-01')]:
        try:
            c = Client(dc, timeout=180)
            inv = c.get_stations(latitude=LAT, longitude=LON, maxradius=0.05,
                                 level='channel', starttime=t0, endtime=t1)
        except Exception as e:
            print(f'[{dc}] {label}: query failed ({str(e)[:60]})'); continue
        rows = []
        for net in inv:
            for sta in net:
                d, _, _ = gps2dist_azimuth(LAT, LON, sta.latitude, sta.longitude)
                dep = sorted({ch.depth for ch in sta})
                rows.append((d, net.code, sta.code, dep,
                             sorted({ch.code[:2] for ch in sta}),
                             str(sta.start_date)[:10], str(sta.end_date)[:10],
                             sta.site.name[:40]))
        rows.sort()
        print(f'\n[{dc}] {label}: {len(rows)} stations within 5.5 km')
        for d, n, s, dep, ch, sd, ed, nm in rows:
            deep = max(dep) if dep else 0
            flag = '  <-- DOWNHOLE' if deep > 50 else ''
            print(f'   {n}.{s:6s} {d/1000:5.2f} km  depths {dep}  {",".join(ch)}'
                  f'  {sd}..{ed}  {nm}{flag}')
