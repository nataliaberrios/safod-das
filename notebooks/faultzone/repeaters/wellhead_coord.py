"""Which SAFOD coordinate is the wellhead? Settle it from survey metadata.

Two are in use, ~1.1 km apart:
  Ellsworth's R script : 35.982,    -120.544     ("distance from MH030")
  ReapetingEvents.ipynb: 35.974204, -120.552141  (= the SF.MH029 station location)

Both feed the 45-degree fault-parallel rotation and every distance, so the choice
propagates into every along-strike coordinate.

The PGSI 2005 survey gives the downhole receiver string in UTM NAD27 with known
well depths. Extrapolating the string to zero depth gives the collar independently
of either guess.
"""
import numpy as np, pandas as pd
try:
    import utm as _utm            # pure python, no PROJ database needed
except ImportError:
    _utm = None
from pyproj import Transformer    # preferred, but its database is often unset here


def utm_to_ll(x, y, zone=10, north=True, epsg=26710):
    """UTM -> lat/lon, pyproj first and the utm package as fallback.

    pyproj honours the NAD27 datum (EPSG:26710); the utm package assumes WGS84.
    The NAD27-WGS84 shift in central California is roughly 100 m, which is
    immaterial here: the two candidate coordinates are 1.1 km apart, so 100 m
    cannot change which one is closer. The datum used is reported either way."""
    try:
        tr = Transformer.from_crs(f'EPSG:{epsg}', 'EPSG:4326', always_xy=True)
        lon, lat = tr.transform(x, y)
        if np.isfinite(lat) and abs(lat) < 90:
            return lat, lon, f'pyproj EPSG:{epsg}'
    except Exception:
        pass
    if _utm is not None:
        lat, lon = _utm.to_latlon(x, y, zone, northern=north)
        return lat, lon, 'utm package (WGS84, ~100 m datum offset)'
    raise RuntimeError('no usable UTM transform')

P = ('/oak/stanford/groups/ettore88/data/SAFOD/Paulsson_array_data/'
     'PGSIarray_rec_coords_pos1.txt')
d = pd.read_csv(P, sep=r'\s+')
print(f'PGSI pos1: {len(d)} levels, WELL_DEP {d.WELL_DEP.min():.2f}'
      f'-{d.WELL_DEP.max():.2f} m')
top = d.iloc[0]
print(f'shallowest: {top.WELL_DEP:.2f} m at UTM E {top.REC_X:.2f} N {top.REC_Y:.2f}'
      f' elev {top.REC_ELEV:.2f}')

CAND = {'Ellsworth 35.982,-120.544 (MH030)': (35.982, -120.544),
        'notebook 35.974204,-120.552141 (MH029)': (35.974204, -120.552141)}

def report(lat, lon, label):
    print(f'\n{label}: lat {lat:.6f}, lon {lon:.6f}')
    for k, (la, lo) in CAND.items():
        dn = (la - lat) * 111.19
        de = (lo - lon) * 111.19 * np.cos(np.radians(36))
        print(f'    {1000*np.hypot(dn, de):8.0f} m from {k}')

lat, lon, how = utm_to_ll(top.REC_X, top.REC_Y)
report(lat, lon, f'shallowest receiver [{how}]')

# extrapolate the top six levels to zero depth -> the collar
x0 = np.polyval(np.polyfit(d.WELL_DEP[:6], d.REC_X[:6], 1), 0.0)
y0 = np.polyval(np.polyfit(d.WELL_DEP[:6], d.REC_Y[:6], 1), 0.0)
lat0, lon0, how0 = utm_to_ll(x0, y0)
report(lat0, lon0, f'EXTRAPOLATED WELL COLLAR [{how0}]')

print('\nVERDICT')
best = min(CAND.items(), key=lambda kv: np.hypot(
    (kv[1][0]-lat0)*111.19, (kv[1][1]-lon0)*111.19*np.cos(np.radians(36))))
print(f'  closest to the surveyed collar: {best[0]}')
print('  Use that for the fault-parallel rotation. Note Ellsworth calls his')
print('  reference "MH030", which is a survey monument, not the SF.MH029 sensor.')
