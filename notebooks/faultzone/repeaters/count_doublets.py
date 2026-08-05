"""
The yes/no: are there repeating-earthquake DOUBLETS inside the DAS coverage?

Coda-wave interferometry needs TWO events from the same repeating sequence
recorded on the SAME instrument. One event matched to a known sequence is not
enough. So the question that decides the whole repeater pivot is simply: how
many candidate doublets fall on days the fiber was recording?

Important limitation, stated up front because it bounds everything below: true
repeater identification requires waveform cross-correlation (CC > 0.9 or so).
This script uses CATALOG COLOCATION as a proxy, and catalog absolute locations
at Parkfield carry ~0.2-0.5 km horizontal error while real repeaters are <20 m
across. So every number here is an UPPER BOUND on candidate pairs, not a
repeater count. What it can settle is the yes/no: whether there is enough
co-located, similar-magnitude seismicity on covered days to be worth running
waveform CC on at all.

Context from the literature review: Schaff & Waldhauser report >40% of Parkfield
seismicity at CC>0.9 is repeating, and >50% of events in the region sit in
repeating clusters -- so a substantial fraction of these candidates should
survive waveform screening.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CAT = os.path.join(HERE, 'parkfield_catalog_2024_2025.csv')
DAYS = os.path.join(HERE, 'das_days_cemented.txt')

# Repeaters are <20 m across; catalog error is ~0.2-0.5 km. Cluster at a scale
# set by the catalog, not by the physics, and say so.
R_KM = 0.5          # horizontal + vertical separation to call a pair co-located
DM = 0.4            # magnitude similarity (repeaters have near-identical size)

# SAFOD target repeating sequences -- the patch the borehole was sited on
SAFOD_LAT, SAFOD_LON, SAFOD_Z = 35.974, -120.552, 2.7


def main():
    cat = pd.read_csv(CAT)
    cat['time'] = pd.to_datetime(cat['time'], format='mixed', utc=True)
    cat['day'] = cat['time'].dt.strftime('%Y-%m-%d')
    days = set(open(DAYS).read().split())
    print(f'catalog: {len(cat)} events, {cat.time.min():%Y-%m-%d} to '
          f'{cat.time.max():%Y-%m-%d}')
    print(f'DAS coverage: {len(days)} days\n')

    cat['covered'] = cat['day'].isin(days)
    cov = cat[cat['covered']].copy().reset_index(drop=True)
    print(f'events on covered days: {len(cov)} / {len(cat)} '
          f'({100*len(cov)/len(cat):.0f}%)')
    print(f'  magnitude range {cov.mag.min():.2f} to {cov.mag.max():.2f}')
    print(f'  depth range {cov.depth.min():.2f} to {cov.depth.max():.2f} km\n')

    # pairwise co-location among covered events
    lat = cov.latitude.values; lon = cov.longitude.values
    dep = cov.depth.values; mag = cov.mag.values
    klat = 111.19
    klon = 111.19 * np.cos(np.radians(np.nanmean(lat)))

    pairs = []
    for i in range(len(cov)):
        for j in range(i + 1, len(cov)):
            dx = (lon[i] - lon[j]) * klon
            dy = (lat[i] - lat[j]) * klat
            dz = dep[i] - dep[j]
            d = np.sqrt(dx * dx + dy * dy + dz * dz)
            if d <= R_KM and abs(mag[i] - mag[j]) <= DM:
                pairs.append((i, j, d))
    print(f'candidate co-located pairs on covered days: {len(pairs)}')

    # collapse pairs into clusters (a sequence may repeat more than twice)
    parent = list(range(len(cov)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i, j, _ in pairs:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)

    groups = {}
    for i in range(len(cov)):
        groups.setdefault(find(i), []).append(i)
    seqs = [g for g in groups.values() if len(g) >= 2]
    seqs.sort(key=lambda g: -len(g))
    print(f'candidate sequences with >=2 covered events: {len(seqs)}')
    print(f'total events in those sequences: {sum(len(g) for g in seqs)}\n')

    dsaf = np.sqrt(((lon - SAFOD_LON) * klon) ** 2 +
                   ((lat - SAFOD_LAT) * klat) ** 2)
    print(f'{"n":>3} {"depth km":>9} {"mag":>6} {"km from SAFOD":>14}  span')
    for g in seqs[:20]:
        t = cov.time.iloc[g].sort_values()
        gap = (t.iloc[-1] - t.iloc[0]).days
        print(f'{len(g):3d} {np.mean(dep[g]):9.2f} {np.mean(mag[g]):6.2f} '
              f'{np.mean(dsaf[g]):14.2f}  {t.iloc[0]:%Y-%m-%d} -> '
              f'{t.iloc[-1]:%Y-%m-%d} ({gap} d)')

    near = [g for g in seqs if np.mean(dsaf[g]) < 5.0]
    print(f'\ncandidate sequences within 5 km of SAFOD: {len(near)}')
    deep = [g for g in seqs if 1.0 < np.mean(dep[g]) < 12.0]
    print(f'candidate sequences at seismogenic depth (1-12 km): {len(deep)}')

    out = os.path.join(HERE, 'doublet_candidates.csv')
    rows = []
    for k, g in enumerate(seqs):
        for i in g:
            rows.append(dict(seq=k, time=cov.time.iloc[i], mag=mag[i],
                             depth=dep[i], lat=lat[i], lon=lon[i],
                             km_from_safod=dsaf[i], id=cov.id.iloc[i]))
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f'\nwrote {out}')
    print('\nNOTE: these are catalog-colocation candidates, an upper bound. '
          'Waveform\ncross-correlation on the DAS or HRSN records is required to '
          'confirm any of them\nas true repeaters.')


if __name__ == '__main__':
    main()
