"""
Phase A1 -- which DDRT events did the DAS actually record?

Every later count depends on this, and both previous attempts got it wrong in
different directions.

  * The prior screen (awd_clean/repeating/repeater_catalog_report.md) audited
    against SAFOD_2024_2025.csv only. That manifest ends 2025-07-25, so anything
    later is invisible to it -- including 2025-07-27, one of the events in a pair
    HRSN later confirmed at CC 0.991.

  * My own build_event_list.py used DAY-level coverage: an event counted as covered
    if the DAS recorded *something* that calendar day. A day can be listed while the
    specific minute is a gap, so that over-counts.

This does it at file-interval resolution against BOTH manifests. An event is
covered only if some file's [t0, t1) actually contains the window it needs.

The window matters and is not just the origin instant. For the correlation to work
the record must span pre-event noise through coda:

    origin - PRE_S  ...  origin + POST_S

with PRE_S/POST_S matching cache_all/ so the audit describes the same thing the
extractor will produce. A partial overlap is reported separately from full
coverage, because a partially covered event can still be usable for a short
direct-arrival window even when it fails for coda.

Reads only the shared SAFOD manifests and awd_clean/repeating/DDRT_DAS.txt.
Writes only into faultzone/repeaters/.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DDRT = os.path.abspath(os.path.join(
    HERE, '..', '..', 'awd_clean', 'repeating', 'DDRT_DAS.txt'))
OAK = '/oak/stanford/groups/ettore88/data/SAFOD'
MANIFESTS = [
    os.path.join(OAK, 'SAFODAS1-harddrive-transfer', 'SAFOD_2024_2025.csv'),
    os.path.join(OAK, 'SAFOD-Harvest-2026-01-28', 'SAFOD_vertical_2026_01_28.csv'),
]
CACHE_PRE_S, CACHE_POST_S = 5.0, 20.0     # matches extract_all.py
SHORT_PRE_S, SHORT_POST_S = 2.0, 6.0      # minimum for a direct-arrival window

SAFOD_LAT, SAFOD_LON = 35.982, -120.544
STRIKE_DEG = 45.0                          # Ellsworth's rotation


def load_manifests():
    """Both manifests as [t0, t1) intervals, via the loader that already works.

    My first attempt re-parsed the manifests by hand with header=None and five
    column names. The files actually carry a HEADER and ten columns, so every
    field was mis-assigned; and rows with nSamples = -1 are invalid placeholders
    which must be dropped before the timestamps are parsed. Both mistakes
    silently produced 0/329 coverage.

    g1_coda_snr.load_manifest already handles the header, the -1 rows, and the
    union of BOTH manifests, and caches to a pickle. Use it rather than repeat it.
    """
    import sys
    sys.path.insert(0, HERE)
    from g1_coda_snr import load_manifest
    db = load_manifest()
    print(f'  {len(db)} valid intervals, '
          f'{db.t0.min():%Y-%m-%d} to {db.t1.max():%Y-%m-%d}')
    return db[['t0', 't1']].sort_values('t0').reset_index(drop=True)


def load_ddrt():
    cols = ['yr', 'mo', 'day', 'hr', 'mi', 'sec', 'lat', 'lon', 'depth', 'mag']
    df = pd.read_csv(DDRT, sep=r'\s+', header=None, usecols=range(10),
                     names=cols, skip_blank_lines=True).dropna()
    df['time'] = pd.to_datetime(dict(
        year=df.yr.astype(int), month=df.mo.astype(int), day=df.day.astype(int),
        hour=df.hr.astype(int), minute=df.mi.astype(int),
        second=df.sec.astype(float).round().astype(int).clip(0, 59)), utc=True)
    # Ellsworth's fault-parallel rotation, from 'DAS repeater search.src.txt'
    klat = 111.19
    klon = 111.19 * np.cos(np.radians(SAFOD_LAT))
    n = (df.lat - SAFOD_LAT) * klat
    e = (df.lon - SAFOD_LON) * klon
    s = np.radians(STRIKE_DEG)
    df['r'] = n * np.cos(s) + e * np.sin(s)      # fault-normal
    df['t'] = -n * np.sin(s) + e * np.cos(s)     # along-strike
    df['km_safod'] = np.sqrt(n ** 2 + e ** 2)
    return df.sort_values('time').reset_index(drop=True)


def covered(db, origin, pre, post):
    """Is [origin-pre, origin+post] inside the union of manifest intervals?

    Returns (full, partial, fraction_covered). Intervals are merged on the fly so
    an event spanning a file boundary still counts as fully covered.
    """
    a = origin - pd.Timedelta(seconds=pre)
    b = origin + pd.Timedelta(seconds=post)
    sel = db[(db.t0 < b) & (db.t1 > a)]
    if sel.empty:
        return False, False, 0.0
    iv = sorted(zip(sel.t0, sel.t1))
    merged, cs, ce = [], iv[0][0], iv[0][1]
    for s0, s1 in iv[1:]:
        if s0 <= ce:
            ce = max(ce, s1)
        else:
            merged.append((cs, ce)); cs, ce = s0, s1
    merged.append((cs, ce))
    got = sum(max(pd.Timedelta(0), min(e, b) - max(s, a)).total_seconds()
              for s, e in merged)
    need = (b - a).total_seconds()
    frac = got / need
    return frac >= 0.999, frac > 0.0, frac


def main():
    print('MANIFESTS')
    db = load_manifests()
    print(f'  union: {len(db)} intervals\n')

    ev = load_ddrt()
    print(f'DDRT: {len(ev)} events, {ev.time.min():%Y-%m-%d} to '
          f'{ev.time.max():%Y-%m-%d}, M{ev.mag.min():.2f}-{ev.mag.max():.2f}')
    print(f'  (Ellsworth reports 329 in the email -- '
          f'{"MATCH" if len(ev) == 329 else "MISMATCH, investigate"})\n')

    rows = []
    for _, e in ev.iterrows():
        f_full, f_part, f_frac = covered(db, e.time, CACHE_PRE_S, CACHE_POST_S)
        s_full, _, s_frac = covered(db, e.time, SHORT_PRE_S, SHORT_POST_S)
        rows.append(dict(time=e.time, mag=e.mag, depth=e.depth, lat=e.lat,
                         lon=e.lon, r=e.r, t=e.t, km_safod=e.km_safod,
                         cov_full=f_full, cov_partial=f_part, cov_frac=f_frac,
                         cov_short=s_full, cov_short_frac=s_frac))
    E = pd.DataFrame(rows)
    E['tag'] = E.time.dt.strftime('ev_%Y%m%dT%H%M%S')
    E.to_csv(os.path.join(HERE, 'phaseA_events.csv'), index=False)

    print('COVERAGE')
    print(f'  full {CACHE_PRE_S:.0f}s pre / {CACHE_POST_S:.0f}s post : '
          f'{int(E.cov_full.sum()):3d} / {len(E)}')
    print(f'  short window ({SHORT_PRE_S:.0f}/{SHORT_POST_S:.0f}s), usable for '
          f'direct arrivals only: {int(E.cov_short.sum()):3d} / {len(E)}')
    print(f'  any overlap at all                                    : '
          f'{int(E.cov_partial.sum()):3d} / {len(E)}')
    print(f'  no data                                               : '
          f'{int((~E.cov_partial).sum()):3d} / {len(E)}')

    print('\n  by year:')
    for y, g in E.groupby(E.time.dt.year):
        print(f'    {y}: {int(g.cov_full.sum()):3d}/{len(g):3d} full, '
              f'{int(g.cov_short.sum()):3d} short')

    # ---- pairs -------------------------------------------------------
    cov = E[E.cov_full].reset_index(drop=True)
    n = len(cov)
    print(f'\nPAIRS from {n} fully covered events: {n*(n-1)//2}')
    iu = np.triu_indices(n, 1)
    P = pd.DataFrame(dict(i=iu[0], j=iu[1]))
    klon = 111.19 * np.cos(np.radians(SAFOD_LAT))
    dn = (cov.lat.values[iu[1]] - cov.lat.values[iu[0]]) * 111.19
    de = (cov.lon.values[iu[1]] - cov.lon.values[iu[0]]) * klon
    P['h_km'] = np.sqrt(dn ** 2 + de ** 2)
    P['dz_km'] = np.abs(cov.depth.values[iu[1]] - cov.depth.values[iu[0]])
    P['sep_km'] = np.sqrt(P.h_km ** 2 + P.dz_km ** 2)
    P['dmag'] = np.abs(cov.mag.values[iu[1]] - cov.mag.values[iu[0]])
    P['mag_max'] = np.maximum(cov.mag.values[iu[0]], cov.mag.values[iu[1]])
    P['days'] = np.abs((cov.time.values[iu[1]] - cov.time.values[iu[0]])
                       / np.timedelta64(1, 'D'))
    P['t_i'] = cov.time.values[iu[0]]; P['t_j'] = cov.time.values[iu[1]]
    P['m_i'] = cov.mag.values[iu[0]];  P['m_j'] = cov.mag.values[iu[1]]
    P['tag_i'] = cov.tag.values[iu[0]]; P['tag_j'] = cov.tag.values[iu[1]]

    print('\nSCREENS (Ellsworth: close in magnitude and location, prefer larger)')
    for lbl, m in [
            ('prior screen (sep<1 km, dM<0.5)', (P.h_km <= 1.0) &
             (P.dz_km <= 1.0) & (P.dmag <= 0.5)),
            ('  ... and >7 d apart', (P.h_km <= 1.0) & (P.dz_km <= 1.0) &
             (P.dmag <= 0.5) & (P.days > 7)),
            ('tight (sep<50 m, dM<0.15, >30 d)', (P.sep_km < 0.05) &
             (P.dmag < 0.15) & (P.days > 30)),
            ('  ... and M>=1.0 (Ellsworth prefers larger)', (P.sep_km < 0.05) &
             (P.dmag < 0.15) & (P.days > 30) & (P.mag_max >= 1.0)),
            ('very tight (sep<20 m ~ source radius at M1)', (P.sep_km < 0.02) &
             (P.dmag < 0.15) & (P.days > 30))]:
        print(f'  {lbl:<46} {int(m.sum()):5d}')

    P.to_csv(os.path.join(HERE, 'phaseA_pairs.csv'), index=False)
    print(f'\nwrote phaseA_events.csv ({len(E)} events) and '
          f'phaseA_pairs.csv ({len(P)} pairs)')

    # ---- reconcile with what we already believe -----------------------
    print('\nRECONCILIATION with the 8 HRSN-confirmed pairs')
    hp = os.path.join(HERE, 'hrsn_extended.csv')
    if os.path.exists(hp):
        h = pd.read_csv(hp).dropna(subset=['hrsn'])
        h = h[h.hrsn > 0.90]
        key = lambda a, b: '|'.join(sorted([pd.Timestamp(a).strftime('%Y%m%dT%H%M'),
                                            pd.Timestamp(b).strftime('%Y%m%dT%H%M')]))
        have = {key(r.t_i, r.t_j): r.hrsn for _, r in h.iterrows()}
        mine = {key(r.t_i, r.t_j) for _, r in P.iterrows()}
        found = sum(1 for k in have if k in mine)
        print(f'  {found}/{len(have)} appear among fully covered pairs here')
        for k, cc in sorted(have.items(), key=lambda x: -x[1]):
            print(f'    HRSN {cc:.3f}  {k}  '
                  f'{"present" if k in mine else "NOT fully covered by this audit"}')


if __name__ == '__main__':
    main()
