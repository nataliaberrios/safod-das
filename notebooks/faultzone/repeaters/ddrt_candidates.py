"""
Rebuild the repeater candidate list on the DOUBLE-DIFFERENCE relocated catalog.

Why redo this. G2 failed (0/5 sequences, best CC 0.893) on a candidate list I
built from the standard USGS catalog with single-link clustering at 0.5 km. Two
defects, both mine:

  * Standard-catalog absolute locations at Parkfield carry ~0.2-0.5 km error,
    while true repeaters are <20 m across. Clustering at 0.5 km therefore groups
    genuinely different patches by construction, and low waveform similarity is
    the expected result -- it was telling me the clustering was wrong, not that
    repeaters are absent.
  * Single-link clustering chains: sequence 0 spanned 1 km of depth despite the
    0.5 km threshold.

William Ellsworth (co-author of Poupinet, Ellsworth & Frechet 1984, the paper
that introduced doublet monitoring) supplied the right catalog: NCEDC
Double-Difference (DDRT), 2024-05-01 to 2026-04-01, within 15 km of SAFOD at
35.982, -120.544. DD RELATIVE locations are accurate to tens of metres, which is
the scale repeaters actually live at.

His guidance, from the email: focus on spots showing multiple dates at the same
location, and prefer larger events. His own script rotates into fault-parallel
coordinates on a 45 degree strike; that rotation is reproduced here.

Changes from my first attempt, each addressing a named defect:
  * DD locations instead of standard catalog
  * COMPLETE linkage instead of single -- no chaining
  * 150 m cluster diameter instead of 500 m
  * magnitude spread < 0.3 required
  * pairs must be > 30 days apart to be useful for CWI

Reads awd_clean/repeating/ (read-only); writes only to faultzone/repeaters/.
"""
import os
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DDRT = ('/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/'
        'repeating/DDRT_DAS.txt')
# Union of BOTH cemented-fiber manifests. Using only SAFODAS1 undercounted
# coverage by 127 days and stopped six months early at 2025-07; the Harvest
# manifest (SAFOD_vertical_2026_01_28.csv) extends it to 2026-01-27 on the same
# cable, so pairs spanning the two are still valid CWI pairs.
DAYS = os.path.join(HERE, 'das_days_all.txt')

SAFOD_LAT, SAFOD_LON = 35.982, -120.544     # MH030, Ellsworth's reference
STRIKE_DEG = 45.0
# Selectable. 150 m yields only 2 candidates, which cannot test a gate that
# requires 3 to pass. 250 m yields 8. Widening the search to make the gate
# answerable is not the same as lowering the CC threshold -- that stays at 0.90.
CLUSTER_M = float(os.environ.get('CLUSTER_M', '150'))
MAG_SPREAD = 0.3
MIN_DAYS = 30.0


def main():
    cols = ['yr', 'mo', 'day', 'hr', 'mi', 'sec', 'lat', 'lon', 'depth', 'mag']
    df = pd.read_csv(DDRT, sep=r'\s+', header=None, usecols=range(10),
                     names=cols, skip_blank_lines=True).dropna()
    df['time'] = pd.to_datetime(dict(
        year=df.yr.astype(int), month=df.mo.astype(int), day=df.day.astype(int),
        hour=df.hr.astype(int), minute=df.mi.astype(int),
        second=df.sec.astype(float).round().astype(int).clip(0, 59)), utc=True)
    print(f'DDRT catalog: {len(df)} events, '
          f'{df.time.min():%Y-%m-%d} to {df.time.max():%Y-%m-%d}')
    print(f'  depth {df.depth.min():.2f}-{df.depth.max():.2f} km, '
          f'mag {df.mag.min():.2f}-{df.mag.max():.2f}')

    # local cartesian, then Ellsworth's fault-parallel / fault-normal rotation
    klat = 111.19
    klon = 111.19 * np.cos(np.radians(SAFOD_LAT))
    n = (df.lat - SAFOD_LAT) * klat
    e = (df.lon - SAFOD_LON) * klon
    s = np.radians(STRIKE_DEG)
    df['r'] = n * np.cos(s) + e * np.sin(s)     # fault-normal
    df['t'] = -n * np.sin(s) + e * np.cos(s)    # along-strike
    df['km_safod'] = np.sqrt(n ** 2 + e ** 2)

    days = set(open(DAYS).read().split())
    df['covered'] = df.time.dt.strftime('%Y-%m-%d').isin(days)
    cov = df[df.covered].reset_index(drop=True)
    print(f'\non cemented-fiber covered days: {len(cov)}/{len(df)}')
    if len(cov) < 3:
        print('too few covered events to cluster')
        return

    # complete linkage: every member within CLUSTER_M of every other member,
    # so a cluster cannot chain across a kilometre the way single link did
    XYZ = np.column_stack([cov.t, cov.r, cov.depth]) * 1000.0     # metres
    Z = linkage(pdist(XYZ), method='complete')
    cov['cl'] = fcluster(Z, t=CLUSTER_M, criterion='distance')

    seqs = []
    for c, g in cov.groupby('cl'):
        if len(g) < 2:
            continue
        if g.mag.max() - g.mag.min() > MAG_SPREAD:
            continue
        span = (g.time.max() - g.time.min()).total_seconds() / 86400
        if span < MIN_DAYS:
            continue
        d = pdist(np.column_stack([g.t, g.r, g.depth]) * 1000.0)
        seqs.append(dict(cl=c, n=len(g), span=span, diam=float(d.max()),
                         depth=float(g.depth.mean()), mag=float(g.mag.mean()),
                         km=float(g.km_safod.mean()), g=g))
    seqs.sort(key=lambda s: (-s['n'], -s['mag']))

    print(f'\nsequences: >=2 events, diameter <{CLUSTER_M:.0f} m, '
          f'mag spread <{MAG_SPREAD}, span >{MIN_DAYS:.0f} d\n')
    print(f'{"cl":>5}{"n":>4}{"diam_m":>9}{"depth":>8}{"mag":>7}'
          f'{"km_safod":>10}{"span_d":>9}   dates')
    for s in seqs:
        ds = ', '.join(f'{t:%Y-%m-%d}' for t in sorted(s['g'].time))
        print(f'{s["cl"]:>5}{s["n"]:>4}{s["diam"]:9.0f}{s["depth"]:8.2f}'
              f'{s["mag"]:7.2f}{s["km"]:10.2f}{s["span"]:9.0f}   {ds[:60]}')
    print(f'\ntotal: {len(seqs)} sequences, '
          f'{sum(s["n"] for s in seqs)} events')

    if seqs:
        rows = []
        for k, s in enumerate(seqs):
            for _, e_ in s['g'].iterrows():
                rows.append(dict(seq=k, cl=s['cl'], time=e_.time, mag=e_.mag,
                                 depth=e_.depth, lat=e_.lat, lon=e_.lon,
                                 km_from_safod=e_.km_safod, t_strike=e_.t,
                                 r_normal=e_.r, diam_m=s['diam']))
        out = os.path.join(HERE, os.environ.get('OUT_CSV', 'ddrt_candidates.csv'))
        pd.DataFrame(rows).to_csv(out, index=False)
        print(f'wrote {out}')

    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    sc = ax[0].scatter(df.t, -df.depth, s=8 * (df.mag + 1),
                       c=df.time.astype('int64'), cmap='turbo', alpha=0.75)
    ax[0].set(xlabel='distance along strike (km)', ylabel='depth (km)',
              title='A  DDRT cross-section, coloured by date')
    cb = plt.colorbar(sc, ax=ax[0])
    cb.set_ticks([df.time.astype('int64').min(), df.time.astype('int64').max()])
    cb.set_ticklabels([f'{df.time.min():%Y-%m}', f'{df.time.max():%Y-%m}'])

    ax[1].scatter(df.t, -df.depth, s=6, c='0.85')
    ax[1].scatter(cov.t, -cov.depth, s=8, c='0.6')
    for k, s in enumerate(seqs):
        ax[1].scatter(s['g'].t, -s['g'].depth, s=40 * (s['g'].mag + 1),
                      marker='o', facecolors='none', edgecolors='C3', lw=1.4)
        ax[1].annotate(f'{s["n"]}', (s['g'].t.mean(), -s['g'].depth.mean()),
                       fontsize=7, color='C3',
                       xytext=(3, 3), textcoords='offset points')
    ax[1].set(xlabel='distance along strike (km)', ylabel='depth (km)',
              title=f'B  {len(seqs)} candidate sequences '
                    f'(complete linkage, <{CLUSTER_M:.0f} m)')
    for a in ax:
        a.grid(alpha=0.3)
    fig.suptitle('Repeater candidates from the double-difference relocated catalog',
                 fontsize=12)
    fig.tight_layout()
    p = os.path.join(HERE, os.environ.get('OUT_CSV', 'ddrt_candidates.csv').replace('.csv', '.png'))
    fig.savefig(p, dpi=140)
    print(f'wrote {p}')


if __name__ == '__main__':
    main()
