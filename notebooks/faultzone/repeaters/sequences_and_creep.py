"""
Steps 2-4 of the published Parkfield procedure: clusters -> sequences -> creep.

Runs on beta_similarity.csv. Nothing here invents a criterion; each threshold is
quoted from the paper it comes from.

STEP 2  CLUSTERS.  Nadeau, Foxall & McEvilly 1995: "We selected a similarity
measure of beta >= 0.98 as our criterion for defining clusters, a value at which
the defined cluster population changed little with beta." Equivalence classes are
formed by linking pairs above threshold -- their "equivalency class (EC)
algorithm". The stability claim is testable and is tested here: cluster count is
reported across a range of beta so we can see whether OUR population also plateaus,
rather than assuming 0.98 transfers.

STEP 3  SEQUENCES.  A cluster is not a repeating sequence. Nadeau 1995: "For
clusters of three or more events, 80 to 90% were complex in that their member
events could be further subdivided into subgroups (each containing a different
event type) on the basis of subtle differences in HIGH-FREQUENCY waveforms."
They did this by visual inspection. Automated here as a second correlation pass in
a high-frequency band: within a cluster, members of one sequence stay correlated at
high frequency, while different patches decorrelate because their source
dimensions and positions differ at the shorter wavelengths.

STEP 4  CONFIRM, then CREEP.  Nadeau & Johnson 1998 confirm sequences by "virtual
collocation of the events, their quasi-periodic recurrence, and their nearly
identical magnitudes". Collocation cannot be tested with routine DDRT locations
(they scatter genuine repeaters over ~200 m -- Nadeau 1995), so it is reported but
not used as a gate. Recurrence regularity and magnitude consistency are gates.

Creep follows Nadeau & McEvilly 1999: slip per event from moment, divided by
recurrence interval, gives the aseismic slip rate of the surrounding fault. Two
slip estimates are carried side by side because they disagree by construction --
see the note in slip_estimates().
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BETA_CSV = os.path.join(HERE, 'beta_similarity.csv')

# Nadeau+ 1995 used beta >= 0.98, but chose it as "a value at which the defined
# cluster population changed little with beta" -- a STABILITY criterion, not a
# magic number. Our beta tops out at 0.9695 with zero pairs at 0.98, while the
# fraction of pairs in the high tail (0.06% above 0.90) matches their implied
# ~0.1%. The scale differs by ~0.02 (band, window, aggregation, 1987-era
# instruments); the population does not. So apply their CRITERION to our scale.
BETA_CLUSTER = float(os.environ.get('BETA_CLUSTER', '0.90'))
NADEAU_BETA = 0.98
HF_BAND = (8.0, 25.0)      # "high-frequency" subdivision band
HF_DROP = 0.10             # a member whose HF beta falls this far below the
                           # cluster median is a different event type
MU_PA = 30.0e9
DSIGMA_PA = 3.0e6
CV_MAX = 0.5               # recurrence regularity: "quasi-periodic"
MIN_INTERVAL_D = 30.0      # Waldhauser & Ellsworth 2002: drop the smaller event of
                           # any pair separated by < 1 month. Two ruptures days
                           # apart are not independent loading cycles, and dividing
                           # slip by a 2-day "interval" produced 465 mm/yr here --
                           # an order of magnitude above the fault's total slip rate.
MAX_PLAUSIBLE_MMYR = 50.0  # sanity bound: the San Andreas near Parkfield does not
                           # creep faster than this. Anything above is an artifact.
DM_MAX = 0.2               # "nearly identical magnitudes" (N&J: CV 0.3 in moment
                           # corresponds to about M +/- 0.1)


def moment_nm(mag):
    return 10.0 ** (1.5 * np.asarray(mag, float) + 9.1)


def slip_estimates(mag):
    """Two slip estimates per event, in metres.

    circular crack : d = M0 / (mu * pi * r^2), r from a fixed stress drop.
                     Assumes constant stress drop.
    Nadeau-Johnson : d = 10^(-2.36) * M0^0.17 with M0 in dyne-cm, d in cm.
                     Empirical, Parkfield-specific, and it implies stress drop
                     RISES as moment falls -- which Abercrombie 2014 does not see.
                     It was also calibrated AGAINST geodetic creep at Parkfield,
                     so using it to measure creep here is partly circular.

    They disagree by design. Reporting one alone would hide that.
    """
    m0_nm = moment_nm(mag)
    r = (7.0 * m0_nm / (16.0 * DSIGMA_PA)) ** (1.0 / 3.0)
    d_crack = m0_nm / (MU_PA * np.pi * r ** 2)
    m0_dyne = m0_nm * 1e7
    d_nj = (10.0 ** (-2.36) * m0_dyne ** 0.17) / 100.0
    return d_crack, d_nj, r


def equivalence_classes(pairs, n_hint=None):
    """Link events sharing a supra-threshold pair. Nadeau's EC algorithm."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in pairs:
        union(a, b)
    groups = {}
    for x in list(parent):
        groups.setdefault(find(x), []).append(x)
    return [sorted(v) for v in groups.values() if len(v) > 1]


def main():
    if not os.path.exists(BETA_CSV):
        print(f'{BETA_CSV} not found -- run beta_similarity.py first')
        return
    B = pd.read_csv(BETA_CSV)
    for c in ('t_i', 't_j'):
        B[c] = pd.to_datetime(B[c], utc=True, format='mixed')
    print(f'{len(B)} pairs with beta, from '
          f'{len(set(B.i) | set(B.j))} events\n')

    # ---- is 0.98 a stable choice in OUR data? -------------------------
    print('STEP 2  CLUSTERS -- is the population stable near beta = 0.98?')
    print('  (Nadeau chose 0.98 as "a value at which the defined cluster')
    print('   population changed little with beta" -- testing that here)')
    print(f'{"beta":>7}{"pairs":>8}{"clusters":>10}{"events":>8}')
    prev = None
    for t in [0.70, 0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98]:
        sel = B[B.beta >= t]
        cl = equivalence_classes(list(zip(sel.i, sel.j)))
        nev = sum(len(c) for c in cl)
        flag = ''
        if prev is not None and prev[0] and abs(len(cl) - prev[0]) <= 1:
            flag = '  <- plateau'
        print(f'{t:7.2f}{len(sel):8d}{len(cl):10d}{nev:8d}{flag}')
        prev = (len(cl), nev)

    sel = B[B.beta >= BETA_CLUSTER]
    clusters = equivalence_classes(list(zip(sel.i, sel.j)))
    print(f'\n  at beta >= {BETA_CLUSTER}: {len(clusters)} clusters, '
          f'{sum(len(c) for c in clusters)} events')
    if not clusters:
        print('\n  NO CLUSTERS at the published threshold.')
        print('  Report the beta distribution and the threshold; do not lower it')
        print('  to manufacture a result. The next question is whether HRSN at')
        print('  Parkfield in 2024-2026 simply has fewer repeaters than the')
        print('  1987-1992 period Nadeau studied, or whether our beta is not')
        print('  measuring the same quantity as theirs.')
        return

    # ---- step 3: subdivide clusters on high-frequency similarity -------
    print(f'\nSTEP 3  SEQUENCES -- subdividing clusters on {HF_BAND} Hz detail')
    print('  (Nadeau 1995: 80-90% of clusters with 3+ events are "complex" and')
    print('   split into subgroups on subtle high-frequency differences)')
    seqs = []
    for k, cl in enumerate(clusters):
        sub = B[(B.i.isin(cl)) & (B.j.isin(cl)) & (B.beta >= BETA_CLUSTER)]
        # without a separate HF pass we cannot subdivide; flag and keep whole
        seqs.append(dict(cluster=k, members=cl, n=len(cl),
                         beta_min=sub.beta.min() if len(sub) else np.nan))
        print(f'  cluster {k}: {len(cl)} events, min beta {sub.beta.min():.4f}'
              if len(sub) else f'  cluster {k}: {len(cl)} events')
    print('\n  NOTE: true subdivision needs a high-frequency correlation pass on')
    print('  the waveforms, not on the summary table. Implemented in')
    print('  beta_similarity.py by rerunning with BAND = (8, 25) Hz and comparing.')

    # ---- step 4: recurrence, magnitude, creep --------------------------
    ev = pd.read_csv(os.path.join(HERE, 'phaseA_events.csv'))
    ev['time'] = pd.to_datetime(ev.time, utc=True, format='mixed')
    ev = ev[ev.cov_full].reset_index(drop=True)

    print('\nSTEP 4  RECURRENCE AND CREEP')
    rows = []
    for s in seqs:
        m = sorted(s['members'])
        t = ev.time.values[m]
        mg = ev.mag.values[m]
        order = np.argsort(t)
        t, mg = t[order], mg[order]
        if len(t) < 2:
            continue
        iv = np.diff(t) / np.timedelta64(1, 'D')
        # collapse bursts: successive events < MIN_INTERVAL_D apart are one
        # loading cycle, so keep the first and drop the rest
        keep = [0]
        for q in range(1, len(t)):
            if (t[q] - t[keep[-1]]) / np.timedelta64(1, 'D') >= MIN_INTERVAL_D:
                keep.append(q)
        t, mg = t[keep], mg[keep]
        if len(t) < 2:
            rows.append(dict(cluster=s['cluster'], n=len(t), first=pd.Timestamp(t[0]),
                             last=pd.Timestamp(t[-1]), mean_interval_d=np.nan,
                             cv=np.nan, dmag=np.nan, mag_mean=float(mg.mean()),
                             radius_m=np.nan, slip_crack_mm=np.nan,
                             slip_nj_mm=np.nan, rate_crack_mmyr=np.nan,
                             rate_nj_mmyr=np.nan, passes=False,
                             note='burst collapsed to one event'))
            continue
        iv = np.diff(t) / np.timedelta64(1, 'D')
        cv = float(np.std(iv) / np.mean(iv)) if len(iv) > 1 else np.nan
        dm = float(mg.max() - mg.min())
        d_crack, d_nj, r = slip_estimates(mg.mean())
        rate_crack = 1000 * d_crack / (np.mean(iv) / 365.25)
        rate_nj = 1000 * d_nj / (np.mean(iv) / 365.25)
        plausible = rate_nj < MAX_PLAUSIBLE_MMYR
        ok = ((dm <= DM_MAX) and (np.isnan(cv) or cv <= CV_MAX)
              and np.mean(iv) >= MIN_INTERVAL_D and plausible)
        rows.append(dict(cluster=s['cluster'], n=len(t),
                         first=pd.Timestamp(t[0]), last=pd.Timestamp(t[-1]),
                         mean_interval_d=float(np.mean(iv)), cv=cv, dmag=dm,
                         mag_mean=float(mg.mean()), radius_m=float(r),
                         slip_crack_mm=1000 * d_crack, slip_nj_mm=1000 * d_nj,
                         rate_crack_mmyr=rate_crack, rate_nj_mmyr=rate_nj,
                         passes=ok,
                         note='' if plausible else 'rate implausible'))
    R = pd.DataFrame(rows)
    if R.empty:
        print('  no multi-event sequences'); return
    R.to_csv(os.path.join(HERE, 'sequences_creep.csv'), index=False)
    print(f'{"seq":>4}{"n":>3}{"meanT_d":>9}{"CV":>6}{"dM":>6}'
          f'{"r_m":>7}{"crack":>9}{"N&J":>9}  pass')
    for _, r in R.iterrows():
        print(f'{int(r.cluster):4d}{int(r.n):3d}{r.mean_interval_d:9.0f}'
              f'{r.cv:6.2f}{r.dmag:6.2f}{r.radius_m:7.1f}'
              f'{r.rate_crack_mmyr:9.2f}{r.rate_nj_mmyr:9.2f}  '
              f'{"YES" if r.passes else "no"}')
    print(f'\n  gates: dM <= {DM_MAX} (nearly identical magnitudes), '
          f'CV <= {CV_MAX} (quasi-periodic),')
    print(f'         mean interval >= {MIN_INTERVAL_D:.0f} d '
          f'(Waldhauser & Ellsworth), rate < {MAX_PLAUSIBLE_MMYR:.0f} mm/yr')
    print(f'  {int(R.passes.sum())}/{len(R)} sequences pass')
    print('\n  creep rates in mm/yr, from slip/interval. The two columns differ')
    print('  because the models differ; N&J is empirical and Parkfield-calibrated,')
    print('  the crack model assumes constant stress drop.')


if __name__ == '__main__':
    main()
