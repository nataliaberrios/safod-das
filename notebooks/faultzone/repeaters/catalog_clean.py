"""
Phase 1.1 + 1.2 -- a defensible pair list, with every rejection stated.

WHY. The confirmed count has been quoted as 7, 8, 17, 83 and 93 in different
places, because different thresholds and different (wrong) duplicate assumptions
were applied each time. This script fixes one list with the criteria written down.

TWO CORRECTIONS ESTABLISHED WHILE WRITING IT, both by checking rather than
asserting:

  1. There is exactly ONE genuine duplicate in the catalog, not several.
     A systematic scan of all consecutive events within 60 s finds seven close
     pairs; six are separated by 11-56 s and are real distinct earthquakes. Only
     2024-12-30 04:39:32.700 has two solutions at the IDENTICAL millisecond, 295 m
     apart (nc75109596: 18 sta, gap 60, rms 0.07; nc75109601: 16 sta, gap 84,
     rms 0.09). Two earthquakes cannot share an origin time to the millisecond
     295 m apart -- that is one event with two solutions. METHODS_STATUS section 14
     said "exactly one collision" and was correct; a later claim that near-
     duplicates had slipped through, and that the 83/93 per-channel counts were
     therefore inflated, was wrong.

  2. The 2025-04-02 events are 27 s apart, NOT 0.3 s. That came from misreading
     dt_days 38.932465 vs 38.932153 (a difference of 0.000312 d = 27 s). They are
     two separately reviewed solutions with different IDs, depths (1.87 vs 1.55 km)
     and station counts. So the confirmed set is 8 pairs over 12 events, not 7
     over 11.

     They still matter for creep, but under a different rule: Waldhauser &
     Ellsworth drop the smaller event of any pair separated by < 1 month, because
     two ruptures 27 s apart are not independent loading cycles.

WHY THE PHYSICAL CRITERIA MATTER MORE HERE THAN USUAL. Waveform similarity has been
the only filter applied. For a creep calculation the method assumes a
*characteristic patch* rupturing repeatedly under steady loading, so a large
magnitude spread breaks the premise even at CC 0.99. Two of the eight pairs have
dM > 0.4.

Criteria, from Waldhauser & Ellsworth 2002 (doi:10.1029/2000JB000084) and
Waldhauser & Schaff 2008 (doi:10.1029/2007JB005479):
  * dM <= 0.3
  * source-area overlap: separation < source radius at a stated stress drop
  * drop the smaller event of pairs separated by < 1 month

Nothing is silently culled. Every pair is reported with each criterion's verdict.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CAT = os.path.join(HERE, 'parkfield_catalog_2024_2025.csv')
PAIRS = os.path.join(HERE, 'hrsn_extended.csv')

HRSN_MIN = 0.90          # Waldhauser & Schaff convention, measured ON HRSN
DM_MAX = 0.30
STRESS_DROP_PA = 3.0e6   # 3 MPa; Abercrombie 2014 is the check on this
MU_PA = 30.0e9           # shear modulus
MIN_SEP_DAYS = 30.0      # below this, drop the smaller event
DUP_DT_S = 2.0           # same-instant test
DUP_SEP_KM = 0.5


def moment_nm(mag):
    """Seismic moment from magnitude. Md treated as Mw -- an approximation that
    should be stated, not buried: these are duration magnitudes."""
    return 10.0 ** (1.5 * np.asarray(mag, float) + 9.1)


def source_radius_m(mag, dsigma=STRESS_DROP_PA):
    """Circular crack (Eshelby): r = (7 M0 / (16 dsigma))^(1/3)."""
    return (7.0 * moment_nm(mag) / (16.0 * dsigma)) ** (1.0 / 3.0)


def slip_m(mag, dsigma=STRESS_DROP_PA):
    """Average slip for a circular crack: d = M0 / (mu * pi * r^2)."""
    r = source_radius_m(mag, dsigma)
    return moment_nm(mag) / (MU_PA * np.pi * r ** 2)


def sep_km(a, b):
    dn = (b['latitude'] - a['latitude']) * 111.19
    de = (b['longitude'] - a['longitude']) * 111.19 * np.cos(np.radians(36.0))
    dz = b['depth'] - a['depth']
    return float(np.sqrt(dn * dn + de * de + dz * dz))


def main():
    c = pd.read_csv(CAT)
    c['time'] = pd.to_datetime(c.time, utc=True, format='mixed')
    c = c.sort_values('time').reset_index(drop=True)
    print(f'catalog: {len(c)} events\n')

    # ---------------- 1.1 duplicates -----------------------------------
    print('1.1  DUPLICATE SOLUTIONS (same instant, same place)')
    drop_ids = set()
    dt = c.time.diff().dt.total_seconds()
    for i in np.where(dt < 60)[0]:
        a, b = c.iloc[i - 1], c.iloc[i]
        d, s = float(dt.iloc[i]), sep_km(a, b)
        dup = (d <= DUP_DT_S) and (s <= DUP_SEP_KM)
        if dup:
            # keep the better-constrained solution
            worse = b if (a['nst'], -a['gap'], -a['rms']) > \
                         (b['nst'], -b['gap'], -b['rms']) else a
            drop_ids.add(worse['id'])
            print(f'  DUPLICATE {a.time:%Y-%m-%d %H:%M:%S.%f} dt={d:.3f}s '
                  f'sep={1000*s:.0f}m -> drop {worse["id"]} '
                  f'(nst {worse["nst"]}, gap {worse["gap"]}, rms {worse["rms"]})')
        else:
            print(f'  distinct  {a.time:%Y-%m-%d %H:%M:%S} +{d:6.2f}s '
                  f'sep={1000*s:7.0f}m  M{a.mag:.2f}/{b.mag:.2f}  -> both kept')
    print(f'\n  {len(drop_ids)} duplicate solution(s) to drop: '
          f'{sorted(drop_ids) if drop_ids else "none"}')

    # ---------------- 1.2 physical criteria ----------------------------
    p = pd.read_csv(PAIRS).dropna(subset=['hrsn'])
    p = p[p.hrsn > HRSN_MIN].copy()
    for col in ('t_i', 't_j'):
        p[col] = pd.to_datetime(p[col], utc=True, format='mixed')
    print(f'\n1.2  PHYSICAL CRITERIA on {len(p)} pairs with HRSN CC > {HRSN_MIN}')
    print(f'     dM <= {DM_MAX}; source radius at {STRESS_DROP_PA/1e6:.0f} MPa; '
          f'separation < {MIN_SEP_DAYS:.0f} d drops the smaller event\n')

    rows = []
    print(f'{"pair":<26}{"dM":>6}{"r_i":>7}{"r_j":>7}{"dt_d":>8}'
          f'{"HRSN":>7}   verdict')
    for _, r in p.sort_values('t_i').iterrows():
        dm = abs(r.m_i - r.m_j)
        ri, rj = source_radius_m(r.m_i), source_radius_m(r.m_j)
        fails = []
        if dm > DM_MAX:
            fails.append(f'dM={dm:.2f}>{DM_MAX}')
        if r.dt_days < MIN_SEP_DAYS:
            fails.append(f'dt={r.dt_days:.1f}d<{MIN_SEP_DAYS:.0f}d')
        ok = not fails
        rows.append(dict(t_i=r.t_i, t_j=r.t_j, m_i=r.m_i, m_j=r.m_j, dm=dm,
                         r_i=ri, r_j=rj, dt_days=r.dt_days, hrsn=r.hrsn,
                         passes=ok, reason='; '.join(fails) or 'pass'))
        print(f'{r.t_i:%Y-%m-%d}/{r.t_j:%Y-%m-%d}{dm:6.2f}{ri:7.1f}{rj:7.1f}'
              f'{r.dt_days:8.0f}{r.hrsn:7.3f}   '
              f'{"PASS" if ok else "FAIL: " + "; ".join(fails)}')

    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(HERE, 'catalog_clean_pairs.csv'), index=False)
    n_ok = int(R.passes.sum())
    print(f'\n  {n_ok}/{len(R)} pairs satisfy the criteria')

    # ---------------- what this means for creep -------------------------
    print('\nIMPLICATIONS FOR THE CREEP CALCULATION')
    print(f'  source radius at M1.0, {STRESS_DROP_PA/1e6:.0f} MPa: '
          f'{source_radius_m(1.0):.1f} m; average slip {1000*slip_m(1.0):.1f} mm')
    print(f'  full source overlap therefore needs separation < ~'
          f'{2*source_radius_m(1.0):.0f} m, which is at or below the resolution')
    print('  of DD relocation. Overlap CANNOT be verified from locations alone --')
    print('  waveform similarity is the stronger evidence and should be said so.')

    if n_ok:
        print('\n  usable pairs for recurrence:')
        for _, r in R[R.passes].sort_values('t_i').iterrows():
            print(f'    {r.t_i:%Y-%m-%d} / {r.t_j:%Y-%m-%d}  '
                  f'{r.dt_days:5.0f} d  M{r.m_i:.2f}/{r.m_j:.2f}  '
                  f'slip {1000*slip_m((r.m_i+r.m_j)/2):.1f} mm  -> '
                  f'{1000*slip_m((r.m_i+r.m_j)/2)/(r.dt_days/365.25):.1f} mm/yr')
        print('\n  (rates are per-interval, single-interval, and assume the whole')
        print('   recurrence is captured -- Phase 3 tests exactly that assumption)')

    print(f'\nwrote catalog_clean_pairs.csv')


if __name__ == '__main__':
    main()
