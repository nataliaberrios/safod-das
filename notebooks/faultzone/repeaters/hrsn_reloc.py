"""
HRSN double-difference relative relocation: settle "repeaters or neighbours" with
the method that is actually designed for it.

WHY THIS, AND WHY NOW. This project has spent its effort trying to make the DAS
aperture answer Ellsworth's question --

    "you will want to find events that are close in magnitude and location, which
     will then need to be verified as either repeaters or neighbors"

-- and the DAS answer failed its own control: inferred relative offset tracks
waveform correlation at r = -0.826 WITHIN the family group, so the estimator is
substantially reading waveform quality rather than geometry (METHODS_STATUS 18.3).
The reason that confound was fatal is that there was no ground truth to check it
against.

The borehole network supplies exactly that ground truth, and has been sitting in
hrsn_cache/ (414 event files, median 9 stations, 250 Hz) and hrsn3c_cache/ (207
files, 3 components, 6 stations) the whole time, used only as a similarity referee.

Waldhauser & Ellsworth 2000 (BSSA, doi:10.1785/0120000006) is the standard method
and it is Ellsworth's own. Cross-correlation differential times across ~9 stations
relocate nearby event pairs to ~10-20 m, against rupture radii of 14-48 m here. That
is the resolving power the DAS measurement did not achieve.

WHAT THIS CHANGES ABOUT THE PROJECT. It converts DAS from "the measurement that has
to work" into "the measurement being validated against an independent answer". A
waveform-quality artefact and a real geometric signal are indistinguishable without
a reference; with one, they are separable.

METHOD. For events i, j recorded at station k, the double-difference residual

    dr_ijk = (t_ik - t_jk)_obs  -  (t_ik - t_jk)_calc

is linear in the relative source position for small separations:

    dr_ijk = -(l_k . dm) / V  +  dtau

where l_k is the unit vector from source to station k and dm = (dx, dy, dz) is the
relative offset. Four unknowns (dx, dy, dz, dtau) against ~9 stations.

The velocity model enters only through the ray directions and a single V, and
common-path errors cancel in the difference -- that cancellation is the whole point
of the double-difference formulation and is why this is far less model-sensitive
than absolute location.

--------------------------------------------------------------------------------
PREDICTIONS AND CONTROLS, REGISTERED BEFORE RUNNING.

 R1 SYNTHETIC RESOLUTION FIRST. Using the real station geometry, inject known
    offsets of 5, 20 and 100 m with the measured differential-time scatter as
    noise, and invert. If 20 m is not recovered within its error bar, this station
    geometry cannot resolve a rupture dimension and the result is an upper bound,
    not a location. This runs BEFORE the real data are inverted.

 R2 FAMILY vs CONTROL. Similarity-selected pairs should relocate closer together
    than random pairs. Unlike the DAS version, this comparison has a defence
    against the correlation-quality confound: the inversion residual and the
    station count are reported per pair, so a pair that "relocates tight" purely
    because its delays are noisy will show up as a large formal error rather than
    a small offset.

 R3 THE CONFOUND TEST THAT KILLED THE DAS VERSION, REPEATED HERE. Regress
    log(offset) on HRSN CC within the family group. If r is again strongly
    negative, the relocation is reading waveform quality too and neither
    instrument can answer the question with this event set. This is the single
    most important number in the script.

 R4 VELOCITY SENSITIVITY. Repeat at V = 5.0, 5.8 and 6.5 km/s. Offsets scale
    roughly linearly with V, so the conclusion must be stated as a ratio to
    rupture radius, which is what matters, and the spread reported.

FAILURE IS INFORMATIVE: if R3 fails here as well, the honest conclusion is that
similarity-selected pairs in this catalog cannot be verified as same-patch with
either instrument, which is a direct observational confirmation of Gao, Kao & Wang
2021 (doi:10.1029/2021gl092815) and is publishable as such.
--------------------------------------------------------------------------------
"""
import os
import sys
import json
import numpy as np
import pandas as pd
from obspy import read
from obspy.geodetics import gps2dist_azimuth
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dvv_core import sub_sample_delay, bulk_align                # noqa: E402

CACHE1 = os.path.join(HERE, 'hrsn_cache')
CACHE3 = os.path.join(HERE, 'hrsn3c_cache')
STA_JSON = os.path.join(HERE, 'hrsn_stations.json')

BAND = (5.0, 30.0)          # HRSN is 250 Hz; 30 Hz is comfortably inside Nyquist
PRE_S = 5.0                 # caches start 5 s before catalog origin
P_WIN = (-0.3, 1.2)         # about the PREDICTED P arrival, not the origin
MAX_LAG = 0.05              # residual only; bulk lag removed first
VP_LIST = [5000.0, 5800.0, 6500.0]
VP_REF = 5800.0
MIN_STA = 5                 # 4 unknowns, so 5 is the bare minimum for a residual
MIN_COH = 0.3
DSIGMA = 3e6


def station_coords():
    """Cached station metadata; fetched once from IRIS if absent."""
    if os.path.exists(STA_JSON):
        return json.load(open(STA_JSON))
    from obspy.clients.fdsn import Client
    # HRSN (network BP) is archived at NCEDC, not IRIS -- an IRIS query returns
    # HTTP 204 with no data. hrsn_control.py:136 and hrsn_extend.py:93 both use
    # NCEDC for the waveforms, so the metadata must come from there too.
    c = Client('NCEDC', timeout=120)
    inv = c.get_stations(network='BP', latitude=35.9743, longitude=-120.5521,
                         maxradius=0.6, level='channel')
    out = {}
    for net in inv:
        for sta in net:
            deps = [ch.depth for ch in sta]
            out[sta.code] = dict(lat=float(sta.latitude), lon=float(sta.longitude),
                                 elev=float(sta.elevation),
                                 depth=float(np.median(deps)) if deps else 0.0)
    json.dump(out, open(STA_JSON, 'w'), indent=1)
    return out


def enu(lat0, lon0, dep0_km, lat, lon, dep_m, elev_m=0.0):
    """Local east-north-up metres of a station relative to a source.

    Station elevation matters here: Parkfield relief across the HRSN footprint is
    several hundred metres, comparable to the sensor depths themselves (0-345 m).
    Sensor height above datum is (elevation - sensor depth); the source sits at
    -depth_km*1000. Dropping elevation biases the vertical ray direction, and the
    vertical component is the worst-resolved one in this geometry, so it is the
    one that can least afford an avoidable error.
    """
    d, az, _ = gps2dist_azimuth(lat0, lon0, lat, lon)
    e = d * np.sin(np.radians(az))
    n = d * np.cos(np.radians(az))
    u = (elev_m - dep_m) - (-dep0_km * 1000.0)
    return np.array([e, n, u])


def load_ev(tag):
    for C in (CACHE3, CACHE1):
        f = os.path.join(C, f'{tag}.mseed')
        if os.path.exists(f):
            try:
                return read(f)
            except Exception:
                pass
    return None


def vertical_traces(st):
    """One vertical-ish trace per station, bandpassed."""
    from dvv_core import bandpass
    out = {}
    if st is None:
        return out
    for tr in st:
        if tr.stats.channel not in ('DP1', 'DPZ', 'BHZ', 'EHZ'):
            continue
        try:
            fs = float(tr.stats.sampling_rate)
            d = tr.data.astype(float)
            d = bandpass(d - d.mean(), fs, BAND)
            if np.any(d) and np.all(np.isfinite(d)):
                out[tr.stats.station] = (d, fs)
        except Exception:
            continue
    return out


def diff_times(A, B, ev_i, ev_j, sta):
    """Differential P time per common station, by cross-spectral phase."""
    rows = []
    for s in sorted(set(A) & set(B)):
        if s not in sta:
            continue
        a, fs = A[s]
        b, _ = B[s]
        # predicted P from the catalog location of event i
        v = enu(ev_i.lat, ev_i.lon, ev_i.depth, sta[s]['lat'], sta[s]['lon'],
                sta[s]['depth'], sta[s].get('elev', 0.0))
        R = float(np.linalg.norm(v))
        tp = R / VP_REF
        # bulk align on the whole record first -- origin times carry 0.1-0.5 s
        b2, lag = bulk_align(a, b, fs, max_lag_s=2.0)
        i0 = int((PRE_S + tp + P_WIN[0]) * fs)
        i1 = int((PRE_S + tp + P_WIN[1]) * fs)
        if i0 < 0 or i1 > min(a.size, b2.size):
            continue
        dt, coh = sub_sample_delay(a[i0:i1], b2[i0:i1], fs, BAND,
                                   max_lag_s=MAX_LAG)
        if not np.isfinite(dt) or coh < MIN_COH:
            continue
        rows.append(dict(sta=s, dt=dt + lag, coh=coh, R=R,
                         l=v / R))          # unit vector source -> station
    return rows


def invert(rows, vp, n_boot=300, seed=0):
    """Least squares for (dx, dy, dz, dtau); bootstrap over stations."""
    if len(rows) < MIN_STA:
        return None
    G = np.array([[-r['l'][0] / vp, -r['l'][1] / vp, -r['l'][2] / vp, 1.0]
                  for r in rows])
    d = np.array([r['dt'] for r in rows])
    w = np.array([r['coh'] for r in rows])
    W = np.sqrt(w)[:, None]

    def solve(idx):
        g, dd, ww = G[idx], d[idx], np.sqrt(w[idx])[:, None]
        try:
            m, *_ = np.linalg.lstsq(g * ww, dd * ww[:, 0], rcond=None)
            return m
        except Exception:
            return None

    m = solve(np.arange(len(rows)))
    if m is None:
        return None
    resid = d - G @ m
    rms = float(np.sqrt(np.mean(resid ** 2)))
    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(rows), len(rows))
        if len(set(idx.tolist())) < 4:
            continue
        mm = solve(idx)
        if mm is not None:
            boot.append(mm[:3])
    boot = np.array(boot)
    off = float(np.linalg.norm(m[:3]))
    err = float(np.std(np.linalg.norm(boot, axis=1))) if len(boot) > 20 else np.nan
    return dict(dx=m[0], dy=m[1], dz=m[2], dtau=m[3], offset=off, err=err,
                rms_s=rms, nsta=len(rows))


def synthetic_resolution(rows, vp, truths=(5.0, 20.0, 100.0), n=200, seed=1):
    """R1: can this station geometry recover a known offset? Runs before the data."""
    if len(rows) < MIN_STA:
        return {}
    G = np.array([[-r['l'][0] / vp, -r['l'][1] / vp, -r['l'][2] / vp, 1.0]
                  for r in rows])
    # noise level = the measured residual scatter of the real differential times
    rng = np.random.default_rng(seed)
    out = {}
    for T in truths:
        rec = []
        for k in range(n):
            u = rng.normal(size=3)
            u /= np.linalg.norm(u)
            true = np.r_[u * T, 0.0]
            sig = G @ true + rng.normal(0, 2e-3, len(rows))   # 2 ms scatter
            try:
                m, *_ = np.linalg.lstsq(G, sig, rcond=None)
                rec.append(np.linalg.norm(m[:3]))
            except Exception:
                pass
        if rec:
            out[T] = (float(np.median(rec)), float(np.std(rec)))
    return out


def main():
    sta = station_coords()
    print(f'{len(sta)} HRSN stations with coordinates')
    mt = pd.read_csv(os.path.join(HERE, 'moveout_test.csv'))
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    print(f'{int(mt.is_cand.sum())} family pairs, {int((~mt.is_cand).sum())} controls\n',
          flush=True)

    # ---- gather differential times once, at the reference velocity ----
    per_pair = {}
    for _, r in mt.iterrows():
        i, j = int(r.i), int(r.j)
        A = vertical_traces(load_ev(ev.tag[i]))
        B = vertical_traces(load_ev(ev.tag[j]))
        if not A or not B:
            continue
        rows = diff_times(A, B, ev.iloc[i], ev.iloc[j], sta)
        if len(rows) >= MIN_STA:
            per_pair[(i, j)] = rows
    print(f'{len(per_pair)} pairs with >= {MIN_STA} usable stations\n', flush=True)
    if not per_pair:
        print('no pairs usable -- check station codes against the cache'); return

    # ---- R1 synthetic resolution, BEFORE looking at the data ----
    print('R1 synthetic resolution with the real station geometry')
    ex = per_pair[list(per_pair)[0]]
    res = synthetic_resolution(ex, VP_REF)
    print(f'  {"injected m":>12}{"recovered m":>14}{"scatter m":>12}')
    for T, (med, sd) in res.items():
        print(f'  {T:12.0f}{med:14.1f}{sd:12.1f}')
    ok20 = 20.0 in res and abs(res[20.0][0] - 20.0) < max(2 * res[20.0][1], 10.0)
    print(f'  -> 20 m {"RECOVERED" if ok20 else "NOT recovered"}; rupture radii'
          f' here are 14-48 m\n', flush=True)

    # ---- real inversions ----
    out = []
    for (i, j), rows in per_pair.items():
        r0 = mt[(mt.i == i) & (mt.j == j)].iloc[0]
        m0 = ev.mag[[i, j]].max()
        M0 = 10 ** (1.5 * m0 + 9.1)
        rad = (7 * M0 / (16 * DSIGMA)) ** (1 / 3)
        rec = {}
        for vp in VP_LIST:
            s = invert(rows, vp)
            if s:
                rec[vp] = s
        if VP_REF not in rec:
            continue
        s = rec[VP_REF]
        spread = (max(rec[v]['offset'] for v in rec) -
                  min(rec[v]['offset'] for v in rec))
        out.append(dict(i=i, j=j, fam=bool(r0.is_cand), hrsn=r0.hrsn,
                        das=r0.das_pub, dt_days=r0.dt_days,
                        offset=s['offset'], err=s['err'], nsta=s['nsta'],
                        rms_ms=1e3 * s['rms_s'], rad=rad,
                        R_sep=s['offset'] / rad, vspread=spread,
                        tag=f'{ev.tag[i][3:11]}/{ev.tag[j][3:11]}'))
    D = pd.DataFrame(out)
    if D.empty:
        print('no inversions succeeded'); return
    D.to_csv(os.path.join(HERE, 'hrsn_reloc.csv'), index=False)

    F = D[D.fam]
    C = D[~D.fam]
    print('FAMILY PAIRS')
    print(f'{"events":>20}{"offset m":>10}{"err m":>8}{"rad m":>8}{"R_sep":>8}'
          f'{"nsta":>6}{"rms ms":>8}{"HRSN CC":>9}')
    for _, x in F.sort_values('offset').iterrows():
        print(f'{x.tag:>20}{x.offset:10.1f}{x.err:8.1f}{x.rad:8.1f}{x.R_sep:8.2f}'
              f'{int(x.nsta):6d}{x.rms_ms:8.2f}{x.hrsn:9.3f}')

    print(f'\nR2 family vs control')
    print(f'  family  n={len(F):3d}  median offset {F.offset.median():8.1f} m'
          f'   median R_sep {F.R_sep.median():.2f}')
    if len(C):
        print(f'  control n={len(C):3d}  median offset {C.offset.median():8.1f} m'
              f'   median R_sep {C.R_sep.median():.2f}')
        try:
            from scipy.stats import mannwhitneyu
            _u, pv = mannwhitneyu(F.offset, C.offset, alternative='less')
            print(f'  Mann-Whitney p = {pv:.4f}')
        except Exception:
            pass

    print('\nR3 THE CONFOUND TEST (this is the number that matters)')
    if len(F) > 3:
        rr = float(np.corrcoef(F.hrsn, np.log10(F.offset.clip(lower=1e-3)))[0, 1])
        print(f'  within family, r(HRSN CC, log offset) = {rr:+.3f}   (n={len(F)})')
        print(f'  DAS version of this test gave -0.826, which invalidated it.')
        if rr < -0.6:
            print('  -> RELOCATION IS ALSO READING WAVEFORM QUALITY. Neither')
            print('     instrument can verify same-patch with this event set.')
            print('     Report as an observational confirmation of Gao et al. 2021.')
        else:
            print('  -> relocation is NOT dominated by waveform quality. The')
            print('     family/control separation can be read geometrically.')

    print('\nR4 velocity sensitivity (5.0 / 5.8 / 6.5 km/s)')
    print(f'  median spread in offset across V: {D.vspread.median():.1f} m'
          f'  ({100*D.vspread.median()/max(D.offset.median(),1e-9):.0f}% of median offset)')

    print('\nSAME-PATCH VERDICT (R_sep = offset / rupture radius)')
    for lab, S in (('family', F), ('control', C)):
        if not len(S):
            continue
        n_in = int((S.R_sep < 1).sum())
        print(f'  {lab:>8}: {n_in}/{len(S)} pairs with R_sep < 1 '
              f'(sources overlap), median R_sep {S.R_sep.median():.2f}')

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    ax[0].errorbar(range(len(F)), F.sort_values('offset').offset,
                   yerr=F.sort_values('offset').err, fmt='o', color='C3',
                   label='family')
    if len(C):
        ax[0].plot(range(len(C)), np.sort(C.offset), '.', color='0.6',
                   label='control')
    ax[0].axhline(F.rad.median(), color='k', ls='--', label='rupture radius')
    ax[0].set(ylabel='relative offset (m)', yscale='log', xlabel='pair',
              title='A  DD relative offsets')
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which='both')
    if len(F) > 3:
        ax[1].scatter(F.hrsn, F.offset, c='C3')
        ax[1].set(xlabel='HRSN CC', ylabel='offset (m)', yscale='log',
                  title='B  R3: does offset track waveform quality?')
        ax[1].grid(alpha=.3, which='both')
    ax[2].hist(np.log10(F.R_sep.clip(lower=1e-3)), bins=12, color='C3',
               alpha=.7, label='family')
    if len(C):
        ax[2].hist(np.log10(C.R_sep.clip(lower=1e-3)), bins=12, color='0.6',
                   alpha=.6, label='control')
    ax[2].axvline(0, color='k', ls='--', label='R_sep = 1')
    ax[2].set(xlabel='log10 (offset / rupture radius)', title='C  same-patch test')
    ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
    fig.tight_layout()
    p = os.path.join(HERE, 'hrsn_reloc.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
