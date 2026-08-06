"""
Step 1 of the published Parkfield procedure: measure Nadeau's similarity beta.

WHAT BETA IS, from the source rather than paraphrase.

  Nadeau, Foxall & McEvilly 1995 (Science 267:503):
    "The similarity measure, beta, ... is based on a network-wide characterization
     of maximum cross-correlation coefficient values for P and S waves between
     pairs of earthquakes"
    "We selected a similarity measure of beta >= 0.98 as our criterion for defining
     clusters, a value at which the defined cluster population changed little"
    "The gap in the range 0.6 < beta < 0.9, 200 m < offset < 500 m generally
     separates highly similar clustered from nonclustered behavior"

  Nadeau & Johnson 1998 (BSSA 88:790):
    "correlation coefficients ... generally exceed a value of 0.98, and by virtual
     collocation of the events, their quasi-periodic recurrence, and their nearly
     identical magnitudes"

So beta = per-station maximum cross-correlation, computed separately on the P and
the S window, then aggregated across the network. Implemented here as the median
over stations of the per-station value, with the mean and station count also
reported so the aggregation choice can be revisited.

WHY THIS AND NOT A LOCATION SCREEN. Nadeau & Johnson state plainly that waveform
similarity is "more discriminating for broadband borehole recordings than typical
selection methods based on event size, recurrence times, or locations." This
project has now twice led with a location screen, and Nadeau 1995 explains why that
fails: for the beta > 0.9 population, ROUTINE hypocentre locations scatter over up
to 200 m, while waveform-based relative relocation puts the same events within
10-20 m. Separations computed from the routine DDRT catalog are therefore
uninformative about co-location -- measured here at 167-771 m for pairs that
correlate at 0.92-0.99.

WHAT THIS SCRIPT DOES NOT DO. beta >= 0.98 defines a CLUSTER, not a sequence.
Nadeau 1995: "For clusters of three or more events, 80 to 90% were complex in that
their member events could be further subdivided into subgroups (each containing a
different event type) on the basis of subtle differences in high-frequency
waveforms." Subdivision is step 2, in a separate script.

Instrument: HRSN (network BP), the borehole array the Parkfield repeater catalogs
were built on, so beta here is directly comparable to the published values.
"""
import os
import sys
import numpy as np
import pandas as pd
from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from scipy.signal import butter, sosfiltfilt, correlate as sp_correlate

HERE = os.path.dirname(os.path.abspath(__file__))
WF = os.path.join(HERE, 'hrsn3c_cache')     # three-component cache
CHANS = 'DP1,DP2,DP3'

# WHY THREE COMPONENTS. Nadeau's beta is "a network-wide characterization of
# maximum cross-correlation coefficient values for P and S waves". A single
# component discards two thirds of the particle motion and caps the achievable
# correlation: the DP1-only run reached max beta = 0.9737 and put ZERO pairs at
# the published 0.98 threshold, while reproducing the correct RANKING (the eight
# pairs confirmed earlier ranked 1,2,3,4,7,32,51 of 21461). That pattern -- right
# order, depressed values -- is what an information deficit looks like, not an
# absence of repeaters.
os.makedirs(WF, exist_ok=True)

BAND = (1.5, 15.0)        # Schaff & Waldhauser 2005 NCSN convention
PRE_S, POST_S = 5.0, 20.0
VP, VS = 5.5, 3.2         # km/s, for windowing only -- a lag search follows
P_WIN = (-0.3, 1.2)       # s about the predicted P
S_WIN = (-0.3, 2.0)       # s about the predicted S
BULK_LAG_S = 2.0          # one alignment per station on the full record
RESID_LAG_S = 0.10        # residual search inside the phase window only

# WHY TWO LAG SCALES. Catalog origin times carry 0.1-0.5 s of error, so some lag
# search is mandatory -- correlating at zero lag drives identical waveforms to
# zero. But searching +/-1 s INSIDE a 1.5 s window maximises over ~30 independent
# lags of a ~45-DOF correlation, and the expected maximum of that is ~0.4. The
# first run of this script did exactly that and produced a null with median
# beta = 0.427, which swamped Nadeau's 0.6-0.9 gap. Aligning once on the full
# record removes the origin-time error, after which a +/-0.1 s residual search is
# enough and the null falls back toward zero.
MIN_STA = 3               # a beta from fewer stations is not "network-wide"


def fetch(client, t0, tag, stations):
    f = os.path.join(WF, f'{tag}.mseed')
    if os.path.exists(f):
        try:
            return read(f)
        except Exception:
            os.remove(f)
    try:
        st = client.get_waveforms('BP', stations, '*', CHANS,
                                  UTCDateTime(t0) - PRE_S,
                                  UTCDateTime(t0) + POST_S)
    except Exception:
        return None
    if len(st) == 0:
        return None
    st.write(f, format='MSEED')
    return st


def bulk_shift(a, b, fs, max_lag_s=BULK_LAG_S):
    """Integer-sample shift aligning b to a on the whole record."""
    n = min(a.size, b.size)
    aa, bb = a[:n] - a[:n].mean(), b[:n] - b[:n].mean()
    if not np.any(aa) or not np.any(bb):
        return 0
    cc = sp_correlate(bb, aa, 'full', method='fft')   # O(n log n)
    mid, pad = n - 1, int(max_lag_s * fs)
    seg = cc[max(mid - pad, 0):mid + pad + 1]
    if seg.size == 0:
        return 0
    return int(np.argmax(np.abs(seg))) + max(mid - pad, 0) - mid


def phase_windows(st, ev, sta_xy):
    """Per-station normalised P and S windows, FFT-ready.

    Windows are placed from predicted travel times using constant VP/VS. That is
    crude, but a +/-1 s lag search follows and the window is 1.5-2.3 s long, so the
    phase cannot fall outside it for these distances. Precision comes from the lag
    search, not from the prediction.
    """
    out = {}
    if st is None:
        return out
    for tr in st:
        s = f'{tr.stats.station}.{tr.stats.channel}'
        base = tr.stats.station
        if base not in sta_xy:
            continue
        slat, slon, sdep = sta_xy[base]
        d, _, _ = gps2dist_azimuth(ev['lat'], ev['lon'], slat, slon)
        hyp = np.sqrt((d / 1000.0) ** 2 + (ev['depth'] + sdep / 1000.0) ** 2)
        tp, ts = hyp / VP, hyp / VS
        fs = tr.stats.sampling_rate
        try:
            x = tr.data.astype(float)
            x -= x.mean()
            sos = butter(4, list(BAND), btype='band', fs=fs, output='sos')
            x = sosfiltfilt(sos, x)
        except Exception:
            continue
        w = {}
        for name, (t_pred, win) in (('P', (tp, P_WIN)), ('S', (ts, S_WIN))):
            i0 = int((PRE_S + t_pred + win[0]) * fs)
            i1 = int((PRE_S + t_pred + win[1]) * fs)
            if i0 < 0 or i1 > x.size or i1 - i0 < 32:
                continue
            w[name] = (i0, i1)          # bounds only; cut later, after shifting
        if w:
            w['_full'] = (x, fs)
            out[s] = w
    return out


def cut(x, i0, i1, shift=0):
    """Extract and normalise [i0,i1) with the index bounds shifted.

    The bulk lag is applied to the INDICES, not by rolling an already-cut window.
    Rolling wraps: a 0.5 s origin-time error is 125 samples at 250 Hz against a
    375-sample window, so np.roll would fold a third of the window in from the
    opposite end and quietly corrupt the correlation.
    """
    a, b = i0 + shift, i1 + shift
    if a < 0 or b > x.size or b - a < 32:
        return None
    seg = x[a:b] - x[a:b].mean()
    n = np.sqrt(np.sum(seg ** 2))
    return seg / n if (np.isfinite(n) and n > 0) else None


def maxcc(a, b, fs, lag_s=RESID_LAG_S):
    """Maximum normalised cross-correlation over a small residual lag."""
    n = min(a.size, b.size)
    a, b = a[:n], b[:n]
    cc = sp_correlate(b, a, 'full', method='fft')
    mid, pad = n - 1, max(int(lag_s * fs), 1)
    seg = cc[max(mid - pad, 0):mid + pad + 1]
    return float(np.max(np.abs(seg))) if seg.size else np.nan


def beta(A, B):
    """Network-wide beta: median over stations of the per-station max-CC,
    computed on P and S separately then combined per station."""
    percomp = {}
    for s in sorted(set(A) & set(B)):
        # align this station's pair once on the full record, then window
        k = 0
        if '_full' in A[s] and '_full' in B[s]:
            fa, fs0 = A[s]['_full']; fb, _ = B[s]['_full']
            k = bulk_shift(fa, fb, fs0)
        xa, fs = A[s]['_full']; xb, _ = B[s]['_full']
        ph = []
        for name in ('P', 'S'):
            if name in A[s] and name in B[s]:
                ia0, ia1 = A[s][name]; ib0, _ = B[s][name]
                L = min(ia1 - ia0, B[s][name][1] - ib0)
                wa = cut(xa, ia0, ia0 + L)
                wb = cut(xb, ib0, ib0 + L, shift=k)
                if wa is None or wb is None:
                    continue
                c = maxcc(wa, wb, fs)
                if np.isfinite(c):
                    ph.append(c)
        if ph:
            percomp.setdefault(s.split('.')[0], []).append(float(np.mean(ph)))
    # average the components of a station, then take the network median
    vals = [float(np.mean(v)) for v in percomp.values()]
    if len(vals) < MIN_STA:
        return np.nan, np.nan, len(vals)
    return float(np.median(vals)), float(np.mean(vals)), len(vals)


def main():
    ev = pd.read_csv(os.path.join(HERE, 'phaseA_events.csv'))
    ev['time'] = pd.to_datetime(ev.time, utc=True, format='mixed')
    ev = ev[ev.cov_full].reset_index(drop=True)
    print(f'{len(ev)} DAS-covered DDRT events\n', flush=True)

    client = Client('NCEDC', timeout=180)
    inv = client.get_stations(network='BP', channel='DP1', level='channel',
                              starttime='2024-01-01', endtime='2026-06-01')
    sta_xy = {s.code: (s.latitude, s.longitude,
                       float(np.median([c.depth for c in s])))
              for n in inv for s in n}
    stations = ','.join(sorted(sta_xy))
    print(f'HRSN stations: {len(sta_xy)}\n', flush=True)

    W, ok = {}, []
    for k, e in ev.iterrows():
        st = fetch(client, e.time, e.tag, stations)
        w = phase_windows(st, e, sta_xy)
        if len(w) >= MIN_STA:
            W[k] = w; ok.append(k)
        if (k + 1) % 25 == 0:
            print(f'  {k+1}/{len(ev)} fetched, {len(ok)} usable', flush=True)
    print(f'\n{len(ok)} events with >= {MIN_STA} stations\n', flush=True)

    rows = []
    for ii, i in enumerate(ok):
        for j in ok[ii + 1:]:
            if (abs((ev.time[j] - ev.time[i]).total_seconds()) < 2.0):
                continue      # catalog double-listing, not two earthquakes
            b_med, b_mean, nst = beta(W[i], W[j])
            if not np.isfinite(b_med):
                continue
            rows.append(dict(i=i, j=j, t_i=ev.time[i], t_j=ev.time[j],
                             m_i=ev.mag[i], m_j=ev.mag[j],
                             beta=b_med, beta_mean=b_mean, nsta=nst,
                             days=abs((ev.time[j] - ev.time[i]).total_seconds())
                             / 86400,
                             dmag=abs(ev.mag[j] - ev.mag[i]),
                             sep_km=np.sqrt(
                                 ((ev.lat[j]-ev.lat[i])*111.19)**2 +
                                 ((ev.lon[j]-ev.lon[i])*111.19 *
                                  np.cos(np.radians(36)))**2 +
                                 (ev.depth[j]-ev.depth[i])**2)))
        if (ii + 1) % 25 == 0:
            print(f'  beta: {ii+1}/{len(ok)} rows done, {len(rows)} pairs',
                  flush=True)
    B = pd.DataFrame(rows)
    B.to_csv(os.path.join(HERE, 'beta_similarity.csv'), index=False)
    print(f'\n{len(B)} pairs with beta\n')

    print('BETA DISTRIBUTION')
    for q in [50, 90, 99, 99.9]:
        print(f'  {q:5.1f}th pct: {np.percentile(B.beta, q):.4f}')
    print(f'  max: {B.beta.max():.4f}')
    for t in [0.6, 0.7, 0.8, 0.9, 0.95, 0.98]:
        print(f'  beta > {t:.2f}: {int((B.beta > t).sum()):6d} pairs '
              f'({100*(B.beta > t).mean():6.3f}%)')

    print('\nIS NADEAU\'S GAP PRESENT IN OUR DATA?')
    print('  (1995: "the gap in the range 0.6 < beta < 0.9 ... generally separates')
    print('   highly similar clustered from nonclustered behavior")')
    h, edges = np.histogram(B.beta, bins=np.arange(0, 1.02, 0.02))
    for k in range(len(h)):
        lo = edges[k]
        if lo < 0.3:
            continue
        bar = '#' * int(60 * h[k] / max(h.max(), 1))
        print(f'  {lo:.2f}-{edges[k+1]:.2f} {h[k]:6d} {bar}')

    print('\nPAIRS AT NADEAU\'S CLUSTER THRESHOLD (beta >= 0.98)')
    c = B[B.beta >= 0.98].sort_values('beta', ascending=False)
    if len(c):
        print(f'{"events":<26}{"beta":>7}{"nsta":>6}{"days":>8}{"dM":>6}'
              f'{"sep_m":>8}')
        for _, r in c.head(40).iterrows():
            print(f'{r.t_i:%Y-%m-%d}/{r.t_j:%Y-%m-%d}{r.beta:7.4f}'
                  f'{int(r.nsta):6d}{r.days:8.0f}{r.dmag:6.2f}'
                  f'{1000*r.sep_km:8.0f}')
    print(f'\n  {len(c)} pairs at beta >= 0.98')
    print(f'  involving {len(set(c.i) | set(c.j))} distinct events')


if __name__ == '__main__':
    main()
