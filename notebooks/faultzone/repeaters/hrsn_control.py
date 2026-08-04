"""
CONTROL: correlate the same event pairs on HRSN borehole seismometers.

The question this settles. My DAS all-pairs search gives a median CC of 0.081 and
nothing above 0.85 at long time separation, over 14 months at Parkfield -- where
Schaff & Waldhauser report >40% of seismicity is repeating. Either my DAS
processing is destroying the discrimination, or these events genuinely are not
repeaters. Those have opposite remedies and I cannot tell them apart from DAS
alone.

HRSN (network BP) is the 13-station borehole array the Parkfield repeater
catalogs were actually built on. It records the same earthquakes. So:

  * if HRSN gives CC ~0.95 for a pair where DAS gives 0.79
        -> the DAS processing is the problem, and the fault is mine
  * if HRSN also gives ~0.79 for that pair
        -> those events are not repeaters, and the DAS is telling the truth

Random pairs are included alongside the candidates so the HRSN null is measured
on the same footing rather than assumed.

Processing is deliberately matched to the DAS chain: same band, same window
relative to catalog origin, same lag-searched normalised cross-correlation.
"""
import os
import numpy as np
import pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from scipy.signal import butter, sosfiltfilt, correlate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
WF = os.path.join(HERE, 'hrsn_cache')
os.makedirs(WF, exist_ok=True)

BAND = (5.0, 20.0)
PRE_S, POST_S = 5.0, 20.0
WIN = (-1.0, 10.0)          # same window as the DAS correlation
MAX_LAG_S = 2.0
N_RANDOM = 40               # random pairs, to measure the HRSN null
STATIONS = 'CCRB,EADB,FROB,GHIB,JCNB,LCCB,MMNB,RMNB,SCYB,SMNB,VARB,VCAB'


def fetch(client, t0, tag):
    """Download and cache one event's HRSN waveforms."""
    f = os.path.join(WF, f'{tag}.mseed')
    if os.path.exists(f):
        try:
            from obspy import read
            return read(f)
        except Exception:
            os.remove(f)
    try:
        st = client.get_waveforms('BP', STATIONS, '*', 'DP1',
                                  UTCDateTime(t0) - PRE_S,
                                  UTCDateTime(t0) + POST_S)
    except Exception as e:
        print(f'   {tag}: fetch failed ({str(e)[:60]})', flush=True)
        return None
    if len(st) == 0:
        return None
    st.write(f, format='MSEED')
    return st


def traces(st, t0):
    """Per-station normalised trace over the analysis window."""
    out = {}
    if st is None:
        return out
    for tr in st:
        try:
            x = tr.copy()
            x.detrend('demean')
            fs = x.stats.sampling_rate
            sos = butter(4, list(BAND), btype='band', fs=fs, output='sos')
            d = sosfiltfilt(sos, x.data.astype(np.float64))
            i0 = int((PRE_S + WIN[0]) * fs)
            i1 = int((PRE_S + WIN[1]) * fs)
            if i1 > d.size:
                continue
            seg = d[max(i0, 0):i1]
            seg = seg - seg.mean()
            n = np.sqrt(np.sum(seg ** 2))
            if n > 0 and np.isfinite(n):
                out[tr.stats.station] = (seg / n, fs)
        except Exception:
            continue
    return out


def cc_pair(A, B):
    """Median lag-searched CC across stations common to both events."""
    common = sorted(set(A) & set(B))
    vals = []
    for s in common:
        a, fs = A[s]
        b, _ = B[s]
        n = min(a.size, b.size)
        a, b = a[:n], b[:n]
        cc = correlate(b, a, mode='full', method='fft')
        mid = n - 1
        pad = int(MAX_LAG_S * fs)
        seg = cc[max(mid - pad, 0):mid + pad + 1]
        if seg.size == 0:
            continue
        vals.append(seg[int(np.argmax(np.abs(seg)))])
    if not vals:
        return np.nan, 0
    return float(np.median(vals)), len(vals)


def main():
    d = np.load(os.path.join(HERE, 'correlate_all.npz'), allow_pickle=True)
    C = d['C']
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    ev['time'] = pd.to_datetime(ev['time'], utc=True, format='mixed')
    N = len(ev)
    iu = np.triu_indices(N, 1)
    cc_das = C[iu]
    dt_d = np.array([abs((ev.time[j] - ev.time[i]).total_seconds()) / 86400
                     for i, j in zip(*iu)])

    # candidates: highest DAS CC with a real time baseline, plus random controls
    cand = [k for k in np.argsort(cc_das)[::-1] if dt_d[k] > 20][:10]
    rng = np.random.default_rng(0)
    rand = rng.choice(np.where(dt_d > 20)[0], size=N_RANDOM, replace=False)
    sel = list(dict.fromkeys(list(cand) + list(rand)))
    need = sorted({i for k in sel for i in (iu[0][k], iu[1][k])})
    print(f'{len(sel)} pairs ({len(cand)} candidates + random controls), '
          f'{len(need)} events to fetch\n', flush=True)

    client = Client('NCEDC', timeout=120)
    store = {}
    for i in need:
        e = ev.iloc[i]
        st = fetch(client, e['time'], e['tag'])
        T = traces(st, e['time'])
        store[i] = T
        print(f'  {e["tag"]}  M{e["mag"]:.2f}  {len(T)} stations', flush=True)

    rows = []
    for k in sel:
        i, j = iu[0][k], iu[1][k]
        h, nst = cc_pair(store.get(i, {}), store.get(j, {}))
        rows.append(dict(i=i, j=j, das=cc_das[k], hrsn=h, nsta=nst,
                         dt_days=dt_d[k], is_cand=k in cand,
                         t_i=ev.time[i], t_j=ev.time[j],
                         m_i=ev.mag[i], m_j=ev.mag[j]))
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(HERE, 'hrsn_control.csv'), index=False)

    ok = R.dropna(subset=['hrsn'])
    print(f'\n{len(ok)}/{len(R)} pairs with HRSN data\n')
    print('CANDIDATES (highest DAS CC):')
    print(f'{"DAS":>7}{"HRSN":>8}{"nsta":>6}{"dt_d":>8}   events')
    for _, r in ok[ok.is_cand].sort_values('das', ascending=False).iterrows():
        print(f'{r.das:7.3f}{r.hrsn:8.3f}{int(r.nsta):6d}{r.dt_days:8.1f}   '
              f'{r.t_i:%Y-%m-%d} M{r.m_i:.2f} / {r.t_j:%Y-%m-%d} M{r.m_j:.2f}')

    ctrl = ok[~ok.is_cand]
    print(f'\nRANDOM CONTROLS ({len(ctrl)} pairs):')
    print(f'  DAS  median {ctrl.das.median():.3f}  max {ctrl.das.max():.3f}')
    print(f'  HRSN median {ctrl.hrsn.median():.3f}  max {ctrl.hrsn.max():.3f}')

    print('\nVERDICT:')
    if len(ok[ok.is_cand]):
        cd, ch = ok[ok.is_cand].das.max(), ok[ok.is_cand].hrsn.max()
        print(f'  best candidate: DAS {cd:.3f} vs HRSN {ch:.3f}')
        if ch > cd + 0.10:
            print('  -> HRSN resolves these pairs BETTER than DAS. The DAS '
                  'processing is\n     losing discrimination; the fault is in my '
                  'chain, not the data.')
        elif abs(ch - cd) <= 0.10:
            print('  -> HRSN and DAS agree. These events are genuinely not '
                  'repeaters,\n     and the DAS is reporting the truth.')
        else:
            print('  -> DAS exceeds HRSN, which is unexpected; check the '
                  'station set and band.')

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].scatter(ctrl.das, ctrl.hrsn, s=18, c='0.6', label='random pairs')
    ax[0].scatter(ok[ok.is_cand].das, ok[ok.is_cand].hrsn, s=55, c='C3',
                  label='DAS candidates')
    lim = [-0.2, 1.0]
    ax[0].plot(lim, lim, 'k--', lw=1)
    ax[0].set(xlim=lim, ylim=lim, xlabel='DAS CC', ylabel='HRSN CC',
              title='A  Same pairs, two instruments')
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)

    ax[1].hist(ctrl.das, bins=20, alpha=0.6, label='DAS, random')
    ax[1].hist(ctrl.hrsn, bins=20, alpha=0.6, label='HRSN, random')
    ax[1].set(xlabel='CC', ylabel='pairs', title='B  Null on each instrument')
    ax[1].legend(fontsize=8)
    fig.suptitle('Control: is low DAS correlation the processing or the events?',
                 fontsize=12)
    fig.tight_layout()
    p = os.path.join(HERE, 'hrsn_control.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
