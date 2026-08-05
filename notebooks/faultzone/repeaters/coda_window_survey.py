"""
Where does usable coda actually exist? Measure it instead of guessing again.

WHY THIS EXISTS. G3 (dvv_hrsn.py) returned dv/v values of 1-3% that looked like a
seasonal signal and were not one. Coda CC came back at 0.01-0.63 (median 0.22),
which is far below what stretching needs, so the estimator was fitting noise. The
random-pair control confirmed it: pairs sharing NO source returned |dv/v| ~0.8%,
so the ~1.3% from real repeaters was indistinguishable from the noise output.

The cause was my lapse windows. I used 3-8, 5-12 and 8-18 s after catalog origin
for M0-2 events at 3-13 km hypocentral distance. S arrives roughly 2-5 s after
origin at those distances, so the earliest window can contain the S arrival itself
and the latest is ambient noise. G1's "coda SNR > 3 to 10 s lapse" was measured on
the LARGEST events (M1.6-1.9) and amplitude-above-noise is a far weaker condition
than two events' coda correlating with each other.

I have now guessed windows twice. This script measures three things as functions of
lapse time so the window is chosen from data:

  1. CODA ENVELOPE DECAY vs pre-event noise -> where is there signal at all?
  2. PAIR CODA CC vs lapse -> where do two occurrences actually correlate?
  3. dv/v STABILITY vs lapse, with the random-pair null measured in every window
     -> where does the estimator separate repeaters from noise?

The deliverable is a usable lapse range per pair, or a demonstration that none
exists at HRSN for events this small -- which would itself be the answer, and would
push the measurement toward the direct-arrival timing observable instead of coda.

Reads only hrsn_cache/; no downloads.
"""
import os
import sys
import numpy as np
import pandas as pd
from obspy import read
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dvv_core import bandpass, stretch_dvv, bulk_align          # noqa: E402

WF = os.path.join(HERE, 'hrsn_cache')
BAND = (5.0, 20.0)
PRE_S = 5.0
# fine grid of short windows, starting much earlier than G3 used
STARTS = np.arange(0.5, 12.1, 0.5)
WIN_LEN = 2.0
SNR_MIN = 3.0
CC_MIN = 0.70            # what stretching actually needs


def load(tag):
    f = os.path.join(WF, f'{tag}.mseed')
    if not os.path.exists(f):
        return None
    try:
        return read(f)
    except Exception:
        return None


def traces(st):
    out = {}
    if st is None:
        return out
    for tr in st:
        fs = tr.stats.sampling_rate
        d = bandpass(tr.data.astype(float) - float(np.mean(tr.data)), fs, BAND)
        out[tr.stats.station] = (d, fs)
    return out


def envelope_snr(d, fs):
    """RMS in each window divided by pre-event noise RMS."""
    n0, n1 = int(0.3 * fs), int((PRE_S - 0.5) * fs)
    noise = np.sqrt(np.mean(d[n0:n1] ** 2))
    if not np.isfinite(noise) or noise <= 0:
        return np.full(STARTS.size, np.nan)
    out = np.full(STARTS.size, np.nan)
    for k, s0 in enumerate(STARTS):
        i0 = int((PRE_S + s0) * fs)
        i1 = int((PRE_S + s0 + WIN_LEN) * fs)
        if i1 <= d.size:
            out[k] = np.sqrt(np.mean(d[i0:i1] ** 2)) / noise
    return out


def pair_profile(A, B):
    """Per-lapse-window coda CC and dv/v, medianed over common stations."""
    common = sorted(set(A) & set(B))
    cc = np.full((len(common), STARTS.size), np.nan)
    dv = np.full((len(common), STARTS.size), np.nan)
    sn = np.full((len(common), STARTS.size), np.nan)
    for si, s in enumerate(common):
        a, fs = A[s]; b, _ = B[s]
        # ALIGN FIRST. Without this the catalog origin-time error (0.1-0.5 s)
        # drives every windowed correlation to zero -- see bulk_align's docstring.
        b, _lag = bulk_align(a, b, fs, max_lag_s=2.0)
        sn[si] = np.fmin(envelope_snr(a, fs), envelope_snr(b, fs))
        for k, s0 in enumerate(STARTS):
            i0 = int((PRE_S + s0) * fs)
            i1 = int((PRE_S + s0 + WIN_LEN) * fs)
            if i1 > min(a.size, b.size):
                continue
            aw, bw = a[i0:i1], b[i0:i1]
            na = np.sqrt(np.sum(aw ** 2)); nb = np.sqrt(np.sum(bw ** 2))
            if na <= 0 or nb <= 0:
                continue
            # residual lag search inside the window, for the same reason
            r = np.correlate(bw / nb, aw / na, 'full')
            m0 = aw.size - 1
            pd_ = int(0.3 * fs)
            sg = r[max(m0 - pd_, 0):m0 + pd_ + 1]
            cc[si, k] = (float(sg[int(np.argmax(np.abs(sg)))])
                         if sg.size else np.nan)
            d_, c_, _, _ = stretch_dvv(aw, bw, fs, eps_max=0.04, n_eps=121)
            dv[si, k] = d_
    with np.errstate(invalid='ignore'):
        return (np.nanmedian(cc, axis=0), np.nanmedian(dv, axis=0),
                np.nanmedian(sn, axis=0), len(common))


def main():
    R = pd.read_csv(os.path.join(HERE, 'hrsn_control.csv'))
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    R['t_i'] = pd.to_datetime(R['t_i'], utc=True, format='mixed')
    R['t_j'] = pd.to_datetime(R['t_j'], utc=True, format='mixed')
    cand = R[R.is_cand & (R.hrsn > 0.78)]
    ctrl = R[~R.is_cand].head(25)
    print(f'{len(cand)} repeater pairs, {len(ctrl)} control pairs')
    print(f'{STARTS.size} windows of {WIN_LEN} s from {STARTS[0]} to '
          f'{STARTS[-1]} s after origin\n', flush=True)

    CC, DV, SN, lab = [], [], [], []
    CCc, DVc = [], []
    for kind, D in [('rep', cand), ('ctrl', ctrl)]:
        for _, r in D.iterrows():
            sa, sb = load(ev.tag[int(r.i)]), load(ev.tag[int(r.j)])
            if sa is None or sb is None:
                continue
            cc, dv, sn, n = pair_profile(traces(sa), traces(sb))
            if kind == 'rep':
                CC.append(cc); DV.append(dv); SN.append(sn)
                lab.append(f'{r.t_i:%Y-%m-%d}/{r.t_j:%Y-%m-%d}')
                print(f'  rep  {lab[-1]}  {n} sta  '
                      f'peak coda CC {np.nanmax(cc):.3f} at '
                      f'{STARTS[int(np.nanargmax(cc))]:.1f} s', flush=True)
            else:
                CCc.append(cc); DVc.append(dv)
    if not CC:
        print('nothing loaded'); return
    CC = np.array(CC); DV = np.array(DV); SN = np.array(SN)
    CCc = np.array(CCc); DVc = np.array(DVc)

    cc_m = np.nanmedian(CC, axis=0)
    sn_m = np.nanmedian(SN, axis=0)
    ccc_m = np.nanmedian(CCc, axis=0)

    print('\nLAPSE-TIME PROFILE (median over pairs)')
    print(f'{"lapse":>7}{"SNR":>8}{"codaCC":>9}{"ctrlCC":>9}'
          f'{"|dvv|%":>9}{"ctrl|dvv|%":>12}   usable?')
    for k, s0 in enumerate(STARTS):
        dvr = 100 * np.nanmedian(np.abs(DV[:, k]))
        dvc = 100 * np.nanmedian(np.abs(DVc[:, k])) if DVc.size else np.nan
        use = 'YES' if (sn_m[k] > SNR_MIN and cc_m[k] > CC_MIN) else ''
        print(f'{s0:7.1f}{sn_m[k]:8.2f}{cc_m[k]:9.3f}{ccc_m[k]:9.3f}'
              f'{dvr:9.3f}{dvc:12.3f}   {use}')

    good = (sn_m > SNR_MIN) & (cc_m > CC_MIN)
    print('\nVERDICT')
    print(f'  requirement: coda SNR > {SNR_MIN} AND pair coda CC > {CC_MIN}')
    if good.any():
        lo, hi = STARTS[good].min(), STARTS[good].max() + WIN_LEN
        print(f'  -> USABLE lapse range at HRSN: {lo:.1f} to {hi:.1f} s '
              f'after origin')
        print(f'     ({good.sum()} of {STARTS.size} windows qualify)')
        print(f'     G3 used 3-8, 5-12, 8-18 s -- '
              f'{"overlapping" if hi > 3 else "entirely outside"} this range')
    else:
        b = int(np.nanargmax(cc_m))
        print(f'  -> NO window satisfies both. Best is {STARTS[b]:.1f} s with '
              f'coda CC {cc_m[b]:.3f}, SNR {sn_m[b]:.2f}.')
        print(f'     Coda-wave interferometry is NOT available at HRSN for events '
              f'this small.')
        print(f'     The measurement must move to DIRECT-ARRIVAL differential '
              f'timing, which\n     is the plan\'s Observable A and does not need '
              f'coda. DAS may still work:\n     900 channels stack down the noise '
              f'in a way 8 stations cannot.')

    sep = np.nanmedian(np.abs(DV), axis=0) - np.nanmedian(np.abs(DVc), axis=0)
    if np.any(np.isfinite(sep)):
        b = int(np.nanargmax(sep))
        print(f'\n  best repeater-vs-null separation at lapse {STARTS[b]:.1f} s: '
              f'{100*sep[b]:.3f}% (repeaters minus null)')
        print('  note: separation is meaningless where coda CC is below '
              f'{CC_MIN} -- check the table')

    np.savez(os.path.join(HERE, 'coda_window_survey.npz'),
             starts=STARTS, CC=CC, DV=DV, SN=SN, CCc=CCc, DVc=DVc,
             labels=np.array(lab), win_len=WIN_LEN)

    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    a = ax[0]
    for i in range(CC.shape[0]):
        a.plot(STARTS + WIN_LEN / 2, SN[i], color='C0', alpha=0.3, lw=0.9)
    a.plot(STARTS + WIN_LEN / 2, sn_m, 'C0-', lw=2.2, label='median')
    a.axhline(SNR_MIN, color='C3', ls='--', label=f'SNR = {SNR_MIN}')
    a.set(yscale='log', xlabel='lapse time (s after origin)',
          ylabel='coda RMS / pre-event noise RMS',
          title='A  Is there coda at all?')
    a.legend(fontsize=8); a.grid(alpha=0.3)

    a = ax[1]
    for i in range(CC.shape[0]):
        a.plot(STARTS + WIN_LEN / 2, CC[i], color='C2', alpha=0.35, lw=0.9)
    a.plot(STARTS + WIN_LEN / 2, cc_m, 'C2-', lw=2.2, label='repeaters (median)')
    if CCc.size:
        a.plot(STARTS + WIN_LEN / 2, ccc_m, '0.5', lw=2, ls='--',
               label='random pairs')
    a.axhline(CC_MIN, color='C3', ls='--', label=f'CWI needs {CC_MIN}')
    a.axhspan(3, 18, color='C1', alpha=0.12)
    a.text(10.5, 0.9, 'windows G3 used', fontsize=7, color='C1')
    a.set(xlabel='lapse time (s after origin)', ylabel='pair coda CC',
          title='B  Do the two occurrences correlate?')
    a.legend(fontsize=8); a.grid(alpha=0.3)

    a = ax[2]
    a.plot(STARTS + WIN_LEN / 2, 100 * np.nanmedian(np.abs(DV), axis=0),
           'C2-', lw=2, label='repeaters |dv/v|')
    if DVc.size:
        a.plot(STARTS + WIN_LEN / 2, 100 * np.nanmedian(np.abs(DVc), axis=0),
               '0.5', ls='--', lw=2, label='random pairs |dv/v| (noise floor)')
    a.set(xlabel='lapse time (s after origin)', ylabel='|dv/v| (%)',
          title='C  Does the estimator separate\nsignal from noise anywhere?')
    a.legend(fontsize=8); a.grid(alpha=0.3)

    fig.suptitle('Where does usable coda exist? Measuring the lapse window instead '
                 'of assuming it', fontsize=12)
    fig.tight_layout()
    p = os.path.join(HERE, 'coda_window_survey.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
