"""
Does moveout correction before channel-stacking recover the DAS/HRSN CC deficit?

THE CONTRADICTION THIS SETTLES. correlate_all.py:61 asserts the earthquake arrival
is FLAT across the 900 channels, so the earthquake is the common mode. That is
backwards. For a VERTICAL fiber, a near-vertically incident wave propagates ALONG
the fiber, so the apparent velocity equals the formation velocity -- the steepest
moveout the geometry can produce, ~919/3000 = 0.31 s end to end. Broadside
incidence is the case that gives a flat arrival. My own nano_velocity_scan.py
measured 2975 m/s by slant-stacking this very fiber, and Lellouch et al. (2019,
doi:10.1785/0120190176) DETECT earthquakes on it from the moveout slope. A flat
arrival would make both impossible.

WHY IT MATTERS QUANTITATIVELY. Stacking N traces whose arrivals are spread over a
delay T is convolution with a boxcar of width T: coherent signal survives as
|sinc(pi f T)| while incoherent noise falls as 1/sqrt(N). At 12.5 Hz with
T = 0.31 s, sinc = 0.043 and 1/sqrt(700) = 0.038, so the 700-channel stack buys
1.1 dB over a SINGLE channel. Correct the moveout and the same stack buys
sqrt(700) = 28.5 dB. The predicted recovery is ~27 dB, which is the entire
measured DAS/HRSN deficit (7-17 dB, mean 13 dB).

A SECOND SYMPTOM, SAME CAUSE. The DAS random-pair null sits at +0.124 while HRSN's
sits at -0.004 under identical processing. For independent signals CC must be
zero-mean, so +0.124 means every DAS window shares a common component. Boxcar
smearing of width T imposes nearly the same narrowband ringing on every event, so
every pair correlates positively. If moveout is the cause, correcting it must pull
the null back toward zero. That is an independent prediction, not a restatement.

FALSIFIABLE, AND SET UP TO FAIL LOUDLY. The 7 pairs confirmed on HRSN (CC
0.917-0.995) are ground truth. If alignment moves DAS from 0.63-0.79 toward HRSN
parity, the diagnosis holds and the deficit was mine. If DAS does not move, the
diagnosis is wrong, the deficit is a real coherence ceiling, and I need to say so.
The slowness scan INCLUDES p = 0, so the data decides between flat and dipping
rather than my argument deciding.

Alignment cannot manufacture similarity: it is a rigid per-channel time shift set
by one scalar per event, with no reference to the other event of the pair.

Design: four variants isolate the two factors independently, so a change can be
attributed to alignment or to the window rather than to both at once.
  flat/origin window   -- reproduces the published 0.63-0.79
  aligned/origin window -- alignment alone
  flat/P window        -- window alone
  aligned/P window     -- both (the literature-standard configuration)
"""
import os
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt, correlate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache_all')

CH_LO, CH_HI = 100, 800
DX = 1.0210                  # channel spacing, m (900 ch over 919 m)
BAND = (5.0, 20.0)
PRE_S = 5.0                  # cached window starts 5 s before catalog origin

# slowness scan. +/-6.7e-4 s/m is |V| > 1500 m/s, covering body waves and the
# tube wave, both signs (the wave arrives from below, but the channel->depth
# orientation is not assumed here). p = 0 is included so "flat" is a candidate
# the data can select rather than an assumption.
P_MAX = 6.7e-4
N_P = 81
SEM_WIN = 0.5                # semblance window, s
T_SEARCH = (5.0, 9.0)        # search for P in this part of the cached window
V_BODY = 2000.0              # |V| above this counts as a body wave, not tube

WIN_ORIGIN = (-1.0, 10.0)    # relative to catalog origin (the old choice)
WIN_P = (-0.2, 1.8)          # relative to the picked P (Waldhauser/Schaff style)
MAX_LAG_S = 2.0


def prep(tag):
    """Load one event, bandpass, drop dead channels, L2-normalise per channel."""
    f = os.path.join(CACHE, f'{tag}.npz')
    if not os.path.exists(f):
        return None, None
    d = np.load(f)
    X, fs = d['X'], float(d['fs'])
    sos = butter(4, list(BAND), btype='band', fs=fs, output='sos')
    A = sosfiltfilt(sos, X[CH_LO:CH_HI].astype(np.float64), axis=-1)
    A -= A.mean(axis=1, keepdims=True)
    rms = np.sqrt((A ** 2).mean(axis=1))
    good = np.isfinite(rms) & (rms > 0)
    # drop channels whose RMS is a wild outlier: dead fibre or a coupling glitch
    if good.sum() > 20:
        med = np.median(rms[good])
        good &= (rms > 0.05 * med) & (rms < 20.0 * med)
    A = A[good] / rms[good, None]          # per-channel L2 (Lellouch step 1)
    z = (np.arange(CH_LO, CH_HI)[good] - CH_LO) * DX
    return (A, z, fs), int(good.sum())


def slant_scan(A, z, fs):
    """2-D (slowness, time) semblance scan. Returns best p, best t, semblance."""
    n = A.shape[1]
    F = np.fft.rfft(A, axis=-1)
    freq = np.fft.rfftfreq(n, 1.0 / fs)
    ps = np.linspace(-P_MAX, P_MAX, N_P)
    w = int(SEM_WIN * fs)
    i0, i1 = int(T_SEARCH[0] * fs), int(T_SEARCH[1] * fs)
    N = A.shape[0]

    best = (0.0, 0.0, -1.0)
    S = np.zeros((N_P, i1 - i0))
    for k, p in enumerate(ps):
        # exact sub-sample alignment by FFT phase shift; integer shifts at 100 Hz
        # would leave 72 deg of phase error at 20 Hz
        ph = np.exp(2j * np.pi * np.outer(p * z, freq))
        B = np.fft.irfft(F * ph, n=n, axis=-1)
        num = np.cumsum(B.sum(axis=0) ** 2)
        den = np.cumsum((B ** 2).sum(axis=0))
        sem = ((num[w:] - num[:-w]) / (N * (den[w:] - den[:-w] + 1e-300)))
        seg = sem[max(i0 - 0, 0):i1]
        if seg.size < i1 - i0:
            seg = np.pad(seg, (0, i1 - i0 - seg.size))
        S[k] = seg
        j = int(np.argmax(seg))
        if seg[j] > best[2]:
            best = (p, (i0 + j) / fs, float(seg[j]))
    return best[0], best[1], best[2], S, ps


def align(A, z, fs, p):
    n = A.shape[1]
    F = np.fft.rfft(A, axis=-1)
    freq = np.fft.rfftfreq(n, 1.0 / fs)
    return np.fft.irfft(F * np.exp(2j * np.pi * np.outer(p * z, freq)),
                        n=n, axis=-1)


def stack(A, fs, i0, i1):
    if i0 < 0 or i1 > A.shape[1]:
        return None
    s = A[:, i0:i1].mean(axis=0)
    s = s - s.mean()
    nn = np.sqrt(np.sum(s ** 2))
    return s / nn if (np.isfinite(nn) and nn > 0) else None


def cc(a, b, fs):
    if a is None or b is None:
        return np.nan
    n = min(a.size, b.size)
    a, b = a[:n], b[:n]
    c = correlate(b, a, mode='full', method='fft')
    mid, pad = n - 1, int(MAX_LAG_S * fs)
    seg = c[max(mid - pad, 0):mid + pad + 1]
    return float(seg[int(np.argmax(np.abs(seg)))]) if seg.size else np.nan


def main():
    R = pd.read_csv(os.path.join(HERE, 'hrsn_control.csv'))
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    need = sorted(set(R.i) | set(R.j))
    print(f'{len(R)} pairs from the HRSN control, {len(need)} events\n', flush=True)

    E, diag = {}, None
    for c, i in enumerate(need):
        tag = ev.tag[i]
        got, ngood = prep(tag)
        if got is None:
            print(f'  {tag}: not cached', flush=True)
            continue
        A, z, fs = got
        p, tp, sem, S, ps = slant_scan(A, z, fs)
        V = np.inf if p == 0 else 1.0 / p
        B = align(A, z, fs, p)
        o0 = int((PRE_S + WIN_ORIGIN[0]) * fs)
        o1 = int((PRE_S + WIN_ORIGIN[1]) * fs)
        q0 = int((tp + WIN_P[0]) * fs)
        q1 = int((tp + WIN_P[1]) * fs)
        E[i] = dict(fs=fs, p=p, V=V, tp=tp, sem=sem, ngood=ngood,
                    flat_o=stack(A, fs, o0, o1), algn_o=stack(B, fs, o0, o1),
                    flat_p=stack(A, fs, q0, q1), algn_p=stack(B, fs, q0, q1))
        print(f'  {tag}  ch={ngood:3d}  V_app={V:+9.0f} m/s  '
              f'tP={tp - PRE_S:+5.2f}s  semblance={sem:.3f}', flush=True)
        if diag is None and i in set(R[R.is_cand].i):
            diag = dict(tag=tag, A=A, B=B, z=z, fs=fs, S=S, ps=ps, tp=tp, p=p)

    vs = np.array([e['V'] for e in E.values() if np.isfinite(e['V'])])
    sems = np.array([e['sem'] for e in E.values()])
    body = np.abs(vs) > V_BODY
    print(f'\nMOVEOUT: is the arrival flat or dipping?')
    print(f'  events with a finite fitted slowness: {vs.size}/{len(E)}  '
          f'(p = 0 exactly would mean flat)')
    print(f'  |V_app| > {V_BODY:.0f} m/s (body wave): {int(body.sum())}/{vs.size}')
    if body.any():
        print(f'  body-wave |V_app|: median {np.median(np.abs(vs[body])):.0f} m/s, '
              f'IQR {np.percentile(np.abs(vs[body]), 25):.0f}-'
              f'{np.percentile(np.abs(vs[body]), 75):.0f}')
        T = 919.0 / np.median(np.abs(vs[body]))
        f0 = np.mean(BAND)
        sinc = abs(np.sinc(f0 * T))
        print(f'  => moveout across 919 m: {T*1000:.0f} ms; at {f0:.1f} Hz the '
              f'un-corrected stack retains {100*sinc:.1f}% of signal amplitude')
        print(f'  => predicted gain from alignment: '
              f'{20*np.log10(1.0/max(sinc, 1e-9)):.1f} dB')
    print(f'  semblance: median {np.median(sems):.3f}, max {sems.max():.3f}')

    rows = []
    for _, r in R.iterrows():
        a, b = E.get(int(r.i)), E.get(int(r.j))
        if a is None or b is None:
            continue
        fs = a['fs']
        rows.append(dict(
            i=int(r.i), j=int(r.j), is_cand=bool(r.is_cand),
            dt_days=r.dt_days, hrsn=r.hrsn, das_pub=r.das,
            flat_o=cc(a['flat_o'], b['flat_o'], fs),
            algn_o=cc(a['algn_o'], b['algn_o'], fs),
            flat_p=cc(a['flat_p'], b['flat_p'], fs),
            algn_p=cc(a['algn_p'], b['algn_p'], fs),
            dV=abs(a['V'] - b['V']) if np.isfinite(a['V'] + b['V']) else np.nan))
    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(HERE, 'moveout_test.csv'), index=False)

    cd = D[D.is_cand].sort_values('hrsn', ascending=False)
    ct = D[~D.is_cand]
    print(f'\nCONFIRMED REPEATERS ({len(cd)} pairs, HRSN is ground truth)')
    print(f'{"HRSN":>7}{"flat/org":>10}{"algn/org":>10}{"flat/P":>9}'
          f'{"algn/P":>9}{"dt_d":>8}')
    for _, r in cd.iterrows():
        print(f'{r.hrsn:7.3f}{r.flat_o:10.3f}{r.algn_o:10.3f}'
              f'{r.flat_p:9.3f}{r.algn_p:9.3f}{r.dt_days:8.0f}')

    print(f'\nNULL ({len(ct)} random pairs) -- must move toward 0.000 if boxcar '
          f'smearing caused the +0.124 bias')
    for c_ in ['flat_o', 'algn_o', 'flat_p', 'algn_p']:
        print(f'  {c_:>8}: median {ct[c_].median():+.3f}   '
              f'max {ct[c_].max():.3f}   sd {ct[c_].std():.3f}')

    print('\nSEPARATION (repeaters vs their own null -- the figure of merit):')
    for c_ in ['flat_o', 'algn_o', 'flat_p', 'algn_p']:
        mad = 1.4826 * np.median(np.abs(ct[c_] - ct[c_].median()))
        z_ = (cd[c_].median() - ct[c_].median()) / mad if mad > 0 else np.nan
        print(f'  {c_:>8}: repeater median {cd[c_].median():.3f}, '
              f'null {ct[c_].median():+.3f} +/- {mad:.3f}  ->  {z_:.1f} MAD')

    print('\nVERDICT:')
    g = cd.algn_p.median() - cd.flat_o.median()
    print(f'  best-variant repeater median {cd.algn_p.median():.3f} vs '
          f'published flat {cd.flat_o.median():.3f}  (change {g:+.3f})')
    if g > 0.10:
        print('  -> Moveout correction recovers the deficit. The flat stack was '
              'smearing\n     the arrival; the deficit was processing, not physics.')
    elif g > 0.02:
        print('  -> Partial recovery. Alignment helps but does not close the gap; '
              'a\n     residual coherence ceiling remains to be explained.')
    else:
        print('  -> NO recovery. The sinc/moveout diagnosis is WRONG. The DAS '
              'deficit is\n     not un-corrected moveout and must be explained '
              'some other way.')

    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.26)

    if diag:
        fs = diag['fs']
        t = np.arange(diag['A'].shape[1]) / fs - PRE_S
        sl = slice(int((diag['tp'] - 0.6) * fs), int((diag['tp'] + 1.4) * fs))
        ax = fig.add_subplot(gs[0, 0])
        v = np.percentile(np.abs(diag['A'][:, sl]), 99)
        ax.imshow(diag['A'][:, sl], aspect='auto', cmap='seismic',
                  vmin=-v, vmax=v, extent=[t[sl][0], t[sl][-1],
                                           diag['z'][-1], diag['z'][0]])
        ax.plot(diag['tp'] - PRE_S + diag['p'] * diag['z'], diag['z'], 'k--', lw=1.4)
        ax.set(xlabel='time from origin (s)', ylabel='depth along fibre (m)',
               title=f'A  Raw gather, fitted moveout\n{diag["tag"]}')
        ax.grid(False)

        ax = fig.add_subplot(gs[0, 1])
        tt = np.arange(diag['S'].shape[1]) / fs + T_SEARCH[0] - PRE_S
        ax.imshow(diag['S'], aspect='auto', cmap='magma',
                  extent=[tt[0], tt[-1], diag['ps'][-1] * 1e4, diag['ps'][0] * 1e4])
        ax.axhline(0, color='c', lw=1, ls=':')
        ax.plot(diag['tp'] - PRE_S, diag['p'] * 1e4, 'co', ms=7, mfc='none', mew=2)
        ax.set(xlabel='time from origin (s)',
               ylabel=r'slowness ($10^{-4}$ s/m)',
               title='B  Semblance scan\n(dotted line = flat, p=0)')
        ax.grid(False)

        ax = fig.add_subplot(gs[0, 2])
        q = slice(int((diag['tp'] - 0.2) * fs), int((diag['tp'] + 1.8) * fs))
        tq = np.arange(q.stop - q.start) / fs
        f_ = diag['A'][:, q].mean(axis=0)
        b_ = diag['B'][:, q].mean(axis=0)
        ax.plot(tq, f_ / np.abs(b_).max(), lw=1, c='0.55', label='flat stack')
        ax.plot(tq, b_ / np.abs(b_).max(), lw=1.2, c='C0', label='aligned stack')
        ax.set(xlabel='time from P (s)', ylabel='amplitude (common scale)',
               title=f'C  Same channels, two stacks\n'
                     f'aligned/flat amplitude = '
                     f'{np.abs(b_).max()/max(np.abs(f_).max(),1e-30):.1f}x')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 0])
    xs = np.arange(len(cd))
    for c_, lb, col in [('flat_o', 'flat / origin win', '0.6'),
                        ('algn_o', 'aligned / origin win', 'C0'),
                        ('algn_p', 'aligned / P win', 'C2')]:
        ax.plot(xs, cd[c_].values, 'o-', ms=5, c=col, label=lb, lw=1.2)
    ax.plot(xs, cd.hrsn.values, 'k*--', ms=11, label='HRSN (ground truth)')
    ax.set(xlabel='confirmed repeater pair', ylabel='CC', ylim=(0, 1.05),
           title='D  Does alignment close the gap to HRSN?')
    ax.legend(fontsize=8, loc='lower left')
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 1])
    for c_, lb, col in [('flat_o', 'flat / origin', '0.6'),
                        ('algn_p', 'aligned / P', 'C2')]:
        ax.hist(ct[c_], bins=18, alpha=0.55, color=col, label=f'null, {lb}')
    for c_, col in [('flat_o', '0.3'), ('algn_p', 'C2')]:
        for v_ in cd[c_]:
            ax.axvline(v_, color=col, lw=1.1, alpha=0.85)
    ax.set(xlabel='CC', ylabel='pairs',
           title='E  Null (bars) vs repeaters (lines)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 2])
    vv = np.array([e['V'] for e in E.values()])
    vv = vv[np.isfinite(vv)]
    ax.hist(np.abs(vv), bins=np.logspace(3, 4.3, 26), color='C0', alpha=0.75)
    ax.axvline(2975, color='C3', lw=2,
               label='2975 m/s (my slant-stack,\nNano fibre)')
    ax.axvline(V_BODY, color='k', ls=':', lw=1.2, label='tube/body cut')
    ax.set(xscale='log', xlabel=r'fitted $|V_{app}|$ (m/s)', ylabel='events',
           title='F  Fitted apparent velocity\nper event')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    fig.suptitle('Moveout correction before channel-stacking: does it recover the '
                 'DAS/HRSN correlation deficit?', fontsize=13)
    p_ = os.path.join(HERE, 'moveout_test.png')
    fig.savefig(p_, dpi=140, bbox_inches='tight')
    print(f'\nwrote {p_}')


if __name__ == '__main__':
    main()
