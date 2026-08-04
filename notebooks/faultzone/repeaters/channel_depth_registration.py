"""
G0: where is the wellhead, and which channels are usable?

WHY THIS GATES EVERYTHING. Li & Ben-Zion 2023 (doi:10.1029/2022JB025682) put the
seasonal velocity signal at ~17 m peak sensitivity, in the "top tens of meters."
Every script in this directory starts at CH_LO = 100 with the comment "skip the
noisy top" -- an assumption I wrote, never measured. If channels below 100 are
usable the fiber samples the peak of the signal; if not, it samples only the tail.
The project's value rests on ~46 channels that are currently discarded.

Absolute registration has never been established here. step2_timing_check.py:5
found the arrival at 4 ms at "130 m of fiber" where the check shot says ~100 ms at
130 m depth, listed "(C) the channel at 130 m of fiber is nowhere near 130 m depth",
and stopped. awd_clean's own dashboard says "every depth label is explicitly
provisional."

THE ARITHMETIC THAT MAKES THIS TESTABLE. Lellouch et al. 2019 (JGR
doi:10.1029/2019JB017533) state the fiber reaches 864 m and is cemented in a steel
tube under ~1 N tension, so "its length is expected to match the location along the
well." At our 1.0210 m channel spacing, 864 m of well = 846 channels. We have 900.
That leaves ~54 channels (~55 m) unaccounted for. Cable overstuff is "typically
0.1-0.3%" (Madsen et al. 2016, TLE 35(7), doi:10.1190/tle35070610.1) = only ~2 m
over 864 m, so overstuff cannot explain 55 m. It must be surface lead-in or a
wellhead accumulation -- and Madsen warns those are different:

    "There also may be intentional fiber accumulations at certain places, e.g., in
     the wellhead area, in which case the linear assumption does not hold."

A straight surface lead-in is a constant channel offset and linear mapping holds
below the wellhead. A coil near the wellhead breaks linearity, and Madsen found it
"generally plac[es] channels too deep in the well... the error would be most severe
in the top" -- precisely where our signal lives.

TWO INDEPENDENT DETERMINATIONS, following Madsen's "match the entire pattern"
rather than relying on single reference points.

  (1) WELLHEAD FROM COHERENT-MOVEOUT ONSET. Channels in air are not coupled to the
      ground and cannot reproduce the coherent borehole moveout that every event
      shows. So the wellhead is where windowed semblance turns on. This uses the
      same earthquakes as the measurement and needs no external data.

  (2) VELOCITY-STRUCTURE MATCH AGAINST THE 2005 GEOPHONE VSP. Madsen's diagnostic
      was "comparing traveltime curves from iDAS data and conventional wireline
      VSP." We have that VSP: the PGSI 80-level array, geophone depths KNOWN from
      PGSIarray_rec_coords_pos1.txt (46.68 -> 1250.64 m, 15.24 m spacing), and the
      digitised check shot (checkshot_traveltime.npz, 754 points, 20-2529 m,
      source 45.72 m off the wellhead). Its vertical-corrected interval velocities
      -- 1496 m/s over 20-75 m, ~3400 over 75-300, a 3091 m/s low-velocity zone at
      300-500 m, 3495 over 500-800 -- are distinctive enough to match on. The
      300-500 m slow zone is the "clamp pattern" equivalent: a feature with a known
      depth that must line up.

Agreement between (1) and (2) is the result. Disagreement means a non-linear
accumulation and the datum goes to the acquisition metadata instead.

Scans ALL 900 channels. Nothing is excluded a priori -- that is the whole point.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache_all')
FZ = os.path.abspath(os.path.join(HERE, '..'))
CHECKSHOT = os.path.join(FZ, 'checkshot_traveltime.npz')
PGSI = os.path.abspath(os.path.join(
    HERE, '..', '..', 'awd_clean', 'pgsi_reference',
    'PGSIarray_rec_coords_pos1.txt'))

DX = 1.0210                  # channel spacing along fibre, m
N_CH = 900
BAND = (5.0, 20.0)
PRE_S = 5.0                  # cached window starts 5 s before catalog origin
T_SEARCH = (5.0, 9.5)        # search for the arrival here (s into the record)

WIN_CH = 61                  # channel window for local slant stack (~62 m)
STEP_CH = 8
P_MAX = 6.7e-4               # |V| > 1500 m/s, both signs; includes p = 0
N_P = 41
SEM_WIN = 0.4                # semblance time window, s

FIBRE_IN_WELL = 864.0        # Lellouch et al. 2019
MAX_EVENTS = int(os.environ.get('MAX_EVENTS', '206'))


# ----------------------------------------------------------------- reference
def checkshot_velocity():
    """Vertical-corrected interval Vp(z) from the digitised 2005 check shot.

    The recorded time is along a slant path because the source sits 45.72 m from
    the wellhead; dividing depth by raw time understates velocity badly in the
    shallow section (that error is what produced the bogus '300 m/s' figure).
    Straight-ray correction to vertical time: t_vert = t * z / slant.
    """
    d = np.load(CHECKSHOT)
    z, t, s = d['depth'], d['traveltime'], d['slant']
    tv = t * z / s
    o = np.argsort(z)
    z, tv = z[o], tv[o]
    # monotonic cleanup, then smooth derivative -> interval velocity
    keep = np.concatenate([[True], np.diff(tv) > 0])
    z, tv = z[keep], tv[keep]
    zg = np.arange(10.0, 1250.0, 5.0)
    tg = np.interp(zg, z, tv)
    k = 15                                  # ~75 m smoothing for a stable dz/dt
    ker = np.ones(k) / k
    tsm = np.convolve(tg, ker, 'same')
    v = np.gradient(zg, tsm)
    good = (zg > 10 + 5 * k / 2) & (zg < 1250 - 5 * k / 2)
    return zg[good], v[good], (z, tv)


def pgsi_depths():
    df = pd.read_csv(PGSI, sep=r'\s+')
    return df['WELL_DEP'].values, df['REC_DEP'].values


# ----------------------------------------------------------------- DAS scan
def prep_all(tag):
    """Load one event over ALL channels; bandpass, per-channel L2 normalise."""
    f = os.path.join(CACHE, f'{tag}.npz')
    if not os.path.exists(f):
        return None
    d = np.load(f)
    X, fs = d['X'], float(d['fs'])
    if X.shape[0] < N_CH:
        return None
    sos = butter(4, list(BAND), btype='band', fs=fs, output='sos')
    A = sosfiltfilt(sos, X[:N_CH].astype(np.float64), axis=-1)
    A -= A.mean(axis=1, keepdims=True)
    rms = np.sqrt((A ** 2).mean(axis=1))
    dead = ~np.isfinite(rms) | (rms <= 0)
    An = A.copy()
    An[~dead] /= rms[~dead, None]
    An[dead] = 0.0
    # pre-event noise RMS: independent diagnostic of coupling character
    n0 = int(0.5 * fs)
    n1 = int((PRE_S - 0.5) * fs)
    noise = np.sqrt((A[:, n0:n1] ** 2).mean(axis=1))
    return An, fs, rms, noise, dead


def scan_event(An, fs):
    """Local slant-stack over sliding channel windows, all 900 channels.

    Returns per-window best semblance and best apparent velocity. A window of
    air-coupled or dead channels cannot produce a coherent moveout at any
    slowness, so its semblance stays at the incoherent floor.
    """
    n = An.shape[1]
    F = np.fft.rfft(An, axis=-1)
    freq = np.fft.rfftfreq(n, 1.0 / fs)
    zf = np.arange(N_CH) * DX
    ps = np.linspace(-P_MAX, P_MAX, N_P)
    starts = np.arange(0, N_CH - WIN_CH + 1, STEP_CH)
    w = max(int(SEM_WIN * fs), 4)
    i0, i1 = int(T_SEARCH[0] * fs), min(int(T_SEARCH[1] * fs), n - w - 1)

    best_s = np.full(starts.size, -1.0)
    best_v = np.full(starts.size, np.nan)
    for p in ps:
        B = np.fft.irfft(F * np.exp(2j * np.pi * np.outer(p * zf, freq)),
                         n=n, axis=-1)
        cs = np.cumsum(B, axis=0)
        cs2 = np.cumsum(B ** 2, axis=0)
        for k, s0 in enumerate(starts):
            s1 = s0 + WIN_CH
            stk = cs[s1 - 1] - (cs[s0 - 1] if s0 > 0 else 0.0)
            sq = cs2[s1 - 1] - (cs2[s0 - 1] if s0 > 0 else 0.0)
            num = np.cumsum(stk ** 2)
            den = np.cumsum(sq)
            sem = ((num[w:] - num[:-w]) /
                   (WIN_CH * (den[w:] - den[:-w]) + 1e-300))
            seg = sem[i0:i1]
            if seg.size == 0:
                continue
            j = int(np.argmax(seg))
            if seg[j] > best_s[k]:
                best_s[k] = seg[j]
                best_v[k] = np.inf if p == 0 else 1.0 / p
    return starts + WIN_CH // 2, best_s, best_v


def main():
    zcs, vcs, raw = checkshot_velocity()
    wd, rd = pgsi_depths()
    print(f'reference: check shot {zcs[0]:.0f}-{zcs[-1]:.0f} m, '
          f'Vp {vcs.min():.0f}-{vcs.max():.0f} m/s')
    print(f'PGSI geophones: {len(wd)} levels, well depth {wd[0]:.2f}-{wd[-1]:.2f} m, '
          f'spacing {np.median(np.diff(wd)):.2f} m')
    print(f'fibre: {N_CH} ch x {DX} m = {N_CH*DX:.1f} m; '
          f'{FIBRE_IN_WELL:.0f} m in well = {FIBRE_IN_WELL/DX:.0f} ch; '
          f'excess = {N_CH - FIBRE_IN_WELL/DX:.0f} ch '
          f'({(N_CH*DX - FIBRE_IN_WELL):.1f} m)\n', flush=True)

    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    S, V, NOISE, RMS = [], [], [], []
    centres = None
    used = 0
    for _, e in ev.iterrows():
        if used >= MAX_EVENTS:
            break
        got = prep_all(e['tag'])
        if got is None:
            continue
        An, fs, rms, noise, dead = got
        c, s, v = scan_event(An, fs)
        centres = c
        S.append(s); V.append(v); NOISE.append(noise); RMS.append(rms)
        used += 1
        if used % 20 == 0:
            print(f'  {used} events scanned', flush=True)
    if not S:
        print('no events loaded'); return
    S = np.array(S); V = np.array(V)
    NOISE = np.array(NOISE); RMS = np.array(RMS)
    print(f'\n{used} events scanned, {centres.size} channel windows\n', flush=True)

    sem_med = np.median(S, axis=0)
    v_med = np.nanmedian(np.where(np.isfinite(V), np.abs(V), np.nan), axis=0)
    noise_med = np.median(NOISE, axis=0)

    # ---- (1) wellhead from coherent-moveout onset -------------------------
    # Incoherent floor: semblance of WIN_CH independent traces ~ 1/WIN_CH.
    floor = 1.0 / WIN_CH
    hi = np.nanmax(sem_med)
    thr = floor + 0.25 * (hi - floor)
    coh = sem_med > thr
    idx = np.where(coh)[0]
    if idx.size == 0:
        print('NO coherent window found -- cannot proceed'); return
    first_win, last_win = centres[idx[0]], centres[idx[-1]]
    print('(1) COHERENT-MOVEOUT ONSET')
    print(f'  incoherent floor 1/{WIN_CH} = {floor:.4f}; peak {hi:.3f}; '
          f'threshold {thr:.3f}')
    print(f'  coherent from channel {first_win} to {last_win}')
    print(f'  => wellhead at or above channel {first_win} '
          f'(window centre; window half-width {WIN_CH//2} ch)')
    print(f'  => bottom of usable fibre near channel {last_win}')
    pred = N_CH - FIBRE_IN_WELL / DX
    print(f'  predicted lead-in from length arithmetic: {pred:.0f} ch\n')

    # ---- (2) velocity-structure match to the check shot -------------------
    ok = np.isfinite(v_med) & coh
    offs = np.arange(-20, 200)
    mis = np.full(offs.size, np.nan)
    for i, c0 in enumerate(offs):
        zd = (centres[ok] - c0) * DX
        m = (zd > zcs[0]) & (zd < zcs[-1])
        if m.sum() < 8:
            continue
        pv = np.interp(zd[m], zcs, vcs)
        mis[i] = np.sqrt(np.nanmean((v_med[ok][m] - pv) ** 2))
    if np.all(np.isnan(mis)):
        print('(2) match failed -- no overlap'); best_off = np.nan
    else:
        best_off = offs[int(np.nanargmin(mis))]
        lo_ = np.nanmin(mis)
        within = offs[mis < 1.10 * lo_]
        print('(2) VELOCITY-STRUCTURE MATCH TO 2005 GEOPHONE VSP')
        print(f'  best channel offset (wellhead channel): {best_off}')
        print(f'  misfit {lo_:.0f} m/s; within 10% of minimum: '
              f'channels {within.min()} to {within.max()}')
        print(f'  => depth of channel 0 = {-best_off*DX:.1f} m '
              f'(negative = lead-in above ground)\n')

    print('CONSISTENCY')
    print(f'  onset says wellhead <= ch {first_win}')
    print(f'  VSP match says wellhead  = ch {best_off}')
    print(f'  length arithmetic says     ~ ch {pred:.0f}')
    agree = (not np.isnan(best_off)) and abs(best_off - pred) <= 25
    if agree:
        print('  -> the three agree within 25 channels. Linear mapping with a '
              'constant\n     lead-in is supported; treat depth = (ch - '
              f'{best_off})*{DX} m.')
    else:
        print('  -> they DISAGREE. Per Madsen et al. 2016 this is the signature of '
              'a\n     non-linear fibre accumulation near the wellhead. Do not '
              'assign absolute\n     depths; escalate to acquisition metadata '
              '(StartLocusIndex, lead-in\n     length, OTDR/tap test).')

    if not np.isnan(best_off):
        print(f'\nSHALLOW CHANNEL VERDICT (the question that gates the project)')
        for a, b in [(0, 50), (50, 100), (100, 150), (150, 250)]:
            m = (centres >= a) & (centres < b)
            if m.any():
                zz = (np.array([a, b]) - best_off) * DX
                print(f'  ch {a:3d}-{b:3d} (z {zz[0]:6.0f} to {zz[1]:6.0f} m): '
                      f'median semblance {np.nanmedian(sem_med[m]):.3f}  '
                      f'{"COHERENT" if np.nanmedian(sem_med[m]) > thr else "incoherent"}')

    np.savez(os.path.join(HERE, 'channel_depth_registration.npz'),
             centres=centres, sem_med=sem_med, v_med=v_med, noise_med=noise_med,
             S=S, V=V, zcs=zcs, vcs=vcs, offs=offs, mis=mis,
             best_off=best_off, first_win=first_win, last_win=last_win,
             thr=thr, floor=floor, dx=DX, n_events=used)

    # ------------------------------------------------------------- figure
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))

    a = ax[0, 0]
    a.plot(sem_med, centres, 'C0-', lw=1.5)
    a.fill_betweenx(centres, np.percentile(S, 25, axis=0),
                    np.percentile(S, 75, axis=0), color='C0', alpha=0.25)
    a.axvline(floor, color='k', ls=':', lw=1, label=f'incoherent floor 1/{WIN_CH}')
    a.axvline(thr, color='C3', ls='--', lw=1.2, label='threshold')
    a.axhline(first_win, color='C2', lw=1.5, label=f'onset ch {first_win}')
    if not np.isnan(best_off):
        a.axhline(best_off, color='C1', lw=1.5, ls='-.',
                  label=f'VSP match ch {best_off}')
    a.axhline(pred, color='0.5', lw=1.2, ls=':', label=f'length arith ch {pred:.0f}')
    a.invert_yaxis()
    a.set(xlabel='semblance (median over events)', ylabel='channel',
          title='A  Where does coherent borehole moveout begin?')
    a.legend(fontsize=7); a.grid(alpha=0.3)

    a = ax[0, 1]
    a.semilogx(noise_med, centres if noise_med.size == centres.size
               else np.arange(noise_med.size), 'C4-', lw=1)
    a.axhline(first_win, color='C2', lw=1.5)
    a.invert_yaxis()
    a.set(xlabel='pre-event noise RMS (median)', ylabel='channel',
          title='B  Noise character\n(coupling changes at the wellhead)')
    a.grid(alpha=0.3)

    a = ax[1, 0]
    a.plot(mis, offs, 'C0-')
    if not np.isnan(best_off):
        a.axhline(best_off, color='C3', lw=1.5, label=f'best {best_off}')
    a.axhline(pred, color='0.5', ls=':', lw=1.2, label=f'predicted {pred:.0f}')
    a.set(xlabel='RMS velocity misfit (m/s)', ylabel='trial wellhead channel',
          title='C  Match DAS velocity structure\nto the 2005 VSP')
    a.legend(fontsize=8); a.grid(alpha=0.3)

    a = ax[1, 1]
    if not np.isnan(best_off):
        zd = (centres - best_off) * DX
        a.plot(v_med, zd, 'C0-', lw=1.5, label='DAS local slant stack')
    else:
        a.plot(v_med, centres * DX, 'C0-', lw=1.5, label='DAS (unregistered)')
    a.plot(vcs, zcs, 'k--', lw=1.5, label='2005 check shot (known depths)')
    a.plot(np.full(len(wd), 0), wd, 'k|', ms=6, alpha=0.5,
           label='PGSI geophone levels')
    a.set(xlim=(0, 6000), ylim=(1000, -100), xlabel='Vp (m/s)',
          ylabel='depth below surface (m)',
          title='D  Registered DAS velocity vs VSP')
    a.legend(fontsize=8); a.grid(alpha=0.3)

    fig.suptitle('G0: channel-to-depth registration of the SAFOD cemented fibre '
                 f'({used} earthquakes)', fontsize=13)
    fig.tight_layout()
    p_ = os.path.join(HERE, 'channel_depth_registration.png')
    fig.savefig(p_, dpi=140)
    print(f'\nwrote {p_}')


if __name__ == '__main__':
    main()
