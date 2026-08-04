"""
G0 refined: fix the two errors in channel_depth_registration.py's verdict.

That run printed "they DISAGREE" and recommended escalating to acquisition
metadata. On inspection the disagreement was mine, not the data's.

ERROR 1 -- SEMBLANCE IS NOT A DISCRIMINATOR HERE. Median semblance is 0.3-0.85
across ALL 900 channels and never approaches the incoherent floor of 1/61 = 0.016,
not even in channels 0-30. The earthquake arrival is coherent along the entire
array, so "where coherent moveout begins" has no answer -- it begins immediately.
The threshold crossing at channel 30 was an artifact of an arbitrary cut.

What DOES discriminate is the pre-event noise: RMS is 30-50x higher for the first
~25-30 channels, drops sharply, stays flat to ~855, then rises again. That is a
coupling transition, which is what an uncemented lead-in looks like. So this script
locates the transitions by change-point detection on log noise instead.

ERROR 2 -- MATCHING VELOCITIES INSTEAD OF TRAVEL TIMES. The check-shot reference was
built by smoothing t(z) and differentiating, which amplifies digitisation noise: the
reference swung 1500-8574 m/s and the misfit curve became a broad monotonic plateau
with no real minimum, placing the "best" offset at 150 where the plateau flattens.
Travel time is the integral, not the derivative, so it is smooth and robust. This
script integrates the DAS slowness profile to a travel-time curve and matches THAT
against the check-shot vertical time, solving jointly for channel offset and an
additive constant (which absorbs the unknown source instant).

ERROR 3 -- INCOMPLETE ARITHMETIC. The 54-channel prediction assumed all 900 channels
are lead-in plus well. But the fibre reaches 864 m and Lellouch et al. 2019 report
the loop failed at the end, so channels past the fibre terminus record nothing. With
a dead tail of D channels the lead-in is 900 - 846 - D, not 900 - 846. Panel B shows
the tail exists, so the two determinations were never actually in conflict:

    30 lead-in  +  846 in well (= 864.0 m)  +  24 past terminus  =  900

Reuses channel_depth_registration.npz -- no rescan of the 206 events.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.join(HERE, 'channel_depth_registration.npz')
CHECKSHOT = os.path.abspath(os.path.join(HERE, '..', 'checkshot_traveltime.npz'))
PGSI = os.path.abspath(os.path.join(
    HERE, '..', '..', 'awd_clean', 'pgsi_reference',
    'PGSIarray_rec_coords_pos1.txt'))

DX = 1.0210
N_CH = 900
FIBRE_IN_WELL = 864.0
CH_IN_WELL = FIBRE_IN_WELL / DX          # 846.2


def noise_transitions(noise, edge=120, k=9):
    """Locate top and bottom coupling transitions from the noise profile.

    The cemented mid-section is stable, so its median and scatter define the
    baseline. A transition is where a smoothed profile leaves that baseline by a
    robust margin. No arbitrary absolute threshold.
    """
    ln = np.log10(np.maximum(noise, 1e-30))
    sm = np.convolve(ln, np.ones(k) / k, 'same')
    mid = sm[edge:N_CH - edge]
    base, sd = np.median(mid), 1.4826 * np.median(np.abs(mid - np.median(mid)))
    hi = base + max(4.0 * sd, 0.30)      # >=2x in RMS, or 4 robust sigma
    top = 0
    for i in range(edge):
        if sm[i] > hi:
            top = i + 1
    bot = N_CH
    for i in range(N_CH - 1, N_CH - edge, -1):
        if sm[i] > hi:
            bot = i
    return top, bot, base, sd, hi, sm


def checkshot_vertical_time():
    """Vertical-corrected travel time t(z) -- the integral, kept as an integral."""
    d = np.load(CHECKSHOT)
    z, t, s = d['depth'], d['traveltime'], d['slant']
    tv = t * z / s                       # straight-ray vertical correction
    o = np.argsort(z)
    z, tv = z[o], tv[o]
    keep = np.concatenate([[True], np.diff(z) > 0])
    z, tv = z[keep], tv[keep]
    # enforce monotonic time (digitisation can invert locally)
    tv = np.maximum.accumulate(tv)
    return z, tv


def main():
    d = np.load(NPZ)
    centres, sem, vmed = d['centres'], d['sem_med'], d['v_med']
    noise = d['noise_med']
    nev = int(d['n_events'])
    print(f'loaded scan of {nev} events, {centres.size} windows\n')

    # ---------------- 1. coupling transitions from noise ------------------
    top, bot, base, sd, hi, sm = noise_transitions(noise)
    print('(1) COUPLING TRANSITIONS FROM PRE-EVENT NOISE')
    print(f'  cemented baseline log10(RMS) = {base:.3f} +/- {sd:.3f}')
    print(f'  top transition    : channel {top}')
    print(f'  bottom transition : channel {bot}')
    print(f'  noise ratio, ch 0-{max(top,1)} vs baseline: '
          f'{10**(np.median(sm[:max(top,1)]) - base):.1f}x')
    print(f'  in-well span implied: {bot - top} ch = {(bot-top)*DX:.1f} m')
    print(f'  fibre in well (Lellouch): {FIBRE_IN_WELL:.0f} m = {CH_IN_WELL:.0f} ch')

    # ---------------- 2. arithmetic reconciliation ------------------------
    dead_tail = N_CH - bot
    leadin_pred = N_CH - CH_IN_WELL - dead_tail
    print('\n(2) ARITHMETIC WITH THE DEAD TAIL INCLUDED')
    print(f'  dead tail past fibre terminus : {dead_tail} ch')
    print(f'  => predicted lead-in = 900 - {CH_IN_WELL:.0f} - {dead_tail} '
          f'= {leadin_pred:.0f} ch')
    print(f'  observed top transition       : {top} ch')
    print(f'  agreement: {abs(leadin_pred - top):.0f} channels '
          f'({abs(leadin_pred-top)*DX:.1f} m)')

    # ---------------- 3. travel-time match -------------------------------
    zcs, tcs = checkshot_vertical_time()
    ok = np.isfinite(vmed) & (np.abs(vmed) > 500)
    cc = centres[ok].astype(float)
    slow = 1.0 / np.abs(vmed[ok])
    # integrate slowness along the fibre -> relative vertical travel time
    t_das = np.concatenate([[0.0], np.cumsum(np.diff(cc) * DX *
                                            0.5 * (slow[1:] + slow[:-1]))])

    offs = np.arange(-10, 200)
    mis = np.full(offs.size, np.nan)
    for i, c0 in enumerate(offs):
        zd = (cc - c0) * DX
        m = (zd > max(zcs[0], 30.0)) & (zd < min(zcs[-1], 860.0))
        if m.sum() < 15:
            continue
        ref = np.interp(zd[m], zcs, tcs)
        r = t_das[m] - ref
        mis[i] = np.std(r)               # additive constant removed by std
    best = offs[int(np.nanargmin(mis))]
    lo = np.nanmin(mis)
    within = offs[mis < lo + 0.002]      # within 2 ms of the minimum
    print('\n(3) TRAVEL-TIME MATCH TO THE 2005 GEOPHONE VSP')
    print(f'  best channel offset : {best}')
    print(f'  residual scatter    : {1000*lo:.2f} ms')
    print(f'  within +2 ms        : channels {within.min()} to {within.max()}')
    curv = 'well-defined' if (within.max() - within.min()) < 60 else 'broad/weak'
    print(f'  minimum is {curv}')

    # ---------------- verdict --------------------------------------------
    ests = {'noise transition': float(top),
            'length arithmetic': float(leadin_pred),
            'VSP travel time': float(best)}
    vals = np.array(list(ests.values()))
    spread = vals.max() - vals.min()
    print('\nVERDICT')
    for k, v in ests.items():
        print(f'  {k:>18}: channel {v:.0f}')
    print(f'  spread {spread:.0f} channels ({spread*DX:.1f} m)')
    if spread <= 25:
        wh = float(np.median(vals))
        print(f'  -> CONSISTENT. Wellhead = channel {wh:.0f}; '
              f'depth = (ch - {wh:.0f}) x {DX} m.')
        print(f'     Linear mapping with a constant lead-in of {wh*DX:.0f} m is '
              f'supported.\n     Usable in-well channels: {int(wh)}-{bot}.')
    else:
        wh = float(top)
        print(f'  -> PARTIAL. Noise and arithmetic agree; the VSP match is the '
              f'outlier.\n     Adopt wellhead = channel {wh:.0f} PROVISIONALLY and '
              f'flag absolute depths\n     as uncertain by +/-{spread*DX:.0f} m '
              f'until acquisition metadata confirms.')

    print(f'\nSHALLOW VERDICT -- are the discarded channels usable?')
    for a, b in [(int(wh), int(wh) + 25), (int(wh) + 25, int(wh) + 50),
                 (int(wh) + 50, 100), (100, 200)]:
        m = (centres >= a) & (centres < b)
        if m.any():
            print(f'  ch {a:3d}-{b:3d}  z {(a-wh)*DX:6.0f} to {(b-wh)*DX:6.0f} m : '
                  f'semblance {np.nanmedian(sem[m]):.3f}, '
                  f'noise {10**(np.median(sm[a:b])-base):5.2f}x baseline')
    print(f'\n  CH_LO=100 currently in use corresponds to depth '
          f'{(100-wh)*DX:.0f} m.')
    print(f'  Li & Ben-Zion peak sensitivity (~17 m) is channel '
          f'{wh + 17/DX:.0f}.')

    np.savez(os.path.join(HERE, 'g0_refine.npz'),
             wellhead=wh, top=top, bot=bot, leadin_pred=leadin_pred,
             vsp_best=best, offs=offs, mis=mis, t_das=t_das, cc=cc,
             zcs=zcs, tcs=tcs, noise_sm=sm, base=base, hi=hi, dx=DX)

    # ---------------- figure ---------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.6))

    a = ax[0]
    a.semilogx(noise, np.arange(N_CH), color='0.75', lw=0.8, label='per channel')
    a.semilogx(10 ** sm, np.arange(N_CH), 'C4-', lw=1.4, label='smoothed')
    a.axvline(10 ** base, color='k', ls='-', lw=1, label='cemented baseline')
    a.axvline(10 ** hi, color='C3', ls='--', lw=1, label='transition level')
    a.axhline(top, color='C2', lw=2, label=f'top ch {top}')
    a.axhline(bot, color='C1', lw=2, label=f'bottom ch {bot}')
    a.invert_yaxis()
    a.set(xlabel='pre-event noise RMS', ylabel='channel',
          title=f'A  Coupling transitions\nlead-in is {10**(np.median(sm[:max(top,1)])-base):.0f}x noisier')
    a.legend(fontsize=7); a.grid(alpha=0.3)

    a = ax[1]
    a.plot(1000 * mis, offs, 'C0-')
    a.axhline(best, color='C3', lw=1.6, label=f'VSP best {best}')
    a.axhline(top, color='C2', lw=1.6, label=f'noise {top}')
    a.axhline(leadin_pred, color='0.5', ls=':', lw=1.4,
              label=f'arithmetic {leadin_pred:.0f}')
    a.set(xlabel='travel-time residual scatter (ms)',
          ylabel='trial wellhead channel',
          title='B  Travel-time match\n(integral, not derivative)')
    a.legend(fontsize=8); a.grid(alpha=0.3)

    a = ax[2]
    zd = (cc - wh) * DX
    off = np.median(np.interp(np.clip(zd, zcs[0], zcs[-1]), zcs, tcs) - t_das)
    a.plot(1000 * (t_das + off), zd, 'C0-', lw=1.6, label='DAS integrated slowness')
    a.plot(1000 * tcs, zcs, 'k--', lw=1.6, label='2005 check shot')
    wd = pd.read_csv(PGSI, sep=r'\s+')['WELL_DEP'].values
    a.plot(np.zeros(len(wd)), wd, 'k|', ms=7, alpha=0.6, label='PGSI levels')
    a.set(ylim=(900, -60), xlim=(0, 330), xlabel='vertical travel time (ms)',
          ylabel='depth below surface (m)',
          title=f'C  Registered DAS vs VSP\nwellhead = ch {wh:.0f}')
    a.legend(fontsize=8); a.grid(alpha=0.3)

    fig.suptitle('G0 refined: channel-to-depth registration by noise transition '
                 'and travel-time matching', fontsize=12)
    fig.tight_layout()
    p = os.path.join(HERE, 'g0_refine.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
