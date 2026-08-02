"""
Locate the physical end / break of the Nano (cemented) fiber from ambient noise.

Why not RMS: past a fiber break the interrogator still reports a noise floor --
that noise is optical and electronic, not strain -- so a flat RMS profile cannot
tell live fiber from no fiber. nano_far_end_check.py hit exactly that wall.

The discriminator is spatial coherence. Live fiber sensing ambient ground motion
sees wavelengths of tens to hundreds of metres, so neighbouring channels are
correlated. Past a break there is no backscatter, the phase is random, and
channels decorrelate.

One catch: the gauge length is 16.4588 m = 13 channels, so channels closer than
13 apart share raw samples and correlate even on pure instrument noise. The test
therefore has to be read at separations *beyond* the gauge length -- k >= 20 is
where a real answer lives.

Lellouch et al. 2019 (10.1029/2019JB017533) put this fiber's end at 864 m with a
failed end loop capping usable data at 800 m. At 1.26606 m/channel that predicts
a transition at ch 632 (800 m) or ch 682 (864 m).

Runs on ambient files from before the AWD survey window -- no weight drops, so
what is left is ground noise and instrument noise, which is the contrast we want.
"""
import sys
from pathlib import Path
SEARCH_DIRS = [
    '/home/groups/edunham/nberrios/safod_das/DAS-utilities/python',
    '/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python',
]
for _p in SEARCH_DIRS:
    if Path(_p).exists() and _p not in sys.path:
        sys.path.insert(0, _p)

import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from DASutils import readFile_protobuf

NANO_DIR = '/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/Nano/'
FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'

DX = 1.26606202          # channel_spacing from the .pb acquisition_stats
GL_M = 16.4588           # gauge_length, same source
GL_CH = int(round(GL_M / DX))

# Band matters more than anything else here. Sintela delivers strain *rate*, so
# the spectrum is the ground-motion spectrum multiplied by omega^2 -- a first
# pass at 2-90 Hz came back with 40-90 Hz power ~100x the 2-10 Hz power and
# coherence flat at ~0.1 everywhere, because that band is photon noise, which is
# spatially incoherent on live and dead fiber alike. Coherent ground motion lives
# at low frequency, where wavelengths (V/f ~ 100-1000 m) exceed the separations
# being tested.
BANDS = {'low 1-8 Hz': (1.0, 8.0), 'high 30-90 Hz': (30.0, 90.0)}
PROBE_BAND = 'low 1-8 Hz'

# Keeping the common mode (median=False) is what makes the low band coherent,
# but common-mode *interrogator* noise -- laser phase drift, shared electronics --
# is also coherent across every channel, live or dead, and would mimic a live
# fiber past a break. Re-running with the across-channel median removed
# separates them: genuine local ground motion survives, common mode does not.
MEDIAN_MODES = [False, True]
LAGS = [1, 5, 13, 20, 40, 80]      # channel separations to correlate across
DUR_S = 240.0                      # ambient window per file

# published landmarks, in channels
LANDMARKS = {'Lellouch 800 m cap': 800.0 / DX, 'fiber end 864 m': 864.0 / DX}


def ambient_files(n=3):
    """Files from before the survey window, so no weight drops are present."""
    allf = sorted(glob.glob(NANO_DIR + '*.pb'))
    pre = [f for f in allf if '2026-06-16_0' in os.path.basename(f)]
    return pre[20:20 + n] if len(pre) >= 20 + n else allf[:n]


def lagged_corr(x, lag):
    """Zero-time-lag Pearson correlation between channel i and channel i+lag."""
    a, b = x[:-lag], x[lag:]
    num = np.sum(a * b, axis=1)
    den = np.sqrt(np.sum(a * a, axis=1) * np.sum(b * b, axis=1))
    out = np.full(x.shape[0], np.nan)
    out[:-lag] = num / np.where(den > 0, den, np.nan)
    return out


def main():
    files = ambient_files()
    print(f'gauge length {GL_M} m = {GL_CH} channels -> read lags k >= {GL_CH+1}')
    print('ambient files:')
    for f in files:
        print('   ', os.path.basename(f))

    # median=False matters: the default subtracts the across-channel median at
    # every sample, which strips the common-mode ground motion this test needs.
    results, rms = {}, {}
    for bname, (fmin, fmax) in BANDS.items():
        for med in MEDIAN_MODES:
            if bname != PROBE_BAND and med:
                continue                      # only need both modes on the probe band
            key = bname + (' (median removed)' if med else '')
            corr_stack, rms_stack = [], []
            for f in files:
                DAS, info = readFile_protobuf([f], fmin=fmin, fmax=fmax,
                                              desampling=False, median=med)
                fs = info['fs']
                seg = DAS[:, :int(DUR_S * fs)]
                seg = seg - seg.mean(axis=1, keepdims=True)
                corr_stack.append(np.array([lagged_corr(seg, k) for k in LAGS]))
                rms_stack.append(np.sqrt(np.mean(seg ** 2, axis=1)))
                del DAS, seg
            results[key] = np.nanmedian(np.array(corr_stack), axis=0)
            rms[key] = np.nanmedian(np.array(rms_stack), axis=0)
            print(f'{key}: median k=20 coherence '
                  f'{np.nanmedian(results[key][LAGS.index(20)]):.3f}')

    # The median-removed profile is the one that can actually be trusted to
    # report live fiber, so probe that.
    corr = results[PROBE_BAND + ' (median removed)']       # (n_lags, n_ch)
    band = rms[PROBE_BAND] / rms['high 30-90 Hz']         # red-vs-white ratio
    common = results[PROBE_BAND][LAGS.index(20)]          # with common mode kept
    n_ch = corr.shape[1]
    dist = np.arange(n_ch) * DX

    # Transition: use the first lag beyond the gauge length, smooth it, and find
    # where it drops below halfway between its live-zone and tail levels.
    probe = corr[LAGS.index(20)]
    sm = np.convolve(np.nan_to_num(probe), np.ones(15) / 15, mode='same')
    live = np.nanmedian(sm[50:250])
    tail = np.nanmedian(sm[-60:])

    print(f'\nk=20 coherence, common mode removed: '
          f'shallow (60-320 m) {live:.3f}, tail {tail:.3f}')

    # Verdict before any transition-picking. This test only means something if
    # local (non-common-mode) coherence is actually present somewhere to lose.
    # It is not: at 1-8 Hz the array is short against the wavelength so the
    # wavefield is nearly all common mode, and removing it leaves incoherent
    # local noise everywhere past ~300 m. At 30-90 Hz there is no coherent
    # ground motion to begin with. So there is no band in which live and dead
    # fiber look different, and this method cannot locate the fiber end here.
    inconclusive = live < 0.3
    if inconclusive:
        print('\nVERDICT: INCONCLUSIVE. Local coherence never rises far enough '
              'above the tail\n  for its collapse to mark anything. Do not read '
              'a fiber end off this.\n  The published 800 m failure depth can be '
              'neither confirmed nor refuted\n  from ambient noise -- an OTDR '
              'trace is required.')

    half = 0.5 * (live + tail)
    below = np.where(sm < half)[0]
    below = below[below > 250]
    trans = -1 if inconclusive else (below[0] if below.size else -1)
    if trans > 0:
        print(f'transition at ch {trans} = {trans*DX:.1f} m '
              f'(+/- {GL_CH} ch = {GL_M:.1f} m from gauge length)')
        for name, c in LANDMARKS.items():
            print(f'   vs {name}: ch {c:.0f} ({c*DX:.1f} m), '
                  f'offset {(trans-c)*DX:+.1f} m')
    else:
        print('no sustained transition found -- coherence never collapses')

    print('\n  ch   fiber_m' + ''.join(f'   k={k:<3d}' for k in LAGS) + '   lo/hi')
    for c in range(0, n_ch, 40):
        print(f'{c:5d} {c*DX:8.1f}' +
              ''.join(f' {corr[i, c]:7.3f}' for i in range(len(LAGS))) +
              f' {band[c]:7.2f}')

    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for i, k in enumerate(LAGS):
        style = '--' if k <= GL_CH else '-'
        ax[0].plot(dist, corr[i], style, lw=1.0, label=f'k={k}'
                   + (' (within GL)' if k <= GL_CH else ''))
    ax[0].plot(dist, common, lw=1.6, color='k', alpha=0.5,
               label='k=20, common mode kept')
    ax[0].axhline(0, color='0.6', lw=0.8)
    ax[0].set_ylabel('channel-to-channel correlation')
    ax[0].legend(fontsize=8, ncol=2)
    ax[0].set_title('Nano (cemented) fiber: ambient-noise coherence vs distance\n'
                    'black = common mode kept, colours = common mode removed',
                    fontsize=10)
    ax[1].semilogy(dist, band, lw=0.9, color='C4')
    ax[1].set_ylabel('RMS ratio  1-8 Hz / 30-90 Hz')
    ax[1].set_xlabel('distance along fiber (m)')
    for a in ax:
        a.grid(alpha=0.3)
        for name, c in LANDMARKS.items():
            a.axvline(c * DX, ls=':', c='C3', lw=1.2)
        if trans > 0:
            a.axvline(trans * DX, ls='--', c='C2', lw=1.4)
    ax[0].text(0.99, 0.02, 'red dotted = published landmarks, green dashed = measured',
               transform=ax[0].transAxes, ha='right', fontsize=8, color='0.3')
    fig.tight_layout()
    out = os.path.join(FIG_DIR, 'nano_fiber_end_coherence.png')
    fig.savefig(out, dpi=140)
    print('\nSaved', out)


if __name__ == '__main__':
    main()
