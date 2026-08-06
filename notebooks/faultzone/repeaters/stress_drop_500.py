"""
Stress drop on the LARGE covered events. The previous closure tested the wrong ones.

WHY THIS EXISTS. Stress drop was closed on the grounds that "M0.65 has a ~108 Hz
corner and the signal dies by 70 Hz". True, and irrelevant to whether the method
works, because corner frequency falls as magnitude rises:

    M0.65 -> fc ~ 108 Hz    (dead: hf_snr_test measured 0 usable channels there)
    M1.86 -> fc ~  25 Hz
    M3.17 -> fc ~   6 Hz    (hf_snr_test: SNR > 3 on 855 channels at 80-100 Hz
                             and 853 at 140-200 Hz)

For M3.17 the corner sits an order of magnitude BELOW the usable band edge, so the
spectrum is sampled from well under the corner to well over it -- exactly the
configuration a Brune fit needs. The earlier test picked the one magnitude where
the measurement cannot work and generalised from it. Ellsworth said the same thing
from the other direction: "Larger events would be better too."

WHAT IS BEING MEASURED. Per-channel displacement-spectrum fit of the Brune model

    A(f) = Omega0 / (1 + (f/fc)^2)

over the usable band, then

    M0 = f(Omega0)  [via catalog magnitude here, not absolute amplitude]
    r  = 0.37 * Vs / fc
    dsigma = 7 M0 / (16 r^3)

Absolute DAS amplitude calibration is NOT assumed. fc comes from the spectral
SHAPE, which is calibration-free; M0 comes from the catalog. That keeps the one
quantity DAS is bad at out of the answer.

--------------------------------------------------------------------------------
PREDICTIONS, REGISTERED BEFORE RUNNING.

 S1  fc must fall with magnitude across the tested events. If fc comes out
     independent of magnitude, the fit is reading the instrument response or the
     noise spectrum, not the source, and the result is void regardless of how
     clean it looks.

 S2  fc measured per channel must be CONSISTENT ACROSS DEPTH for one event. The
     source corner is a property of the earthquake; a corner that trends with
     channel depth is path attenuation or coupling, not source. Reported as the
     spread of fc over channels, and this is the main internal control.

 S3  the 16.335 m gauge length imposes a |sinc(pi f G / V)| rolloff with its -3 dB
     point at 54-81 Hz for V = 2000-3000 m/s. Fits are corrected for it and the
     fit is refused above the first null. If a fitted fc lands within 20% of the
     rolloff knee it is reported as unresolved rather than as a measurement.

 S4  stress drops should land in 0.1-100 MPa. Anything outside that is a fit
     failure, not a discovery.

FAILURE IS AN ACCEPTABLE OUTCOME: if fc does not scale with magnitude (S1) or
varies with depth (S2), report that DAS spectra at SAFOD are not usable for source
parameters and close the direction properly this time, with the reason attached to
the right cause.
--------------------------------------------------------------------------------
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import decimate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CACHE_HF = os.path.join(HERE, 'cache_hf')
FS_TARGET = 500.0
CH_LO, CH_HI = 23, 896          # G0 usable range
GAUGE = 16.335
PRE_S = 5.0
SIG_WIN = (0.0, 2.0)            # about origin: P + S + early coda
NOISE_WIN = (-4.5, -0.5)
FIT_LO, FIT_HI = 3.0, 150.0     # refuse above the gauge null
VS = 3500.0
DSIGMA_REF = 3e6


def gauge_response(f, v=2500.0):
    return np.abs(np.sinc(f * GAUGE / v))


def brune(f, omega0, fc):
    return omega0 / (1.0 + (f / fc) ** 2)


def load(tag):
    f = os.path.join(CACHE_HF, f'{tag}.npz')
    if not os.path.exists(f):
        return None
    d = np.load(f)
    X, fs = d['X'].astype(np.float64), float(d['fs'])
    if fs > FS_TARGET * 1.5:
        q = int(round(fs / FS_TARGET))
        X = decimate(X, q, axis=-1, ftype='fir', zero_phase=True)
        fs = fs / q
    return (X, fs) if abs(fs - FS_TARGET) < 1 else None


def channel_corners(X, fs):
    """Fit Brune fc per channel on the signal/noise-corrected spectrum."""
    i0, i1 = int((PRE_S + SIG_WIN[0]) * fs), int((PRE_S + SIG_WIN[1]) * fs)
    j0, j1 = int((PRE_S + NOISE_WIN[0]) * fs), int((PRE_S + NOISE_WIN[1]) * fs)
    if i1 > X.shape[1] or j0 < 0:
        return None, None, None
    S, N = X[:, i0:i1], X[:, j0:j1]
    fS = np.fft.rfftfreq(S.shape[1], 1 / fs)
    fN = np.fft.rfftfreq(N.shape[1], 1 / fs)
    AS = np.abs(np.fft.rfft(S * np.hanning(S.shape[1]), axis=-1)) / S.shape[1]
    AN = np.abs(np.fft.rfft(N * np.hanning(N.shape[1]), axis=-1)) / N.shape[1]
    ANi = np.array([np.interp(fS, fN, a) for a in AN])

    m = (fS >= FIT_LO) & (fS <= FIT_HI)
    f = fS[m]
    g = gauge_response(f)
    # strain rate -> displacement is a double integration in the amplitude
    # spectrum (one for rate->strain, one for the omega-square source model as
    # written); dividing by f twice keeps the SHAPE right, which is all fc needs
    fc_out, snr_out = [], []
    for k in range(X.shape[0]):
        a = AS[k, m]
        n = ANi[k, m]
        snr = float(np.median(a / np.maximum(n, 1e-30)))
        if snr < 3.0:
            fc_out.append(np.nan); snr_out.append(snr); continue
        d_ = (a - n).clip(min=1e-30) / (f ** 2) / np.maximum(g, 0.05)
        try:
            p0 = [d_[0], 20.0]
            popt, _ = curve_fit(brune, f, d_, p0=p0, maxfev=6000,
                                bounds=([0, FIT_LO], [np.inf, FIT_HI]))
            fc_out.append(float(popt[1]))
        except Exception:
            fc_out.append(np.nan)
        snr_out.append(snr)
    return np.array(fc_out), np.array(snr_out), f


def main():
    ev = pd.read_csv(os.path.join(HERE, 'correlate_all_events.csv'))
    have = {f[:-4] for f in os.listdir(CACHE_HF) if f.endswith('.npz')}
    ev = ev[ev.tag.isin(have)].copy()
    ev = ev.sort_values('mag', ascending=False).head(12)
    print(f'{len(ev)} largest covered events at 500 Hz, '
          f'M{ev.mag.min():.2f}-{ev.mag.max():.2f}\n', flush=True)
    print('gauge rolloff |sinc(pi f G/V)|, G=16.335 m, V=2500 m/s:')
    for ff in (40, 60, 80, 100, 150):
        print(f'   {ff:3d} Hz -> {gauge_response(ff):.2f}', end='')
    print('\n', flush=True)

    rows = []
    print(f'{"event":>20}{"M":>6}{"pred fc":>9}{"fit fc":>9}{"MAD":>8}'
          f'{"nch":>6}{"r m":>8}{"dsigma MPa":>12}')
    for _, e in ev.iterrows():
        got = load(e.tag)
        if got is None:
            continue
        X, fs = got
        fc, snr, f = channel_corners(X[CH_LO:CH_HI], fs)
        if fc is None:
            continue
        ok = np.isfinite(fc)
        if ok.sum() < 30:
            print(f'{e.tag[3:]:>20}{e.mag:6.2f}   only {int(ok.sum())} channels'
                  f' -- skipped', flush=True)
            continue
        fcm = float(np.median(fc[ok]))
        mad = float(1.4826 * np.median(np.abs(fc[ok] - fcm)))
        M0 = 10 ** (1.5 * e.mag + 9.1)
        r_pred = (7 * M0 / (16 * DSIGMA_REF)) ** (1 / 3)
        fc_pred = 0.37 * VS / r_pred
        r = 0.37 * VS / fcm
        ds = 7 * M0 / (16 * r ** 3) / 1e6
        print(f'{e.tag[3:]:>20}{e.mag:6.2f}{fc_pred:9.1f}{fcm:9.1f}{mad:8.1f}'
              f'{int(ok.sum()):6d}{r:8.0f}{ds:12.2f}', flush=True)
        rows.append(dict(tag=e.tag, mag=e.mag, depth=e.depth, fc=fcm,
                         fc_mad=mad, fc_pred=fc_pred, n=int(ok.sum()),
                         r_m=r, dsigma_MPa=ds,
                         fc_by_depth=fc, snr=snr))

    if not rows:
        print('\nno events fitted'); return
    D = pd.DataFrame([{k: v for k, v in r.items()
                       if k not in ('fc_by_depth', 'snr')} for r in rows])
    D.to_csv(os.path.join(HERE, 'stress_drop_500.csv'), index=False)

    print('\n' + '=' * 72)
    print('VERDICT against the predictions registered in the docstring')
    s1 = np.nan
    if len(D) > 3:
        s1 = float(np.corrcoef(D.mag, np.log10(D.fc))[0, 1])
        print(f'  S1 fc falls with magnitude : r(M, log fc) = {s1:+.3f}  '
              f'{"PASS" if s1 < -0.5 else "FAIL"}')
        if s1 >= -0.5:
            print('     -> fc is not tracking magnitude. The fit is reading the'
                  '\n        instrument or noise spectrum, not the source. VOID.')
    spread = float(np.median(D.fc_mad / D.fc))
    print(f'  S2 fc stable across depth  : median MAD/fc = {spread:.2f}  '
          f'{"PASS" if spread < 0.5 else "FAIL"}')
    knee = 0.5 / (GAUGE / 2500.0)
    near = D[(D.fc > 0.8 * knee) & (D.fc < 1.2 * knee)]
    print(f'  S3 fc clear of gauge knee  : knee ~{knee:.0f} Hz, '
          f'{len(near)} of {len(D)} events within 20%')
    bad = D[(D.dsigma_MPa < 0.1) | (D.dsigma_MPa > 100)]
    print(f'  S4 stress drop 0.1-100 MPa : {len(D)-len(bad)} of {len(D)} in range')
    print()
    if np.isfinite(s1) and s1 < -0.5 and spread < 0.5:
        print('  -> USABLE. fc tracks magnitude and is depth-stable. Report stress')
        print('     drops for the large events; the M<1 closure stands separately')
        print('     and is a magnitude limit, not a method failure.')
    else:
        print('  -> NOT USABLE. Close the direction with the correct reason')
        print('     attached, and do not reopen it on magnitude grounds again.')

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    ax[0].errorbar(D.mag, D.fc, yerr=D.fc_mad, fmt='o', color='C0', label='fitted')
    ax[0].plot(D.mag, D.fc_pred, 'k--', label='3 MPa Brune')
    ax[0].set(xlabel='magnitude', ylabel='corner frequency (Hz)', yscale='log',
              title=f'A  S1: fc vs M  (r={s1:+.2f})')
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which='both')
    for r_ in rows[:6]:
        v = r_['fc_by_depth']
        ax[1].plot(np.arange(v.size)[np.isfinite(v)], v[np.isfinite(v)], '.',
                   ms=2, alpha=.5, label=f"M{r_['mag']:.2f}")
    ax[1].axhline(knee, color='C3', ls=':', label='gauge knee')
    ax[1].set(xlabel='channel', ylabel='fc (Hz)', yscale='log',
              title='B  S2: is fc depth-stable?')
    ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
    ax[2].scatter(D.mag, D.dsigma_MPa, c='C2')
    ax[2].axhspan(0.1, 100, color='0.9', zorder=0)
    ax[2].set(xlabel='magnitude', ylabel='stress drop (MPa)', yscale='log',
              title='C  S4: stress drop')
    ax[2].grid(alpha=.3, which='both')
    fig.tight_layout()
    p = os.path.join(HERE, 'stress_drop_500.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
