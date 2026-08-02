"""
Spectral content of the AWD direct-P arrival on both SAFOD fibers.

Answers three things the depth/velocity work left open:
  1. Which frequencies actually carry the weight-drop signal, and how that band
     narrows with depth -- i.e. why SNR dies near 530 m.
  2. Whether 20-50 Hz is the right working band or just an inherited habit.
  3. How the cemented and wireline fibers compare spectrally.

Two details that matter for reading these plots:

* The signal window follows the arrival, it is not fixed. Direct P moves out at
  ~2975 m/s (nano_velocity_scan.py), so a fixed window would compare the
  arrival at shallow depth against empty record at depth and manufacture a
  spectral trend that isn't there.

* Units are reconciled before any cross-fiber panel. Sintela delivers strain
  RATE; OptaSense delivers strain, and DASutils only differentiates it when
  readFile_HDF gets diff=True, which paired_stack_job.py does not pass. Left
  alone the two differ by a factor of omega -- a slope across the band, exactly
  the thing a spectral comparison is trying to measure. Differentiating the
  stored Deep stack here is equivalent and costs nothing.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, get_window

FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
NPZ = os.path.join(FIG_DIR, 'epoch_stacks_paired.npz')

DX_NANO, DX_DEEP = 1.26606202, 2.0419
PRE_S = 0.5
V_P = 2975.0                 # measured, nano_velocity_scan.py
SIG_S = 0.10                 # window following the arrival; long enough for the
                             # direct-P wavelet, short enough to exclude coda
DEPTHS_M = [150, 250, 350, 450, 550]
SPEC_DEPTH_M = 250

# A single 0.1 s window gives a 2-DOF periodogram -- roughly 100% scatter per
# frequency bin, which is enough to make one depth's spectrum cross another's at
# random. Averaging the power spectra over the channels within +/- this distance
# buys ~2*BAND_CH DOF. The arrival is moveout-corrected per channel before the
# average, so this stacks signal rather than smearing it.
AVG_HALFWIN_M = 20.0


def stack(d, key):
    n = d['n_common']
    good = n > 0
    w = n[good].astype(float)
    return np.tensordot(w, d[key][good], axes=(0, 0)) / w.sum(), int(w.sum())


def windowed_spectra(sec, dx, fs, i0):
    """Amplitude spectra of noise (pre-drop) and signal (following the arrival)."""
    n_ch = sec.shape[0]
    nsig = int(SIG_S * fs)
    z = np.arange(n_ch) * dx
    taper = get_window('hann', nsig)
    freq = np.fft.rfftfreq(nsig, 1 / fs)
    sig = np.full((n_ch, freq.size), np.nan)
    noi = np.full((n_ch, freq.size), np.nan)
    for c in range(n_ch):
        a = i0 + int(z[c] / V_P * fs)
        if a + nsig > sec.shape[1] or i0 - nsig < 0:
            continue
        # keep power, not amplitude -- power is what averages correctly
        sig[c] = np.abs(np.fft.rfft(sec[c, a:a + nsig] * taper)) ** 2
        noi[c] = np.abs(np.fft.rfft(sec[c, i0 - nsig:i0] * taper)) ** 2
    return freq, sig, noi, z


def band_average(psd, z, zq, half=AVG_HALFWIN_M):
    """Mean power spectrum over channels within +/- half of depth zq."""
    m = np.abs(z - zq) <= half
    if not m.any():
        return None
    return np.nanmean(psd[m], axis=0)


def main():
    d = np.load(NPZ)
    fs = float(d['fs'])
    i0 = int(PRE_S * fs)
    nano, ndrops = stack(d, 'nano_stacks')
    deep, _ = stack(d, 'deep_stacks')
    print(f'{ndrops} drops stacked; nano {nano.shape}, deep {deep.shape}')

    # strain -> strain rate, so the two fibers are in the same units
    deep_rate = np.gradient(deep, 1.0 / fs, axis=-1)

    f, sig_n, noi_n, z_n = windowed_spectra(nano, DX_NANO, fs, i0)
    f, sig_d, noi_d, z_d = windowed_spectra(deep_rate, DX_DEEP, fs, i0)

    # channel-averaged SNR profile: smooth both numerator and denominator first
    zq_grid = np.arange(60, 700, 10.0)
    snr_map = np.array([
        (band_average(sig_n, z_n, zq) / band_average(noi_n, z_n, zq))
        for zq in zq_grid])

    band = (f >= 5) & (f <= 120)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    # A: signal vs noise spectra at a few depths
    for zq in DEPTHS_M:
        ps = band_average(sig_n, z_n, zq)
        if ps is None:
            continue
        ax[0, 0].loglog(f[band], np.sqrt(ps[band]), lw=1.4, label=f'{zq} m')
    ax[0, 0].loglog(f[band], np.sqrt(np.nanmean(noi_n, axis=0))[band],
                    'k--', lw=1.3, label='noise (all channels)')
    ax[0, 0].set_xlabel('frequency (Hz)')
    ax[0, 0].set_ylabel(r'amplitude ($\mu\varepsilon$/s)')
    ax[0, 0].set_title('A  Cemented fiber: direct-P spectrum vs depth')
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(alpha=0.3, which='both')

    # B: SNR as a function of frequency and depth -- the usable band
    snr_db = 10 * np.log10(snr_map[:, band])
    im = ax[0, 1].pcolormesh(f[band], zq_grid, snr_db, cmap='magma',
                             vmin=0, vmax=25, shading='auto')
    ax[0, 1].contour(f[band], zq_grid, snr_db, levels=[9.54],
                     colors='c', linewidths=1.6)
    for fq in (20, 50):
        ax[0, 1].axvline(fq, color='w', ls=':', lw=1.0)
    ax[0, 1].invert_yaxis()
    ax[0, 1].set_xscale('log')
    ax[0, 1].set_xlabel('frequency (Hz)')
    ax[0, 1].set_ylabel('distance along fiber (m)')
    ax[0, 1].set_title('B  SNR (dB); cyan = SNR 3, shaded = 20-50 Hz')
    plt.colorbar(im, ax=ax[0, 1], label='SNR (dB)')

    # C: spectrogram of the stacked drop at one depth
    c = int(SPEC_DEPTH_M / DX_NANO)
    # short window so the direct-P wavelet is resolved, and referenced to the
    # pre-drop level per frequency so the arrival is visible rather than the
    # overall spectral shape
    ff, tt, Sxx = spectrogram(nano[c], fs=fs, nperseg=64, noverlap=60)
    pre = tt < PRE_S - 0.05
    ref = np.median(Sxx[:, pre], axis=1, keepdims=True)
    im = ax[1, 0].pcolormesh(tt - PRE_S, ff, 10 * np.log10(Sxx / (ref + 1e-30)),
                             cmap='magma', vmin=0, vmax=30, shading='auto')
    ax[1, 0].axhline(20, c='c', ls=':', lw=1)
    ax[1, 0].axhline(50, c='c', ls=':', lw=1)
    ax[1, 0].set_ylim(0, 150)
    ax[1, 0].set_xlim(-0.2, 1.0)
    ax[1, 0].set_xlabel('time after drop (s)')
    ax[1, 0].set_ylabel('frequency (Hz)')
    ax[1, 0].set_title(f'C  Spectrogram at {SPEC_DEPTH_M} m, dB above pre-drop')
    plt.colorbar(im, ax=ax[1, 0], label='dB above pre-drop level')

    # D: cemented vs wireline, units reconciled
    zq = 250
    for psd, zz, col, lab in [(sig_n, z_n, 'C0', 'cemented'),
                              (sig_d, z_d, 'C1', 'wireline')]:
        ps = band_average(psd, zz, zq)
        if ps is not None:
            ax[1, 1].loglog(f[band], np.sqrt(ps[band]), col, lw=1.5,
                            label=f'{lab} (signal)')
    for psd, zz, col, lab in [(noi_n, z_n, 'C0', 'cemented'),
                              (noi_d, z_d, 'C1', 'wireline')]:
        ps = band_average(psd, zz, zq)
        if ps is not None:
            ax[1, 1].loglog(f[band], np.sqrt(ps[band]), col, lw=1.0, alpha=0.45,
                            label=f'{lab} (noise)')
    ax[1, 1].set_xlabel('frequency (Hz)')
    ax[1, 1].set_ylabel(r'amplitude ($\mu\varepsilon$/s)')
    ax[1, 1].set_title(f'D  Both fibers at {zq} m, both as strain rate')
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].grid(alpha=0.3, which='both')

    fig.suptitle(f'SAFOD AWD spectral content -- {ndrops} stacked weight drops',
                 fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, 'awd_spectra.png')
    fig.savefig(out, dpi=140)

    # where does the band actually live?
    print('\n depth   peak_f   SNR>3 band      band-limited SNR (20-50 Hz)')
    for zq in DEPTHS_M:
        ps, pn = band_average(sig_n, z_n, zq), band_average(noi_n, z_n, zq)
        if ps is None:
            continue
        s = np.sqrt(ps / pn)
        pk = f[band][np.nanargmax(s[band])]
        ok = f[band][s[band] > 3]
        rng = f'{ok.min():5.1f}-{ok.max():5.1f} Hz' if ok.size else '    none    '
        b = (f >= 20) & (f <= 50)
        print(f'{zq:5d} m {pk:7.1f} Hz  {rng}   {np.nanmedian(s[b]):6.2f}')

    np.savez(os.path.join(FIG_DIR, 'awd_spectra.npz'),
             freq=f, sig_nano=sig_n, noise_nano=noi_n, z_nano=z_n,
             sig_deep=sig_d, noise_deep=noi_d, z_deep=z_d,
             snr_map=snr_map, snr_depth=zq_grid,
             n_drops=ndrops, v_p=V_P, sig_win_s=SIG_S,
             avg_halfwin_m=AVG_HALFWIN_M)
    print('\nSaved', out)


if __name__ == '__main__':
    main()
