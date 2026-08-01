"""
Far-end liveness check for the Nano fiber.

Question: is the full ~864 m of Nano fiber actually sensing, or does the
record go dead somewhere before the last channel (unbonded tail, turnaround,
gauge-length taper, bad splice)?

Method: for each channel, compare energy in a signal window after the weight
drop against energy in the pre-drop noise window, averaged over epochs. A
live channel has signal/noise > 1 and a noise floor comparable to its
neighbors. A dead tail shows up as SNR collapsing to 1 and/or the noise RMS
dropping off a cliff (no strain being recorded at all).

Reads the stacks already on disk -- no re-reading of .pb files.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
NPZ = os.path.join(FIG_DIR, 'epoch_stacks_paired.npz')

PRE_S, POST_S = 0.5, 3.0
SIG_S = 0.6          # signal window length after the drop

d = np.load(NPZ)
nano = d['nano_stacks']            # (n_ep, n_ch, n_t)
fs = float(d['fs'])
dx = float(d['dx_nano'])
n_common = d['n_common']
n_ep, n_ch, n_t = nano.shape
print(f'nano stacks {nano.shape}  fs={fs}  dx={dx:.4f} m  '
      f'-> fiber length {n_ch*dx:.1f} m')

i_drop = int(PRE_S * fs)
noise = nano[:, :, :i_drop]                              # 0.5 s before drop
sig = nano[:, :, i_drop:i_drop + int(SIG_S * fs)]        # 0.6 s after drop

# keep only epochs that actually stacked drops
good = n_common > 0
print(f'{good.sum()}/{n_ep} epochs with common drops')
noise, sig = noise[good], sig[good]

noise_rms = np.sqrt(np.mean(noise ** 2, axis=2))    # (n_ep_good, n_ch)
sig_rms = np.sqrt(np.mean(sig ** 2, axis=2))
snr = sig_rms / np.where(noise_rms > 0, noise_rms, np.nan)

noise_med = np.nanmedian(noise_rms, axis=0)
sig_med = np.nanmedian(sig_rms, axis=0)
snr_med = np.nanmedian(snr, axis=0)

depth = np.arange(n_ch) * dx

# Where does the signal die? SNR decays smoothly into a ~1.1 asymptote, so a
# bare threshold picks up noise wobble hundreds of channels into the dead zone.
# Take the last channel *before* SNR stays under threshold for good instead.
THRESH = 1.5
live = snr_med > THRESH
last_live = -1
for c in range(n_ch - 1, -1, -1):
    if live[c] and live[max(0, c - 20):c + 1].mean() > 0.5:
        last_live = c
        break
print(f'\nsignal detection limit (median SNR > {THRESH}, sustained): '
      f'ch {last_live} = {last_live*dx:.1f} m along fiber')

# Is the tail dead fiber, or live fiber with no signal? Compare the noise floor
# past the detection limit against the mid-fiber noise floor. Dead/unbonded
# fiber loses the noise too; live fiber keeps recording ambient strain.
if 0 < last_live < n_ch - 20:
    tail = np.nanmedian(noise_med[last_live + 1:])
    body = np.nanmedian(noise_med[:last_live + 1])
    print(f'noise floor  body {body:.4g}  tail {tail:.4g}  '
          f'(ratio {tail/body:.2f}) -> '
          f'{"tail still live" if 0.3 < tail/body < 3 else "tail looks dead"}')

# Hairpin test: if the fiber doubles back, channels either side of the midpoint
# sample the same depths and their SNR profiles should mirror each other.
mid = n_ch // 2
ks = np.arange(1, min(mid, n_ch - mid))
print('\nhairpin mirror test (SNR down-leg vs up-leg about ch %d):' % mid)
for k in [15, 60, 150, 300]:
    if mid + k < n_ch:
        print(f'  k={k:4d}  ch{mid-k:4d} {snr_med[mid-k]:7.2f}   '
              f'ch{mid+k:4d} {snr_med[mid+k]:6.2f}')

print('\n  ch   fiber_m   noise_rms   sig_rms    SNR')
for c in list(range(0, 60, 10)) + list(range(n_ch - 80, n_ch, 10)) + [n_ch - 1]:
    print(f'{c:5d} {c*dx:8.1f} {noise_med[c]:11.4g} {sig_med[c]:9.4g} '
          f'{snr_med[c]:7.2f}')

fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
ax[0].semilogy(depth, noise_med, lw=0.9, color='0.35')
ax[0].set_ylabel('pre-drop noise RMS')
ax[0].set_title('Nano fiber: is the far end live?')
ax[1].semilogy(depth, sig_med, lw=0.9, color='C0')
ax[1].set_ylabel(f'signal RMS (0-{SIG_S}s)')
ax[2].plot(depth, snr_med, lw=0.9, color='C3')
ax[2].axhline(1.5, ls='--', c='0.5', lw=0.8)
ax[2].set_ylabel('median SNR')
ax[2].set_xlabel('distance along fiber (m)')
for a in ax:
    a.grid(alpha=0.3)
    if last_live >= 0:
        a.axvline(last_live * dx, ls=':', c='C2', lw=1.2)
fig.tight_layout()
out = os.path.join(FIG_DIR, 'nano_far_end_liveness.png')
fig.savefig(out, dpi=140)
print('\nSaved', out)
