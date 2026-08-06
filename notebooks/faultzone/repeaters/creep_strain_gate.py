"""
Can the cemented fibre see fault creep? Feasibility gate for quasi-static strain.

THE ARGUMENT. A fibre cemented into a borehole IS a strainmeter. The SAFOD cemented
fibre spans 0-864 m depth about 1.8 km from the surface trace of an actively
creeping San Andreas. Parkfield already has borehole strainmeters that record creep
events routinely; every one of them is a POINT instrument. This is 900 gauges over
864 m of depth, so it measures a strain PROFILE.

That matters because the depth distribution of the strain field is what constrains
the depth extent of the slip. The Madden & Rowe creep workshop asks exactly that --
"Does the depth of the creep matter?" -- and notes the deepest instrument in use
anywhere is a 20 m inclinometer.

WHY THIS IS ONLY NOW POSSIBLE TO TRY. The archive turns out not to be high-passed:

    MinimumFrequency = 0.0      MaximumFrequency = 250.0
    RawData is int32 optical phase with a large DC offset (mean ~1.6e6)

Optical phase is proportional to STRAIN, not strain rate. DASutils.readFile_HDF
differentiates it when diff=True, which is what every script in this project has
used. Reading with diff=False keeps the quasi-static band, and no integration is
needed -- so there is no drift accumulated by integrating, only the optical drift
already present.

EXPECTED SIGNAL. Parkfield creep events are described in the workshop as pulses of
a few mm over about a day. For slip u at distance r the strain scale is ~u/r, so
1 mm at 1.8 km is ~5e-7 = 0.5 microstrain, and a few mm is a few microstrain.
LF-DAS resolves tens of nanostrain over hours in industrial use, so the signal
should sit 1-2 orders of magnitude above the floor IF drift is controlled.

DRIFT IS THE WHOLE PROBLEM, AND THE DEPTH ARRAY IS THE DEFENCE. Laser and
interrogator drift is common-mode: it moves all channels together. Real deformation
has a depth pattern. So the observable is the strain DIFFERENCE between depth
intervals, with the across-channel median removed. This is the same logic that
makes the interval dv/v measurement immune to origin time, applied to strain.

--------------------------------------------------------------------------------
GATE CRITERIA, FIXED BEFORE RUNNING. This is a feasibility test, not a creep
measurement, and it is allowed to fail.

 C1 THE RECORD MUST BE STABLE ENOUGH TO SEE ANYTHING. After common-mode removal,
    the residual strain noise over 1 hour must be below 100 nanostrain. If it is
    not, a 0.5 microstrain creep signal is unreachable and the direction closes.

 C2 DRIFT MUST BE COMMON-MODE, NOT PER-CHANNEL. Report the fraction of total
    variance removed by subtracting the across-channel median. If common-mode
    removal takes out most of the drift, the array defence works. If each channel
    drifts independently, it does not, and no amount of averaging helps.

 C3 THE RESIDUAL MUST BE DEPTH-COHERENT, NOT WHITE ACROSS CHANNELS. Neighbouring
    channels share rock; if the residual decorrelates from one channel to the next
    it is instrument noise rather than deformation. Measured as the correlation
    length in channels, which must exceed the 16-channel gauge length to mean
    anything.

 C4 NO INTERPRETATION WITHOUT AN INDEPENDENT EVENT. Any transient found here is
    NOT a creep event until it coincides with one on a Parkfield creepmeter or
    borehole strainmeter. This project has already produced one result that
    survived every internal check and died to an external control; strain
    transients in DAS have too many instrumental explanations to be trusted alone.

WHAT FAILURE MEANS. If C1-C3 fail the honest statement is that the cemented fibre
as archived cannot do quasi-static strain, which is worth writing down because the
obvious reading of "MinimumFrequency = 0.0" is that it can.
--------------------------------------------------------------------------------

Reads a bounded number of consecutive files, so it does not hammer Lustre.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from g1_coda_snr import load_manifest                            # noqa: E402

for p in ['/home/groups/edunham/nberrios/safod_das/DAS-utilities/python',
          '/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python']:
    if os.path.exists(p) and p not in sys.path:
        sys.path.insert(0, p)
import DASutils                                                  # noqa: E402
import h5py                                                      # noqa: E402

CH_LO, CH_HI = 23, 896          # G0 usable range
DX = 1.0209523
GAUGE_CH = 16
N_FILES = int(os.environ.get('CREEP_NFILES', 360))     # 360 x 60 s = 6 h
START = os.environ.get('CREEP_START', '2024-06-01')
DECIM = 60                      # samples kept per 60 s file -> 1 Hz


def read_strain(fn):
    """One 60 s file as STRAIN (phase), decimated to 1 Hz. No differentiation."""
    try:
        with h5py.File(fn, 'r') as h:
            g = [k for k in h['Acquisition'].keys() if 'Raw' in k][0]
            d = h[f'Acquisition/{g}/RawData']
            X = np.asarray(d[:, :], dtype=np.float64)
            if X.shape[0] < X.shape[1]:
                X = X.T                       # -> (channel, time)
            else:
                X = X.T if X.shape[1] == 900 else X
            t0 = np.asarray(h[f'Acquisition/{g}/RawDataTime'][:2], dtype=np.int64)
    except Exception as e:
        return None, None, str(e)[:60]
    if X.shape[0] != 900:
        X = X.T
    n = X.shape[1] // DECIM * DECIM
    if n == 0:
        return None, None, 'short'
    # block mean = anti-aliased decimation to 1 Hz; we want the LOW band, so a
    # boxcar average is the right operation and cheap
    Y = X[:, :n].reshape(X.shape[0], -1, DECIM).mean(axis=2)
    return Y, t0[0], None


def main():
    db = load_manifest()
    db = db.sort_values('t0').reset_index(drop=True)
    sel = db[db['t0'] >= pd.Timestamp(START, tz='UTC')].head(N_FILES)
    sel = sel[sel['fn'].map(os.path.exists)]
    print(f'{len(sel)} consecutive files from {START} '
          f'({len(sel)/60:.1f} h at 60 s each)', flush=True)
    if len(sel) < 30:
        print('too few files'); return

    cols, times, bad = [], [], 0
    for k, (_, r) in enumerate(sel.iterrows()):
        Y, t0, err = read_strain(r['fn'])
        if Y is None:
            bad += 1; continue
        cols.append(Y[CH_LO:CH_HI])
        times.append(r['t0'])
        if k % 60 == 0:
            print(f'  {k}/{len(sel)}  {r["t0"]}', flush=True)
    if not cols:
        print('nothing read'); return
    S = np.concatenate(cols, axis=1)            # (channel, time) at 1 Hz
    nch, nt = S.shape
    print(f'\nstrain matrix {nch} channels x {nt} s '
          f'({nt/3600:.2f} h), {bad} files unreadable', flush=True)

    # phase offset per channel is arbitrary -- only changes matter
    S = S - S[:, :1]
    raw_var = float(np.var(S))

    # C2: common-mode removal
    cm = np.median(S, axis=0)
    R = S - cm[None, :]
    res_var = float(np.var(R))
    removed = 1.0 - res_var / max(raw_var, 1e-30)
    print(f'\nC2 common-mode: {100*removed:.1f}% of variance removed by the '
          f'across-channel median')

    # convert to strain: phase counts -> strain needs the interrogator constant,
    # which is not in the header. Report in NATIVE COUNTS and as a relative
    # quantity, and state the conversion is missing rather than inventing one.
    print(f'   raw std {np.std(S):.4g} counts, residual std {np.std(R):.4g} counts')

    # C1: 1-hour residual stability
    hr = min(3600, nt)
    seg = R[:, :hr]
    drift_hr = float(np.median(np.abs(seg[:, -1] - seg[:, 0])))
    noise_hr = float(np.median(np.std(np.diff(seg, axis=1), axis=1)))
    print(f'\nC1 1-hour stability (counts): median |end-start| = {drift_hr:.4g}, '
          f'sample-to-sample noise = {noise_hr:.4g}')
    print(f'   ratio drift/noise = {drift_hr/max(noise_hr,1e-30):.1f}')

    # C3: depth coherence of the residual
    Rz = R - R.mean(axis=1, keepdims=True)
    Rz /= np.maximum(Rz.std(axis=1, keepdims=True), 1e-30)
    lags = np.arange(1, min(120, nch // 3))
    cvals = [float(np.mean(np.sum(Rz[:-L] * Rz[L:], axis=1) / Rz.shape[1]))
             for L in lags]
    cvals = np.array(cvals)
    below = np.flatnonzero(cvals < 1 / np.e)
    Lc = int(lags[below[0]]) if below.size else int(lags[-1])
    print(f'\nC3 depth coherence: 1/e correlation length = {Lc} channels '
          f'({Lc*DX:.1f} m); gauge length is {GAUGE_CH} channels')
    print(f'   {"PASS" if Lc > GAUGE_CH else "FAIL"} -- '
          f'{"structure beyond the gauge" if Lc > GAUGE_CH else "no structure beyond the gauge; this is instrument noise"}')

    np.savez_compressed(os.path.join(HERE, 'creep_strain_gate.npz'),
                        S=S.astype(np.float32), R=R.astype(np.float32),
                        cm=cm.astype(np.float32),
                        times=np.array([str(t) for t in times]),
                        removed=removed, Lc=Lc)

    print('\n' + '=' * 70)
    print('VERDICT against criteria fixed before running')
    c2 = removed > 0.5
    c3 = Lc > GAUGE_CH
    print(f'  C2 drift is common-mode  : {100*removed:5.1f}% removed   '
          f'{"PASS" if c2 else "FAIL"}')
    print(f'  C3 residual depth-coherent: Lc = {Lc} ch          '
          f'{"PASS" if c3 else "FAIL"}')
    print(f'  C1 absolute strain floor  : NOT EVALUABLE -- the phase-to-strain')
    print(f'     constant is not in the HDF5 header. Everything above is in')
    print(f'     native counts. Get the interrogator constant from OptaSense or')
    print(f'     calibrate against the active-source data before quoting strain.')
    print(f'  C4 external validation    : NOT DONE. No transient here is a creep')
    print(f'     event until it coincides with a Parkfield creepmeter.')
    print()
    if c2 and c3:
        print('  -> FEASIBLE SO FAR. The array defence against drift works and the')
        print('     residual has depth structure. Next: a longer window, the strain')
        print('     calibration, and a creepmeter comparison.')
    else:
        print('  -> NOT FEASIBLE as processed. Say so plainly: the obvious reading')
        print('     of MinimumFrequency = 0.0 is that quasi-static strain is')
        print('     available, and that reading would be wrong.')

    fig, ax = plt.subplots(2, 2, figsize=(14, 8))
    tt = np.arange(nt) / 3600.0
    ax[0, 0].plot(tt, cm, 'C0', lw=.8)
    ax[0, 0].set(xlabel='hours', ylabel='counts', title='A  common mode (drift)')
    ax[0, 0].grid(alpha=.3)
    v = np.percentile(np.abs(R), 98)
    im = ax[0, 1].imshow(R, aspect='auto', vmin=-v, vmax=v, cmap='RdBu_r',
                         extent=[0, tt[-1], (CH_HI - 1) * DX, CH_LO * DX])
    ax[0, 1].set(xlabel='hours', ylabel='depth (m)',
                 title='B  residual strain after common-mode removal')
    plt.colorbar(im, ax=ax[0, 1], fraction=.046)
    ax[1, 0].plot(lags, cvals, 'C2')
    ax[1, 0].axhline(1 / np.e, color='k', ls='--', label='1/e')
    ax[1, 0].axvline(GAUGE_CH, color='C3', ls=':', label='gauge length')
    ax[1, 0].set(xlabel='channel lag', ylabel='correlation',
                 title=f'C  C3: depth coherence (Lc = {Lc} ch)')
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=.3)
    for frac, lab in ((0.15, 'shallow'), (0.5, 'mid'), (0.85, 'deep')):
        k = int(frac * nch)
        ax[1, 1].plot(tt, R[k], lw=.7, label=f'{lab} ({(CH_LO+k)*DX:.0f} m)')
    ax[1, 1].set(xlabel='hours', ylabel='counts',
                 title='D  residual at three depths')
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=.3)
    fig.tight_layout()
    p = os.path.join(HERE, 'creep_strain_gate.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
