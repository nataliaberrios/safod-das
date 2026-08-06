"""
Deep-cable quasi-static strain: a fibre that CROSSES the creeping San Andreas.

WHY THE DEEP CABLE AND NOT THE CEMENTED ONE. The cemented fibre (0-864 m) sits about
1.8 km from the fault trace and never reaches it, so creep only shows up as a faint
far-field stretch. The deep fibre is 3200 loci at 2.04 m = 6534 m of glass in a
~3.4 km hole, i.e. a HAIRPIN that goes down past the SAFOD damage zone
(3150-3414 m, SDZ at 3192 m) and comes back up. A fibre crossing a slipping fault
records the slip directly, as a step localised at the crossing, rather than as a
distant elastic strain.

Every creep instrument in the Madden & Rowe workshop is a surface measurement, and
the deepest instrument named anywhere in it is a 20 m borehole inclinometer.

WHAT THE DATA IS. LF_DAS_1Hz_serial: 1249 files, 2026-03-28 to 2026-05-01, each
2050 channels x 1800 s at 1 Hz. About 26 days continuous, ~18 GB total, so the whole
record fits in one job.

    per-channel mean |value| 1.1e8   -> ABSOLUTE phase, i.e. strain, not rate
    30-min drift, median      5.2e7   -> ~30x the within-file spread of 1.8e6

So the signal band is present and drift is the whole problem.

THREE DEFENCES AGAINST DRIFT, IN INCREASING ORDER OF STRENGTH.

 1. Common-mode removal. Laser and interrogator drift moves all channels together;
    deformation has a depth pattern. Subtract the across-channel median.

 2. Depth coherence. Neighbouring channels share rock. A residual that decorrelates
    channel-to-channel is instrument noise regardless of how large it is.

 3. THE HAIRPIN, which is the one this geometry uniquely allows. The fibre passes
    each depth TWICE, once descending and once ascending. Real ground strain at a
    given depth must appear on BOTH legs with the same sign and size. An optical or
    interrogator artefact is a function of distance ALONG THE FIBRE and will not
    respect the fold. This is a genuine physical control, not a statistical one,
    and it is the strongest evidence available here.

--------------------------------------------------------------------------------
GATES, FIXED BEFORE RUNNING.

 D1 FIND THE FOLD FROM THE DATA, NOT FROM METADATA. The turnaround channel is
    located by maximising the correlation between the descending and ascending
    legs of the channel-wise strain pattern. awd_clean records deep-cable depths
    as "explicitly provisional", so the fold must be measured. If no fold position
    produces a clear correlation peak, the hairpin control is unavailable and
    everything below is weaker.

 D2 COMMON-MODE MUST DOMINATE THE DRIFT. Report the fraction of variance removed
    by subtracting the across-channel median. If per-channel drift survives, the
    array defence fails.

 D3 RESIDUAL MUST BE DEPTH-COHERENT beyond the gauge length (10.21 m gauge /
    2.04 m spacing = 5 channels).

 D4 LEG AGREEMENT. Correlation between the two legs at matched depth, after
    common-mode removal. This is the number that decides whether any transient is
    ground motion.

 D5 NO INTERPRETATION WITHOUT AN EXTERNAL EVENT. A transient is not a creep event
    until it coincides with one on a Parkfield creepmeter or borehole strainmeter.

NOT CALIBRATED, AND SAID PLAINLY: the phase-to-strain constant is not in these
files. Everything is reported in native units. Do not quote microstrain from this
script.
--------------------------------------------------------------------------------
"""
import os
import sys
import glob
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
LF = '/oak/stanford/groups/ettore88/data/SAFOD/SAFOD-deep/LF_DAS_1Hz_serial'
DX_DEEP = 2.0419047
GAUGE_CH = 5                      # 10.21 m gauge / 2.04 m spacing
DECIM = int(os.environ.get('DEEP_DECIM', 60))     # 1 Hz -> 1/min
NFILE = int(os.environ.get('DEEP_NFILE', 0))      # 0 = all


def load_all():
    fs = sorted(glob.glob(os.path.join(LF, '*.h5')))
    if NFILE:
        fs = fs[:NFILE]
    print(f'{len(fs)} LF files', flush=True)
    cols, stamps = [], []
    for k, f in enumerate(fs):
        try:
            with h5py.File(f, 'r') as h:
                D = np.asarray(h['Data'][:], dtype=np.float64)
        except Exception:
            continue
        if D.ndim != 2:
            continue
        n = D.shape[1] // DECIM * DECIM
        if n == 0:
            continue
        cols.append(D[:, :n].reshape(D.shape[0], -1, DECIM).mean(axis=2))
        stamps.append(os.path.basename(f)[6:-3])
        if k % 200 == 0:
            print(f'  {k}/{len(fs)}', flush=True)
    if not cols:
        return None, None
    nch = min(c.shape[0] for c in cols)
    return np.concatenate([c[:nch] for c in cols], axis=1), stamps


def find_fold(profile):
    """D1: locate the hairpin turnaround by leg-to-leg correlation."""
    n = profile.size
    best, bestc = None, -2.0
    curve = []
    for fold in range(int(0.25 * n), int(0.75 * n)):
        m = min(fold, n - fold)
        if m < 100:
            curve.append(np.nan); continue
        a = profile[fold - m:fold]
        b = profile[fold:fold + m][::-1]        # ascending leg, depth-reversed
        a = a - a.mean(); b = b - b.mean()
        d = np.sqrt((a ** 2).sum() * (b ** 2).sum())
        c = float((a * b).sum() / d) if d > 0 else np.nan
        curve.append(c)
        if np.isfinite(c) and c > bestc:
            bestc, best = c, fold
    return best, bestc, np.array(curve), int(0.25 * n)


def main():
    S, stamps = load_all()
    if S is None:
        print('no data'); return
    nch, nt = S.shape
    print(f'\nstrain matrix {nch} channels x {nt} samples '
          f'({nt*DECIM/86400:.2f} days at {DECIM}s)', flush=True)

    S = S - S[:, :1]                       # phase offsets are arbitrary
    raw_var = float(np.var(S))

    # D2 common-mode
    cm = np.median(S, axis=0)
    R = S - cm[None, :]
    removed = 1.0 - float(np.var(R)) / max(raw_var, 1e-30)
    print(f'\nD2 common-mode removes {100*removed:.1f}% of variance')
    print(f'   raw std {np.std(S):.4g}, residual std {np.std(R):.4g} (native units)')

    # D1 fold, from the time-averaged residual profile
    prof = np.nanstd(R, axis=1)
    fold, foldc, curve, off = find_fold(prof)
    print(f'\nD1 hairpin fold: channel {fold} (correlation {foldc:+.3f})')
    print(f'   implies {fold*DX_DEEP:.0f} m of fibre to the turnaround, '
          f'{nch*DX_DEEP:.0f} m total')
    print(f'   {"PASS" if foldc > 0.3 else "FAIL"} -- '
          f'{"fold located" if foldc > 0.3 else "no clear fold; hairpin control unavailable"}')

    # D3 depth coherence
    Rz = R - R.mean(axis=1, keepdims=True)
    Rz /= np.maximum(Rz.std(axis=1, keepdims=True), 1e-30)
    lags = np.arange(1, min(200, nch // 4))
    cv = np.array([float(np.mean(np.sum(Rz[:-L] * Rz[L:], axis=1) / Rz.shape[1]))
                   for L in lags])
    bel = np.flatnonzero(cv < 1 / np.e)
    Lc = int(lags[bel[0]]) if bel.size else int(lags[-1])
    print(f'\nD3 depth coherence: 1/e length {Lc} channels ({Lc*DX_DEEP:.0f} m); '
          f'gauge is {GAUGE_CH} channels')
    print(f'   {"PASS" if Lc > GAUGE_CH else "FAIL"}')

    # D4 leg agreement through time
    legc = np.nan
    if fold and foldc > 0.3:
        m = min(fold, nch - fold)
        A = R[fold - m:fold]
        B = R[fold:fold + m][::-1]
        a = A - A.mean(axis=1, keepdims=True)
        b = B - B.mean(axis=1, keepdims=True)
        num = (a * b).sum(axis=1)
        den = np.sqrt((a ** 2).sum(axis=1) * (b ** 2).sum(axis=1))
        per = num / np.maximum(den, 1e-30)
        legc = float(np.nanmedian(per))
        print(f'\nD4 leg agreement at matched depth: median r = {legc:+.3f} '
              f'over {m} depth pairs')
        print(f'   {"PASS" if legc > 0.3 else "FAIL"} -- '
              f'{"both legs see the same thing: consistent with ground strain" if legc > 0.3 else "legs disagree: signal is a function of fibre distance, i.e. instrumental"}')

    np.savez_compressed(os.path.join(HERE, 'deep_creep_strain.npz'),
                        R=R.astype(np.float32), cm=cm.astype(np.float32),
                        prof=prof, fold=fold if fold else -1, foldc=foldc,
                        Lc=Lc, legc=legc, removed=removed, decim=DECIM,
                        stamps=np.array(stamps))

    print('\n' + '=' * 70)
    print('VERDICT')
    ok = (foldc > 0.3) and (Lc > GAUGE_CH) and (np.isfinite(legc) and legc > 0.3)
    print(f'  D1 fold located          : {foldc:+.3f}  {"PASS" if foldc>0.3 else "FAIL"}')
    print(f'  D2 common-mode dominant  : {100*removed:.1f}%')
    print(f'  D3 depth-coherent        : Lc={Lc} ch  {"PASS" if Lc>GAUGE_CH else "FAIL"}')
    print(f'  D4 both legs agree       : {legc:+.3f}  '
          f'{"PASS" if np.isfinite(legc) and legc>0.3 else "FAIL"}')
    print(f'  D5 external validation   : NOT DONE -- no transient is a creep event')
    print(f'     until it matches a Parkfield creepmeter.')
    print(f'  CALIBRATION              : phase-to-strain constant absent from the')
    print(f'     files. Native units only. Do not quote microstrain.')
    print()
    print('  -> ' + ('DEEP FIBRE CARRIES DEPTH-COHERENT DEFORMATION SIGNAL that '
                     'survives\n     the hairpin test. Next: creepmeter comparison '
                     'and calibration.'
                     if ok else
                     'NOT ESTABLISHED. Report which gate failed and what that '
                     'bounds.'))

    fig, ax = plt.subplots(2, 2, figsize=(14, 8))
    tt = np.arange(nt) * DECIM / 86400.0
    ax[0, 0].plot(tt, cm, lw=.7)
    ax[0, 0].set(xlabel='days', ylabel='native units', title='A  common mode (drift)')
    ax[0, 0].grid(alpha=.3)
    v = np.percentile(np.abs(R), 98)
    im = ax[0, 1].imshow(R, aspect='auto', vmin=-v, vmax=v, cmap='RdBu_r',
                         extent=[0, tt[-1], nch, 0])
    if fold:
        ax[0, 1].axhline(fold, color='k', lw=1, ls='--')
    ax[0, 1].set(xlabel='days', ylabel='channel', title='B  residual (dashed = fold)')
    plt.colorbar(im, ax=ax[0, 1], fraction=.046)
    ax[1, 0].plot(np.arange(off, off + curve.size), curve, 'C2')
    if fold:
        ax[1, 0].axvline(fold, color='C3', ls='--', label=f'fold {fold}')
        ax[1, 0].legend(fontsize=8)
    ax[1, 0].set(xlabel='candidate fold channel', ylabel='leg correlation',
                 title=f'C  D1: hairpin fold ({foldc:+.2f})')
    ax[1, 0].grid(alpha=.3)
    ax[1, 1].plot(prof, np.arange(nch), lw=.7)
    if fold:
        ax[1, 1].axhline(fold, color='k', ls='--')
    ax[1, 1].invert_yaxis()
    ax[1, 1].set(xlabel='residual std (native)', ylabel='channel',
                 title='D  strain amplitude along fibre')
    ax[1, 1].grid(alpha=.3)
    fig.tight_layout()
    p = os.path.join(HERE, 'deep_creep_strain.png')
    fig.savefig(p, dpi=140)
    print(f'\nwrote {p}')


if __name__ == '__main__':
    main()
