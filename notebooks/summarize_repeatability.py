"""
Pull the headline repeatability numbers out of the saved AWD products, so the
claims on any poster/figure come from the arrays rather than from reading a PNG.

Reports, per fiber and band, over the depth interval where the direct-P is
actually above noise: waveform correlation between repeats, arrival-time lag
scatter, and NRMS. Also states the depth at which each metric fails.
"""
import os
import glob
import numpy as np

FIG_DIR = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026'
DX_NANO = 1.26606202
DX_DEEP = 2.0419

# Lellouch et al. 2019: the CEMENTED fiber ends at 864 m with a failed end loop,
# usable to 800 m. These bounds apply to the cemented ("shallow"/Nano) fiber only
# -- the wireline fiber is a separate installation reaching ~3.5 km, so applying
# them there would be meaningless.
FIBER_END_M = 864.0
FAILURE_M = 800.0
# Wireline turnaround located from data (deep_hairpin_turnaround_burst24.png):
# array ch 1702 = 3475 m along fiber. Its depth registration is unverified.
DEEP_TURNAROUND_M = 3475.0


def describe(path):
    d = np.load(path, allow_pickle=True)
    print(f'\n=== {os.path.basename(path)} ===')
    for k in d.files:
        v = d[k]
        if v.ndim == 0:
            print(f'  {k:32s} scalar = {v}')
        else:
            finite = np.isfinite(v) if v.dtype.kind == 'f' else np.ones(v.shape, bool)
            print(f'  {k:32s} shape {str(v.shape):16s} dtype {v.dtype} '
                  f'finite {finite.sum()}/{v.size}')
    return d


def profile_stats(d, dx, label, cemented):
    """Find any per-channel metric arrays and report where they hold up.

    These products store their own depth axis (`scan_depth` / `analysis_depth_m`)
    because the scan channels are subsampled across the fiber -- 126 scan points
    span ~900 m, not 126*dx. Use it; synthesising arange(n)*dx picks the wrong
    channels entirely.
    """
    depth = None
    for key in ('scan_depth', 'analysis_depth_m'):
        if key in d.files:
            depth = d[key]
            print(f'  depth axis from `{key}`: {depth.min():.1f} to '
                  f'{depth.max():.1f} m over {depth.size} scan channels')
            break
    if depth is None:
        print(f'  ({label}: no stored depth axis, skipping)')
        return
    n_ch_guess = depth.size
    if cemented and depth.max() > FAILURE_M:
        print(f'  !! {(depth > FIBER_END_M).sum()} scan channels lie past the '
              f'{FIBER_END_M:.0f} m fiber end, '
              f'{(depth > FAILURE_M).sum()} past the {FAILURE_M:.0f} m '
              f'failure depth (Lellouch et al. 2019)')
    elif not cemented:
        print(f'  (wireline fiber: Lellouch\'s bounds do not apply here. '
              f'Turnaround located\n   from data at {DEEP_TURNAROUND_M:.0f} m '
              f'along fiber; depth axis is ch*dx with zero\n   offset and is '
              f'NOT calibrated, since StartLocusIndex=1800.)')

    # Zone edges differ by fiber: the cemented one has a published failure depth,
    # the wireline one only has the data-derived turnaround.
    zones = ([(130, 400, 'signal 130-400 m'),
              (400, FAILURE_M, 'fading 400-800 m'),
              (FAILURE_M, 1e9, 'past failure depth')] if cemented else
             [(130, 400, 'signal 130-400 m'),
              (400, 1000, 'fading 400-1000 m'),
              (1000, 1e9, 'below 1000 m')])

    for k in d.files:
        v = d[k]
        if not (v.ndim == 1 and v.size == n_ch_guess and v.dtype.kind == 'f'):
            continue
        m = np.isfinite(v)
        if m.sum() < 10:
            continue
        # Quote the good zone and the dead zone separately -- lumping them is
        # exactly how noise gets reported as a measurement.
        for lo, hi, tag in zones:
            win = m & (depth > lo) & (depth < hi)
            if win.sum() < 5:
                continue
            print(f'  {k:24s} {tag:20s} n={win.sum():4d}  '
                  f'median {np.median(v[win]):9.4f}  '
                  f'p10 {np.percentile(v[win],10):9.4f}  '
                  f'p90 {np.percentile(v[win],90):9.4f}')


for path in sorted(glob.glob(os.path.join(FIG_DIR, '*repeatability*.npz'))
                   + glob.glob(os.path.join(FIG_DIR, '*dvv*.npz'))):
    d = describe(path)
    name = os.path.basename(path)
    cemented = 'deep' not in name          # 'shallow' files are the cemented fiber
    dx = DX_NANO if cemented else DX_DEEP
    profile_stats(d, dx, name, cemented)
