"""
Step 1a: extract the ground-truth traveltime curve from the PGSI check-shot gather.

`check shot top deployment.png` is a picked check-shot gather for the SAFOD main
hole: HOFFSET 45.72 m, WELL_DEP 6.68 to 2530.3 m. The first-break curve is already
drawn on it by the original processing, so digitising the raster gives t(z) without
re-picking waveforms.

This route was chosen because route 1b failed. Re-picking `1078.clean.txt` produced
incoherent first breaks under either channel ordering (`pgsi_ordering_check.png`),
with no coherent energy anywhere near the 3000-5000 m/s guides, and an arc peaking in
the middle of the array -- the signature of a source at depth, not a surface shot. So
that file is not delivering check-shot first breaks under the assumed 200 ms pre-shot
delay, and the picked raster is the more trustworthy product.

Axis calibration, measured from the image rather than assumed:
  * Depth: 12 tick marks at x = 101.5 ... 858.5 px, evenly spaced (69 px). Labels are
    NOT evenly spaced in metres (268.6, then mostly 228.6 = 750 ft, with two 213 m
    steps), so the x-axis is trace index and the composite gather has depth gaps
    between array positions. Depth is therefore interpolated between labelled ticks,
    not fitted with a single linear scale.
  * Time: 19 horizontal gridlines, row 130.0 = 50 ms through row 735.5 = 950 ms,
    giving 1.4864 ms/px. The plot's top frame at row 95 back-solves to t = -2 ms,
    i.e. 0 within a pixel -- an independent check that the calibration is right.

Writes t(z) plus an overlay figure so the extraction can be checked by eye rather
than trusted.
"""
import os
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PNG = ('/home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/'
       'pgsi_reference/check shot top deployment.png')
OUT = '/home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone'
FIG = os.path.join(OUT, 'figures')

# measured tick positions (px) and their WELL_DEP labels (m)
TICK_X = np.array([101.5, 170.5, 238.5, 307.5, 376.5, 445.5,
                   514.5, 583.5, 651.5, 720.5, 789.5, 858.5])
TICK_Z = np.array([6.68, 275.28, 503.88, 732.48, 960.77, 1189.4,
                   1402.7, 1631.1, 1859.7, 2088.3, 2301.7, 2530.3])

ROW_50MS, MS_PER_PX = 130.0, 900.0 / 605.5
TOP_ROW, BOT_ROW = 96, 745
GRIDLINES = np.array([95.0, 130.0, 164.0, 197.5, 231.0, 264.5, 298.0, 331.5, 365.0,
                      399.0, 433.5, 467.0, 500.5, 534.0, 568.0, 601.0, 635.0,
                      668.5, 702.0, 735.5])
HOFFSET = 45.72

DARK = 110          # grey level counting as "ink"

# A "first dark pixel" rule fails here: every trace is drawn as a thin dark
# vertical line all the way up the plot, so it triggers immediately and the pick
# jumps to the top of the frame (81% of points then needed monotonic repair).
# The first-break band is instead distinguished by HORIZONTAL continuity -- many
# adjacent columns going black together. Thin vertical lines score low on that,
# the arrival band scores high.
HWIN = 15           # horizontal window (px) for the ink-density measure
DENS = 0.55         # fraction of that window that must be ink


def t_of_row(r):
    return 50.0 + (r - ROW_50MS) * MS_PER_PX


def main():
    a = np.array(Image.open(PNG).convert('L')).astype(float)
    print(f'image {a.shape}')

    # Gridlines are dark across the full width and would be picked as the first
    # break in every column; blank them and interpolate over the gap.
    mask = np.zeros(a.shape[0], bool)
    for g in GRIDLINES:
        mask[int(round(g)) - 1:int(round(g)) + 2] = True
    print(f'masking {mask.sum()} gridline rows')

    dark = (a < DARK).astype(float)
    dark[mask, :] = np.nan                      # gridlines excluded, not counted
    # horizontal ink density
    k = np.ones(HWIN) / HWIN
    dens = np.apply_along_axis(
        lambda r: np.convolve(np.nan_to_num(r), k, 'same'), 1, dark)
    dens[mask, :] = 0.0

    xs = np.arange(int(TICK_X[0]), int(TICK_X[-1]) + 1)
    rows = np.full(xs.size, np.nan)
    for i, x in enumerate(xs):
        col = dens[TOP_ROW:BOT_ROW, x]
        hit = np.where(col > DENS)[0]
        if hit.size:
            rows[i] = TOP_ROW + hit[0]

    ok = np.isfinite(rows)
    print(f'{ok.sum()}/{xs.size} columns picked')

    z = np.interp(xs, TICK_X, TICK_Z)
    t = t_of_row(rows) * 1e-3                      # seconds

    # A wavefront cannot arrive earlier deeper. Enforce it, then report how much
    # had to be corrected -- a large correction would mean the extraction is bad.
    t_raw = t.copy()
    tm = np.maximum.accumulate(np.where(ok, t, -np.inf))
    n_fix = int(np.sum((tm > t_raw + 1e-9) & ok))
    print(f'monotonic correction applied to {n_fix}/{int(ok.sum())} points '
          f'({100*n_fix/max(ok.sum(),1):.1f}%)')

    zz, tt = z[ok], tm[ok]
    slant = np.sqrt(zz ** 2 + HOFFSET ** 2)
    print(f'\nt = {tt[0]*1e3:.1f} ms at {zz[0]:.1f} m  ->  '
          f'{tt[-1]*1e3:.1f} ms at {zz[-1]:.1f} m')
    for lo, hi in [(0, 500), (500, 1000), (1000, 1500), (1500, 2530)]:
        m = (zz >= lo) & (zz < hi)
        if m.sum() > 5:
            A = np.vstack([tt[m], np.ones(m.sum())]).T
            v = np.linalg.lstsq(A, slant[m], rcond=None)[0][0]
            print(f'  {lo:5d}-{hi:5d} m : Vp = {v:6.0f} m/s   (n={m.sum()})')

    np.savez(os.path.join(OUT, 'checkshot_traveltime.npz'),
             depth=zz, traveltime=tt, slant=slant, hoffset=HOFFSET,
             source='check shot top deployment.png (digitised)')

    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    ax[0].imshow(a, cmap='gray', vmin=0, vmax=255)
    ax[0].plot(xs[ok], rows[ok], 'r-', lw=1.2, label='digitised first break')
    ax[0].set(xlim=(80, 900), ylim=(760, 80), title='A  Overlay on the original raster')
    ax[0].legend(loc='lower left', fontsize=8)
    ax[0].axis('off')

    ax[1].plot(tt * 1e3, zz, 'k-', lw=1.6)
    ax[1].invert_yaxis()
    ax[1].set(xlabel='traveltime (ms)', ylabel='WELL_DEP (m)',
              title=f'B  Ground-truth t(z), HOFFSET {HOFFSET} m')
    ax[1].grid(alpha=0.3)
    fig.suptitle('SAFOD check-shot traveltime curve, digitised from the picked gather',
                 fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, 'checkshot_traveltime.png')
    fig.savefig(out, dpi=140)
    print('\nSaved', out)


if __name__ == '__main__':
    main()
