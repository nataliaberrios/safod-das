#!/usr/bin/env python3
"""Does tapering the F-K wedge remove the ringing that manufactures moveout?

THE PROBLEM.  `ambient_fk_transfer_test.fk_filter` applies a binary indicator in
the 2-D Fourier domain:

    mask = (af>=5)&(af<=20)&(v>=2500)&(v<=4500)&(ak>0)

That is a brick wall, not a filter design.  Multiplication by a sharp-edged
region in (f,k) is convolution in (t,x) with the region's inverse transform,
which for a boxcar has long sinc tails oriented along the wedge's own
velocities.  The operator therefore stamps a dipping texture onto ANY input.
`ambient_fk_white_noise_v1` measured the consequence: the narrow 2.8-3.8 km/s
fan returns |score| 0.644 at 3200 m/s from pure Gaussian noise, against 0.703
from real data.

THE FIX TESTED HERE.  Replace the indicator with a raised-cosine (Tukey) taper
on both the frequency edges and the velocity edges, so the mask goes smoothly
0 -> 1 -> 0 over declared transition widths.  A smoother spectral window has a
more compact kernel and therefore less ringing.

WHAT IS MEASURED.  Three things, in increasing order of what matters:

  1. Kernel compactness.  The impulse response of the mask itself -- no data
     involved.  Fraction of kernel energy inside the central lobe, and the decay
     of the tails.  This is the direct, data-free measure of ringing.
  2. White-noise floor.  Gaussian noise through the complete operator.  If
     tapering works, the score it manufactures from noise should fall.
  3. Real-data score, and the margin between the two.  A taper is only useful if
     it lowers the noise floor by MORE than it lowers the real score.

WHAT THIS CANNOT FIX, stated up front.  Tapering reduces ringing; it does not
make a velocity-selective filter into evidence for that velocity.  Any fan, hard
or tapered, passes only energy consistent with its own velocity band, so its
output necessarily resembles that moveout.  The decisive gate remains the
pre-filter channel scramble -- real spectra, real nonstationarity, spatial order
destroyed -- which the production fan fails.  This script therefore reports a
margin, not a verdict.

Output: ambient_fk_taper_test.{npz,png,txt}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, detrend, sosfiltfilt

HERE = Path(__file__).resolve().parent
STEM = HERE / "ambient_fk_taper_test"

FMIN, FMAX = 5.0, 20.0
VMIN, VMAX = 2500.0, 4500.0
V_REF = 3200.0
SEED = 20260814


def ramp(x):
    """Raised-cosine 0->1 on x in [0,1], clipped outside."""
    x = np.clip(x, 0.0, 1.0)
    return 0.5 * (1.0 - np.cos(np.pi * x))


def build_mask(nch, nt, fs, dx, f_taper_hz, v_taper_frac, branch="negative"):
    """Signed velocity wedge, optionally raised-cosine tapered on both edges.

    f_taper_hz = 0 and v_taper_frac = 0 reproduces the production brick wall
    exactly, so the hard case is this same code path with zero widths.
    """
    f = np.fft.fftfreq(nt, 1.0 / fs)
    k = np.fft.fftfreq(nch, dx)
    K, F = np.meshgrid(k, f, indexing="ij")
    af, ak = np.abs(F), np.abs(K)
    v = af / np.maximum(ak, 1e-12)

    if f_taper_hz <= 0:
        wf = ((af >= FMIN) & (af <= FMAX)).astype(float)
    else:
        wf = np.minimum(ramp((af - FMIN) / f_taper_hz + 1.0) *
                        ramp((FMAX - af) / f_taper_hz + 1.0), 1.0)
        wf = np.where((af >= FMIN - f_taper_hz) & (af <= FMAX + f_taper_hz), wf, 0.0)

    if v_taper_frac <= 0:
        wv = ((v >= VMIN) & (v <= VMAX)).astype(float)
    else:
        lo, hi = VMIN * (1 - v_taper_frac), VMAX * (1 + v_taper_frac)
        wv = np.minimum(ramp((v - lo) / (VMIN - lo)), ramp((hi - v) / (hi - VMAX)))
        wv = np.clip(wv, 0.0, 1.0)
        wv = np.where((v >= lo) & (v <= hi), wv, 0.0)

    m = wf * wv * (ak > 0)
    sign = (F * K) < 0 if branch == "negative" else (F * K) > 0
    return m * sign


def kernel_metrics(mask):
    """Impulse response of the mask: how concentrated is it?"""
    ker = np.fft.fftshift(np.real(np.fft.ifft2(mask)))
    e = ker ** 2
    tot = e.sum()
    c0, c1 = ker.shape[0] // 2, ker.shape[1] // 2
    out = {}
    for half in (2, 4, 8, 16):
        blk = e[c0 - half:c0 + half + 1, c1 - half:c1 + half + 1].sum()
        out["energy_within_%d" % half] = float(blk / tot)
    # tail level: median |kernel| outside a 32-sample box, relative to the peak
    m = np.ones_like(e, dtype=bool)
    m[c0 - 32:c0 + 33, c1 - 32:c1 + 33] = False
    out["tail_median_rel_peak"] = float(np.median(np.abs(ker[m])) / np.abs(ker).max())
    return ker, out


def preprocess(x, fs):
    x = detrend(x, axis=1, type="linear")
    sos = butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos")
    x = sosfiltfilt(sos, x, axis=1)
    n = max(3, int(5.0 * fs))
    m = uniform_filter1d(np.abs(x), size=n, axis=1, mode="nearest")
    floor = np.percentile(m, 5, axis=1, keepdims=True) * 0.1 + 1e-12
    return x / np.maximum(m, floor)


def correlate_top(x, fs, dx, max_lag=0.35):
    """Channel-0 virtual source against receivers every 50 m."""
    nch, nt = x.shape
    targets = [int(round(50.0 * j / dx)) for j in range(1, 15)]
    targets = [t for t in targets if t < nch]
    ml = int(round(max_lag * fs))
    lags = np.arange(-ml, ml + 1) / fs
    nfft = 1 << int(np.ceil(np.log2(2 * nt - 1)))
    A = np.fft.rfft(x[0], n=nfft)
    B = np.fft.rfft(x[targets], n=nfft, axis=1)
    C = np.fft.irfft(np.conj(A)[None, :] * B, n=nfft, axis=1)
    C = np.concatenate((C[:, -ml:], C[:, :ml + 1]), axis=1)
    den = np.sqrt(np.sum(x[0] ** 2) * np.sum(x[targets] ** 2, axis=1))
    ok = den > 0
    C[ok] /= den[ok][:, None]
    return lags, C, np.asarray(targets) * dx


def score_at(top, lags, dist, v):
    return float(np.median([r[np.argmin(np.abs(lags - d / v))] for r, d in zip(top, dist)]))


def apply_mask(x, mask):
    return np.fft.ifft2(np.fft.fft2(x) * mask).real


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    rng = np.random.default_rng(SEED)
    fs, dx = 250.0, 2.0419046878814697      # post-2x-decimation sampling, as production
    nch, nt = 400, 3000                     # 400 channels, 12 s

    configs = [
        ("hard (production)", 0.0, 0.0),
        ("taper f=1 Hz, v=5%", 1.0, 0.05),
        ("taper f=2 Hz, v=10%", 2.0, 0.10),
        ("taper f=4 Hz, v=20%", 4.0, 0.20),
    ]

    say("F-K wedge taper test -- does smoothing the mask remove the ringing?")
    say("mask: 5-20 Hz, 2500-4500 m/s, F*K<0; grid %d ch x %d samples, fs %.0f Hz, dx %.3f m"
        % (nch, nt, fs, dx))
    say("")
    say("--- 1. kernel compactness (no data involved) ---")
    say("%-22s %10s %10s %10s %12s" % ("config", "E<=2", "E<=8", "E<=16", "tail/peak"))
    masks, kers = {}, {}
    for name, ft, vt in configs:
        m = build_mask(nch, nt, fs, dx, ft, vt)
        ker, met = kernel_metrics(m)
        masks[name] = m; kers[name] = ker
        say("%-22s %10.4f %10.4f %10.4f %12.2e"
            % (name, met["energy_within_2"], met["energy_within_8"],
               met["energy_within_16"], met["tail_median_rel_peak"]))

    say("")
    say("--- 2. white-noise floor and 3. real-data score ---")
    say("    white noise: 8 Gaussian realizations through the complete operator")

    white = {name: [] for name, _, _ in configs}
    for _ in range(8):
        raw = rng.standard_normal((nch, nt))
        xp = preprocess(raw, fs)
        for name, _, _ in configs:
            y = apply_mask(xp, masks[name])
            lags, top, dist = correlate_top(y, fs, dx)
            white[name].append(abs(score_at(top, lags, dist, V_REF)))

    # Real data, built through an IDENTICAL path to the noise: raw 500 Hz /
    # 1.021 m excerpt -> preprocess at full rate -> 2x decimate -> mask, which is
    # what production does. Duration-matched to the noise realizations (8 x 12 s)
    # so the two floors are directly comparable.
    real = {name: [] for name, _, _ in configs}
    import pandas as pd, h5py
    CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/"
               "SAFOD_2024_2025.csv")
    db = pd.read_csv(CSV, sep=r"\s+"); db = db[db.nSamples > 0].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    day = db[db.t.dt.strftime("%Y-%m-%d") == "2024-12-20"].sort_values("t").reset_index(drop=True)
    picked = 0
    for row in day.iloc[::180].itertuples(index=False):
        if picked >= 8:
            break
        f = Path(str(row.file).replace("/data/SAFODAS1-harddrive-transfer",
                                       "/data/SAFOD/SAFODAS1-harddrive-transfer"))
        if not f.is_file():
            continue
        with h5py.File(f, "r") as h:
            g = h["Acquisition/Raw[0]"]
            d = g["RawData"][:12000, :800].astype(np.float64).T   # 24 s, 800 ch
            fs_raw = float(g.attrs.get("OutputDataRate", 500.0))
        d = np.diff(d, axis=1, prepend=d[:, :1])                  # phase -> strain rate
        xp_full = preprocess(d, fs_raw)
        xr = xp_full[::2, ::2][:nch, :nt]                         # match noise grid
        if xr.shape != (nch, nt):
            continue
        for name, _, _ in configs:
            y = apply_mask(xr, masks[name])
            lags, top, dist = correlate_top(y, fs, dx)
            real[name].append(abs(score_at(top, lags, dist, V_REF)))
        picked += 1
    say("    real data: %d duration-matched 12 s excerpts from 2024-12-20" % picked)
    say("")
    say("%-22s %13s %12s %13s %10s" % ("config", "white median", "white 95th",
                                        "real median", "real/white"))
    base = float(np.median(white["hard (production)"]))
    margins = {}
    for name, _, _ in configs:
        w = np.array(white[name])
        r = np.array(real[name]) if real.get(name) else np.array([np.nan])
        marg = np.median(r) / np.median(w) if np.median(w) > 0 else np.nan
        margins[name] = marg
        say("%-22s %13.4f %12.4f %13.4f %10.2f"
            % (name, np.median(w), np.percentile(w, 95), np.median(r), marg))
    say("")
    say("  'real/white' is the quantity that decides whether the fan adds information.")
    say("  A taper is only useful if it raises this ratio -- i.e. suppresses the")
    say("  manufactured floor by more than it suppresses the real score.")

    say("")
    say("--- interpretation ---")
    best = min(configs[1:], key=lambda c: np.median(white[c[0]]))
    red = np.median(white[best[0]]) / base
    say("  Tapering changes the manufactured white-noise score by %.2fx at best"
        " (%s)." % (red, best[0]))
    say("  real/white ratio: hard %.2f -> best taper %.2f"
        % (margins["hard (production)"], margins[best[0]]))
    if red < 0.5:
        say("  The taper substantially suppresses the ringing floor.")
    elif red < 0.9:
        say("  The taper reduces the ringing floor modestly.")
    else:
        say("  The taper does NOT meaningfully reduce the manufactured score.")
    say("")
    say("  Either way this does not license the F-K result. A velocity fan passes")
    say("  only energy consistent with its own band, so its output resembles that")
    say("  moveout whatever the edge shape. The gate is the pre-filter channel")
    say("  scramble on real data, which the production fan fails.")

    fig, ax = plt.subplots(2, len(configs), figsize=(4 * len(configs), 7),
                           constrained_layout=True)
    for j, (name, _, _) in enumerate(configs):
        ax[0, j].imshow(np.fft.fftshift(masks[name]), aspect="auto", cmap="magma",
                        origin="lower")
        ax[0, j].set_title(name, fontsize=9); ax[0, j].set_xticks([]); ax[0, j].set_yticks([])
        if j == 0: ax[0, j].set_ylabel("mask in (f,k)")
        c0, c1 = kers[name].shape[0] // 2, kers[name].shape[1] // 2
        cut = kers[name][c0 - 40:c0 + 41, c1 - 60:c1 + 61]
        ax[1, j].imshow(cut, aspect="auto", cmap="RdBu_r", origin="lower",
                        vmin=-np.abs(cut).max() * 0.15, vmax=np.abs(cut).max() * 0.15)
        ax[1, j].set_xticks([]); ax[1, j].set_yticks([])
        ax[1, j].set_xlabel("white |score| %.3f" % np.median(white[name]), fontsize=9)
        if j == 0: ax[1, j].set_ylabel("kernel (t,x), clipped")
    fig.suptitle("F-K wedge: mask shape (top) and its impulse response (bottom).\n"
                 "Tails in the kernel are the ringing that stamps moveout onto any input.",
                 fontsize=10)
    fig.savefig(str(STEM) + ".png", dpi=190)

    np.savez(str(STEM) + ".npz",
             configs=np.array([c[0] for c in configs]),
             white_scores=np.array([white[c[0]] for c in configs]),
             kernels=np.array([kers[c[0]] for c in configs]))
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
