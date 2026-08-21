"""Virtual-source gathers on the June 2026 AWD data, by correlation and by deconvolution.

Why this is not the ambient interferometry already in this repo
--------------------------------------------------------------
`stack_daily.py`, `sanity/sanity_cc.py` and
`ambient_transfer_test.normalized_corr_pairs` all cross-correlate a virtual
source against the array, but all three run on 2024-25 *ambient noise*, where
you correlate because you do not know when or where the source was.

Here the source is known. The point of correlating is not to find a Green's
function -- it is to **remove the source**. Correlating channel A with channel B
cancels the source term common to both and leaves the response between A and B,
as though a source sat at A. That redatums the surface weight drop down into the
borehole and takes the near-surface path with it (Bakulin & Calvert 2006, "the
virtual source method").

That matters here specifically. `docs/paper1/STATUS.md` already concludes the
source is the limiter -- sigma_alpha ~ 0.30 ms common-mode timing and a 39%
amplitude CV. Anything that cancels the source term attacks the measured
bottleneck directly.

Correlation and deconvolution, and why both
-------------------------------------------
Classic Bakulin-Calvert sums over a *range* of surface source positions to
satisfy stationary phase. The AWD did not move (`README.md`: fixed for the
survey, ~15 m from the Nano wellhead), so that construction is not available and
the correlation gather here should be read as a redatumed section, not as a
clean Green's function.

Deconvolution interferometry does not need a source aperture. Dividing B by A in
the frequency domain cancels the source spectrum whatever its shape, which is
the Snieder & Safak (2006) construction later applied to borehole arrays by
Nakata & Snieder. For *monitoring* that is the useful one: you do not need a
correct Green's function, you need a repeatable waveform whose changes track the
medium.

Both are computed so they can be compared. Deconvolution is the one to trust for
source cancellation; correlation is the one with the better SNR.

A side benefit worth stating: both operations cancel any gain applied equally to
the whole record, so the read-time Tukey taper documented in `PREPROCESSING.md`
cannot bias these gathers.

What this deliberately does NOT do
----------------------------------
No F-K filter, no velocity wedge, no semblance scan, no permutation null. This
is a plain gather -- lag against distance -- so the wavefield can be looked at
before any selection is applied to it. Several virtual-source channels are used
so the result cannot be an artifact of always sitting at channel 0.

Input
-----
`canonical_epoch_stacks_paired_deep_all.npz` (859 paired drops in 46 burst
stacks, built by `paired_stack_job_deep_all.py`).

Outputs
-------
figures/awd_2026/plain_look/vs_fig01_correlation_gathers.png
figures/awd_2026/plain_look/vs_fig02_deconvolution_gathers.png
figures/awd_2026/plain_look/vs_fig03_moveout_and_wiggles.png
figures/awd_2026/plain_look/awd_virtual_source.npz
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt

HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
OUT_DIR = Path(
    "/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026/plain_look"
) / "virtual_source"

PRE_S = 0.5                     # the stacks are cut -0.5 to +3.0 s
WATER_LEVEL = 0.01              # fraction of mean |A|^2, standard spectral-division floor

# Per-fiber configuration. Deep needs a much longer lag window: over its 2.8 km
# aperture at ~1545 m/s the differential travel time reaches ~1.8 s, so the
# 0.35 s window that suits Nano's 600 m would truncate the gather entirely.
#
# Deep is restricted to the outbound leg. The fiber reverses at channel 1700 (safod_geometry.py), and
# past that "distance along fiber" stops being monotonic in depth, which would
# make a redatumed moveout meaningless. 200-3000 m is channels 98-1469, safely
# inside the outbound limb, and is the range validated for the Deep guided mode.
# GEOMETRY GATE, added 2026-08-20. A virtual source only redatums the surface
# AWD source *downhole* if the source channel is actually downhole. Both fibre
# configurations previously violated that:
#
#   Nano  source  50 m -> channel  39, ABOVE the wellhead at channel 73
#                                       (nano_find_wellhead.txt, entry at 92 m)
#   Deep  source 400 m -> channel 196, SURFACE LEAD-IN (safod_geometry.py,
#                                       first in-hole channel is 211 = 431 m)
#
# Both apertures also started above fibre entry (Nano at 0 m, Deep at 200 m).
# Sources are now inside the ground on Nano, and inside the near-vertical
# outbound section on Deep (channels 211-949 = 431-1938 m), so that along-fibre
# distance is a depth. `check_geometry()` re-derives this at run time and
# refuses to proceed rather than trusting these literals.
NANO_WELLHEAD_CH = 73           # nano_find_wellhead.txt: amplitude AND coherence
DEEP_FIRST_IN_HOLE_CH = 211     # safod_geometry.py
DEEP_LAST_NEAR_VERTICAL_CH = 949

FIBERS = {
    "nano": dict(
        stack_key="nano_stacks", dx_key="dx_nano",
        band=(20.0, 50.0),          # Nano working band
        aperture_m=(100.0, 600.0),  # starts below the 92 m wellhead; signal dies
                                    # well before the 926 m fiber end
        max_lag_s=0.35,
        sources_m=[150.0, 250.0, 350.0, 450.0],
        v_ref=2975.0,
        label="Nano (cemented)",
    ),
    "deep": dict(
        stack_key="deep_stacks", dx_key="dx_deep",
        band=(15.0, 30.0),          # validated Deep guided-mode band
        aperture_m=(440.0, 3000.0), # starts below fibre entry at 431 m. Above
                                    # ~1938 m the hole deviates past 5 deg, so
                                    # along-fibre distance stops being depth.
        max_lag_s=2.00,
        sources_m=[500.0, 900.0, 1300.0, 1700.0],
        v_ref=1544.6,               # frozen outbound trajectory
        label="Deep (wireline, outbound leg)",
    ),
}


def check_geometry(fiber, dx, sources_m, aperture_m):
    """Refuse to build a virtual source on fibre that is not in the ground."""
    bad = []
    if fiber == "nano":
        lo = NANO_WELLHEAD_CH
        for m in list(sources_m) + [aperture_m[0]]:
            c = int(round(m / dx))
            if c < lo:
                bad.append("%.0f m -> ch %d is ABOVE the Nano wellhead (ch %d)"
                           % (m, c, lo))
    else:
        import safod_geometry as GEO
        for m in list(sources_m) + [aperture_m[0]]:
            c = int(round(m / dx))
            if c < DEEP_FIRST_IN_HOLE_CH:
                bad.append("%.0f m -> ch %d is SURFACE LEAD-IN: %s"
                           % (m, c, GEO.describe(c)))
        for m in sources_m:
            c = int(round(m / dx))
            if c > DEEP_LAST_NEAR_VERTICAL_CH:
                bad.append("%.0f m -> ch %d is DEVIATED, along-fibre distance is "
                           "not depth: %s" % (m, c, GEO.describe(c)))
    if bad:
        raise SystemExit("geometry check failed for %s:\n  %s"
                         % (fiber, "\n  ".join(bad)))
    print("geometry OK (%s): sources %s m, aperture from %.0f m"
          % (fiber, [int(x) for x in sources_m], aperture_m[0]))

# Set by main() from the chosen fiber.
BAND = FIBERS["nano"]["band"]
APERTURE_M = FIBERS["nano"]["aperture_m"]
MAX_LAG_S = FIBERS["nano"]["max_lag_s"]
SOURCE_DEPTHS_M = FIBERS["nano"]["sources_m"]


def bandpass(x, fs, band):
    sos = butter(4, list(band), btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, np.asarray(x, float), axis=-1)


def weighted_stack(stacks, counts):
    """Drop-count-weighted mean over bursts, as every other script here does."""
    good = counts > 0
    w = counts[good].astype(float)
    return np.tensordot(w, stacks[good], axes=(0, 0)) / w.sum(), int(w.sum())


def _spectra(section, source_trace, nfft):
    B = np.fft.rfft(section, n=nfft, axis=-1)
    A = np.fft.rfft(source_trace, n=nfft)
    return A, B


def correlate_gather(section, source_trace, fs, max_lag=MAX_LAG_S):
    """B correlated with A: cancels the source term, keeps its spectrum."""
    n = section.shape[-1]
    nfft = 1 << int(np.ceil(np.log2(2 * n - 1)))
    A, B = _spectra(section, source_trace, nfft)
    cc = np.fft.irfft(B * np.conj(A)[None, :], n=nfft, axis=-1)
    return _center(cc, fs, max_lag)


def deconvolve_gather(section, source_trace, fs, max_lag=MAX_LAG_S,
                      water=WATER_LEVEL):
    """B deconvolved by A: cancels the source term *and* its spectrum.

    Water-level regularisation, the standard guard against dividing by the
    near-zeros of |A|. eps is a fraction of the mean power of A, so the floor
    scales with the trace rather than being an absolute number."""
    n = section.shape[-1]
    nfft = 1 << int(np.ceil(np.log2(2 * n - 1)))
    A, B = _spectra(section, source_trace, nfft)
    power = np.abs(A) ** 2
    eps = water * float(np.mean(power))
    dec = np.fft.irfft(B * np.conj(A)[None, :] / (power + eps)[None, :],
                       n=nfft, axis=-1)
    return _center(dec, fs, max_lag)


def _center(x, fs, max_lag):
    """Roll zero lag to the middle and trim to +/- max_lag."""
    ml = int(round(max_lag * fs))
    out = np.concatenate((x[:, -ml:], x[:, :ml + 1]), axis=-1)
    lags = np.arange(-ml, ml + 1) / fs
    return lags, out


def norm_rows(g):
    p = np.max(np.abs(g), axis=1, keepdims=True)
    return g / np.where(p > 0, p, 1.0)


def draw(ax, lags, gather, z, title, v_ref=None, z_src=None):
    v = np.percentile(np.abs(gather), 99.0)
    im = ax.pcolormesh(lags, z, gather, cmap="seismic", vmin=-v, vmax=v,
                       shading="auto")
    if v_ref is not None and z_src is not None:
        # the moveout a wave leaving the virtual source would follow
        for sign in (+1, -1):
            ax.plot(sign * np.abs(z - z_src) / v_ref, z, "k--", lw=0.9, alpha=0.6)
    if z_src is not None:
        ax.axhline(z_src, color="lime", lw=1.0, alpha=0.8)
    ax.set_xlabel("lag (s)")
    ax.set_ylabel("distance along fiber (m)")
    ax.set_title(title, fontsize=10)
    ax.invert_yaxis()
    return im


def main():
    global BAND, APERTURE_M, MAX_LAG_S, SOURCE_DEPTHS_M
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fiber", choices=tuple(FIBERS), default="nano")
    fiber = parser.parse_args().fiber
    cfg = FIBERS[fiber]
    BAND = cfg["band"]
    APERTURE_M = cfg["aperture_m"]
    MAX_LAG_S = cfg["max_lag_s"]
    SOURCE_DEPTHS_M = cfg["sources_m"]
    suffix = "" if fiber == "nano" else f"_{fiber}"

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    d = np.load(STACKS)
    fs = float(d["fs"])
    dx = float(d[cfg["dx_key"]])
    # Refuse to build a virtual source on fibre that is not in the ground.
    check_geometry(fiber, dx, SOURCE_DEPTHS_M, APERTURE_M)
    counts = d["n_common"]
    nano, n_drops = weighted_stack(d[cfg["stack_key"]], counts)
    print(f"{cfg['label']} stack {nano.shape}, fs={fs}, dx={dx:.6f}, "
          f"{n_drops} drops over {int((counts > 0).sum())} bursts", flush=True)
    print(f"band {BAND[0]:g}-{BAND[1]:g} Hz, max lag {MAX_LAG_S:g} s", flush=True)

    c0, c1 = int(APERTURE_M[0] / dx), min(int(APERTURE_M[1] / dx), nano.shape[0])
    sec = bandpass(nano[c0:c1], fs, BAND)
    z = (np.arange(c0, c1)) * dx
    src_ch = [min(int(round(zs / dx)) - c0, sec.shape[0] - 1)
              for zs in SOURCE_DEPTHS_M]
    print(f"aperture {z[0]:.0f}-{z[-1]:.0f} m ({sec.shape[0]} channels); "
          f"virtual sources at {[f'{z[c]:.0f} m' for c in src_ch]}", flush=True)

    corr, deco = {}, {}
    for c in src_ch:
        lags, corr[c] = correlate_gather(sec, sec[c], fs, MAX_LAG_S)
        _, deco[c] = deconvolve_gather(sec, sec[c], fs, MAX_LAG_S)

    # v_ref only draws a reference moveout; it is not fitted here.
    v_ref = cfg["v_ref"]

    for name, data, fname, note in [
        ("cross-correlation", corr, f"vs_fig01_correlation_gathers{suffix}.png",
         "source term cancelled, source spectrum retained"),
        ("deconvolution", deco, f"vs_fig02_deconvolution_gathers{suffix}.png",
         "source term and spectrum both cancelled"),
    ]:
        fig, axes = plt.subplots(1, len(src_ch), figsize=(4.2 * len(src_ch), 6.5),
                                 sharey=True)
        for a, c in zip(np.atleast_1d(axes), src_ch):
            im = draw(a, lags, norm_rows(data[c]), z,
                      f"virtual source at {z[c]:.0f} m",
                      v_ref=v_ref, z_src=z[c])
            plt.colorbar(im, ax=a, label="amplitude / trace peak")
        fig.suptitle(
            f"{cfg[chr(39)+chr(39)] if False else cfg['label']} AWD virtual-source gathers by {name} -- {note}\n"
            f"{n_drops} drops, {BAND[0]:g}-{BAND[1]:g} Hz, traces normalised; "
            f"dashed = {v_ref:g} m/s reference moveout, green = source channel",
            fontsize=12)
        fig.tight_layout()
        fig.savefig(OUT_DIR / fname, dpi=140)
        plt.close(fig)

    # ---- a plain side-by-side, plus wiggles, at one source position -------
    c = src_ch[1]
    fig, ax = plt.subplots(1, 3, figsize=(17, 6.5))
    im = draw(ax[0], lags, norm_rows(corr[c]), z, "correlation", v_ref, z[c])
    plt.colorbar(im, ax=ax[0], label="amplitude / trace peak")
    im = draw(ax[1], lags, norm_rows(deco[c]), z, "deconvolution", v_ref, z[c])
    plt.colorbar(im, ax=ax[1], label="amplitude / trace peak")

    step = max(1, sec.shape[0] // 40)
    for k in range(0, sec.shape[0], step):
        tr = deco[c][k]
        pk = np.max(np.abs(tr))
        if pk > 0:
            ax[2].plot(lags, 2.2 * step * dx * tr / pk + z[k], "k", lw=0.5)
    ax[2].axhline(z[c], color="lime", lw=1.0)
    for sign in (+1, -1):
        ax[2].plot(sign * np.abs(z - z[c]) / v_ref, z, "r--", lw=1.0, alpha=0.7)
    ax[2].invert_yaxis()
    ax[2].set_xlim(lags[0], lags[-1])
    ax[2].set_xlabel("lag (s)")
    ax[2].set_ylabel("distance along fiber (m)")
    ax[2].set_title("deconvolution, wiggles", fontsize=10)

    fig.suptitle(f"AWD virtual source at {z[c]:.0f} m: correlation vs deconvolution "
                 f"-- {n_drops} drops, {BAND[0]:g}-{BAND[1]:g} Hz", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"vs_fig03_moveout_and_wiggles{suffix}.png", dpi=140)
    plt.close(fig)

    np.savez(OUT_DIR / f"awd_virtual_source{suffix}.npz",
             lags=lags, z=z, fs=fs, dx=dx, band=np.asarray(BAND),
             n_drops=n_drops, water_level=WATER_LEVEL,
             source_channels=np.asarray(src_ch),
             source_depths=np.asarray([z[c] for c in src_ch]),
             correlation=np.stack([corr[c] for c in src_ch]).astype(np.float32),
             deconvolution=np.stack([deco[c] for c in src_ch]).astype(np.float32))

    # plain numbers, so the figures are not the only record.
    #
    # A per-channel max|amplitude| pick is not a usable moveout estimator on
    # these gathers: the largest sample on a trace is as often noise or the
    # near-zero-lag peak as it is the arrival, so the fitted slope is garbage.
    # Slant-stack instead -- sum along each trial moveout and take the trial
    # that maximises coherent energy. Standard, and it uses the whole gather
    # rather than one sample per trace.
    print("\n--- slant-stack apparent speed of the redatumed arrival ---")
    v_grid = np.arange(1000.0, 6000.0, 25.0)
    speeds = {}
    for c in src_ch:
        g = deco[c]
        offset = z - z[c]
        use = offset > 60.0          # causal side, clear of the source itself
        if use.sum() < 20:
            continue
        gg, off = g[use], offset[use]
        best = (np.nan, -1.0)
        curve = []
        for v in v_grid:
            shift = np.round(off / v * fs).astype(int)
            idx = shift + (len(lags) // 2)
            ok = (idx >= 0) & (idx < gg.shape[1])
            if ok.sum() < 20:
                curve.append(0.0)
                continue
            amp = gg[np.arange(gg.shape[0])[ok], idx[ok]]
            # semblance: coherent energy over total energy along the trajectory
            s = float(amp.sum() ** 2 / (ok.sum() * np.sum(amp ** 2))) \
                if np.any(amp) else 0.0
            curve.append(s)
            if s > best[1]:
                best = (v, s)
        speeds[c] = (best[0], best[1], np.asarray(curve))
        print(f"  source {z[c]:5.0f} m: {best[0]:6.0f} m/s "
              f"(semblance {best[1]:.3f}, {int(use.sum())} channels)")
    if speeds:
        vals = np.array([speeds[c][0] for c in speeds])
        print(f"  across {len(vals)} virtual sources: median {np.median(vals):.0f} m/s, "
              f"spread {vals.min():.0f}-{vals.max():.0f}")
        np.savez(OUT_DIR / f"awd_virtual_source_speeds{suffix}.npz",
                 v_grid=v_grid,
                 source_depths=np.asarray([z[c] for c in speeds]),
                 best_v=vals,
                 best_semblance=np.asarray([speeds[c][1] for c in speeds]),
                 curves=np.stack([speeds[c][2] for c in speeds]))
    print("\nwrote", OUT_DIR)


if __name__ == "__main__":
    main()
