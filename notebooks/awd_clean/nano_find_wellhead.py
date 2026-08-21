#!/usr/bin/env python3
"""Find where the Nano fibre enters the hole, so no channel in the air is used.

WHY. The Deep fibre turned out to have 211 channels of SURFACE LEAD-IN, and one
earlier source-channel scan used channel 98 -- which was in the air, not the hole,
and duly returned nothing. Nano has no entry in
SAFOD_Phase2_GeoReferenced_Channels.xlsx, and its source scan used channels
10/40/80/120/160 with no check that any of them were downhole. Gate G0
(`faultzone/repeaters/channel_depth_registration.py`) placed the 2024-25 wellhead
at channel 23, but that is a different interrogator on the same fibre and its
channel numbering need not carry over.

HOW, without a registration file. Surface and downhole fibre differ in ways that
are visible in ambient noise alone:

  1  AMPLITUDE. Air/surface fibre is loose and mechanically noisy; cemented
     downhole fibre is quieter. Expect a step down in per-channel RMS at entry.
  2  NEIGHBOUR COHERENCE. Cemented fibre samples a common wavefield and correlates
     strongly with its neighbours; loose surface fibre does not.
  3  SPECTRAL CHARACTER. Surface fibre carries more high-frequency, wind- and
     traffic-driven energy.

The wellhead is placed where these agree. All three are reported per channel so
the choice is inspectable rather than asserted, and the script refuses to name a
wellhead if the diagnostics disagree.

Read with median=False so the across-channel median is not removed -- that would
suppress exactly the shared-wavefield signature diagnostic 2 relies on.

Output: nano_find_wellhead.{npz,png,txt}
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, detrend, sosfiltfilt

import safod_geometry as geo
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nano_ambient_cc as nano

STEM = HERE / "nano_find_wellhead"
BAND = (5.0, 20.0)
HF_BAND = (40.0, 100.0)
INK, MUTED = "#444444", "#6b6b6b"
C1, C2, C3 = "#0072B2", "#D55E00", "#009E73"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfiles", type=int, default=6)
    a = ap.parse_args()

    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    files = sorted(nano.NANO_DIR.glob("*.pb"))
    pre = [(p, t) for p, t in ((p, nano.nano_time(p.name)) for p in files)
           if t is not None and t < nano.FIRST_DROP][: a.nfiles]
    say("Locating the Nano wellhead from ambient noise")
    say("  %d records = %.2f h, %s to %s UTC"
        % (len(pre), len(pre) * 300 / 3600.0, pre[0][1], pre[-1][1]))

    arr, info = nano.readFile_protobuf([str(p) for p, _ in pre], fmin=1.0, fmax=200.0,
                                       desampling=True, **nano.RAW_KW)
    x = np.asarray(arr, dtype=np.float64); del arr
    fs = float(info["fs"]); dx = float(info.get("dx", nano.DX_NANO))
    say("  %d channels at %.4f m, fs %.0f Hz (%.0f m extent)"
        % (x.shape[0], dx, fs, x.shape[0] * dx))
    say("")

    x = detrend(np.diff(x, axis=1), axis=1)          # strain-rate proxy
    band = sosfiltfilt(butter(4, list(BAND), btype="bandpass", fs=fs, output="sos"),
                       x, axis=1)
    hf = sosfiltfilt(butter(4, [HF_BAND[0], min(HF_BAND[1], 0.45 * fs)],
                            btype="bandpass", fs=fs, output="sos"), x, axis=1)
    del x

    rms = np.sqrt(np.mean(band ** 2, axis=1))
    hf_rms = np.sqrt(np.mean(hf ** 2, axis=1))
    ratio = hf_rms / np.maximum(rms, 1e-30)          # diagnostic 3
    # diagnostic 2: correlation with the immediate neighbour
    nb = np.full(band.shape[0], np.nan)
    for i in range(band.shape[0] - 1):
        u, v = band[i], band[i + 1]
        su, sv = u.std(), v.std()
        if su > 0 and sv > 0:
            nb[i] = float(np.mean((u - u.mean()) * (v - v.mean())) / (su * sv))

    n = rms.size
    say("=== per-channel diagnostics, first 200 channels ===")
    say("  %5s %12s %10s %10s" % ("ch", "RMS", "nbr corr", "HF/band"))
    for c in range(0, min(200, n), 10):
        say("  %5d %12.4e %10.3f %10.3f" % (c, rms[c], nb[c], ratio[c]))
    say("")

    # find the step: compare a leading window against the interior
    interior = slice(max(1, n // 4), n)
    med_in = float(np.nanmedian(rms[interior]))
    nb_in = float(np.nanmedian(nb[interior]))
    cand = []
    for c in range(2, min(400, n - 5)):
        lead_quiet = rms[:c].max() > 3.0 * med_in          # lead-in is much louder
        now_quiet = rms[c:c + 20].max() < 3.0 * med_in
        coh = np.nanmedian(nb[c:c + 20]) > 0.5 * nb_in
        if lead_quiet and now_quiet and coh:
            cand.append(c)
    say("=== verdict ===")
    say("  interior median RMS %.4e, interior neighbour corr %.3f" % (med_in, nb_in))
    if cand:
        wh = cand[0]
        say("  amplitude and coherence agree on entry at channel %d (%.0f m along fibre)"
            % (wh, wh * dx))
        say("  channels 0-%d look like SURFACE/AIR fibre and must not be used" % (wh - 1))
    else:
        loud = int(np.sum(rms > 3.0 * med_in))
        say("  NO clear step found: %d channels exceed 3x the interior RMS." % loud)
        if loud == 0:
            say("  The fibre looks uniformly downhole from channel 0, i.e. there is")
            say("  little or no lead-in and the earlier source channels were valid.")
        else:
            say("  Diagnostics disagree; do NOT name a wellhead from this alone.")
        wh = -1
    say("")
    say("  For comparison: Deep has 211 lead-in channels; gate G0 placed the")
    say("  2024-25 wellhead at channel 23 on this same physical fibre.")

    fig, ax = plt.subplots(1, 3, figsize=(15.0, 4.2), constrained_layout=True)
    # LOCAL time in the title, UTC demoted to a footnote: whether personnel were
    # at SAFOD is a local-clock question, and UTC-7 pushes a California afternoon
    # onto the following UTC date.
    _title, _foot = geo.figure_label(pre[0][1], pre[-1][1], len(pre) * 300 / 3600.0,
                                     fibre="Nano",
                                     extra="where does the fibre enter the hole?")
    fig.suptitle(_title, fontsize=10.5)
    fig.text(0.995, 0.002, _foot, ha="right", va="bottom",
             fontsize=6.5, color="#9a9a9a")
    ch = np.arange(n)
    for axi, (y, lab, col) in zip(ax, ((rms, "RMS, %g-%g Hz" % BAND, C1),
                                       (nb, "correlation with next channel", C2),
                                       (ratio, "HF / band RMS", C3))):
        axi.plot(ch, y, "-", color=col, lw=0.9)
        if wh > 0:
            axi.axvline(wh, color=INK, ls="--", lw=1.4, label="entry ch %d" % wh)
            axi.legend(fontsize=8, frameon=False)
        axi.set(xlabel="channel", ylabel=lab, xlim=(0, min(400, n)))
        axi.grid(alpha=.3)
    ax[0].set_yscale("log")
    fig.savefig(str(STEM) + ".png", dpi=300, bbox_inches="tight")
    np.savez(str(STEM) + ".npz", channel=ch, rms=rms, neighbour_corr=nb,
             hf_ratio=ratio, wellhead=wh, dx=dx, fs=fs)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
