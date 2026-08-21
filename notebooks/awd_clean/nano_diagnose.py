#!/usr/bin/env python3
"""What are the vertical stripes in the Nano correlation gather?

Vertical stripes in a correlation gather mean energy at a FIXED LAG on every
trace: the same waveform at the same time on all channels, with no moveout.
Evenly spaced stripes across the lag axis mean the correlation is periodic in
lag, which means the input is dominated by a NARROWBAND TONE.

Two things are checked here, because the earlier Nano section had a second
problem as well.

1  THE TONE. Per-channel spectra, in the air section and in the hole, looking for
   a narrowband peak inside the 5-20 Hz analysis band. A tone shared across
   channels with no moveout is the textbook cause of the stripe pattern, and if
   it is present in the air section too it is instrumental rather than seismic.

2  THE SOURCE WAS IN THE AIR. The section that produced the stripes used Nano
   source channel 10, and nano_find_wellhead.py subsequently placed the fibre's
   entry into the hole at channel 73. So it correlated an AIR channel against
   downhole channels, which returns whatever the two share -- instrumental
   common mode at fixed lag -- and cannot show moveout even if moveout exists.
   That figure is void; a valid source is channel 73 or below.

Output: nano_diagnose.{png,txt}
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import welch, detrend

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nano_ambient_cc as nano

STEM = HERE / "nano_diagnose"
WELLHEAD = 73
INK, C1, C2, C3 = "#444444", "#0072B2", "#D55E00", "#009E73"

log = []
def say(s):
    print(s, flush=True); log.append(s)

files = sorted(nano.NANO_DIR.glob("*.pb"))
pre = [(p, t) for p, t in ((p, nano.nano_time(p.name)) for p in files)
       if t is not None and t < nano.FIRST_DROP][:2]      # 10 min is plenty
arr, info = nano.readFile_protobuf([str(p) for p, _ in pre], fmin=1.0, fmax=100.0,
                                   desampling=True, **nano.RAW_KW)
x = np.asarray(arr, dtype=np.float32); del arr
fs = float(info["fs"])
from zoneinfo import ZoneInfo
loc = pre[0][1].astimezone(ZoneInfo("America/Los_Angeles"))
say("Nano diagnostics | %s (local) | %.1f min | %d channels at %.0f Hz"
    % (loc.strftime("%a %d %b %Y %H:%M"), len(pre) * 5, x.shape[0], fs))
say("  wellhead is channel %d; channels 0-%d are AIR" % (WELLHEAD, WELLHEAD - 1))
say("")

say("%8s %9s %10s   %s" % ("channel", "zone", "peak Hz", "top 3 peaks (Hz : x median)"))
spectra = {}
for ch in (10, 40, 73, 100, 200, 400, 600):
    if ch >= x.shape[0]:
        continue
    tr = detrend(np.diff(x[ch].astype(np.float64)))
    f, P = welch(tr, fs=fs, nperseg=4096)
    m = (f >= 3) & (f <= 60)
    ff, PP = f[m], P[m]
    med = float(np.median(PP))
    peaks = []
    for i in np.argsort(PP)[::-1]:
        if all(abs(ff[i] - p0) > 1.0 for p0, _ in peaks):
            peaks.append((float(ff[i]), float(PP[i] / med)))
        if len(peaks) == 3:
            break
    spectra[ch] = (ff, PP / med)
    say("%8d %9s %10.2f   %s" % (ch, "AIR" if ch < WELLHEAD else "in hole",
        peaks[0][0], "  ".join("%.1f : %.0fx" % p for p in peaks)))

say("")
strong = [(c, s[1].max()) for c, s in spectra.items()]
worst = max(strong, key=lambda t: t[1])
say("=== reading ===")
if worst[1] > 20:
    say("  A STRONG NARROWBAND TONE is present (peak %.0fx the median in the band)."
        % worst[1])
    say("  A tone shared across channels correlates at a FIXED LAG on every trace,")
    say("  which is exactly the vertical-stripe pattern, and it carries no moveout.")
    air = [v for c, v in strong if c < WELLHEAD]
    hole = [v for c, v in strong if c >= WELLHEAD]
    if air and hole and max(air) > 0.5 * max(hole):
        say("  It is present in the AIR section too, so it is instrumental or")
        say("  cable-borne rather than a seismic arrival.")
else:
    say("  No dominant narrowband tone; the stripes need another explanation.")
say("")
say("  SEPARATELY: the section that produced the stripes used source channel 10,")
say("  which is in the AIR. Correlating an air channel against downhole channels")
say("  returns only what they share and cannot show moveout. That figure is void")
say("  regardless of the tone; a valid Nano source is channel %d or below." % WELLHEAD)

fig, ax = plt.subplots(figsize=(8.6, 4.6), constrained_layout=True)
fig.suptitle("Nano fibre spectra | %s (local) | why the gather is vertical stripes"
             % loc.strftime("%a %d %b %Y %H:%M"), fontsize=11)
for ch, (ff, PP) in sorted(spectra.items()):
    col = C2 if ch < WELLHEAD else C1
    ax.semilogy(ff, PP, lw=1.1, color=col, alpha=.85,
                label="ch %d (%s)" % (ch, "AIR" if ch < WELLHEAD else "in hole"))
ax.axvspan(5, 20, color=C3, alpha=.10, label="5-20 Hz analysis band")
ax.set(xlabel="frequency (Hz)", ylabel="power / band median", xlim=(3, 60))
ax.legend(fontsize=7, ncol=2, frameon=False); ax.grid(alpha=.3, which="both")
fig.savefig(str(STEM) + ".png", dpi=300, bbox_inches="tight")
Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
say("wrote %s.{png,txt}" % STEM.name)
