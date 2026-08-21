#!/usr/bin/env python3
"""Build Ambient_CC_Steps.ipynb -- the step-by-step ambient-noise CC notebook.

Run with the SYSTEM python3, which has nbformat; the `das` env does not. Execute
the result with the `das` kernel. This mirrors the convention already used by
`build_ambient_fk_qc_notebook.py` in this directory.

    python3 build_ambient_cc_steps_notebook.py
    sbatch exec_ambient_cc_steps.sh        # nbconvert --execute with the das kernel

DESIGN OF THE NOTEBOOK. One panel per processing step, starting from raw data, and
each step's NECESSITY shown by ABLATION: the same pipeline with exactly that one
step deleted, on identical input. Assertion is not evidence; a panel that shows
what breaks is.

FIGURE CONVENTIONS, chosen deliberately rather than by default:
  - categorical hues #0072B2 / #D55E00 / #009E73, fixed order, never cycled.
    Validated colourblind-safe: worst adjacent pair dE 11.0 deuteranopia, 25.8
    normal vision, all >= 3:1 contrast on the surface.
  - #444444 is INK (reference curves, the full pipeline) and #8a8a8a is the
    null/threshold neutral. Neither is a categorical slot.
  - signed correlation amplitude uses a DIVERGING map (RdBu_r) with a neutral
    midpoint at zero, never a rainbow.
  - six ablations are shown as SMALL MULTIPLES, two colours per panel, rather
    than six lines on one axes -- six hues cannot be made mutually distinguishable
    under colour-vision deficiency.
  - single axis everywhere; no dual y-axes.
  - recessive grid and spines; values in ink, identity carried by the mark.
"""
from __future__ import annotations

import sys
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent
OUT = HERE / "Ambient_CC_Steps.ipynb"
NOTEBOOK_VERSION = "v1"

PRELUDE = '''
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path.cwd()
STEPS = np.load(HERE / "deep_cc_steps.npz", allow_pickle=True)

# --- figure conventions (see build_ambient_cc_steps_notebook.py docstring) ---
INK, MUTED, NULLC = "#444444", "#6b6b6b", "#8a8a8a"
C1, C2, C3 = "#0072B2", "#D55E00", "#009E73"     # validated categorical, fixed order
mpl.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 300, "figure.facecolor": "white",
    "font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#b0b0b0", "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6,
    "legend.frameon": False, "legend.fontsize": 8,
    "lines.linewidth": 1.6, "lines.solid_capstyle": "round",
})

LAGS   = STEPS["lags"]
OFFS   = STEPS["offsets"]
VGRID  = STEPS["v_grid"]
FS     = float(STEPS["fs"])
DX     = float(STEPS["dx"])
NWIN   = int(STEPS["n_windows"])
ARRIVAL_V = float(STEPS["arrival_v"])
ARRIVAL_P = float(STEPS["arrival_p"])

def waterfall(ax, data, fs, dx, title, vlim_pct=99.0, t0=0.0):
    """Signed amplitude: diverging map, neutral at zero, symmetric limits."""
    lim = float(np.percentile(np.abs(data), vlim_pct)) or 1.0
    im = ax.imshow(data, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                   interpolation="nearest",
                   extent=[t0, t0 + data.shape[1]/fs, data.shape[0]*dx, 0])
    ax.set(xlabel="time (s)", ylabel="distance along fibre (m)", title=title)
    ax.grid(False)
    return im

def gather_plot(ax, gather, lags, offs, title, colour=INK, ref_v=None, scale=40.0):
    g = gather / np.maximum(np.abs(gather).max(axis=1, keepdims=True), 1e-30)
    for o, row in zip(offs, g):
        ax.plot(lags, -o + row*scale, "-", color=colour, lw=0.7)
    if ref_v:
        ax.plot(offs/ref_v, -offs, "--", color=C3, lw=1.4,
                label="%.0f m/s" % ref_v)
        ax.legend(loc="lower right")
    # Derive the limit from the curve being drawn. A hard +-0.35 s put the
    # arrival off the page: 700 m at 1350 m/s lands at 0.52 s.
    _lim = min(float(abs(lags).max()),
               1.3 * float(abs(offs).max()) / ref_v) if ref_v else float(abs(lags).max())
    ax.set(xlim=(-_lim, _lim), xlabel="correlation lag (s)",
           ylabel="offset from virtual source (m)", title=title)
    ax.grid(False)
'''

CELLS: list[tuple[str, str]] = []


def md(text):
    CELLS.append(("md", text))


def code(text):
    CELLS.append(("code", text))


# ----------------------------------------------------------------- title
md(f"""# Ambient-noise cross-correlation, step by step

**SAFOD borehole DAS, June 2026 Deep fibre — {NOTEBOOK_VERSION}**

This notebook walks the ambient-noise cross-correlation pipeline from raw optical
phase to a recovered arrival, one panel per processing step, and shows that each
step is **necessary** by deleting it and displaying what breaks.

Necessity is demonstrated by **ablation**, not asserted. Every comparison runs the
identical pipeline on identical input with exactly one step removed.

### The pipeline, as reported by Lellouch et al. (2019) §4.1

| step | operation | why |
|---|---|---|
| 1 | raw optical phase, `rad · 2π/2¹⁶` | the recorded quantity |
| 2 | time differentiation | phase → strain-rate proxy; removes drift |
| 3 | running-absolute-mean, 0.1 s | stops loud transients dominating the stack |
| 4 | 30 s windows, 15 s overlap | stationarity |
| 5 | cross-correlate vs the virtual source | the interferometry itself |
| 6 | R±10 neighbour sum | signal-to-noise |
| 7 | stack over windows | noise suppression as 1/√N |
| 8 | 5–20 Hz bandpass **then** crop in lag | order matters — see §8 |
| 9 | envelope moveout scan + per-velocity null | the measurement |

### Provenance

Arrays are read from `deep_cc_steps.npz`, produced by `deep_cc_steps.py`, which
verifies its compact correlator against the authoritative engine
(`ambient_lellouch2019_exact_stack.py`) and refuses to write products if they
disagree. Nothing in this notebook is typed in by hand.
""")

code(PRELUDE + '''
print("data window : %s -> %s UTC" % (STEPS["t_first"], STEPS["t_last"]))
print("records     : %d x 60 s = %.2f h" % (int(STEPS["n_files"]), int(STEPS["n_files"])*60/3600))
print("sample rate : %.0f Hz, channel spacing %.4f m" % (FS, DX))
print("windows     : %d (30 s, 15 s overlap)" % NWIN)
print("first weight drop is %.2f h AFTER the last record used" % float(STEPS["drop_gap_hours"]))
print("            -> this is ambient noise, not the active source")
print("recovered   : %.0f m/s at per-velocity p = %.4f" % (ARRIVAL_V, ARRIVAL_P))
''')

# ----------------------------------------------------------------- step 1
md("""## 1 — Raw data

The recorded quantity is **optical phase**, not ground motion. Two features
dominate and both matter for what follows: a large slowly-varying drift, and a
component that is nearly identical on every channel.

Neither is seismic. This is why step 2 exists.""")

code('''
raw = STEPS["raw_snippet"]
fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.0), constrained_layout=True,
                       gridspec_kw={"width_ratios": [1.35, 1]})
im = waterfall(ax[0], raw, FS, DX, "(a) Raw optical phase, 60 s")
fig.colorbar(im, ax=ax[0], label="phase (rad·2π/2¹⁶)", pad=0.01)

t = np.arange(raw.shape[1]) / FS
for i, ch in enumerate([0, raw.shape[0]//3, 2*raw.shape[0]//3]):
    ax[1].plot(t, raw[ch], color=[C1, C2, C3][i], lw=1.0,
               label="channel index %d" % ch)
ax[1].set(xlabel="time (s)", ylabel="phase (rad·2π/2¹⁶)",
          title="(b) Three raw traces: drift dominates")
ax[1].legend()
plt.show()

print("raw dynamic range      : %.3e" % np.ptp(raw))
cm = raw.mean(axis=0)
print("across-channel mean is %.1f%% of a typical trace's amplitude"
      % (100*np.std(cm)/np.std(raw[raw.shape[0]//2])))
''')

# ----------------------------------------------------------------- step 2
md("""## 2 — Time differentiation → strain rate

**Necessity.** Differentiating suppresses the drift and converts phase to a
strain-rate proxy. Panel (c) shows the consequence of *skipping* it: with raw
phase, the low-frequency drift is orders of magnitude larger than the seismic
band, and after correlation the moveout scan carries no interpretable peak.""")

code('''
rate = STEPS["rate_snippet"]
fig, ax = plt.subplots(1, 3, figsize=(15.0, 3.9), constrained_layout=True)
waterfall(ax[0], rate, FS, DX, "(a) After differentiation")

# amplitude spectra, before and after, one mid-fibre channel
ch = rate.shape[0] // 2
for arr, lab, col in ((STEPS["raw_snippet"][ch], "raw phase", C1),
                      (rate[ch], "differentiated", C2)):
    f = np.fft.rfftfreq(arr.size, 1/FS)
    A = np.abs(np.fft.rfft(arr * np.hanning(arr.size)))
    m = (f > 0.05) & (f <= 120)
    ax[1].loglog(f[m], A[m]/A[m].max(), color=col, label=lab)
ax[1].axvspan(5, 20, color=C3, alpha=0.12)
ax[1].set(xlabel="frequency (Hz)", ylabel="normalised amplitude",
          title="(b) Drift removed; 5–20 Hz band shaded")
ax[1].legend()

ax[2].plot(VGRID/1e3, STEPS["curve_full"], color=INK, label="full pipeline")
ax[2].plot(VGRID/1e3, STEPS["curve_no_diff"], color=C2, label="step 2 deleted")
ax[2].axvline(ARRIVAL_V/1e3, color=C3, ls=":", lw=1.4)
ax[2].set(xlabel="trial velocity (km/s)", ylabel="moveout score",
          title="(c) Ablation: no differentiation")
ax[2].legend()
plt.show()

for k, lab in (("curve_full","full"), ("curve_no_diff","no differentiation")):
    c = STEPS[k]; print("%-22s peak %8.4f at %6.0f m/s" % (lab, c.max(), VGRID[c.argmax()]))
''')

# ----------------------------------------------------------------- step 3
md("""## 3 — Running-absolute-mean temporal normalisation

A 0.1 s running-absolute-mean divides each trace by its own local amplitude
envelope (Bensen et al. 2007, ≈ half the maximum period of a 5–20 Hz band).

**Necessity.** Without it the stack is dominated by whichever windows happen to
contain the loudest transients, so most of the record contributes nothing.""")

code('''
normed = STEPS["normed_snippet"]
ch = rate.shape[0] // 2
t = np.arange(rate.shape[1]) / FS

fig, ax = plt.subplots(1, 3, figsize=(15.0, 3.9), constrained_layout=True)
ax[0].plot(t, rate[ch]/np.abs(rate[ch]).max(), color=C1, lw=0.8,
           label="before (strain rate)")
ax[0].set(xlabel="time (s)", ylabel="normalised amplitude",
          title="(a) Before: a few transients dominate")
ax[0].legend()
ax[1].plot(t, normed[ch]/np.abs(normed[ch]).max(), color=C2, lw=0.8,
           label="after (RAM 0.1 s)")
ax[1].set(xlabel="time (s)", title="(b) After: every window contributes")
ax[1].legend()

ax[2].plot(VGRID/1e3, STEPS["curve_full"], color=INK, label="full pipeline")
ax[2].plot(VGRID/1e3, STEPS["curve_no_ram"], color=C2, label="step 3 deleted")
ax[2].axvline(ARRIVAL_V/1e3, color=C3, ls=":", lw=1.4)
ax[2].set(xlabel="trial velocity (km/s)", ylabel="moveout score",
          title="(c) Ablation: no RAM normalisation")
ax[2].legend()
plt.show()

a = np.abs(rate[ch]); b = np.abs(normed[ch])
print("amplitude concentration, 99th percentile / median:")
print("  before RAM %8.1f" % (np.percentile(a,99)/np.median(a)))
print("  after  RAM %8.1f" % (np.percentile(b,99)/np.median(b)))
''')

# ----------------------------------------------------------------- steps 4-7
md("""## 4–7 — Windows, correlation, the R±10 sum, and stacking

Steps 4 and 5 are the interferometry: correlate each 30 s window against the
fixed virtual source. Step 6 sums each receiver with its ±10 neighbouring
channels, as published. Step 7 stacks the windows.

**Necessity of stacking (step 7).** One window is noise. Panel (a) is a single
window; panel (b) is the full stack. **Necessity of the neighbour sum (step 6)**
is panel (c).""")

code('''
fig, ax = plt.subplots(1, 3, figsize=(15.0, 4.6), constrained_layout=True)
gather_plot(ax[0], STEPS["gather_one_window"], LAGS, OFFS,
            "(a) One 30 s window", colour=NULLC)
gather_plot(ax[1], STEPS["gather_full"], LAGS, OFFS,
            "(b) %d windows stacked" % NWIN, colour=INK, ref_v=ARRIVAL_V)

ax[2].plot(VGRID/1e3, STEPS["curve_full"], color=INK, label="R±10 (published)")
ax[2].plot(VGRID/1e3, STEPS["curve_no_neighbour"], color=C2,
           label="single channel (step 6 deleted)")
ax[2].axvline(ARRIVAL_V/1e3, color=C3, ls=":", lw=1.4)
ax[2].set(xlabel="trial velocity (km/s)", ylabel="moveout score",
          title="(c) Ablation: no neighbour sum")
ax[2].legend()
plt.show()

for k, lab in (("curve_one_window","one window"), ("curve_full","full stack"),
               ("curve_no_neighbour","no neighbour sum")):
    c = STEPS[k]; print("%-22s peak %8.4f at %6.0f m/s" % (lab, c.max(), VGRID[c.argmax()]))
''')

md("""### How much stacking is actually required?

The score at the arrival velocity against the number of windows. This is the
quantity that distinguishes a present-but-faint arrival (grows) from an absent one
(flat) — on the 2024–25 main-hole fibre the equivalent curve is flat, which is how
that dataset was shown to contain no arrival.""")

code('''
n = STEPS["stack_n"]; at_arr = STEPS["stack_at_arrival"]; pk = STEPS["stack_peak"]
fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
ax.semilogx(n, at_arr, "o-", color=C1, ms=6, label="score at %.0f m/s" % ARRIVAL_V)
ax.semilogx(n, pk, "s--", color=C2, ms=5, label="peak of the scan")
ax.set(xlabel="windows stacked", ylabel="moveout score",
       title="Stacking behaviour at the arrival velocity")
ax.legend()
plt.show()
for a_, b_, c_ in zip(n, at_arr, pk):
    print("%5d windows: score at arrival %7.4f | scan peak %7.4f" % (a_, b_, c_))
''')

# ----------------------------------------------------------------- step 8
md("""## 8 — Bandpass **then** crop, not the other way round

The 5–20 Hz bandpass is applied to the **full** correlation, and only then is it
cropped in lag.

**How wide the crop must be.** It was `±0.35 s`, and that was wrong: a 700 m
offset at 1350 m/s arrives at 0.52 s, so the far half of the array was cropped
away. Worse, the moveout scan then took a median over whichever offsets still
fitted — the near ones, which sit closest to the zero-lag lobe and read high — so
the score rose as velocity fell for a purely geometric reason. That artefact is
what put the arrival at 1675 m/s. The window is now derived from the aperture and
the slowest velocity scanned (`arrival_velocities.required_lag_s`).

**Necessity of the order.** Cropping first is a hard truncation, which aliases
out-of-band energy into the retained window. Reversing the two operations changes
the answer, so the order is part of the method rather than an implementation
detail.""")

code('''
fig, ax = plt.subplots(1, 3, figsize=(15.0, 4.6), constrained_layout=True)
gather_plot(ax[0], STEPS["gather_full"], LAGS, OFFS,
            "(a) Band, then crop (correct)", colour=INK, ref_v=ARRIVAL_V)
gather_plot(ax[1], STEPS["gather_crop_then_band"], LAGS, OFFS,
            "(b) Crop, then band (wrong)", colour=C2)
ax[2].plot(VGRID/1e3, STEPS["curve_full"], color=INK, label="band → crop")
ax[2].plot(VGRID/1e3, STEPS["curve_crop_then_band"], color=C2, label="crop → band")
ax[2].axvline(ARRIVAL_V/1e3, color=C3, ls=":", lw=1.4)
ax[2].set(xlabel="trial velocity (km/s)", ylabel="moveout score",
          title="(c) The order changes the result")
ax[2].legend()
plt.show()

a_, b_ = STEPS["curve_full"], STEPS["curve_crop_then_band"]
print("band→crop  peak %8.4f at %6.0f m/s" % (a_.max(), VGRID[a_.argmax()]))
print("crop→band  peak %8.4f at %6.0f m/s" % (b_.max(), VGRID[b_.argmax()]))
print("max fractional difference across the scan: %.3f"
      % float(np.max(np.abs(a_-b_)/np.maximum(np.abs(a_),1e-30))))
''')

# ----------------------------------------------------------------- step 9
md("""## 9 — The measurement, and its null

The envelope moveout score is scanned over trial velocity. Significance uses a
**per-velocity** receiver-order permutation null: at each velocity the observation
is compared against permuted realisations **at that same velocity**.

This choice is not cosmetic. A max-over-grid statistic on this data moved from
p = 0.0002 to p = 0.1596 purely by widening the scan limits, with the peak
unchanged — it was tracking the grid, not the data. Per-velocity nulls cannot do
that.""")

code('''
cs = STEPS["curve_full"]; th = STEPS["null_thresh"]; pv = STEPS["p_per_velocity"]
ac = STEPS["acausal"]
clears = cs > th
k = int(np.argmax(cs))

fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.2), constrained_layout=True)
ax[0].plot(VGRID/1e3, cs, color=INK, lw=1.8, label="causal")
ax[0].plot(VGRID/1e3, ac, color=NULLC, lw=1.1, label="acausal")
ax[0].plot(VGRID/1e3, th, "--", color=C2, lw=1.3, label="per-velocity null, 95th")
ax[0].fill_between(VGRID/1e3, cs, th, where=clears, color=C1, alpha=0.22,
                   label="clears its own null")
ax[0].axvline(3.2, color=C3, ls=":", lw=1.4)
ax[0].annotate("3200 m/s\\n(Lellouch)", (3.2, ax[0].get_ylim()[1]*0.92),
               fontsize=8, color=MUTED, ha="left")
ax[0].set(xlabel="trial velocity (km/s)", ylabel="moveout score",
          title="(a) Moveout scan against a per-velocity null")
ax[0].legend(loc="upper right")

ax[1].semilogy(VGRID/1e3, pv, color=C1, lw=1.6)
ax[1].axhline(0.05, color=C2, ls="--", lw=1.2, label="p = 0.05")
ax[1].axvline(ARRIVAL_V/1e3, color=C3, ls=":", lw=1.4,
              label="peak %.0f m/s" % ARRIVAL_V)
ax[1].set(xlabel="trial velocity (km/s)", ylabel="per-velocity p",
          title="(b) Significance versus velocity")
ax[1].legend()
plt.show()

print("peak %.4f at %.0f m/s, per-velocity p = %.4f" % (cs[k], VGRID[k], pv[k]))
print("velocities clearing their own 95th percentile: %d of %d" % (clears.sum(), cs.size))
if clears.any():
    print("  band %.0f-%.0f m/s" % (VGRID[clears].min(), VGRID[clears].max()))
print("causal/acausal at the peak: %.3f" % (cs[k]/ac[k] if ac[k] else float("nan")))
print("per-velocity p at 3200 m/s (Lellouch): %.4f" % float(np.interp(3200, VGRID, pv)))
''')

# ----------------------------------------------------------------- summary
md("""## Ablation summary — every step, and what deleting it costs

Small multiples rather than six overlaid lines: six hues cannot be made mutually
distinguishable under colour-vision deficiency, so each ablation gets its own
panel against the full pipeline in ink.""")

code('''
abl = [("no_diff", "2 — differentiation"), ("no_ram", "3 — RAM normalisation"),
       ("no_neighbour", "6 — R±10 sum"), ("crop_then_band", "8 — band/crop order"),
       ("one_window", "7 — stacking (1 window)")]
fig, ax = plt.subplots(1, len(abl), figsize=(3.3*len(abl), 3.5),
                       constrained_layout=True, sharey=True)
rows = []
for a_, (key, lab) in zip(ax, abl):
    c = STEPS["curve_" + key]
    a_.plot(VGRID/1e3, STEPS["curve_full"], color=INK, lw=1.6, label="full")
    a_.plot(VGRID/1e3, c, color=C2, lw=1.4, label="deleted")
    a_.axvline(ARRIVAL_V/1e3, color=C3, ls=":", lw=1.2)
    a_.set(xlabel="velocity (km/s)", title="step %s" % lab)
    a_.legend(loc="upper right", fontsize=7)
    rows.append((lab, float(c.max()), float(VGRID[c.argmax()])))
ax[0].set_ylabel("moveout score")
plt.show()

full = STEPS["curve_full"]
print("%-26s %9s %10s" % ("step deleted", "peak", "at (m/s)"))
print("%-26s %9.4f %10.0f   <- full pipeline" % ("none", full.max(), VGRID[full.argmax()]))
for lab, pk, at in rows:
    print("%-26s %9.4f %10.0f" % (lab, pk, at))
''')

md("""## What this arrival is, and what it is not

The recovered velocity is **not** the 3,200 m/s P wave of Lellouch et al. (2019) —
the per-velocity p at 3,200 m/s is printed above and is not significant.

The Deep fibre is **wireline-deployed**, not cemented. `MANUSCRIPT_DISCUSSION.md`
§5.1, working from the AWD active-source survey on this same fibre, reports that
the cemented Nano fibre records a mode near **2,950 m/s** while the wireline Deep
fibre records one near **1,547 m/s**, and reads that as a coupling contrast: *"A
cemented fiber is mechanically tied to the formation and senses strain transmitted
through it; a wireline fiber hangs in borehole fluid and couples preferentially to
energy guided by the fluid column."*

So the ambient result here **independently reproduces the active-source
measurement on the same fibre**, from noise alone and with no source. Two further
observations support that reading:

- the arrival's causal/acausal ratio **flips sign at channel 1702**, where the
  fibre reverses direction — a propagating wave must do this; an artefact has no
  reason to;
- three source channels on the outbound limb agree on the velocity independently.

**Open items.** An input-level null (the receiver-order null permutes the finished
gather, so a pre-gather operator sits outside it); the 2026-05 control window, six
weeks before anyone was on site, which tests whether the anthropogenic surface
activity of 2026-06-15 is what excites the mode; and the equivalent run on the
**cemented Nano** fibre, which is the installation that should be able to see a
body wave.""")


def main():
    nb = nbf.v4.new_notebook()
    nb.cells = [nbf.v4.new_markdown_cell(t) if kind == "md"
                else nbf.v4.new_code_cell(t.strip())
                for kind, t in CELLS]
    nb.metadata.update({
        "kernelspec": {"display_name": "das", "language": "python", "name": "das"},
        "language_info": {"name": "python"},
        "ambient_cc_steps_version": NOTEBOOK_VERSION,
    })
    nbf.write(nb, str(OUT))
    n_code = sum(1 for k, _ in CELLS if k == "code")
    print("wrote %s (%d cells, %d code) version %s"
          % (OUT.name, len(CELLS), n_code, NOTEBOOK_VERSION))
    print("execute with:  sbatch exec_ambient_cc_steps.sh")


if __name__ == "__main__":
    sys.exit(main())
