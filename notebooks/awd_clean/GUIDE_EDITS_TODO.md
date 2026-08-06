# Edits owed to AWD_advisor_figure_guide.tex (and STATUS.md), from the 2026-08-05 preprocessing audit

Nothing in the existing analysis was changed. This is the to-do list that came
out of reading the reader chain end to end; see [`PREPROCESSING.md`](PREPROCESSING.md)
for the evidence behind each item. Reviewed against
`AWD_advisor_figure_guide.tex` at v43 (the unsuffixed file is byte-identical to
`_v43.tex`).

**Headline: the accepted/preliminary structure of the guide survives intact.**
Three text edits and one number to recheck. No figure needs recomputing.

---

## Owed edits

### 1. Multiscale record section — add a caveat, or reselect the drop
*Section: "Supplemental diagnostic: What does stacking actually contribute?",
figure `awd_multiscale_record_sections.png`.*

The selected drop (burst 34, drop 10, per `awd_multiscale_record_sections.txt`)
sits 14.9 s into a 300 s Nano file and **54.9 s into a 60 s Deep file**. Across
the displayed Deep window the read-time Tukey gain falls from ≈0.45 to ≈0.08 —
a **5.3× decay with time**.

The caption's per-trace normalisation ("divided by its RMS amplitude in the
common −0.45 to −0.15 s pre-reference interval") removes the constant part of
that gain but **not the ramp**, because the noise window sits at one end of it.
So the apparent amplitude decay with time in the single-drop Deep panel, out to
2.85 s, is partly instrumental.

The burst-mean (col 2) and 859-drop weighted-mean (col 3) columns are fine.

Cheapest fix: keep the existing mechanical selection rule and add "and whose
reader taper window is flat" as a stated criterion, then note it in the caption.
No re-stacking required.

### 2. Colorbar label — Deep is strain, not strain rate
`multiscale_record_sections.py:330,332` label both colorbars "Band-passed
strain-rate amplitude". Nano (Sintela) is strain rate; Deep (OptaSense) is
strain. `readFile_HDF` differentiates only with `diff=True`, which no AWD caller
passes. The image is in noise-RMS units so the picture is not wrong — only the
label is.

### 3. `docs/paper1/STATUS.md` — recheck the 39% source amplitude CV
Not in this guide, but it is the single number most exposed to the taper: a raw
burst-amplitude statistic. Recompute using the per-drop weights in
`figures/awd_2026/plain_look/taper_audit.csv` before it goes in the paper.

---

## Checked and NOT at risk — do not redo these

| Guide item | Why it holds |
|---|---|
| Core 1 geometry | not a measurement |
| Core 2 Nano mode, signed F-K | ridge position; phase-based |
| Core 3 repeatability hierarchy | `nano_hierarchical_repeatability.py:448` sets `tapering=False` |
| Core 4 installation contrast | phase-only beams; caption already says "coherence and signed slowness, not uncalibrated absolute amplitude" |
| Core 5 aperture / injection–recovery | see Figure 11 below |
| Core 6 PGSI registration | timing and moveout |
| Fig 8 Deep split-sample repeatability | phase-only beam, see below |
| Fig 9 Deep channel-permutation null, p=0.002 | phase-only beam, see below |
| Fig 11 frequency-dependent aperture | see below |
| Figs 18–20 tidal-timescale delay + CCA | see below |
| Figs 21–22 Deep target scan | `fixed_score` is semblance over a 60 ms window; a 2× ramp over 3.5 s changes by 0.7% across 60 ms |
| Whole ambient F-K chain | reads `h5py` directly, never calls DASutils |

### Why the Deep beam results are safe
`deep_tube_validation.py:39-42` divides each spectrum by its own modulus before
beamforming — unit modulus, amplitude discarded. And the taper is
**channel-independent**: it is `g(t)`, identical on every channel. A gain that is
common to all channels cannot create or destroy inter-channel phase coherence.

This specifically answers the worry that a taper-induced amplitude gradient along
fiber could beat the channel-order permutation null in Figure 9. It cannot:
amplitude is normalised out before the beam is formed.

### Why Figure 11 is safe
Signal windows follow `t = −0.022 + x/2975`, noise is −0.38 to −0.14 s. The
signal window at 460 m is 0.133 s later than at 100 m, so a decaying ramp would
exaggerate the SNR decline with distance — by **0.23 dB**, against an observed
21.8 → 7.6 dB decline. Not material.

### Why the tidal branch is safe
`nano_tidal_compaction_regression.py:30-31` builds both its response
(`loo_burst_delay_s`) and its coupling proxy (`burst_signal_rms`) from
`nano_burst_repeatability_hierarchical.csv`, written by the one script that sets
`tapering=False`. It never touches the canonical stacks for waveform samples.

---

## The underlying facts, for the caption text

Read-time `tukey(n_file, alpha=0.4)` is flat only between 20% and 80% of each
file. Over the **859 drops the stacks actually see** (both fibers, cut window
inside the file):

| | attenuated | weight below 0.5 | >20% ramp across the 3.5 s window |
|---|---|---|---|
| Nano (300 s files) | 416 (48.4%) | 226 | 227 (26%) |
| Deep (60 s files) | 298 (34.7%) | 136 | 278 (32%) |

Burst-mean Nano weight swings 0.621 → 0.996 across the survey, a 1.60×
non-physical amplitude excursion peaking around bursts 21–27.

Going forward the fix is one kwarg — `tapering=False` — exactly as
`nano_hierarchical_repeatability.py` already does. Existing products can be
corrected by division using `taper_audit.csv` rather than recomputed.
