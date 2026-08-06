# What actually happens to AWD data before you see it

Written 2026-08-05 by reading `DASutils.py` and every AWD script that calls a
reader. Nothing in the existing analysis was changed to produce this. Line
numbers refer to `DAS-utilities/python/DASutils.py` in this checkout.

The point of this file is that the preprocessing is invisible: it happens inside
one function call, with defaults, in a fixed order, and several scripts then
filter again on top of it.

---

## 1. What the reader does when you don't tell it otherwise

`readFile_protobuf` (line 1765) and `readFile_HDF` (line 1560) take the same
keyword defaults:

| kwarg | default | effect |
|---|---|---|
| `detrend` | `True` | linear detrend per trace |
| `tapering` | `True` | multiply by `tukey(n_samples_in_read, alpha=taper)` |
| `taper` | `0.4` | flat only between 20% and 80% of the record |
| `filter` | `True` | Butterworth bandpass at the `fmin`/`fmax` you pass |
| `order` | `4` | filter order |
| `zerophase` | `True` | `sosfiltfilt`, so no phase shift |
| `median` | `True` | `preprocess_medfilt`: subtract the across-channel median at every time sample |
| `diff` | `False` | differentiate — **only applied if `system='OptaSense'` is also passed** (line 1725) |
| `unwrap` | `False` | int32 overflow correction, OptaSense only |
| `desampling` | `True` | resample to 2.5×`fmax`; every AWD caller passes `False` |

Applied in this order, per channel buffer:

```
detrend  →  Tukey taper  →  bandpass(fmin, fmax)  →  [whole array] median removal  →  ×1e6 to microstrain
```

Two consequences that are easy to miss:

- **The taper spans the whole read, not the window you care about.** Callers pass
  one file at a time, so it is 300 s for Nano and 62 s for Deep. See §4.
- **`median=True` is a physics assumption, not a cleanup step.** It removes laser
  and interrogator noise shared by all channels, and it equally removes any real
  arrival that reaches every channel at once. `nano_fiber_end_coherence.py:60`
  already documents the other side of this: keeping the common mode is what makes
  the low band coherent at all.

## 2. Units: the two fibers are not in the same quantity

| Fiber | System | Stored quantity | What the reader returns |
|---|---|---|---|
| Nano | Sintela protobuf | strain **rate** | microstrain/s |
| Deep | OptaSense HDF5 | **strain** | microstrain |

`readFile_HDF` only differentiates when given `diff=True` *and* an explicit
`system='OptaSense'`. No AWD caller passes either, so Deep comes back as strain.

Three scripts reconcile this by hand with `np.gradient`: `awd_spectra.py:93`,
`wireline_tube_look.py:61`, `tube_wave_gate.py:122`. The paired-stack and
record-section paths do not.

A factor of ω is a 90° phase rotation plus a slope across the band. It changes
waveform shape and cross-fiber amplitude. It does **not** move an F-K ridge or
change a moveout velocity.

Note also `multiscale_record_sections.py:330,332` labels both colorbars
"Band-passed strain-rate amplitude". The Deep panel is strain. The image itself
is in noise-RMS units so the picture is not wrong; the label is.

## 3. Per-script inventory

| Script | Fiber | Read band | Reader settings | Second filter applied later |
|---|---|---|---|---|
| `paired_stack_job.py` | both | 1–100 Hz | all defaults | none in-script |
| `paired_stack_job_deep_all.py` | both | 1–100 Hz | all defaults | none in-script |
| `compute_deep_stack.py` | Deep | 1–100 Hz | all defaults | none in-script |
| `multiscale_record_sections.py` | both | 1–100 Hz | all defaults | Nano 30–60, Deep 3–15 Hz |
| `nano_hierarchical_repeatability.py` | Nano | 30–60 Hz | `tapering=False`, `median=True`, `detrend=True`, explicit | none |
| `nano_fiber_end_coherence.py` | Nano | 1–8 / 30–90 Hz | `median` toggled `False`/`True` deliberately | none |
| ambient F-K chain (`ambient_transfer_test.py` etc.) | — | — | **does not use DASutils**; reads `h5py` directly | own 5–20 Hz + 5 s RAM norm |

Everything downstream of `canonical_epoch_stacks_paired_deep_all.npz`
(`measure_repeatability.py`, `null_tests_and_inspection.py`, `deep_target_*.py`,
`mode_scan`, the tidal regressions) inherits the `paired_stack_job_deep_all.py`
row above, then applies its own band on top.

### Cascaded filters

The stacks were read through a 4th-order zero-phase 1–100 Hz Butterworth. Scripts
then apply another 4th-order zero-phase Butterworth at 30–60, 20–50 or 3–15 Hz.
The result is 8th-order in the overlap with sharper skirts than either alone.
Harmless for the science, but "30–60 Hz" in a caption is not the true response.

## 4. The read-time taper — the one that bites

`tukey(n, alpha=0.4)` is flat only from 20% to 80% of the read. Callers read one
file at a time, so:

- **Nano**, 300 s files: flat 60–240 s. The first and last 60 s are attenuated.
- **Deep**, 62 s files: flat 12–50 s. The first and last 12.4 s are attenuated.

Drops are timed by the survey, not by file boundaries, so a drop's weight is set
by where it happened to land. From `awd_manifest.csv` (`taper_audit.py`):

| | drops | attenuated | weight below 0.5 |
|---|---|---|---|
| Nano | 988 | 483 (48.9%) | 271 |
| Deep | 926 | 322 (34.8%) | 163 |

Burst-mean Nano weight swings from **0.640 to 1.000 across the survey** — a 1.56×
smooth, non-monotone amplitude excursion with no physical cause, peaking around
bursts 21–27.

There are two distinct effects:

**A per-drop scale factor.** Survived by timing and by normalised
cross-correlation. Not survived by per-drop or per-burst amplitude, RMS, SNR or
energy, or by anything regressed against them.

**A gain ramp across the 3.5 s cut window.** The taper is not flat over the
`PRE_S=0.5`/`POST_S=3.0` window. On Deep, 3.5 s is 5.6% of a 62 s file, so a drop
inside the ramp gets a time-varying gain: late arrivals are suppressed relative
to early ones. This is the damaging one, because it survives per-trace noise-RMS
normalisation (the noise window sits at one end of the ramp) and mimics
attenuation with time, and therefore with distance along the fiber.

Worked example — the drop in the published multiscale record section
(`awd_multiscale_record_sections.txt`: burst 34, drop 10, Nano file
`..._16.48.28_...`, drop at 16:48:42.92, i.e. 14.9 s into a 300 s file):

- Nano weight ≈ 0.14 at the window start, ≈ 0.20 at the end — a 1.5× ramp.
- Deep offset is 54.9 s into a 62 s file: weight falls from ≈ 0.67 to ≈ 0.26
  across the displayed record — a **2.6× decay with time**, in the panel used to
  look for late, slow, deep arrivals.

## 5. What this does and does not put at risk

**Not at risk.** Anything phase- or timing-based, because the taper is real,
positive and zero-phase, and the median removal is a per-sample constant across
channels:

- apparent velocity, moveout, semblance peak positions;
- F-K and signed F-K ridge locations, and their permutation nulls;
- normalised cross-correlation repeatability (NCC is amplitude-invariant);
- `nano_hierarchical_repeatability.py` entirely — it sets `tapering=False`;
- the whole ambient F-K chain — different code path, never calls DASutils.

**At risk, needs checking.**

1. **Any burst-level amplitude or RMS observable** built from the canonical
   stacks. The 1.56× taper excursion is smooth and multi-hour.
2. **The burst-delay tidal-timescale regression**, which uses signal RMS as a
   source-history/coupling proxy. A smooth non-monotone amplitude excursion of
   unrelated origin is exactly the pattern that branch is fitting. This is the
   single most important item to re-check.
3. **Late-arrival amplitude in Deep**, including the slow-mode candidates and the
   time-gated branch search, wherever the drop sat inside a 62 s file's ramp.
4. **Cross-fiber amplitude and waveform-shape comparison**, from the strain vs
   strain-rate mismatch in §2.

**The cheap fix**, if any of the above needs redoing: pass `tapering=False` to
the reader and taper the 3.5 s cut window instead, or read with `minTime`/
`maxTime` around each drop as `nano_hierarchical_repeatability.py` already does.
`taper_audit.csv` gives the per-drop weight and ramp, so an existing product can
also be corrected after the fact rather than recomputed.
