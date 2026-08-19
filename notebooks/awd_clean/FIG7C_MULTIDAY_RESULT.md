# Figure 7c across five complete days — paper-faithful

Written 2026-08-14. Companion to `AMBIENT_FIG7C_STATUS.md`, which covers
2024-12-20 in full. This file records the extension to four further days and is
the basis for widening that document's claim from one day to the archive.

> **CENSUS CLAIMS WITHDRAWN 2026-08-14, same day.** Every statement in this file
> that rested on `ambient_fk_energy_census.py` is retracted, following an
> independent audit. The census had no geometric baseline: **98.4 % of in-band
> (f,k) cells lie below 1,500 m/s by construction**, so its headline "81–99 % of
> energy below 1,500 m/s" is *at or below* what white noise gives and cannot
> distinguish the SAFOD wavefield from noise. Measured as energy *density* per
> cell the data are 2.6–2.8× **enriched** at body-wave velocities — the opposite
> of what was claimed. Three further defects: the k = 0 column silently dropped
> 22 % of band energy; no spatial taper, so the fan measurement is largely
> Dirichlet leakage from the k ≈ 0 peak; and channels 0–22 carry **97 %** of the
> array's 5–20 Hz energy, so any census run without the `_ch23-896` suffix is
> determined by 23 uncemented lead-in channels. The "downgoing share never
> exceeds 49.9 %" claim is simply **false** — the tree's own outputs include
> 51.6 %.
>
> **The per-day and multi-day results below are unaffected** and were verified
> independently: the chunk re-sum reproduces the stored aggregates to max abs
> difference 0.000e+00. Only the census-derived *motivation* and the
> census-derived *mechanism* are withdrawn. The non-reproduction therefore
> currently has **no validated mechanism**.

## Why more days

`AMBIENT_FIG7C_STATUS.md` tests 2024-12-20 only, matched to the paper's own
one-day design. One day is a thin basis for a claim about an archive, so the
paper-faithful operator was extended to every other complete day. (The original
motivation — that a raw energy census identified some days as richer than others
— is withdrawn; see the banner above. The days were therefore effectively an
unselected set, which is a *stronger* basis for the result than a selected one.)

## What was run

`ambient_lellouch2019_exact_stack.py` configuration 0 (paper baseline: wellhead
source channel 23, RAM 0.1 s, 30 s windows at 15 s overlap over a contiguous
day, literal R±10 sums, simple unshifted stacking, 5–20 Hz applied to the
stacked correlations, no F–K filter), via
`ambient_lellouch2019_exact_stack_days.sbatch`, jobs `39000540` / `39003189` /
`39004321`, output in `ambient_transfer/lellouch2019_exact_stack_days/`.

Only days whose manifest is **exactly** continuous at 60.000 s are usable — the
operator streams across file boundaries and its guard rightly refuses to splice a
discontinuous day. Measured over the manifest, 2024-11-30, 2024-10-28 and
2025-03-04 each carry two timing anomalies (longest continuous runs 1,278, 1,191
and 826 files) and were excluded rather than silently truncated. **2024-11-30 is
the richest day in the census at 4.9 % and is among those excluded**; the richest
usable day is 2024-06-17 at 4.5 %.

## Result

Each day is a complete 1,440-file, 5,759-window stack with its own 10,000-draw
receiver-order familywise null over the declared 1.5–6.0 km/s scan.

The "census fan %" column that previously appeared here has been removed as
invalid (see banner). Peak velocities are **grid ceilings, not velocity
estimates** — see the pedestal note below.

| date | peak | at (m/s) | @3,200 causal | @3,200 acausal | null 95 % | p |
|---|---:|---:|---:|---:|---:|---:|
| 2024-12-20 | 6.131 | 5850 | 2.752 | 2.831 | 6.324 | 0.1470 |
| 2025-02-24 | 2.532 | 5875 | 1.806 | 1.725 | 3.400 | 0.7517 |
| 2024-06-17 | 1.901 | 5925 | 1.185 | 1.054 | 2.196 | 0.7413 |
| 2024-06-26 | 1.297 | 5925 | 0.948 | 1.105 | 1.531 | 0.4206 |
| 2024-05-11 | 1.026 | 1675 | 1.014 | 0.994 | 1.034 | 0.3091 |
| 2024-11-30 † | 6.536 | 5850 | 2.936 | 3.056 | 6.613 | 0.1345 |

† 21.3 h continuous block (1,278 of 1,441 rows, 5,040 windows) via the opt-in
`--continuous-prefix` flag — the last day untestable under the strict continuity
guard. It is the closest any day comes (p = 0.1345) while still failing, still
peaking at 5,850 m/s rather than 3,200, and still with causal/acausal below 1.

**No day reaches significance; the minimum p over all six days is 0.1345.** Fisher's
combination across the five independent days gives χ² = 9.08 on 10 degrees of
freedom, **p = 0.524** — the χ² ≈ df expected from pure noise, with no residual
signal hiding below per-day significance.

Two further points against the arrival being present but weak:

- **The observed score never clears the PER-VELOCITY null, at any velocity.**
  0 of 181 velocities on 2024-12-20, and at 3,200 m/s specifically 2.752 against
  a per-velocity threshold of 2.811. That is the maximally powerful test with no
  multiplicity penalty at all, so the negative does not depend on the familywise
  correction.
- **The peak velocity is an artefact of the statistic, not a measurement.** The
  score samples each trace in a gate at t = offset/v, so as v rises every gate
  slides toward zero lag where these gathers have a dominant broad lobe. Hence
  `corr(trial velocity, score) = +0.976` and the argmax simply lands wherever the
  grid stops: cap the scan at 3,500 / 4,000 / 8,000 / 20,000 m/s and the "peak"
  moves to 3,475 / 3,925 / 7,700 / 18,775. A moveout-free control — every trace
  replaced by the across-receiver median — reproduces the observed curve to
  within 3.5 % at every velocity, including the exact peak location. **The
  p-values remain valid**, because the receiver-order permutation preserves the
  gate geometry and so carries the identical bias (97 % of null curves also peak
  at ≥ 5,500 m/s), but no number in the "at (m/s)" column should be read as a
  velocity.

## The coherent stack — the most sensitive test available

Combining five p values with Fisher's method is not the same experiment as
stacking the data. A weak coherent arrival grows with the number of stacked
windows while incoherent noise averages down, so one long stack is strictly more
sensitive than several short tests however their p values are pooled.
`ambient_lellouch2019_multiday_stack.py` therefore sums the chunk cross spectra —
exactly additive, so the result is identical to having processed one continuous
record — and imports the scoring, bandpass and null from the single-day operator
so the arithmetic is bit-identical.

**The archive is not homogeneous in acquisition rate.** 2024-05-11 was recorded at
**5000 Hz** (`n_fft` 524288, `ram_samples` 501) against 500 Hz on every other day.
Its cross spectra live on a different frequency grid and cannot be coherently
pooled, so it is excluded from the stack; its one-day result above still stands.
This also explains why that day was consistently slow and memory-hungry. Any
future analysis pooling across this archive must group by acquisition rate.

Stacking the four 500 Hz days — 5,760 files, **23,036 windows, 96.0 hours**:

| quantity | value | requirement |
|---|---:|---|
| peak causal score | 1.912 | — |
| at velocity | 5,925 m/s | paper: ~3,200 m/s |
| score at 3,200 m/s, causal | 1.143 | — |
| score at 3,200 m/s, acausal | 1.156 | causal must dominate |
| causal / acausal at 3,200 m/s | **0.99** | **> 1** |
| receiver-order familywise null 95 % | 2.649 | — |
| **p** | **0.9184** | < 0.05 |
| detectability, peak / null 95 % | **0.72** | must reach 1.00 |

Four times the data does not move the statistic toward its threshold — it sits
further from significance than the single best day. That is the behaviour of
noise, not of an arrival too weak to see in 24 hours.

## Day-pair survey — why a p < 0.05 here is not a reproduction

Stacking the two days that individually came closest (2024-11-30 and 2024-12-20)
returns p = 0.039. That is not a reproduction, and the reason it is not is worth
recording because it is the exact shape of the mistake this whole thread exists
to avoid.

Those two days were chosen *because* they had the lowest p values. Running all
ten pairs of the five usable 500 Hz days (`_fig7c_pair_survey.py`, output in
`fig7c_pair_survey.txt`, 2,000-permutation null):

| pair | peak | at (m/s) | p | causal/acausal at 3,200 |
|---|---:|---:|---:|---:|
| 06-17 + 06-26 | 1.643 | 5925 | 0.8251 | 0.93 |
| 06-17 + 11-30 | 2.463 | 5775 | 0.5582 | 1.22 |
| 06-17 + 12-20 | 2.521 | 5775 | 0.5302 | 1.23 |
| 06-17 + 02-24 | 1.980 | 5925 | 0.7466 | 1.19 |
| 06-26 + 11-30 | 1.618 | 4850 | 0.7851 | 0.89 |
| 06-26 + 12-20 | 1.661 | 4800 | 0.8016 | 0.91 |
| 06-26 + 02-24 | 1.379 | 5925 | 0.5607 | 0.86 |
| **11-30 + 12-20** | 6.373 | 5850 | **0.0390** | 0.97 |
| 11-30 + 02-24 | 5.735 | 5900 | 0.3188 | 1.01 |
| **12-20 + 02-24** | 5.937 | 5850 | **0.0165** | 1.01 |

Two of ten fall below 0.05 against 0.5 expected by chance — binomial P(≥2) ≈ 0.086,
and the pairs share days so they are not independent, making the real correction
weaker still. Selecting the best pair of ten and quoting its uncorrected p is not
a result.

**The decisive objection is not multiplicity, it is velocity.** Every pair that
reaches p < 0.05 peaks at **5,850–5,900 m/s** with causal/acausal of 0.97–1.01. A
peak at the top edge of the scan is a *flat-moveout* feature — energy arriving at
every receiver essentially simultaneously across 700 m — which is the signature of
common-mode or instrumental structure, not a propagating body wave. Figure 7c is a
specific claim: a packet near 3,200 m/s with the causal side dominant. At 3,200 m/s
the score is consistently about half the peak and the causal side never dominates.

So this data does contain a statistically detectable coherent component, and it is
emphatically not the arrival the paper reports. `ambient_lellouch2019_multiday_stack.py`
now requires all three conditions — p < 0.05, peak inside the 2,500–4,000 m/s fan,
and causal dominance — before it will call anything a reproduction.

## Conclusion

Figure 7c does not reproduce on the 2024–2025 archive, on the paper's own
operator, across six independent complete days spanning ten months, nor on a
coherent 96-hour stack of the four days that share an acquisition rate. The
minimum p over six days is 0.1345, and the observed score never clears the
per-velocity null at any of 181 velocities — so the negative does not rest on a
multiplicity correction.

**The mechanism is open.** This document previously concluded that the required
downgoing body-wave energy was absent from the input. That explanation came from
`ambient_fk_energy_census.py` and is withdrawn (see banner): its statistic had no
geometric baseline and could not distinguish the SAFOD wavefield from white
noise. A cross-epoch comparison against Lellouch's released 2017 earthquake
records was then attempted as a replacement and is *also* withdrawn — the two
arms were not processed alike, because `extract_all.py` takes DASutils
`median=True` by default and so stripped the common mode from the 2024–25 arm
only. **No validated mechanism for the non-reproduction currently exists.**

One thing the cross-epoch work does still support: gross sensing failure is ruled
out. The 2024–25 array records earthquake wavefronts at 92 % of the 2017 array's
semblance, at lower SNR. That is a normalised within-record measure of a loud
transient, so it is far less sensitive to the processing asymmetry than the
coherence curves were — but it is not entirely immune, and should be redone with
both arms read identically before being relied on.

This is a statement about this archive, this band, and these days. It is not a
criticism of Lellouch et al. (2019), whose Figure 9 companion model this project
reproduces from the released Figure 7d correlograms at r = 0.948.
