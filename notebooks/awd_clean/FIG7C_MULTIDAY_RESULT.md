# Figure 7c across five complete days — paper-faithful

Written 2026-08-14. Companion to `AMBIENT_FIG7C_STATUS.md`, which covers
2024-12-20 in full. This file records the extension to four further days and is
the basis for widening that document's claim from one day to the archive.

## Why more days

`AMBIENT_FIG7C_STATUS.md` tests 2024-12-20 only, matched to the paper's own
one-day design. But the raw energy census (`ambient_fk_energy_census.py`) shows
the 2.5–4 km/s fan share varies by day, 0.3 % to 4.9 %, so the day that had been
tested was not the most favourable one available. If the arrival exists anywhere
in this archive it should appear on the richest day.

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

| date | census fan % | peak | at (m/s) | @3,200 causal | @3,200 acausal | null 95 % | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2024-12-20 | 4.1 | 6.131 | 5850 | 2.752 | 2.831 | 6.324 | 0.1470 |
| 2025-02-24 | 3.4 | 2.532 | 5875 | 1.806 | 1.725 | 3.400 | 0.7517 |
| 2024-06-17 | 4.5 | 1.901 | 5925 | 1.185 | 1.054 | 2.196 | 0.7413 |
| 2024-06-26 | 3.6 | 1.297 | 5925 | 0.948 | 1.105 | 1.531 | 0.4206 |
| 2024-05-11 | 0.3 | 1.026 | 1675 | 1.014 | 0.994 | 1.034 | 0.3091 |
| 2024-11-30 † | **4.9** | 6.536 | 5850 | 2.936 | 3.056 | 6.613 | 0.1345 |

† 21.3 h continuous block (1,278 of 1,441 rows, 5,040 windows) via the opt-in
`--continuous-prefix` flag. **This is the richest day in the whole census** and the
last one untestable under the strict guard; it is now tested and it is the closest
any day comes (p = 0.1345) while still failing, still peaking at 5,850 m/s rather
than 3,200, and still with causal/acausal below 1.

**No day reaches significance; the minimum p over all six days is 0.1345, on the
richest day in the archive.** Fisher's
combination across the five independent days gives χ² = 9.08 on 10 degrees of
freedom, **p = 0.524** — the χ² ≈ df expected from pure noise, with no residual
signal hiding below per-day significance.

Two further points against the arrival being present but weak:

- **The peak is in the wrong place on four of five days**, sitting at 5,850–5,925
  m/s, the top edge of the declared scan, i.e. near-flat moveout rather than the
  paper's ~3,200 m/s.
- **The census does not predict the outcome.** 2024-06-17 at 4.5 % scores
  p = 0.7413, worse than the poorest day; 2024-11-30 at 4.9 % scores p = 0.1345,
  no better than 2024-12-20 at 4.1 %. Whatever occupies the
  body-wave fan in the raw energy budget is not organised into a receiver-ordered
  arrival, so "pick a better day" is not an available fix.

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
operator, on five independent complete days spanning ten months, including the
most energy-favourable day available, nor on a coherent 96-hour stack of the four
days that share an acquisition rate. Taken with the eight-day raw census — 81–99 %
of 5–20 Hz energy below 1,500 m/s and downgoing share never exceeding 49.9 % — the
required downgoing body-wave energy is not present in the input, so no processing
choice recovers it.

This is a statement about this archive, this band, and these days. It is not a
criticism of Lellouch et al. (2019), whose Figure 9 companion model this project
reproduces from the released Figure 7d correlograms at r = 0.948.
