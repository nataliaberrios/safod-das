# Figure 7c across days — the paper-faithful operator on four complete days

Written 2026-08-14. Companion to `AMBIENT_FIG7C_STATUS.md`, which remains the
authoritative single-day record. This file answers one question that document
left open: **the faithful operator had only ever been run on 2024-12-20. Does
another day reproduce Figure 7c?**

It does not.

## Why another day was worth trying

The raw energy census (`ambient_fk_energy_census.py`) measures the share of
5–20 Hz energy inside the 2.5–4 km/s body-wave fan, with no correlation and no
filtering. That share varies by day, from 0.3 % on 2024-05-11 to 4.9 % on
2024-11-30. **2024-12-20, the only day tested, is not the most favourable day
available.** If the arrival were marginal rather than absent, a richer day was
the cheapest way to find it.

## Which days are usable

The operator streams across file boundaries, and its guard refuses to splice a
discontinuous day. Measured directly on the manifest:

| date | files | intervals ≠ 60.000 s | longest continuous run |
|---|---:|---:|---:|
| 2024-12-20 | 1440 | 0 | 1440 |
| 2024-06-17 | 1440 | 0 | 1440 |
| 2024-06-26 | 1440 | 0 | 1440 |
| 2025-02-24 | 1440 | 0 | 1440 |
| 2024-05-11 | 1440 | 0 | 1440 |
| 2024-11-30 | 1441 | 2 | 1278 |
| 2024-10-28 | 1441 | 2 | 1191 |
| 2025-03-04 | 1441 | 2 | 826 |

The three 1441-file days each carry two timing anomalies — on 2024-11-30 the
intervals at file 1277 are 58.806 s then 1.194 s, which sum to 60 s, so it is
jitter rather than missing data. They are **excluded rather than silently
truncated**, and the guard was not weakened to admit them. That leaves the four
continuous days below, run as configuration 0 (paper baseline) only.

## Result

`ambient_lellouch2019_exact_stack_days.sbatch`, jobs `39000540` (96 chunk tasks)
and `39003189` (aggregates), writing to
`ambient_transfer/lellouch2019_exact_stack_days/` so the products cited by
`AMBIENT_FIG7C_STATUS.md` are untouched. Every day: 1,440 files, 5,759 windows,
wellhead source at channel 23, RAM 0.1 s, 5–20 Hz applied to the stacked
correlations, no F-K filter. The p value is familywise over the declared
1.5–6.0 km/s scan using 10,000 receiver-order permutations.

| date | census fan % | peak | at (m/s) | causal/acausal @3200 | null 95 % | peak / null95 | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2024-12-20 (published) | 4.1 | 6.131 | 5850 | 0.97 | 6.324 | 0.97 | 0.1470 |
| 2024-06-17 | 4.5 | 1.901 | 5925 | 1.12 | 2.196 | 0.87 | 0.7413 |
| 2024-06-26 | 3.6 | 1.297 | 5925 | 0.86 | 1.531 | 0.85 | 0.4206 |
| 2025-02-24 | 3.4 | 2.532 | 5875 | 1.05 | 3.400 | 0.74 | 0.7517 |

**No day clears α = 0.05**, and the best of the four is still the originally
published one at p = 0.147. Three findings hold across all four:

1. **The peak is never at 3,200 m/s.** It sits at 5,850–5,925 m/s, the top edge
   of the declared scan, i.e. near-flat moveout — on every day.
2. **The causal side does not consistently dominate.** The ratio at 3,200 m/s
   ranges 0.86–1.12. Figure 7c rests on downgoing energy exceeding upgoing.
3. **The statistic never reaches its own threshold.** peak/null95 spans
   0.74–0.97 and does not approach 1.

## The result that matters most

**Census fan share does not predict the correlation outcome.** 2024-06-17 has the
largest body-wave fan share of any continuous day (4.5 %) and scores *worst*
(p = 0.74); 2024-12-20 has less (4.1 %) and scores best (p = 0.147). So the
failure is not simply "not enough body-wave energy this day, try a richer one" —
picking the day by the very quantity the method needs does not help. That closes
the day-selection avenue rather than leaving it as an open excuse.

## What this closes

Combined with `AMBIENT_FIG7C_STATUS.md`, the reproduction has now been attempted
on the faithful operator across four independent complete days and seven
pre-declared processing branches, with matched white-noise and channel-permutation
nulls, and an eight-day raw-input census. Nothing reaches significance.

**Recommendation: stop attempting the reproduction and write the negative result.**
Further compute on day selection is not justified by anything measured here.

## Not closed

- 2024-05-11 was submitted but its aggregate did not complete; it is the poorest
  day by census (0.3 %) and would not change the conclusion.
- A true multi-day stack. Chunks store summed cross spectra and are therefore
  exactly summable across days, but the aggregate action is single-date, so
  pooling needs a small addition to the script. Given that no individual day
  shows a trend toward the arrival, this is a completeness exercise.
- The three partially continuous days could be run over their longest continuous
  runs (1278, 1191, 826 files) if a reviewer asks.
