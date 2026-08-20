### 4.6 Adaptive f-k

The adaptive filter of Isken et al. (2022) is the one f-k family that makes
no velocity assumption, and therefore the one not addressed by section 4.5.
Before use it was verified that an exponent of zero is a bit-identical no-op,
that Bartlett 50 %-overlap recombination reconstructs the input to a maximum
relative error of 0.000e+00, that it raises the coherence of a synthetic
3,200 m/s plane wave in noise from 0.0523 to 0.1980, and that it does not
manufacture moveout from pure noise (0.0000 to 0.0000).

Five configurations were run on a common 300-record block of 2024-12-20,
with the prediction fixed in advance that applying the filter *without*
prior static-pattern removal would be the worst configuration, because an
adaptive filter enhances the dominant coherent component and here that is the
static pattern rather than an arrival.

| configuration | peak | at (m/s) | p | pedestal corr | causal/acausal at 3,200 | recovered |
|---|---:|---:|---:|---:|---:|:--:|
| afk1 only, no static removal (predicted WORST) | 12.5285 | 5900 | 0.0060 | +0.987 | 1.03 | no |
| median common mode, no AFK (baseline) | 1.9069 | 5900 | 0.8796 | +0.820 | 0.89 | no |
| median common mode + AFK alpha=1 | 1.2109 | 5900 | 0.5199 | +0.766 | 1.00 | no |
| median common mode + AFK alpha=2 | 1.1221 | 5900 | 0.3465 | +0.737 | 1.01 | no |
| median common mode + rank-2 + AFK alpha=1 | 1.0781 | 3300 | 0.9392 | -0.036 | 1.01 | no |

Recovery requires all five predeclared gates: pedestal suppressed
(|corr| < 0.5), peak inside 2,500-4,000 m/s, peak not at a scan edge,
causal side dominant at 3,200 m/s, and p < 0.05.

**No configuration satisfied the gates.** Adaptive f-k does not
recover the arrival either. Given section 4.7 this is expected rather
than surprising: a filter can only enhance coherent energy that is
present, and the wavefield carries no net downgoing component to
enhance. It does, however, close the remaining f-k avenue by
measurement rather than by argument.

**A filter can manufacture significance, and here one did.** The
configuration with no prior static removal reached p = 0.0060 -- which
in isolation reads as a detection -- while carrying the worst
pedestal diagnostic of the five (+0.987, i.e. almost pure pedestal)
and peaking at 5900 m/s, the top of the scan, rather than at 3,200.
The mechanism is specific and worth stating, because it applies to
any filtered ambient result: the receiver-order null permutes the
FINISHED gather, so an operator applied before the gather is formed
sits outside its own null and its amplification of a coherent
contaminant is never tested. The adaptive filter raised the score
6.6-fold (1.91 to 12.53) and the amplified pedestal then cleared a
null that could not see the amplification. Only the predeclared
gates caught it. An input-level null -- built before the operator
runs -- is the correct control for a filtered result, and is what
the F-K QC workflow already requires elsewhere in this project.

**The cleanest statistic this study produced still shows nothing,
and that is the strongest form of the negative.** Stacking the
removals -- median common mode, then rank-2 subspace, then the
adaptive filter -- drives the pedestal diagnostic monotonically to
-0.036 (median common mode + rank-2 + AFK alpha=1), against +0.987 for the unprocessed baseline. That is
effectively zero: the moveout statistic is finally measuring
moveout rather than proximity to the zero-lag lobe, which no
earlier configuration in this project achieved (previous best
-0.381). With the statistic clean, the result is p = 0.9392.

Its peak falls at 3300 m/s, close to the published 3,200 m/s, and
that coincidence should not be read as encouraging. At p = 0.939 the
observed maximum is LOWER than most receiver-order permutations of
the same data. With no pedestal pulling the peak to the scan
ceiling, the peak is free to land anywhere, and it landed there.

The value of this row is what it rules out. The failure to
reproduce is not an artefact of a broken statistic: when the
statistic is repaired, the arrival is still absent.

On the pre-registered prediction: `afk1 only` gives pedestal
corr = +0.987 and p = 0.0060 against the baseline's +0.820 and 0.8796, so the
pre-registered prediction -- that applying the adaptive filter
without prior static removal makes matters worse -- is CONFIRMED.

