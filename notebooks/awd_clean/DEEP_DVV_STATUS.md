# Deep guided-mode sensitivity — status and claim table

Written 2026-08-06. Companion to the frozen
[`DEEP_DVV_PREREGISTRATION.md`](DEEP_DVV_PREREGISTRATION.md). Where numbers here
disagree with an older notebook cell or a figure caption, this file wins.

The question: what is the smallest fractional change in the apparent speed of the
Deep guided mode that can be recovered reliably, and does it beat the Nano
apparent-moveout observable?

---

## 1. Where this stands

| Claim | Status | Evidence |
|---|---|---|
| Deep outbound is more sensitive than Nano | **Established, restricted to the outbound branch** | reliable detection 5×10⁻³ vs Nano 1×10⁻²; null threshold 1.84×10⁻³ vs 3.25×10⁻³ |
| The Deep *installation as a whole* beats Nano | **Not claimable** | return leg reaches 1×10⁻², equal to Nano, not below it — pre-registered Rule 2 |
| The 11.5× lever arm buys 11.5× sensitivity | **Refuted** | it buys ~2×; per-burst timing repeatability is 7.4× worse and consumes the rest |
| Estimator separates source timing from propagation | **Established** | 5 ms constant shift → 4.87 ms into the intercept, 7×10⁻⁷ into ε |
| Recovery is unbiased on the primary branch | **Established** | bias ≤6×10⁻⁵ at every level, outbound |
| Result carried by one aperture or one burst | **Refuted** | no single aperture moves it >0.55× threshold, no burst >0.27× |
| Combining the two legs improves sensitivity | **Partial — improves the noise floor, not the detection level** | scatter 1.36–1.48× better; reliable detection unchanged at 5×10⁻³ |
| Return-leg weakness is SNR-limited | **Unresolved** | 5 of 6 pre-registered criteria met; the 6th rests on a control that turned out non-diagnostic |
| Sensitivity is attributable to the guided mode specifically | **Inherited, not shown here** | rests on the earlier permutation nulls (p = 0.002), not on this experiment |

## 2. Headline numbers

Held-out population (23 bursts, `epoch % 2 == 0`), 15–30 Hz, primary pass.
Nano row from `nano_dvv_injection_recovery.npz`, same estimator code.

| Observable | Null threshold | Reliable detection | Sign-only | Null scatter (MAD) |
|---|---|---|---|---|
| Nano apparent moveout, 46 bursts | 3.25×10⁻³ | 1×10⁻² | — | 1.86×10⁻³ |
| **Deep outbound** | **1.84×10⁻³** | **5×10⁻³** | 2×10⁻³ | 1.20×10⁻³ |
| Deep return | 3.05×10⁻³ | 1×10⁻² | 5×10⁻³ | 1.32×10⁻³ |
| Deep paired, equal weight | 1.82×10⁻³ | 5×10⁻³ | 2×10⁻³ | 8.84×10⁻⁴ |
| Deep paired, inverse-variance | 1.53×10⁻³ | 5×10⁻³ | 2×10⁻³ | 8.12×10⁻⁴ |
| Deep paired, covariance-aware | 1.67×10⁻³ | 5×10⁻³ | 2×10⁻³ | 7.87×10⁻⁴ |

"Reliable detection" is the Nano-comparable metric: smallest tested |ε| where
≥95% of trials both exceed the empirical null threshold and carry the correct
sign. Frozen trajectories: outbound 1544.6 m/s at t₀ = +0.100 s, return
1549.7 m/s at t₀ = +0.346 s — two independent halves agreeing to 0.3%.

## 3. The lever arm does not convert one-for-one

The central quantitative result of this analysis.

| | Nano | Deep outbound |
|---|---|---|
| Regression lever arm | 0.121 s | 1.398 s (11.5×) |
| Null scatter | 1.86×10⁻³ | 1.20×10⁻³ (1.55×) |
| Implied per-burst timing repeatability | 0.225 ms | 1.68 ms (7.4× **worse**) |

Deep's geometric advantage is real and large, but poorer burst-to-burst timing
repeatability consumes about four fifths of it. This is consistent with the
Deep-fibre per-shot repeatability deficit already recorded in
`docs/paper1/STATUS.md` (per-shot CC 0.054).

Do not quote the `lever_arm_s` field in the generated report as the lever arm:
that field is the channel span (1.808 s). The regression runs over aperture
centres, so the operative figure is **1.398 s**.

## 4. Paired-leg combination

Pre-registered as a comparison against outbound, never an automatic replacement.

**Between-leg null error correlation ρ = +0.121.** The two legs record the same
source impacts, so correlated errors were expected; they are almost absent. That
is why combining helps at all.

The gain is exactly what independent averaging predicts:

| | predicted (independent) | observed |
|---|---|---|
| equal weight | 8.92×10⁻⁴ | 8.84×10⁻⁴ |
| inverse-variance | 8.88×10⁻⁴ | 8.12×10⁻⁴ |

Agreement to ~1% for equal weight. **But the reliable-detection level does not
move**: all three combinations stay at 5×10⁻³. Detection probability at 2×10⁻³
rises only from 0.63 (outbound) to 0.74 (inverse-variance), short of 0.95. The
injection grid jumps 2×10⁻³ → 5×10⁻³, so a real improvement of ~1.4× cannot
resolve into a new tested level.

Contrary to expectation, equal weighting does **not** degrade the result. The two
legs have nearly equal ε-variance (1.20 vs 1.32×10⁻³, 10% apart) despite a 2.4×
beam-power difference.

Recommended framing: outbound remains the headline observable, with the paired
inverse-variance estimator reported as a 1.4× noise-floor reduction that does not
cross to the next tested level.

## 5. Return leg — why it is not classified

Of the six pre-registered criteria for calling the return result SNR-limited:

- **Pass** — synthetic geometry/injection validation (worst scale error 0.0040,
  marginally better than outbound's 0.0048)
- **Pass** — monotonic, correct signs at ±10⁻² (detection 1.00 both)
- **Pass** — source jitter absorbed by the intercept, ε = 7×10⁻⁷
- **Weak pass** — scatter worse than outbound by only 1.10×, threshold by 1.66×,
  both well short of the 2.4× beam-power deficit. Directionally consistent,
  smaller than predicted
- **Marginal** — a small one-signed positive bias (+6.8×10⁻⁵ at zero, +2.4%
  scale error at 10⁻²) that outbound does not carry
- **Fails as written** — wrong-observable controls, see below

**The wrong-observable control is not diagnostic.** Every perturbed-trajectory
pass — displaced intercept, both slowness offsets, both aperture lengths, the
split reference — still detects at 10⁻² with probability 1.00. The mechanism is a
flaw in the control, not in the data: the injection shifts the *entire channel
trace* by −ε·T₀(s), so any coherent energy anywhere in the extraction window
carries the same position-dependent shift and is recoverable. The control asks
"is there coherent energy here", not "is the guided mode responsible".

The decisive evidence that the control rather than the analysis is at fault: it
fires identically on outbound and return, at every aperture length and both
reference constructions. A control that behaves the same on the trusted leg and
the leg under test carries no discriminating information.

The one wrong-observable control that does work is the 60–120 Hz band, which
re-derives its own trajectory where no mode was validated: reliable detection
never reached, threshold 8.8×10⁻³ (4.8× worse than primary).

**Consequence:** under a strict reading the return leg cannot be classified,
because criterion 5 is false — and the same criterion would equally condemn
outbound, which is not a defensible scientific reading. Judging criterion 5 on
the band control instead is a post-hoc reinterpretation and must be labelled
exploratory. This is an open decision, not a settled one.

A correctly specified version would inject only inside a time gate around the
trajectory, so the perturbation attaches to the mode rather than the whole trace.
That is a re-run and is exploratory now.

## 6. Defect found and fixed during this analysis

The first complete run produced a corrupted held-out result: `n_recovered`
exceeded the 23 available bursts, and no reliable-detection level was reached at
any injection magnitude on either leg.

Cause: `_inject_population` created `np.random.default_rng(SEED)` separately for
each population. Two generators seeded identically walk the same PCG64 word
stream, so trial ids collided wherever their draw positions aligned — **1383 of
7038 ids duplicated**, about 20%, nothing like birthday chance. The truth table
was a dict keyed on id alone, so the later `allbursts` row silently overwrote the
`heldout` row and ~20% of held-out trials were scored against the wrong injected
value.

The blinded gathers and the recovery estimates were never affected; only the
join was. Re-running `summarize` on the existing data was sufficient.

Fixed three ways: per-population seeding, a group tag making ids unique by
construction, and a composite `(population, leg, band, trial_id)` join key. An
assertion now requires every injection level to contain exactly one trial per
burst — the check that would have caught this immediately.

**A stale, corrupted `deep_dvv_influence.{txt,csv}` was produced before the fix
and has been deleted.** Any figure or number matching thresholds 1.53×10⁻³
(outbound) or 9.94×10⁻³ (return) is from the corrupted run and must be discarded.

## 7. Interpretation ceiling

The recovered quantity is a fractional change in the apparent along-fibre speed
of the selected Deep guided mode. It is not formation Vp or Vs, not fault-zone
stress, pore pressure, permeability, fracture compliance, or tectonic strain, and
it carries no depth resolution.

This experiment measures the sensitivity of the observable *as constructed*.
Attribution of that sensitivity to the guided mode specifically rests on the
prior permutation nulls in `deep_tube_validation` (p = 0.002), not on this test.

## 8. Known caveats

- The return 3–15 Hz trajectory intercept selected at +0.396 s against a search
  ceiling of +0.400 s — a grid-edge hit. Secondary band only; the primary 15–30 Hz
  legs sit well inside. Frozen, so reported rather than re-searched.
- Held-out false-positive rates are granular at 1/23 = 4.3%; the observed 8.7% is
  two trials.
- With n = 23 nulls the 95% threshold is the ~22nd order statistic, so
  differences of 10–20% between estimators are within its sampling noise. The MAD
  is the stable statistic for ranking; the threshold is not.
- The all-46-burst population is trajectory-contaminated and is reported only as
  a check that the held-out threshold is not a small-N artefact. With 46 bursts
  the return leg also reaches 5×10⁻³.

## 9. Files

| File | Role |
|---|---|
| `deep_dvv_injection_recovery.py` | four blinded stages: freeze, inject, recover, summarize |
| `deep_dvv_synthetic_validation.py` | per-leg geometry and injection validation on noiseless synthetics |
| `deep_dvv_influence.py` | leave-one-aperture-out and leave-one-burst-out jackknife |
| `deep_dvv_paired_legs.py` | three paired-leg combinations against the outbound benchmark |
| `deep_dvv_frozen_trajectory.json` | the frozen observable definition |
| `deep_dvv_injection_recovery.sbatch` | `-p serc`, 96 GB, ~16 min |

Blinded gathers (~4.5 GB) are intermediates and live in
`$SCRATCH/deep_dvv_blind/`. The estimator and the injection grid are imported
from `nano_dvv_injection_recovery.py` so that no part of the Nano/Deep difference
can come from a difference in estimator code.
