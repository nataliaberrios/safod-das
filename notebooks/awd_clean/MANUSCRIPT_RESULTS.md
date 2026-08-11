# Manuscript results — wavefield and repeatability

Results §4.1–4.5, following the plan's §21 step 2. The sensitivity results
(§4.6–4.9: Nano, Deep outbound and return, paired legs, Nano–Deep comparison)
live in [`DEEP_DVV_METHODS_DRAFT.md`](DEEP_DVV_METHODS_DRAFT.md) and
[`DEEP_DVV_STATUS.md`](DEEP_DVV_STATUS.md); methods are in
[`MANUSCRIPT_METHODS.md`](MANUSCRIPT_METHODS.md).

Every number is transcribed from a generated output, named inline.

---

## 4.1 Installation dependence of the observable wavefield

The same repeated surface source, recorded on two borehole DAS installations in
the same hole, yields two different coherent modes.

| | Nano | Deep |
|---|---|---|
| Installation | cemented, shallow | wireline, reversing at channel 1702 |
| Coherent mode apparent speed | ~2950 m s⁻¹ | ~1547 m s⁻¹ |
| Band where strongest | 30–60 Hz | 15–30 Hz |
| Aperture | 927 m | 3475 m outbound, 3059 m return |

The two installations are therefore not interchangeable observations of one
wavefield. This is the paper's organising result: **installation controls which
coherent mode is observable**, and — as the sensitivity results show — how useful
that mode is for time-lapse monitoring.

## 4.2 The Nano coherent fast apparent mode

From `nano_mode_identification.py`, a burst-bootstrap signed-slowness semblance
scan over 46 bursts and 859 common drops, with 999 resamples:

| Band (Hz) | Apparent speed (m s⁻¹) | Bootstrap 95% CI | P(positive ordering) |
|---|---:|---|---:|
| 15–25 | 3300 | [3250, 3350] | 1.000 |
| 25–35 | 2900 | [2900, 2900] | 1.000 |
| 35–45 | 2950 | [2950, 2950] | 1.000 |
| 45–60 | 2950 | [2950, 2950] | 1.000 |
| 60–80 | 2950 | [2950, 2950] | 1.000 |

Positive slowness ordering — energy propagating away from the wellhead — is
preferred in every band with bootstrap probability 1.000.

**Frequency-dependent slowness is resolved**, with trend 0.483
[0.414, 0.549] µs m⁻¹ Hz⁻¹. The mode is therefore dispersive, and its apparent
speed varies by ~12% across the scanned bands.

This matters for interpretation and should be stated rather than buried. Resolved
dispersion is compatible with, but not unique to, guided propagation; unresolved
dispersion would have been necessary but not sufficient evidence for direct
body-wave propagation. **This test cannot separate direct P energy from a weakly
dispersive guided or coupled mode.** The Nano apparent speed is consequently
reported as a property of a coherent fast apparent mode and not as a formation
V_P.

## 4.3 The Deep slow guided mode

Trajectories were selected independently on each leg using only the 23
discovery-half bursts (`deep_dvv_frozen_trajectory.json`):

| Leg | Apparent speed (m s⁻¹) | Intercept (s) | Discovery semblance |
|---|---:|---:|---:|
| Outbound | 1544.6 | +0.100 | 0.330 |
| Return | 1549.7 | +0.346 | 0.252 |

**The two legs agree to 0.3% in apparent speed.** They are separate propagation
branches of one fibre selected by independent semblance searches, so the
agreement is a non-trivial consistency check on the mode rather than a fitted
constraint.

The intercepts differ by 0.246 s, as expected: the return limb is reached later
because the mode has travelled down and back.

The frozen global trajectories are slower than the per-aperture local picks
reported by `deep_tube_candidates.csv`, which span 1389–1562 m s⁻¹ (slownesses
0.64–0.72 ms m⁻¹). A single straight trajectory fitted over the whole leg is
therefore not identical to the local slowness at any one aperture — evidence of
curvature or dispersion along the leg. The manuscript should quote the frozen
global values (1544.6, 1549.7 m s⁻¹) for the observable and cite the local range
separately, rather than presenting the range as uncertainty on the global value.

## 4.4 Independent validation of the Deep mode

From `deep_tube_validation.py`: trajectories selected on the 23 odd-indexed
discovery bursts, evaluated on the 23 even-indexed validation bursts, against 499
channel-order permutations (seed 20260802).

| Leg | Band (Hz) | Validation power, median of six candidates | Permutation p |
|---|---|---:|---:|
| Outbound | 3–15 | 0.235 | 0.002 |
| Outbound | 15–30 | 0.218 | 0.002 |
| Return | 3–15 | 0.032 | 0.002 |
| Return | 15–30 | 0.061 | 0.002 |

`p = 0.002` is the floor for 499 permutations: **zero permutations reached the
observed beam power in any of the four tests.**

Directionality is unanimous. In **100% of bursts**, in all four leg/band tests,
positive-slowness beam power exceeds negative-slowness power. At 15–30 Hz the
positive-to-negative power ratio is 46.9 outbound and 13.1 return for the
highest-ranked candidate.

**A precision note that will otherwise cause trouble.** Two different
"validation beam power" numbers circulate for the same test, and both are
correct:

| Quantity | Outbound 15–30 | Return 15–30 | Ratio |
|---|---:|---:|---:|
| Highest-ranked candidate (`deep_tube_candidates.csv`) | 0.353 | 0.147 | 2.4 |
| Median of six candidates (`deep_tube_validation.txt`) | 0.218 | 0.061 | 3.5 |

The pre-registration's leg designation cites the **rank-1** pair, 0.353 against
0.147, a factor of 2.4. Whichever is quoted, it must be labelled, because the
implied outbound-to-return contrast differs by nearly 50% between them.

## 4.5 Source repeatability, from individual drops to burst stacks

From `nano_hierarchical_repeatability.py`, on a fixed phase-neutral 30–60 Hz
moveout beam. Population is 988 Nano-available drops across all 49 bursts —
deliberately not the 46-burst common set, since Nano source repeatability does
not require Deep coverage.

| Metric (16 / 50 / 84%) | 16% | 50% | 84% |
|---|---:|---:|---:|
| Within-burst drop NCC, signal window | 0.551 | **0.889** | 0.989 |
| Within-burst drop NCC, noise window | 0.172 | **0.288** | 0.435 |
| Within-burst drop delay (ms) | −0.572 | +0.009 | +0.588 |
| Within-burst relative amplitude | 0.882 | 1.000 | 1.335 |
| Individual-drop beam SNR (dB) | 1.26 | 7.72 | 24.90 |
| Across-burst template NCC | 0.952 | **0.976** | 0.993 |
| Across-burst template delay (ms) | −0.254 | **−0.048** | +0.311 |

The coherent mode is genuinely present in individual drops: median signal-window
NCC of 0.889 against 0.288 in an equal-duration noise window, and **all 49 of 49
bursts** have median signal NCC exceeding median noise NCC (exact one-sided
sign test, p = 1.78×10⁻¹⁵).

Individual drops nonetheless vary substantially — beam SNR spans 1.26 to
24.90 dB between the 16th and 84th percentiles, a factor of ~20 in power, and
relative amplitude varies by ±25%.

Stacking resolves this. Independent within-burst substacks converge with drop
count:

| Drops per substack | Median NCC to independent reference |
|---:|---:|
| 1 | 0.740 |
| 2 | 0.782 |
| 4 | 0.850 |
| 8 | 0.909 |

and full burst stacks reach 0.976 against a leave-one-burst-out reference. **A
variable impact source becomes a highly repeatable observable through stacking**,
which is the precondition for the sensitivity results that follow.

Two delay quantities must be kept distinct: the median *within-burst drop* delay
is +0.009 ms, whereas the median *across-burst template* delay is −0.048 ms.

---

## Cross-references for the remaining Results sections

| Section | Location |
|---|---|
| 4.6 Nano apparent-velocity sensitivity | `DEEP_DVV_STATUS.md` §2 |
| 4.7 Deep outbound and return sensitivity | `DEEP_DVV_METHODS_DRAFT.md` Results |
| 4.8 Paired-leg precision | `DEEP_DVV_STATUS.md` §4 |
| 4.9 Nano–Deep comparison and the lever arm | `DEEP_DVV_STATUS.md` §3 |
| Small-signal tidal benchmark | `DEEP_DVV_STATUS.md` §5 — placement undecided |
