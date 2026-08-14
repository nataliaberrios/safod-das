# Lellouch Figure 7c on the 2024–2025 archive — authoritative status

Updated 2026-08-14. This document supersedes every earlier Figure 7c verdict in
this repository.

## Question and answer

**Question.** Does the unfiltered one-day ambient-interferometry calculation
reported for Figure 7c of Lellouch et al. (2019) recover an ordered,
approximately 3.2 km/s packet on a complete day from the 2024–2025 SAFOD
main-hole DAS archive?

**Answer for 20 December 2024: no.** The corrected baseline does not contain a
statistically significant receiver-ordered 3.2 km/s moveout. Its scan maximum is
at 5.85 km/s, near the flat-moveout edge of the search, and is below the 95th
percentile of a receiver-order scan-max null (score 6.13 versus 6.32;
familywise p = 0.147). At 3.2 km/s the causal and acausal scores are comparable
(2.75 and 2.83). Most decisively, permuting receiver channels before
preprocessing leaves the gather almost unchanged (flattened waveform
correlation 0.9976).

This is a matched-day negative reproduction. It is not evidence that ambient
interferometry fails for the full archive, and it is not evidence that F–K
filtering is invalid.

## Exact input and baseline

| Quantity | Value |
|---|---|
| Archive index | `/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv` |
| UTC day | 2024-12-20 |
| Raw HDF5 files | 1,440 one-minute records |
| HDF5 dataset | `/Acquisition/Raw[0]/RawData` |
| Per-file shape | 30,000 time samples × 900 channels |
| Sampling | 500 Hz; 1.020952344 m channel spacing; 16.335238 m gauge length |
| Virtual source | channel 23, provisional G0 wellhead estimate |
| Receiver centers | approximately 50–700 m offset in 50 m increments |
| Nearby-receiver operator | same source correlated with R−10 through R+10; literal sum |
| Temporal preprocessing | continuous time derivative; centered 0.1 s running absolute mean |
| Correlation segmentation | 30 s windows every 15 s, including file-boundary starts |
| Total windows | exactly 5,759 unique windows |
| Figure 7c stacking | ordinary unshifted cross-spectrum average and R±10 sum |
| Final display filter | 5–20 Hz on the full centered correlation, then crop to ±0.35 s |
| Explicitly absent | F–K filter, detrend, input bandpass, common-mode subtraction, whitening, per-window correlation normalization, imported velocity alignment |

The 0.1 s RAM duration follows the Bensen et al. (2007) half-maximum-period
recommendation for a 5–20 Hz band; Lellouch et al. do not report their RAM
duration. A 5 s sensitivity is therefore retained.

## Full-day branch results

Every branch contains 1,440 files and 5,759 windows. The p value is familywise
over the predeclared 1.5–6.0 km/s causal velocity scan using 10,000
receiver-order permutations.

| Branch | Best velocity (m/s) | Score at 3,200 | Acausal at 3,200 | Null 95% scan max | p | Flattened correlation to baseline |
|---|---:|---:|---:|---:|---:|---:|
| Paper baseline | 5,850 | 2.752 | 2.831 | 6.324 | 0.1470 | 1.0000 |
| Source channel 0 | 5,925 | 0.959 | 1.152 | 1.858 | 0.7435 | −0.2853 |
| RAM 5 s | 5,850 | 2.747 | 2.857 | 6.317 | 0.1083 | 0.9999 |
| Common mode, median of all 900 channels | 3,650 | 1.000 | 1.046 | 1.115 | 0.9220 | 0.2209 |
| Stabilized Equation 6 sensitivity | 5,850 | 1.357 | 1.371 | 1.510 | 0.3663 | 0.9723 |
| IID broadband white noise | 4,550 | 1.025 | 0.943 | 1.377 | 0.3371 | −0.0739 |
| Measured receivers permuted before preprocessing | 5,850 | 2.698 | 2.490 | 6.068 | 0.9010 | 0.9976 |

No branch clears α = 0.05. The RAM sensitivity shows that 0.1 versus 5 s does
not control the baseline result. Common-mode subtraction removes the dominant
shared structure but does not reveal the target. The declared Equation 6
sensitivity divides the average cross spectrum by the average source-power
spectrum after all windows are combined, with a 10⁻³ water level referenced to
5–20 Hz; it also does not reveal the target.

## Evaluation of the three disputed claims

1. **Common-mode removal is not reported for ambient Figure 7c.** It must not be
   inserted into the paper baseline. Testing it separately was useful because it
   confirms that the uncorrected gather is dominated by a receiver-order-
   insensitive component, but the corrected gather remains nonsignificant.
2. **Equation 6 is present as a theoretical proportionality, but its estimator
   is underspecified.** The post-average, water-level-stabilized implementation
   is one declared sensitivity, not “the exact hidden step.” Its negative result
   does not prove that every possible source-spectrum estimate is equivalent.
3. **The R±10 sum is explicitly reported and was missing from the legacy
   single-pair test.** It is now implemented literally. It is necessary for
   fidelity but insufficient to recover Figure 7c on this day. Twenty-one
   neighboring correlations imply at most an ideal amplitude-SNR gain of
   sqrt(21) ≈ 4.6 for independent noise, not 21×.

Thus the original opinion was partly right about the missing nearby-receiver
sum, wrong to call common-mode subtraction a published requirement, and too
specific about how Equation 6 must be implemented.

## Validation gates

The operator passes independent tests before measured-data interpretation:

- a known 500 m/s increasing-coordinate synthetic has zero median lag error at
  the 0.01 s synthetic sampling interval;
- the summed-receiver and sum-of-correlations R±10 formulations agree to
  3.28×10⁻⁸ relative error;
- distributed window ownership yields exactly 5,759 unique full-day starts;
- the raw white-noise control is broadband and has interchannel correlation
  −0.00179; and
- full-correlation filter-then-crop output agrees exactly with an independently
  constructed reference.

The last gate was added after visual inspection exposed crop-boundary ringing in
an earlier provisional white-noise figure. All seven aggregates were then rerun
with the corrected order.

## Interpretation and next test

The matched day fails because the dominant correlation waveform is insensitive
to receiver ordering, not because the published nearby-receiver sum was omitted.
The result at 5.85 km/s should not be interpreted as a high-velocity arrival: it
lies near the flat-moveout boundary and survives receiver scrambling.

The clean next test is to freeze this exact operator and apply it to independently
chosen complete days, followed by a predeclared multi-day convergence sequence.
That tests whether 20 December 2024 is unrepresentative without retuning the
processing after seeing each result. F–K-assisted correlations remain a separate
extension; they require their own matched white-noise and pre-filter
channel-scramble controls and cannot be used as proof that the unfiltered Figure
7c calculation succeeded.

## Provenance

| Item | Identifier |
|---|---|
| Paper-faithful operator | `awd_clean/ambient_lellouch2019_exact_stack.py` |
| Hourly full-day matrix | SLURM 38988141; 168/168 tasks completed |
| First aggregate audit | SLURM 38988173 |
| Final full-filter-then-crop aggregate audit | SLURM 38993456; 7/7 tasks completed |
| Paper-operator commit | `9c70083` |
| Sensitivity-order commit | `cc07080` |
| Full-lag filter commit | `bb48942` |
| Advisor notebook | `awd_clean/Ambient_FK_QC_workflow.ipynb`, v8 |
| Result directory | `awd_clean/ambient_transfer/lellouch2019_exact_stack/` |

Primary citation: Lellouch et al. (2019), *Journal of Geophysical Research:
Solid Earth*, https://doi.org/10.1029/2019JB017533. RAM guidance: Bensen et al.
(2007), *Geophysical Journal International*,
https://doi.org/10.1111/j.1365-246X.2007.03374.x.
