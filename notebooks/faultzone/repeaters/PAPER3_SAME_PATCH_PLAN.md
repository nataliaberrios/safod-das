# SAFOD Repeating-Earthquake DAS Project Plan

## Working paper concept

### Working title

**Testing the Same-Patch Assumption of Parkfield Repeating Earthquakes with Downhole DAS at SAFOD**

### Central scientific question

Repeating-earthquake families are commonly interpreted as earthquakes that repeatedly rupture approximately the same fault patch. Their recurrence intervals and estimated slip are then used to infer fault creep or loading rate.

This project asks:

> Do established Parkfield repeating-earthquake families produce downhole DAS wavefields that are consistent with repeated rupture of the same source patch?

The principal novelty is not simply detecting repeating earthquakes with DAS. The goal is to use a dense, continuous, downhole sensing aperture to test whether nominal repeaters have the same depth-dependent waveform, moveout, and differential-arrival pattern.

---

# Predefined publication outcomes

The project should be designed so that it produces a rigorous result regardless of whether the same-patch assumption is supported or challenged.

## Outcome 1: Same-patch assumption supported

Downhole DAS constrains members of established repeater families to source regions that overlap within a small fraction of their estimated rupture dimensions.

Possible conclusion:

> The dense downhole DAS aperture supports the same-patch assumption used in recurrence-based estimates of fault creep.

## Outcome 2: Same-patch assumption challenged

Some nominal repeater families exhibit resolvable differential moveout, source separation, or migration relative to their estimated rupture dimensions.

Possible conclusion:

> Some catalog-defined repeating earthquakes do not repeatedly rupture indistinguishable source regions, complicating the conventional same-patch interpretation.

## Outcome 3: DAS observability or methodology result

The DAS data distinguish known family members from matched controls and provide a new family-discrimination or relative-location observable, even if individual source offsets remain below resolution.

Possible conclusion:

> Downhole DAS provides a spatially distributed measure of repeater similarity and relative source location that complements conventional station-based waveform criteria.

The project is successful if it reaches any one of these defensible outcomes.

---

# Hierarchy of scientific results

## Level 1: Detection

Determine whether known HRSN repeater events can be detected coherently on the SAFOD DAS fibers.

Evidence should include:

- coherent energy across fiber position;
- event-specific moveout;
- improved visibility after moveout correction and stacking;
- detections exceeding time-shifted, randomized, and non-event controls;
- consistency across independent subarrays.

This establishes feasibility but is not the final flagship result.

## Level 2: Family-specific waveform similarity

Determine whether members of the same established family have more similar DAS wavefields than unrelated nearby earthquakes.

Compare:

- within-family event pairs;
- nearby but unrelated event pairs;
- similar-magnitude controls;
- randomized event pairs.

Evaluate similarity using:

- conventional seismic stations;
- DAS beam traces;
- channel-resolved or subarray-resolved DAS waveforms.

## Level 3: Same-patch test

After removing a common event-pair timing shift, determine whether nominal repeaters exhibit residual delays that vary systematically along the DAS aperture.

A conceptual model is

    dt_i(s) = a_i + dt_source,i(s) + dt_medium,i(s),

where:

- `s` is fiber position;
- `a_i` is a common origin-time or pick offset;
- `dt_source,i` is the spatial delay pattern caused by source-location differences;
- `dt_medium,i` is the spatial delay pattern caused by propagation changes.

For an approximately same-patch pair, the residual differential moveout should be small after removing the common timing shift.

## Level 4: Relative source separation

Use the DAS differential-delay pattern to determine whether family members are:

- colocated within uncertainty;
- vertically or laterally separated;
- separated by a meaningful fraction of their rupture dimension.

## Level 5: Implications for recurrence-based creep estimates

Compare inferred source overlap with:

- family recurrence intervals;
- recurrence regularity;
- estimated slip per event;
- catalog-based creep or slip-rate estimates;
- estimated rupture dimensions.

The final claim should be that DAS tests an assumption underlying recurrence-based creep estimates, not that DAS directly measures creep.

---

# Phase 1: Build the event and family inventory

## Objective

Construct a complete, auditable table of candidate repeater events and matched control events.

## Required information

For every candidate event, record:

- event identifier;
- origin time;
- latitude, longitude, and catalog depth;
- magnitude or seismic moment;
- established family identifier;
- recurrence interval;
- conventional-station waveform availability;
- distance from SAFOD;
- predicted P- and S-arrival times;
- DAS availability;
- operating fiber or fibers;
- duration and completeness of the DAS record;
- known timing or instrumentation problems;
- preliminary event-quality flag.

## Event-priority criteria

Prioritize events that are:

1. assigned to established Parkfield repeater families;
2. among the largest events in those families;
3. relatively close to SAFOD;
4. recorded during high-quality DAS intervals;
5. accompanied by conventional seismic recordings;
6. distributed across multiple recurrence intervals or families.

## Initial pilot sample

Begin with approximately:

- 10-20 favorable repeater events;
- 3-5 established families;
- at least one matched unrelated control per repeater event;
- both Nano and Deep DAS where available.

Do not begin with every event in the full archive.

## Phase 1 deliverables

- master event table;
- family membership table;
- DAS availability table;
- matched-control table;
- ranked pilot-event list;
- map of SAFOD, event locations, and conventional stations.

---

# Phase 2: Create a standardized event-review workflow

## Objective

Create one transparent processing notebook or script that produces the same diagnostic page for every event.

## 2.1 Raw data review

For each event, display:

- raw Nano record section;
- raw Deep record section where available;
- predicted P and S arrivals;
- conventional-station waveforms;
- data gaps and bad channels;
- sufficiently long pre-event and post-event windows.

Questions:

- Is the event visibly present?
- Is there coherent moveout?
- Can P and S arrivals be separated?
- Are unrelated transients present?
- Is the timing believable?
- Which parts of the fiber are usable?

## 2.2 Basic preprocessing

Freeze a modest preprocessing sequence using a discovery subset.

Possible steps:

- demean;
- detrend;
- remove bad channels;
- apply a fixed band-pass filter;
- taper;
- retain unnormalized data for quantitative amplitude work;
- use trace normalization only for display when necessary.

Avoid selecting a different filter for every event.

## 2.3 Event-specific moveout

Earthquakes arrive from different locations and directions, so one fixed beam velocity will probably not describe every event.

For each event and phase:

- estimate or predict moveout;
- define a P-wave window;
- define an S-wave window where possible;
- align the selected phase across channels;
- form beam traces for multiple independent subarrays;
- save moveout-fit residuals.

Moveout may be obtained from:

- a velocity model;
- a conventional hypocenter;
- a slowness scan on a discovery aperture;
- an empirical reference event.

The method must be documented and frozen before family comparisons.

## 2.4 Event-quality metrics

Save:

- beam signal-to-noise ratio;
- normalized cross-correlation to a reference event;
- number of usable channels;
- moveout-fit residual;
- timing uncertainty;
- coherence across independent subarrays;
- P- and S-phase quality flags;
- Nano and Deep availability.

## Phase 2 deliverables

- one standardized event-review page per event;
- one event-quality table;
- frozen preprocessing parameters;
- frozen phase windows or window-selection rules;
- frozen moveout-estimation procedure.

---

# Phase 3: Demonstrate that DAS recognizes established families

## Objective

Show that known repeater-family members are more similar than carefully matched unrelated events.

## 3.1 Conventional-station similarity

Use one or more conventional stations to verify that:

- established family members have high waveform similarity;
- unrelated controls generally have lower similarity;
- the selected family assignments are reproduced by standard measurements.

This provides the conventional baseline.

## 3.2 DAS beam similarity

For each event pair:

- moveout-correct the selected phase;
- create one beam trace per subarray;
- measure waveform similarity and relative delay;
- compare family pairs with matched controls.

## 3.3 Channel-resolved or local-aperture similarity

Divide the fiber into local apertures.

For each event pair and aperture, calculate:

- maximum normalized cross-correlation;
- differential delay;
- coherence;
- cross-spectral phase;
- delay uncertainty;
- number of usable channels.

Summarize similarity as a function of fiber position.

A convincing family pair should show a coherent pattern over a substantial portion of the aperture, not merely one high beam correlation.

## Primary Phase 3 result

Create the first go/no-go plot:

    within-family DAS similarity   versus   matched-control DAS similarity

If these populations separate, the project has a strong empirical foundation.

## Phase 3 deliverables

- family-versus-control similarity distributions;
- example strong family pair;
- example matched control pair;
- channel-resolved similarity profiles;
- statistical test of within-family versus control separation.

---

# Phase 4: Measure differential delays across the DAS aperture

## Objective

Determine whether nominal repeaters have indistinguishable or systematically different spatial arrival-time patterns.

## Workflow

For each within-family event pair:

1. align the events approximately using conventional stations or DAS beam traces;
2. divide the DAS aperture into local channel groups;
3. measure the relative delay in each local aperture;
4. estimate uncertainty for every delay;
5. remove a common event-pair timing offset;
6. inspect the residual delay as a function of fiber position;
7. repeat on independent subarrays and, where possible, for both P and S phases.

The central observable is

    dt_res(s) = dt(s) - mean(dt),

where the common event-pair timing shift has been removed.

## Candidate delay estimators

Compare a small number of methods on the pilot dataset:

### Time-domain cross-correlation

Advantages:

- simple;
- interpretable;
- familiar;
- useful for initial quality control.

### Cross-spectral phase

Advantages:

- allows sub-sample delay estimates;
- provides frequency-dependent coherence;
- can estimate uncertainty;
- may perform better for highly similar repeater waveforms.

### Multitaper cross-spectral timing

Advantages:

- improved spectral stability;
- useful uncertainty estimates;
- potentially more robust for short windows.

Select the final method using the discovery subset and evaluate it on held-out families or event pairs.

## Phase 4 deliverables

- residual delay profiles for the strongest family pairs;
- matched-control residual profiles;
- uncertainty estimates;
- independent-subarray comparisons;
- comparison of timing estimators;
- frozen final delay estimator.

---

# Phase 5: Establish null distributions and robustness

## Objective

Demonstrate that the apparent same-patch or source-separation signal is not produced by arbitrary waveform similarity, data quality, or spatial overfitting.

## Null A: Nearby unrelated earthquakes

Match controls by:

- hypocentral distance;
- magnitude;
- signal-to-noise ratio;
- recording quality;
- phase availability.

Test whether within-family residual moveout is smaller than that of nearby unrelated pairs.

## Null B: Family-label permutation

Randomly shuffle family labels while preserving appropriate matching constraints.

Test whether the observed within-family similarity and residual-delay coherence exceed random grouping.

## Null C: Spatial permutation

Randomize or disrupt channel order.

Test whether the result depends on the true spatial organization of the DAS array.

## Null D: Time-shifted windows

Repeat the analysis in pre-event or unrelated time windows.

Test the false-detection rate.

## Null E: Independent subarrays

Estimate a spatial delay pattern using one subarray and test it on another.

This reduces overfitting.

## Null F: Phase consistency

Repeat independently for P and S where both are usable.

The phases need not have identical sensitivity, but a source-offset inference should not be contradicted by another well-recorded phase.

## Null G: Reference-event sensitivity

Repeat key measurements using alternative reference members within a family.

The conclusion should not depend on a single chosen reference event.

## Phase 5 deliverables

- null distributions;
- false-positive rates;
- permutation-test results;
- sensitivity to reference choice;
- independent-subarray validation;
- documented criteria for accepting a family pair.

---

# Phase 6: Relative source-location inversion

## Objective

Convert spatially varying differential delays into constraints on relative source position.

Begin only after the empirical delay measurements pass the Phase 5 tests.

## 6.1 Linearized relative-location model

For small source perturbations,

    dt_j ~= (dt_j/dx_s) dx_s + (dt_j/dy_s) dy_s + (dt_j/dz_s) dz_s + dt_0,

where:

- `j` indexes DAS apertures or channels;
- `dx_s, dy_s, dz_s` are relative source offsets;
- `dt_0` is a common origin-time correction.

Use a velocity model to calculate travel-time derivatives or empirical sensitivity kernels.

Estimate:

- relative horizontal offset;
- relative vertical offset;
- origin-time correction;
- posterior covariance or bootstrap uncertainty;
- trade-offs among source-coordinate components.

## 6.2 Combined DAS and conventional-station inversion

Where appropriate, combine:

- DAS differential delays;
- HRSN or other conventional-station differential times;
- catalog constraints.

The DAS aperture may add strong spatial sampling but may not independently constrain every three-dimensional source component.

## 6.3 Escalation to more advanced modeling

Only pursue more advanced modeling if the linearized approach reveals stable spatial information.

Possible later methods:

- three-dimensional ray tracing;
- double-difference relocation;
- empirical calibration with nearby catalog events;
- waveform-based relative location;
- source-specific forward modeling.

Do not begin with full waveform inversion.

## Phase 6 deliverables

- inferred relative offsets;
- uncertainty ellipses or credible regions;
- resolution tests;
- synthetic injection-recovery;
- comparison with catalog relative locations;
- accepted and rejected event-pair table.

---

# Phase 7: Compare source separation with rupture dimension

## Objective

Test whether inferred source separation is small relative to the physical size of the repeating-earthquake rupture.

Estimate rupture dimension from seismic moment or magnitude using clearly stated assumptions.

Define a normalized separation:

    R_sep = (inferred source separation) / (estimated rupture radius)

Interpretation:

- `R_sep << 1`: strongly consistent with repeated rupture of approximately the same patch;
- `R_sep ~ 1`: possible partial overlap or source migration;
- `R_sep > 1`: difficult to reconcile with the same-patch interpretation unless uncertainty is large.

Report uncertainty in both source separation and rupture size.

## Phase 7 deliverables

- normalized source-separation estimates;
- family-by-family overlap classification;
- uncertainty propagation;
- comparison with conventional repeater assumptions.

---

# Phase 8: Connect source overlap to recurrence and creep interpretation

## Objective

Determine whether the degree of inferred source overlap affects the interpretation of repeater recurrence.

Compare source-overlap metrics with:

- recurrence interval;
- recurrence regularity;
- magnitude stability;
- waveform similarity;
- inferred slip per event;
- recurrence-based slip-rate estimates;
- catalog-based local creep estimates.

Possible outcomes include:

- families consistent with the same patch have the most stable recurrence;
- some nominal families migrate while retaining high conventional waveform similarity;
- recurrence-based slip rates remain robust despite small source offsets;
- source migration introduces meaningful uncertainty into creep estimates.

The final wording should remain cautious:

> Downhole DAS tests the source-overlap assumption underlying recurrence-based estimates of fault creep.

Do not claim that DAS directly measures creep unless an independent measurement supports that statement.

---

# Pilot study design

## Pilot objective

Establish whether the DAS data distinguish known repeater-family pairs from matched controls and whether stable spatial differential delays can be measured.

## Pilot dataset

Use:

- 3-5 established families;
- 10-20 favorable repeater events;
- matched unrelated events;
- conventional-station data;
- Nano DAS;
- Deep DAS where available.

## Pilot milestones

### Milestone 1: Event inventory complete

All pilot events have:

- verified timing;
- DAS availability;
- conventional waveforms;
- quality metadata;
- family and control labels.

### Milestone 2: Standardized event pages complete

Every event has the same raw and processed diagnostic products.

### Milestone 3: Family-versus-control plot complete

Within-family DAS similarity is compared with matched controls.

### Milestone 4: Strongest family-pair delay profile complete

A channel- or aperture-resolved differential-delay profile is produced with uncertainty.

### Milestone 5: Independent validation complete

The strongest result survives:

- independent subarrays;
- alternative reference events;
- matched controls;
- spatial and label permutations.

### Milestone 6: Relative-location feasibility established

Synthetic or empirical tests determine whether plausible source offsets are resolvable.

---

# Decision gates

## Gate 1: Event detectability

Proceed when known repeater events produce coherent DAS detections above null windows.

If detections are weak, optimize only a bounded set of:

- frequency bands;
- apertures;
- phase windows;
- stacking strategies.

## Gate 2: Family discrimination

Proceed to relative timing when within-family similarity exceeds matched controls.

If the populations overlap, investigate whether the limitation is:

- timing;
- signal-to-noise ratio;
- moveout mismatch;
- phase mixing;
- insufficient aperture;
- incorrect family assignment.

## Gate 3: Stable differential delays

Proceed to source inversion when delay profiles are:

- coherent across neighboring apertures;
- reproducible across subarrays;
- stable across reference events;
- distinguishable from controls.

## Gate 4: Source-location resolution

Proceed to recurrence and creep interpretation when inferred offsets are resolved relative to uncertainty and rupture dimension.

If offsets remain unresolved, the paper can still report upper bounds on source separation.

---

# Figure plan

## Figure 1: Geometry and dataset

Include:

- SAFOD location;
- repeater and control-event locations;
- Nano and Deep fiber geometry;
- conventional stations;
- family inventory;
- example event records.

## Figure 2: Example established family

Include:

- conventional waveforms;
- raw DAS record sections;
- moveout-corrected DAS beams;
- local-aperture waveform similarity;
- example P and S phases where available.

## Figure 3: Families versus controls

Include:

- within-family similarity distribution;
- matched-control distribution;
- family-label permutation test;
- depth-coherent similarity metric.

## Figure 4: Differential moveout

Include:

- residual delay versus fiber position for a strong family pair;
- matched unrelated pair;
- uncertainty;
- independent-subarray validation;
- reference-event sensitivity.

## Figure 5: Relative source offsets

Include:

- inferred relative source positions;
- uncertainty regions;
- comparison with catalog locations;
- synthetic resolution test;
- normalized separation relative to rupture radius.

## Figure 6: Implications for repeater interpretation

Include:

- source overlap by family;
- overlap versus recurrence regularity;
- overlap versus recurrence-based slip rate;
- supported, ambiguous, and challenged same-patch classifications.

A concise and broadly important result may be suitable for a five-figure GRL manuscript. A fuller methods, validation, and multi-family treatment is more naturally suited to JGR: Solid Earth.

---

# Journal strategy

## Primary target: JGR: Solid Earth

Most appropriate when the paper includes:

- development of the DAS differential-timing method;
- rigorous null tests;
- source-offset inversion;
- several repeater families;
- implications for fault creep or repeater mechanics.

## Aspirational target: Geophysical Research Letters

Appropriate only if the result is exceptionally clear and broadly important, such as:

> Downhole DAS demonstrates that nominal Parkfield repeaters migrate by a significant fraction of their rupture dimensions.

or

> Downhole DAS tightly constrains Parkfield repeaters to repeatedly rupture the same fault patch.

## Alternative target: Bulletin of the Seismological Society of America

Appropriate if the main contribution becomes:

- relative earthquake location with downhole DAS;
- waveform processing;
- array methodology;
- repeater-family discrimination.

## Alternative target: Seismological Research Letters

Appropriate if the strongest result is:

- first characterization of Parkfield repeaters on the SAFOD DAS arrays;
- dataset and observability;
- limited relative-source resolution.

---

# Reproducibility and anti-overfitting rules

## Separate discovery and validation data

Use different:

- events;
- event pairs;
- families;
- channels;
- or subarrays

for parameter selection and final validation.

## Freeze processing decisions

Before final testing, freeze:

- frequency bands;
- event windows;
- phase windows;
- apertures;
- moveout estimation;
- correlation method;
- quality thresholds;
- reference construction;
- null definitions.

## Maintain a decision log

For every methodological decision, record:

- the choice;
- the reason;
- the data used to make it;
- the data reserved for validation;
- rejected alternatives;
- date of the decision.

## Preserve negative tests

A failed analysis should be saved with:

- the exact hypothesis;
- method;
- result;
- explanation of why it was rejected.

This prevents repeated dead ends and strengthens the eventual methods section.

## Require uncertainty estimates

Do not interpret source separation without:

- delay uncertainty;
- model uncertainty;
- velocity-model sensitivity;
- bootstrap or posterior uncertainty;
- resolution tests.

---

# Weekly progress rule

Every week of work must end with at least one paper-producing artifact:

- a final-quality figure panel;
- a completed methods subsection;
- a frozen processing function;
- a validated event or family table;
- a documented null test;
- a completed synthetic resolution test;
- a concise decision memo.

Generating many exploratory plots without freezing a result does not count as progress toward the paper.

---

# Project guardrails

The project will be pursued aggressively, but the desired conclusion will not be predetermined.

The goal is to force the archive through a rigorous analysis that yields one of three publishable conclusions:

1. the same-patch assumption is supported;
2. the same-patch assumption is challenged;
3. downhole DAS provides a new, quantitatively validated repeater-observability or relative-location method.

The project should not depend on proving migration or proving perfect source overlap.

---

# Immediate next deliverable

Create a pilot dataset containing:

- the 3-5 most favorable established repeater families;
- 10-20 high-quality events;
- matched unrelated controls;
- conventional waveforms;
- Nano and Deep DAS availability;
- one standardized review page for every event.

The first decisive scientific plot is:

    within-family DAS similarity   versus   matched-control DAS similarity

The second decisive plot is:

    channel-resolved residual delay   for the strongest family pair

These two results determine the next stage of the project.

---

---

# APPENDIX (added 2026-08-06): what already exists against this plan

Written by Claude when filing this plan into the repo, so that work already done is
not repeated. Numbers here are cross-referenced to `METHODS_STATUS.md`.

## Already satisfied, at least in pilot form

| Plan item | Status | Evidence |
|---|---|---|
| L1 detection | **done** | G1: 785/900 channels, coda SNR > 3 |
| L2 family vs control similarity, conventional stations (3.1) | **done** | HRSN: 10 pairs, CC 0.786-0.995, vs control median -0.004 |
| L2 on DAS beams (3.2) | **done** | `moveout_test.csv`; aligned DAS CC 0.956 vs HRSN 0.954 (parity) |
| Phase 3 primary go/no-go plot | **done, and it separates** | 500 Hz cache, 5-20 Hz: repeaters 0.742 vs null -0.142, null MAD 0.056, d' = 15.8; best configuration d' = 25.3 |
| 2.3 event-specific moveout | **done** | `moveout_test.slant_scan`, per-event slowness scan including p = 0 |
| Phase 4 cross-spectral phase estimator | **built** | `dvv_core.sub_sample_delay` (Poupinet et al. 1984) |
| Phase 6 synthetic injection-recovery | **built and passing at 0.1 %** | `interval_dvv_gate.py` G2 |
| Null D (time-shifted / acausal) | **built** | G3 in the gate |

**Gate 1 and Gate 2 of this plan are therefore already passed.** The project is at
Gate 3 (stable differential delays), not at the beginning.

## Feasibility already computed for Phases 4, 6 and 7

With 500 Hz extraction plus moveout-aligned stacking and a gauge-length block
bootstrap, for a source offset of one rupture radius:

    SNR                        3.3 - 52.2 (median 35.5)
    3-sigma resolvable offset  1.2 - 13.7 m (median 1.9 m)
    rupture radii              14.4 - 48.4 m

So `R_sep` in Phase 7 is resolvable well below 1 for most pairs. Incidence angles are
13.8-65.3 deg, and sensitivity scales as tan(i), so the near-vertical pairs are the
*weakest*, not the strongest -- the opposite of the intuition.

## The one confound this plan should name explicitly

`dt_source` and `dt_medium` in the Level 3 model are **exactly degenerate per pair**:
both enter as a fractional change in apparent slowness. A velocity change of dv/v is
indistinguishable from an offset of `delta = (dv/v) * R / tan(i)`.

They separate only at the population level, because a medium change is common to
family and control pairs while a source offset is not. This is why Phase 3's matched
controls are load-bearing for Phase 4 and cannot be skipped, and why the statistic
must be the family-vs-control contrast rather than any single pair's number.

## Known gaps against this plan

- **No established family labels.** The 10 pairs are HRSN-similarity-confirmed, not
  drawn from a published Parkfield repeater catalog. Phase 1 is genuinely not done.
  Gao, Kao & Wang 2021 (`10.1029/2021gl092815`) show CC alone is not a valid
  repeater criterion, so this gap is the plan's real starting point.
- **Deep DAS unusable for these events.** The deep fiber started 2026-03-28; the
  repeater pairs are 2024-05 to 2025-07. Deep DAS applies only to future events.
- **Nulls B, C, E, G not built** (label permutation, spatial permutation,
  independent subarrays, reference-event sensitivity).
- **S-phase analysis not done** (Null F).
- **No multitaper estimator** — only time-domain CC and cross-spectral phase.
