## Independent unfiltered verification (running)

- **Slurm jobs:** F–K verification `37478467`; exact unfiltered verification `37478473`. Both resumed from checkpoints on compute nodes.
- **Automatic post-processing:** Slurm `37478474` depends on both analysis jobs and will run the two seasonal aggregators after successful completion.
- Exact-resolution path: `run_seasonal_unfiltered.py` using `ambient_transfer_test.py`.
- Purpose: test whether the weak unfiltered result is reproducible independently of the F–K routine and its downsampling.
- Earthquake intervals remain present; 5-s temporal normalization suppresses amplitude leverage without catalog removal.
- Progress log: `ambient_transfer/seasonal_unfiltered_progress.log`.

## Seasonal ambient-noise/F–K validation (running)

- **Started:** 2026-08-03
- **Manifest:** `ambient_transfer/seasonal_day_selection.json`
- **Design:** two reproducibly selected complete days in each meteorological season (DJF, MAM, JJA, SON), excluding the February baseline day from the random draws.
- **Completeness criterion:** at least 1,400 one-minute files and no inter-file gap greater than 65 s.
- **Processing:** each day independently; 5–20 Hz bandpass, 5 s running-absolute-mean temporal normalization, and signed negative/positive/both F–K branches. Earthquake intervals are retained and amplitude-suppressed, not catalog-removed.
- **Outputs:** resumable 10-file chunks in `ambient_transfer/`; aggregation will be run only after day-level products are complete.
- **Launcher:** `run_seasonal_fk.py`; progress log: `ambient_transfer/seasonal_fk_progress.log`.
- **Notebook status:** live seasonal section promoted as v26; independent unfiltered-verification section promoted as v27; Lellouch axis-convention correction promoted as v28; preliminary five-hour comparison promoted as v29; rendered interim comparison figure promoted as v30; Lellouch Figure 7c-style section promoted as v32; original l7c.png reference embedded as v33; direct black-and-white wiggle products promoted as v34. Scientific figures and conclusions remain pending until day-level and across-day checks pass.

# SAFOD AWD v25 validation status

- Generated: `2026-08-02T23:01:00-07:00`
- Nano hierarchical repeatability job: `37253995` (**COMPLETED**)
- Deep split-sample validation audit: `37217153` (**PASS**)
- Deep time-gated branch search: `37249098` (**COMPLETED**)
- Authoritative notebook: `AWD_results_dashboard.ipynb` (**v18**)
- Advisor PDF source: `AWD_advisor_figure_guide_v16.tex` (**updated; compilation unavailable because no TeX engine is installed**)

## Nano hierarchical repeatability

| Product | Result |
|---|---|
| Individual drops | 988 drops in 49 bursts processed |
| Within-burst signal vs noise | Median NCC 0.889 vs 0.288; 49/49 bursts signal > noise |
| Exact burst-level sign test | `p = 1.776e-15` |
| Across-burst template | Median NCC 0.976; median delay −0.048 ms |
| Independent substack convergence | Median NCC 0.740 (1 drop) to 0.909 (8 drops) |

## Deep localized validation and branch search

| Leg | Band | Fixed-p validation power | Permutation p | Bursts +p > −p |
|---|---:|---:|---:|---:|
| outbound | 3–15 Hz | 0.23474 | 0.0020 | 100% |
| outbound | 15–30 Hz | 0.21768 | 0.0020 | 100% |
| return | 3–15 Hz | 0.03151 | 0.0020 | 100% |
| return | 15–30 Hz | 0.06140 | 0.0020 | 100% |

The time-gated signed-branch search tested 120 branch rows and found five
geometric positive/negative intersections. None passed the joint spatial-order
permutation threshold (`p <= 0.05`). The result supports a repeatable slow
signed mode but does not support a fixed conversion/reflection coordinate.

## Interpretation boundary

The Deep QC and validation products establish computational integrity and a
repeatable directed spatial mode. They do not establish the tube-wave phase,
verified measured depth, a permeable fracture, or permeability. The branch
intersections remain provisional fiber coordinates. The Deep coordinate audit and bounded tube-wave forward consistency test are complete. `StartLocusIndex=1800` corresponds to 3675.42 m in interrogator coordinates and hairpin channel 1702 to 3475.31 m from channel 0; neither is surveyed depth. All 24 validated slow-mode candidates are within the 1.3–1.8 km/s fluid prior (median 1.538 km/s; p05–p95 1.389–1.563 km/s). No physical-depth registration or borehole-log comparison is claimed. The next actionable step is to obtain the deployment drawing/surveyed fiber-to-depth transform and borehole fluid, casing, and log files, then replace this envelope with a real tube-wave dispersion model.

## Targeted Deep strand-depth observability scan (v14)

- Script: `deep_target_scan.py`; Sherlock job `37309599` (**COMPLETED**)
- Products: `deep_target_scan.png/.pdf/.csv/.npz/.json/.txt`
- Scope: provisional 3000–3450 m scan; 200 m windows every 50 m; 3–15 and 15–30 Hz; positive and negative signed slowness; outbound and reversed-return legs.
- Conditional mapping: outbound `s = depth`; reversed return `s = 3475.31 − depth` m. This is not independently surveyed.
- Strand tests: positive 15–30 Hz validates at assumed SDZ 3192 m and CDZ 3302 m on both legs; all four explicit channel-order tests give `p = 0.002`.
- Controls: 3–15 Hz is weak; isolated negative-direction p-values do not form a coherent opposing branch.
- Interpretation: target-window observability under a stated registration hypothesis; not a creep rate, casing-deformation amplitude, reflection depth, or permeability result.

## Deep coordinate registration and bounded forward test

- Script: `deep_registration_forward_model.py`
- Products: `deep_registration_forward_model.png`, `.pdf`, `.csv`, `.npz`, `.json`, `.txt`
- Status: completed preliminary consistency test; not a depth or permeability inversion.
- Coordinate evidence: `dx=2.0419 m/channel`, `N=3200`, `StartLocusIndex=1800`, data-derived hairpin `channel=1702`.
- Limitation: no surveyed fiber origin/path, fluid/casing properties, or borehole logs found in the workspace.

## Source-history regression and tidal-timescale candidate

- Script: `nano_tidal_compaction_regression.py`
- Product: `nano_tidal_compaction_regression.png/.pdf/.csv/.npz/.json/.txt`
- Response: 49 burst-level leave-one-burst-out Nano delays.
- Nuisance proxies: cumulative AWD drops and burst signal RMS; no independent ground-level compaction sensor was found.
- Result: phase-free 24-hour harmonic reduces blocked-CV RMSE 0.297 to 0.223 ms; fitted amplitude 0.283 ms; extrema near 5.0 and 17.0 h; circular 3-burst bootstrap p=0.002.
- Interpretation: candidate nonmonotone tidal-timescale component, not yet a tidal-force or fault-opening/closing detection.
- Completed physical-phase control: UTC-phased low-precision Sun/Moon degree-2 potential. It is weaker than the phase-free candidate; no opening/closing sign is assigned. Next: obtain an independent compaction record, fault orientation, stress convention, and precision strain/traction predictor.

## Exploratory canonical-correlation test

- Script: `nano_tidal_cca.py`
- Product: `nano_tidal_cca.png/.pdf/.csv/.npz/.json/.txt`
- X set: phase-free 24.0 h and 12.4206 h sine/cosine bases.
- Y set: source-corrected delay and burst signal RMS.
- Full-sample canonical correlation: 0.533; row-permutation p=0.019; 3-burst block-permutation p=0.118.
- Five contiguous-fold absolute correlations: 0.697, 0.016, 0.378, 0.057, 0.075.
- Interpretation: exploratory association only; the serial null and held-out instability do not confirm a tidal signal.

## Physically phased astronomical tide-potential control

- Script: `nano_physical_tide_cca.py`
- Products: `nano_physical_tide_cca.png/.pdf/.csv/.npz/.json/.txt` and `nano_physical_tide_predictors.csv`
- Predictor: actual UTC burst times with low-precision Sun/Moon degree-2 tide-generating potential at approximate SAFOD coordinates (35.98 N, 120.55 W).
- Response: source-corrected Nano delay and burst signal RMS; cumulative drops remain a source-history proxy because no independent compaction time series was found.
- Result: full-sample canonical r=0.366; row-permutation p=0.116; three-burst block-permutation p=0.263; held-out absolute correlations 0.780, 0.247, 0.230, 0.501, 0.039.
- Interpretation: the phase-free harmonic is not reproduced as a stable externally phased astronomical response. This scalar potential is not fault-normal traction and does not assign opening/closing.

## Deep target-window burst repeatability (v15)

- Script: `deep_target_burst_repeatability.py`; Sherlock job `37312233` (**COMPLETED**)
- Held-out validation: 23 even-numbered bursts, fixed trajectories selected from discovery epochs.
- 15–30 Hz positive-direction fraction: SDZ outbound 21/23 (`p=3.3e-5`), SDZ return 20/23 (`p=2.4e-4`), CDZ outbound 21/23 (`p=3.3e-5`), CDZ return 16/23 (`p=0.047`).
- Similar outbound control coherence means this is a repeatable mode passing through the target interval, not a strand-localized casing-deformation signature.
- It remains a conditional baseline result; the 24-hour record does not measure a creep rate.

## Passive repeating-earthquake candidate catalog (v16)

- Script: `repeating/repeater_catalog_screen.py`
- Catalog: 329 events, 2024-05-01 through 2026-03-30
- Screen: 6,325 candidate pairs; 6,112 on distinct dates and at least 7 days apart
- Products: `repeating/repeater_catalog_screen.png`, `repeater_catalog_screen.csv`, `repeater_candidate_pairs_ranked.csv`, `repeater_candidate_pairs_ranked.md`, and `repeater_catalog_report.md`
- Interpretation: event selection only. Continuous DAS waveforms covering the catalog timestamps are still required for waveform/coda correlation.

## Passive repeater waveform-coverage audit (v17)

- Inventory: `/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv`
- Valid timed records: 392,079 of 392,095
- Catalog events with exact interval coverage: 146 of 329
- Candidate pairs with coverage at both event times: 1,200
- Products: `repeating/repeater_event_waveform_coverage.csv`, `repeater_candidate_pairs_ranked_with_coverage.csv`, and `repeater_waveform_coverage.md`
- Interpretation: coverage audit only; waveform extraction/correlation is the next stage.

## Covered-pair waveform pilot (v18)

- Script: `repeating/repeater_waveform_extract.py`
- Ten highest-ranked covered pairs tested; 7 had usable windows at both events.
- Median channel-wise maximum correlation: 0.139–0.177; no pilot channel exceeded 0.5 correlation.
- Interpretation: preliminary negative triage result; catalog proximity does not guarantee repeatable DAS waveforms.
- Products: `repeating/repeater_pilot_waveform_correlation.csv` and `.md`

## Ambient-noise replication audit (v19)

- Completed a diagnosis of the earlier five-month Lellouch-reproduction attempt.
- The prior branch used a different 2024--2025 acquisition, source channel 200 rather than the published top/50-m geometries, selected short-file stacks rather than matched one-day stacks, and unnormalized FFT correlations.
- The saved branch is therefore an inconclusive reproduction attempt, not a failure of ambient interferometry.
- Detailed report: ambient_replication_audit_v19.md.
- Next: run the controlled two-geometry transfer test with explicit path manifests, normalized/phase-only correlations, complete UTC-day stacks, seven-day uncertainty, and permutation nulls.

## Controlled Lellouch-style transfer pilot (v20)

- Three independent 10-file chunks (30 segments, approximately 30 min) were processed from the 2024-2025 archive.
- Both published source-receiver geometries were implemented with 5-20 Hz filtering, 5-s running-absolute-mean normalization, pairwise RMS-normalized correlations, and actual 1.020952-m spacing.
- The 3.2-km/s score was -0.00266; the weak maximum was 0.00388 near 1.65 km/s; receiver-position permutation p=0.875.
- Status: negative transfer pilot; not a definitive rejection of the published 2017 result.
- Figure: ambient_transfer/ambient_transfer_v20.png; script: ambient_transfer_test.py.
- Next: extend chunked processing to a complete day and seek original 2017 continuous records if available.

## Duration-matched stack and F-K branch design (v21)

- Seven independent 10-file chunks (70 segments, approximately 70 min) were combined using the published geometry and bounded normalization.
- The 3.2-km/s score is -0.00321; weak peak 0.00481 near 1.40 km/s; receiver-position permutation p=0.561.
- F-K filtering is specified as a separate extension: fixed 5-20 Hz, 2.5-4.5 km/s signed-wavenumber wedges, unfiltered and complementary controls, held-out validation.
- Figure: ambient_transfer/ambient_transfer_v21.png.
- Next: implement and run the F-K wedge branch, then extend baseline to a complete available day if needed.

## Exploratory F-K-filtered ambient correlations (v22)

- Positive signed-wavenumber wedge: 5-20 Hz, apparent velocity 2.5-4.5 km/s, applied before CC.
- Ten-segment score at 3.2 km/s: 0.01257; peak 0.03141 near 3.78 km/s.
- Status: exploratory enhancement only; the mask retains the target velocity interval by construction.
- Required controls: negative and two-sided wedges, equal-bandwidth complement, receiver permutations, and held-out chunks.
- Figure: ambient_transfer/fk_transfer_2025-02-20_start0_n10.png; script: ambient_fk_transfer_test.py.

## Signed F-K control result (v23)

- Positive wedge: 3.2-km/s score 0.0126, p=0.996.
- Negative wedge: peak 3.075 km/s, 3.2-km/s score 0.2846, p=0.002 in 10 segments.
- Two-sided wedge: score 0.1197, p=0.002. Equal-width complement: score -0.1156, p=0.778.
- Three-chunk negative branch: peak 3.075 km/s, score 0.2952, p=0.0005.
- Status: preliminary signed F-K recovery; wedge includes target velocity by construction.
- Figure: ambient_transfer/fk_negative_v23.png.
- Next: hold the signed wedge fixed on independent days and test wedge sensitivity before calling this a Lellouch reproduction.

## Segment-duration clarification (v24)

- One HDF5 segment contains 30,000 samples at 500 Hz: 60 seconds.
- Thus 10 segments ≈10 minutes, 30 segments ≈30 minutes, and 70 segments ≈70 minutes.
- This corrects documentation wording only; computed results are unchanged.

## Longer signed F-K stack (v25)

- Frozen negative wedge extended to seven independent 10-file chunks (70 segments, approximately 70 min).
- Peak 3.075 km/s; 3.2-km/s score 0.2947; receiver-position permutation p=0.0002.
- Seven chunk scores at 3.2 km/s are all positive: 0.285, 0.316, 0.285, 0.283, 0.262, 0.306, 0.327.
- Status: stable 70-minute F-K-selected observable; complete-day and held-out-day validation remain.
- Figure: ambient_transfer/fk_negative_v25.png.
## Conditional 500–520 m anomaly test (v37)

- Assumption: the current records use the same cemented main-hole fiber as Lellouch et al.; current channel 0 and 1.020952 m spacing inherit the Lellouch approximate position-as-depth convention.
- Completed diagnostic: 369 ten-minute negative signed F–K products from four dates. Early 50–350 m apparent velocities are 3.13–3.19 km/s; late 450–650 m velocities are 2.49–2.61 km/s.
- Descriptive breakpoint medians are 350–425 m, so the curvature is repeatable but is not uniquely localized at 500–520 m. This is not a lithology or Vp inversion.
- Products: `fk_geology_anomaly_test.py`, `fk_geology_anomaly_test.json`, and `fk_geology_anomaly_test.png`.

## Follow-up validation (running)

- Seasonal F–K: job `37478467`; exact-resolution unfiltered control: `37478473`; aggregate: `37478474`.
- Multi-band test: job `37492825` (dependency on seasonal F–K), using 3–8, 5–12, 8–20, and 15–30 Hz negative signed F–K bands on deterministic blocks from each selected day.
- Full anomaly rerun after seasonal aggregation: job `37493070` (dependency on `37478474`).
- No seasonal or frequency-band conclusion is assigned until these jobs complete and their day-level products are inspected.

