# Curated ambient-transfer products

This directory contains generated products from the ambient-noise branch. Small
review artifacts are tracked when needed; raw waveform files, distributed
checkpoint arrays, and large intermediate products remain on Sherlock.

## Current evidence status (v48)

F–K filtering is retained as a valid directional wavefield operation. The
2.5–4.5 km/s selected correlation observable repeats across held-out dates and
appears with comparable magnitude on both physical lag branches after the
corrected sign convention. Its apparent velocity remains conditional on the
chosen velocity support:

- On seven held-out dates, the production wedge peaks at 3.075 km/s on both
  branches, but the broad 2.0–5.5 km/s wedge peaks at 2.275 and 1.850 km/s.
  Direction-only filtering reduces the held-out 3.2 km/s scores to 0.0004 and
  0.0107. This completed constraint-relaxation test is the direct check against
  imposing a P-wave-like velocity band.
- The exact-resolution unfiltered 11,523-file aggregate is weak at the fixed
  3.2 km/s trajectory (causal p = 0.341; anti-causal p = 0.202).
- The independent 300-file pre-filter F–K energy test gives geometric
  target/reference ratios 0.95229 and 0.95227, with familywise p = 1.000 for
  both channel-permutation and independent-time-shift nulls.
- The corrected five-hour decimation audit gives nearly identical selected
  peaks for direct, explicitly anti-aliased, and full-resolution processing
  (3.075, 3.075, and 3.050 km/s), so the tested direct slicing does not create
  the selected ridge.
- The precision-safe real-noise injection–recovery calibration uses 300 paired
  files, four velocities, both directions, and an off-wedge control. The 3.2
  km/s post-filter threshold of 0.003 is conditional on an uninjected peak
  already near 3.075 km/s and is not a generic 0.3% detection limit.

These results do not invalidate F–K filtering or prove that weak coherent
energy is absent. They mean the selected ridge must not currently be presented
as ambient Green's-function convergence, an identified P wave, or an unbiased
formation velocity.

## Independent validation products

- `fk_prefilter_energy_v1_n300_r20/fk_prefilter_energy_v1_aggregate.json` and
  `_figure.png` — completed pre-filter statistic, decision rule, and two-null
  comparison.
- `signed_lag_v2/seasonal_signed_fk_v2_aggregate.json`, `.png`, and
  `seasonal_signed_fk_v2_completion_audit.json` — corrected physical-lag
  comparison and independent audit across 11,523 one-minute files.
- `fk_injection_recovery_v2_n300/ambient_fk_injection_recovery_v2_aggregate.json`,
  `.png`, and `ambient_fk_injection_recovery_v2_completion_audit.json` —
  precision-safe injection–recovery summary and independent 2,000-bootstrap
  audit. Per-scenario arrays remain untracked on Sherlock.
- `fk_filter_kernel_residual_v1/` — secondary filter-kernel subtraction
  diagnostic; exploratory because its null families disagree.
- `fk_full_pipeline_null_v2_n300_r20/` — superseded selected-operator
  diagnostic. Because filtering and scoring share the velocity corridor, these
  products cannot accept or reject a physical wave.

## Seasonal and sensitivity products

- `fk_mask_sensitivity_v2/ambient_fk_mask_sensitivity_v2.*` — frozen
  development/held-out comparison of narrow, production, broad, and
  direction-only F–K masks; this is the constraint-relaxation result.
- `fk_seasonal_day_sections.png` and `.pdf` — eight-day correlation sections
  plus across-day aggregate, with common axes and reference moveout.
- `fk_seasonal_day_comparison.png` — legacy scatter retained as supplementary
  quality control only.
- `fk_seasonal_aggregate.json` — selected signed F–K metrics and
  receiver-permutation nulls; these describe the filtered observable.
- `seasonal_unfiltered_aggregate.json` — exact-resolution unfiltered control.
- `frequency_band_anomaly_test.json` — frequency-band moveout/breakpoint
  summary.
- `alias_sensitivity_2024-12-20_start0_n300.json` and `.png` — completed
  corrected five-hour direct/anti-aliased/full-resolution comparison.
- `alias_sensitivity_2024-12-20_start0_n30.*` — withdrawn provenance product
  from the earlier implementation that scored before the intended F–K mask;
  do not use it as an F–K validation.

The v1 injection thresholds are withdrawn for float32 sub-count rounding; see
`../ambient_fk_injection_recovery_v1_WITHDRAWN.md`. The authoritative
interpretation and detailed captions are in `../AWD_results_dashboard.ipynb`
(v48) and `../AWD_advisor_figure_guide.tex` (v48).
