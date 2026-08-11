# Curated ambient-transfer products

This directory contains generated products from the ambient-noise branch. Small
review artifacts are tracked when needed; raw waveform files, distributed
checkpoint arrays, and large intermediate products remain on Sherlock.

## Current evidence status (v45)

F–K filtering is retained as a valid directional wavefield operation. The
2.5–4.5 km/s selected correlation observable repeats across held-out dates, but
its physical interpretation is not independently supported by the current
controls:

- The exact-resolution unfiltered 11,523-file aggregate is weak at the fixed
  3.2 km/s trajectory (causal p = 0.341; anti-causal p = 0.202).
- The independent 300-file pre-filter F–K energy test gives geometric
  target/reference ratios 0.95229 and 0.95227, with familywise p = 1.000 for
  both channel-permutation and independent-time-shift nulls.
- The filter-kernel residual diagnostic is exploratory because the two null
  families disagree and its positive peak lies near the upper search edge.

These results do not invalidate F–K filtering or prove that weak coherent
energy is absent. They mean the selected ridge must not currently be presented
as ambient Green's-function convergence, an identified propagating mode, or
formation velocity.

## Independent validation products

- fk_prefilter_energy_v1_n300_r20/fk_prefilter_energy_v1_aggregate.json —
  completed statistics, decision rule, and cross-node reproducibility metadata.
- fk_prefilter_energy_v1_n300_r20/fk_prefilter_energy_v1_aggregate_figure.png —
  publication-facing pre-filter power and two-null comparison.
- fk_filter_kernel_residual_v1/fk_filter_kernel_residual_v1.json and
  fk_filter_kernel_residual_v1_figure.png — secondary filter-kernel subtraction
  diagnostic.
- fk_full_pipeline_null_v2_n300_r20/ — superseded selected-operator diagnostic.
  Because filtering and scoring use the same velocity corridor, these products
  cannot accept or reject a physical wave.

## Seasonal and sensitivity products

- fk_seasonal_day_sections.png — eight-day correlation sections plus across-day
  aggregate, with common axes and reference moveout.
- fk_seasonal_day_sections.pdf — vector version.
- fk_seasonal_day_comparison.png — legacy scatter retained as supplementary
  quality control only.
- fk_seasonal_aggregate.json — selected signed F–K metrics and
  receiver-permutation nulls; these describe the filtered observable.
- seasonal_unfiltered_aggregate.json — exact-resolution unfiltered control.
- frequency_band_anomaly_test.json — frequency-band moveout/breakpoint summary.
- alias_sensitivity_2024-12-20_start0_n30.* — withdrawn provenance product from
  an implementation that scored before the intended F–K mask; do not use it as
  an F–K validation.

The authoritative interpretation and detailed figure captions are in
../AWD_results_dashboard.ipynb (v45) and
../AWD_advisor_figure_guide.tex (v45).
