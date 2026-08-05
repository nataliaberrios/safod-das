# Curated ambient-transfer products

This directory contains generated products from the ambient-noise branch. Only
the small, review-ready summary artifacts listed below are tracked in GitHub;
raw waveform files, checkpoint arrays, and intermediate chunk products remain
on Sherlock.

## Seasonal validation

- `fk_seasonal_day_comparison.png` — eight-day signed F–K versus day comparison.
- `fk_seasonal_day_comparison.pdf` — vector version of the same figure.
- `fk_seasonal_aggregate.json` — signed F–K aggregate metrics and permutation nulls.
- `seasonal_unfiltered_aggregate.json` — exact-resolution unfiltered control metrics.

## Processing-sensitivity diagnostics

- `frequency_band_anomaly_test.json` — frequency-band moveout/breakpoint summary.
- `alias_sensitivity_2024-12-20_start0_n30.png` — direct, anti-aliased, and full-resolution pilot figure.
- `alias_sensitivity_2024-12-20_start0_n30.json` — numerical metrics for that pilot.

The authoritative interpretation and figure captions are in
`../AWD_results_dashboard.ipynb` (v41) and
`../AWD_advisor_figure_guide.tex` (v41).
