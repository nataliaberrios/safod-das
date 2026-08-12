# Ambient F–K injection–recovery v1 — withdrawn

The v1 sensitivity thresholds are invalid and must not be interpreted or
promoted. The synthetic wave was added directly to large-offset `float32` raw
samples. At the observed raw magnitude, the `float32` spacing was larger than
the smallest injected signals, so part of the 0.003–0.03 RMS ladder was rounded
away before preprocessing.

The replacement v2 workflow keeps the production raw-data path unchanged,
applies the linear detrend and 5–20 Hz bandpass separately to real and
synthetic arrays, and adds their `float64` outputs before the nonlinear 5 s
running-absolute-mean normalization. Unit tests require both survival of a
sub-count injection and equivalence to production preprocessing at zero
amplitude. The v2 result is promoted only after all 300 files, eight
direction–velocity scenarios, aggregation, and an independent completion audit
pass.

The retained v1 products are provenance only.
