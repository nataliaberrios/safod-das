# Ambient signed-lag audit (v43)

## Outcome

The legacy ambient correlation routine extracted zero and positive lags
correctly but assembled negative lags from the tail of a truncated FFT array.
The corrected routine extracts negative lags from the tail of the complete
zero-padded circular-correlation array.  This changes the interpretation of
the signed F–K controls:

- The established negative-mask, positive-lag moveout remains reproducible.
- The legacy positive-mask value at positive lag is an opposite-mask leakage
  control, not a test of a physical opposite-propagating branch.
- A physical positive-mask branch must be scored at negative lag.
- Corrected one-file and 30-file pilots recover comparable approximately
  3.1 km/s branches at opposite correlation lags.

## Coordinate and sign convention

The virtual source is Nano channel 0. Receiver coordinate increases along the
fiber from channel 0, and the top-source pairs are `(0, receiver)` at offsets
of approximately 50–700 m.

For a plane wave traveling toward increasing coordinate,

`u(z,t) = g(t - z/v)`,

the production Fourier convention places energy in `F*K < 0`. The receiver is
delayed relative to channel 0, so `conj(source_fft)*receiver_fft` peaks at
positive lag:

`lag = +offset / v`.

For a wave traveling toward decreasing coordinate,

`u(z,t) = g(t + z/v)`,

the energy lies in `F*K > 0`, and the receiver leads channel 0. The physical
correlation branch therefore peaks at

`lag = -offset / v`.

These are coordinate-direction labels. The Nano channel-orientation evidence
is required before translating them into downhole/upgoing language.

## Original implementation error

With trace length `n`, the FFT was padded to `nfft >= 2*n-1`. The circular
correlation was then truncated to `2*n-1` samples before its tail was used as
the negative-lag half. True negative lags wrap to the end of the full `nfft`
array, so this truncation discarded them.

The corrected extraction is conceptually:

```python
circular = irfft(conj(source_fft) * receiver_fft, n=nfft)
correlation = concatenate((circular[-max_lag_samples:],
                           circular[:max_lag_samples + 1]))
```

The associated lag axis is `[-max_lag, ..., 0, ..., +max_lag]`.

## What is and is not affected

| Product or claim | Status after audit |
|---|---|
| Negative mask evaluated at positive lag | Numerically preserved |
| Positive mask evaluated at positive lag | Preserved only as leakage control |
| Negative-lag half of every legacy correlation section | Invalid; recompute |
| Claim that the physical positive branch is absent | Withdrawn |
| Unfiltered and F–K positive-lag comparison | Preserved as a target-moveout comparison |
| Causal/anti-causal amplitude or symmetry comparison | Pending corrected seasonal rerun |
| Exact apparent velocity as a formation property | Still unsupported because the F–K wedge is selected |

## Regression tests

`test_ambient_signed_lags.py` must pass before any corrected production run.
It currently verifies:

1. A receiver delayed by 0.20 s correlates at `+0.20 s`.
2. Reversing the channel pair correlates at `-0.20 s`.
3. A broadband wave traveling toward increasing coordinate is retained by
   `F*K < 0` and follows positive moveout.
4. A broadband wave traveling toward decreasing coordinate is retained by
   `F*K > 0` and follows negative moveout.

All tests pass in the `das` environment.

## Real-data pilot

The first 30 one-minute files on 2024-12-20 were recomputed from raw DAS data
with the corrected lag extraction. Branch-specific receiver-permutation tests
use the physical lag sign for each mask.

| Mask and physical lag | Peak speed | Peak score | Null 95% | p-value | 3.2 km/s opposite-lag leakage |
|---|---:|---:|---:|---:|---:|
| `F*K < 0`, positive lag | 3.075 km/s | 0.253 | 0.140 | 0.0005 | 0.032 |
| `F*K > 0`, negative lag | 3.075 km/s | 0.264 | 0.152 | 0.0005 | 0.029 |

This pilot supports a two-sided correlation pair. It is not yet the seasonal
result and remains conditional on the fixed 2.5–4.5 km/s F–K wedge.

## Versioned full recomputation

Legacy products remain in `ambient_transfer/`. Corrected products are written
only to `ambient_transfer/signed_lag_v2/`.

- Seasonal array: `37717572` (eight independent days)
- Dependent postprocess: `37717613`
- Array launcher: `seasonal_signed_fk_v2.sbatch`
- Chunk processor: `ambient_signed_fk_v2.py`
- Aggregator: `aggregate_seasonal_signed_fk_v2.py`

The postprocess runs only if every array task succeeds. It will produce a
common-scale causal/anti-causal aggregate figure, day-level branch sections,
branch-specific velocity scores, opposite-lag leakage controls, and 5,000-draw
receiver-permutation nulls.
