# Lellouch Figure 7c on the 2024–2025 archive — status

Written 2026-08-14. Where numbers here disagree with an older notebook cell or a
figure caption, this file wins.

**Question.** Can Lellouch et al. (2019) Figure 7c — the top-source ambient
correlation gather showing "a clear wave packet with an apparent velocity of about
3,200 m/s" — be reproduced from the 2024–2025 continuous SAFOD main-hole archive?

---

## 1. Verdict

**No, and the reason is the input, not the processing. The wavefield property the
method requires is absent from this archive: there is no net downgoing energy in
the body-wave velocity fan on any of eight days spanning ten months.**

This is a data-limited negative result, not a refutation of Lellouch et al. and
not a failure of the pipeline. Two independent lines support it, and the method
itself is validated against the paper's own released products.

| Claim | Status | Evidence |
|---|---|---|
| Figure 7c reproduces on 2024-12-20 | **No** | 1,440 files, 4,320 windows, 4 preprocessing configurations, all p ≥ 0.19 (§3) |
| The picker and moveout machinery are correct | **Established** | reproduces the paper's Figure 9 profile from its released 7d traces, r = 0.948 (§2) |
| The 2.5–4 km/s fan carries usable ambient energy | **No** | 3–5 % of 5–20 Hz energy; 81–99 % sits below 1,500 m/s (§4) |
| That energy is net downgoing, as the paper requires | **No** | downgoing share never exceeds 49.9 % on any day; two summer days favour *upgoing* (§4) |
| Failure is caused by a missing processing step | **Refuted** | four configurations spanning the three candidate steps all fail identically (§3) |
| An F–K fan would help | **No — it manufactures the answer** | fails the pre-filter channel-scramble gate; see `Ambient_FK_QC_workflow.ipynb` |

## 2. The method is validated before it is used

`ambient_lellouch_fig7d_profile.py` applies the paper's picker — three adjacent
samples with the largest correlation plus quadratic interpolation (Nakata &
Snieder, 2012) — to Lellouch's own released Figure 7d correlograms in
`lellouch_traces/`:

| depth | 50 m | 150 m | 250 m | 450 m | 550 m | 700 m |
|---|---|---|---|---|---|---|
| v (m/s) | 2416 | 2803 | 3192 | 3305 | 4061 | 4357 |

Fit `v = 2.548·z + 2354`, r = 0.948, mean 3310 m/s. This reproduces the paper's
Figure 9 interferometry model, including its own description of it as following
"a somewhat linear trend", and it recovers the quoted 3,200 m/s as the *depth
average of a gradient* rather than a single ridge velocity. So the picker works
and the target is well posed.

**Two geometries, two jobs.** Figure 7c (top source, receivers every 50 m) is a
phase-identification step yielding one average velocity, which is then reused as
the moveout-correction velocity. The velocity *model* comes from Figure 7d
(constant 50 m offset). Figure 7c is therefore the correct target for
"reproduce the gather", and it is shift-free: the paper applies the 3,200 m/s
shifts only after an average velocity has been obtained by simple stacking.

## 3. The reproduction attempt

`ambient_lellouch_fig7c.py` + `.sbatch`, jobs `38943063`, `38944071`, `38944073`,
`38944075`. Date 2024-12-20, all 1,440 one-minute files, 4,320 windows of 30 s
with 15 s overlap, source channel 0, receivers every 50 m to 700 m, R±10
neighbour stack (simple sum, no shifts), final 5–20 Hz band. **No F–K filter
anywhere.**

Three steps absent from `ambient_transfer_test.py` were added and then ablated
independently, so none is assumed:

- **strain-rate conversion** — the archive stores `RawDataUnit = 'rad * 2PI/2^16'`,
  `RawDescription = 'Diversity Processed Phase'`, i.e. phase/strain, so the paper's
  strain→strain-rate step requires a time derivative. `ambient_transfer_test.preprocess`
  never differentiates; the production chain has been correlating strain.
- **R±10 neighbour stack** — the paper states it "is required to extract a clear
  signal". The cached products correlate single pairs, ~21× less SNR.
- **spectral normalisation** (eq. 6, `|S(ω)|²`) and **common-mode removal** (median
  across channels, as in the legacy CC pipeline and the v19-audited notebook).

| common mode | whitening | peak \|score\| | at | score at 3200 m/s | perm. p | causal/acausal |
|---|---|---|---|---|---|---|
| on | on | 0.0080 | 5900 | −0.0006 | 0.20 | 0.99 |
| off | on | 0.0011 | 4700 | −0.0005 | 0.19 | 0.88 |
| on | off | 0.0073 | 2750 | +0.0023 | 0.34 | 0.97 |
| off | off | 0.0007 | 5500 | −0.0001 | 0.87 | 0.83 |

No configuration produces a 3,200 m/s peak, and **none reaches the paper's
qualitative requirement that the causal side dominate the acausal side** — the
ratio is ≤ 0.99 everywhere. The receiver-order permutation null is never cleared.

## 4. Why — the input census

`ambient_fk_energy_census.py`, jobs `38944180`, `38945793`. Raw wavefield only:
strain-rate proxy, 5–20 Hz, 2-D FFT in 30 s windows, energy binned by apparent
velocity and by the sign of f·k. No correlation, no filtering, no selection.
Twelve files spread across each day at 2 h intervals.

| date | < 1500 m/s | 2.5–4 km/s fan | downgoing | upgoing |
|---|---|---|---|---|
| 2024-05-11 | 98.8 % | 0.3 % | 47.3 % | 52.7 % |
| 2024-06-17 | 80.9 % | 4.5 % | 45.6 % | 54.4 % |
| 2024-06-26 | 84.8 % | 3.6 % | 46.6 % | 53.4 % |
| 2024-10-28 | 83.6 % | 4.0 % | 49.8 % | 50.2 % |
| 2024-11-30 | 81.0 % | 4.9 % | 49.9 % | 50.1 % |
| 2024-12-20 | 82.8 % | 4.1 % | 49.8 % | 50.2 % |
| 2025-02-24 | 85.3 % | 3.4 % | 49.6 % | 50.4 % |
| 2025-03-04 | 87.7 % | 2.9 % | 49.8 % | 50.2 % |

The paper attributes Figure 7c's causal dominance to ambient sources at the
surface sending energy *down* the hole. In this archive that asymmetry does not
exist on any day: the downgoing share never exceeds 49.9 %, and the two summer
days mildly favour upgoing.

**Aperture-resolution check.** 900 channels × 1.0209 m = 918.9 m, so
Δk = 0.00109 cycles/m and the 3,200 m/s fan sits only 1.4 wavenumber bins from
k = 0 at 5 Hz. Restricting to 12–20 Hz, where the fan is 4.3–5.7 bins out and
properly resolved, the fan share *rises slightly* to a 4.56 % mean while the
downgoing share stays at a 49.6 % maximum. The resolution limit therefore does
not explain the result.

`sanity/README.md` states the governing rule: "if there is no body-wave-velocity
ridge in the wavefield, no CC trick will recover a body-wave Green's function."

## 5. What this is not

- **Not a criticism of Lellouch et al. (2019).** Figure 7c is a real result on
  2017 data, and this project's own picker reproduces its companion Figure 9.
- **Not a claim that ambient interferometry fails at SAFOD.** It is a claim about
  this archive, these eight days, and the 5–20 Hz band.
- **Not a licence to reach for the F–K fan.** The fan produces the expected
  geometry from channel-scrambled input; that is why the QC workflow rejects it.

## 5a. Quiet-window selection — tested and closed

The one remaining live avenue was temporal selection: if the body-wave fraction is
diurnal, a quiet subset might carry the arrival that the 24 h average dilutes. The
`fig7c` chunks are 60 files each, so each chunk is exactly one UTC hour and the
test needs no new compute.

Causal/acausal ratio by UTC hour on 2024-12-20 spans **0.93 to 1.08** — scatter
about 1.0 with no structure. Hour of day is independent of the observable, so this
is a fair test rather than a conditioned one. Stacking the six most
causally-dominant hours, selected on the same day, gives ratio 1.02, peak
\|score\| 0.0082 at 5,900 m/s, score at 3,200 m/s of −0.0001; splitting those six
into halves gives 1.02 and 1.03. **Selection buys nothing, which is what happens
when there is nothing to select.**

Note also that the scan's peak sits at ~5,900 m/s — the top of the trial range —
in 20 of 24 hours. That is the flat-moveout end of the scan, not an arrival.

## 6. What would change the answer

Ordered cheapest first. None is required for the negative result to stand, and the
cheapest one has already been tried (§5a).

1. ~~**Quiet-window selection.**~~ Tested, negative — see §5a.
2. **Higher band.** The fan is better resolved and marginally richer above 12 Hz.
   A 12–30 Hz reproduction departs from the paper but is better matched to the
   aperture and the 16.335 m gauge length (2017 used 10 m).
3. **Multi-day stacking.** Eight complete days are cached. Lellouch needed one, and
   stacking cannot create a directional asymmetry that is absent per-day, so this
   is a completeness exercise rather than a fix.
4. **Instrument-transfer comparison.** This acquisition is OptaSense IU,
   16.335 m gauge, 500 Hz output; the 2017 records are a different instrument.
   Quantifying that difference converts the negative result into a measurement.

## 7. Related correction

All 1,162 cached `ambient_transfer/transfer_*_start*_n*.npz` chunks were written
2026-08-03 16:07 → 2026-08-04 18:14, **before** the negative-lag fix committed in
`0931988` on 2026-08-05 (documented in `ambient_signed_lag_audit_v43.md`). Verified
empirically: the old routine's autocorrelation is asymmetric by 0.025 where it must
be 0, and its negative lags correlate with correct ones at −0.0017 while positive
lags correlate at 1.0000. **Zero chunks were regenerated after the fix**, so any
acausal quantity read from that cache is void; positive lags are unaffected. The
F–K null products are *not* affected — `ambient_fk_full_pipeline_null_v2.py`
reprocesses from raw HDF5.

Separately, the aggregated `transfer_seasonal_<date>.npz` products retain only
`top_stack`; `fixed_stack` is dropped at aggregation and so the Figure 7d
observable was never available downstream.

## 8. Provenance

| Script | Product | Jobs |
|---|---|---|
| `ambient_lellouch_fig7c.py`, `.sbatch` | `ambient_transfer/lellouch_fig7c/` | 38943063, 38944071/73/75 |
| `ambient_fk_energy_census.py` | `ambient_fk_energy_census_<date>.{npz,png,txt}` | 38944180, 38945793 |
| `ambient_lellouch_fig7d_profile.py` | `ambient_lellouch_fig7d_profile.{npz,png,txt}` | run interactively |

`38945793_0` (2024-05-11) was OOM-killed at 32 GB after printing its summary; its
numbers above are recovered from `logs/fkcensall_38945793_0.out` and it has no
`.npz`, so it is excluded from the 12–20 Hz re-analysis.
