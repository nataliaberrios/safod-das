# Why the Figure 7c ambient arrival does not reproduce: the low-wavenumber mechanism

Status document. Written 2026-08-19. Where this file and an older notebook cell,
figure caption or script docstring disagree, **this file wins** for the ambient
Figure 7c mechanism question specifically. It does not overrule
`Ambient_FK_QC_workflow.ipynb` on the F-K fan verdict, which still stands.

> **READ THIS FIRST — there are now TWO findings, and this file holds the second
> one.** Later the same day, `interrogator_and_illumination_v2.py` found something
> upstream of everything below: the 2024-25 ambient field carries **no
> downgoing/upgoing asymmetry** in the body-wave fan, while Lellouch's 2017
> records do (|A| = 0.348, p = 0.005 versus 0.040, p = 0.73, same measurement,
> matched spatial ranks). His downgoing P is an *inference from that asymmetry*,
> so the absence of the asymmetry means there is no arrival of that kind present
> to be recovered — independent of any filter, and independent of the static
> pattern described in this file.
>
> The correct ordering of the two findings is therefore:
>
> 1. **Illumination is absent** (`ambient_directional_asymmetry.py`,
>    `interrogator_and_illumination_v2.py`, `illumination_window_scan.py`) — the
>    wavefield does not contain a net downgoing component in this band. This is
>    the primary explanation and it is about the *recording*. The archive-wide
>    scan is now complete and closes the route rather than merely suggesting it:
>    of **240 windows** spanning 2024-05-21 to 2025-05-06, **11 reach p < 0.05
>    against 12.0 expected by chance** (Binomial 95th percentile 18), so no
>    illuminated window was found anywhere. The decision rule was fixed before the
>    run. `illumination_window_scan.txt`.
> 2. **A static fixed-k pattern dominates** (this file) — which explains why the
>    eight velocity-domain methods all failed *in the same way*, and why the
>    pedestal could not be filtered out. This remains correct and useful, but it
>    is the explanation for the *processing* failures, not for the absence of the
>    arrival.
>
> Do not cite this file alone as the reason Figure 7c does not reproduce. See
> `AMBIENT_CC_LITERATURE_REVIEW.md` §2 for why illumination is the binding
> constraint in borehole geometry, and §5 for the overall status.

## 1. The result

The 2024-25 ambient field at 5-20 Hz is dominated by a **static spatial pattern
at fixed wavenumber**, not by a propagating wavefield. Measured on the first 16
non-zero wavenumber cells of a 700 m aperture, with k = 0 itself excluded:

| | 5-12 Hz | 12-20 Hz | ratio |
|---|---|---|---|
| band centre frequency | 8.50 Hz | 16.00 Hz | **1.882** |
| power-weighted \|k\| centroid | 0.00173 cyc/m | 0.00175 cyc/m | **1.011** |
| cells from k = 0 | 1.21 | 1.23 | |
| share of in-band energy | 39.33 % | 39.08 % | |

The discriminator is that apparent velocity is `v = f / k`:

- a **wave** at fixed `v` must have its wavenumber track frequency, so the
  centroid ratio should equal the frequency ratio, **1.882**;
- a **static pattern** at fixed `k` predicts a centroid ratio of **1.000**.

Observed: **1.011**. The frequency nearly doubled and the wavenumber moved by
1.1 %. The feature is at fixed wavenumber.

**Sensitivity — the test had the range to see a wave and did not.** A 3200 m/s
arrival sits at 1.86 cells at 8.5 Hz and 3.50 cells at 16.0 Hz, both inside the
16-cell window measured, so a dominant arrival would have moved the centroid.

Produced by `ambient_fixed_k_test.py` from the `K0_REMOVE=1` census products
`ambient_apparent_velocity_census_{5-12Hz,12-20Hz}_k0rm.npz`.

## 2. Why this explains every failed method

A static spatial pattern **has no velocity**. Every method tried against the
pedestal separates energy *by velocity*:

| method | script |
|---|---|
| F-K velocity fan, brick-wall and raised-cosine | `ambient_fk_taper_test.py` |
| linear Radon / tau-p with a slowness mute sweep | `ambient_radon_slant_stack.py` |
| rank-k coherent-subspace projection, global and windowed | `--svd-rank` in `ambient_lellouch2019_exact_stack.py` |
| phase cross-correlation and phase-weighted stacking | `ambient_pcc_pws.py` |
| median flat-event removal along the offset axis | `ambient_flat_event_removal.py` |
| per-sample median and mean common-mode removal | `--common-mode`, `--common-mode-estimator` |

All eight failed in the same way, and this is why: they are asking a
velocity-domain question of something that does not have a velocity. The
corollary is about **method, not absence** — the removal has to be **spatial**:
estimate and divide out the static per-channel response.
`C2_PERMEABILITY_FOLLOWUP.md` already proposes that calibration and records 5x
headroom from it; this is independent motivation for doing it.

## 3. What this does NOT say

It says the **dominant** low-wavenumber feature is not a wave. It does **not**
say no arrival exists. A weak arrival can hide beneath a feature holding ~39 % of
the in-band energy. Do not cite this file as evidence that the ambient arrival is
absent.

Independently established and unchanged: 6 days give min p = 0.1345; the coherent
96 h stack gives p = 0.9184; **0 of 181 velocities clear the per-velocity null**
(at 3,200 m/s specifically, 2.752 against a per-velocity threshold of 2.811); and
the Figure 7d constant-offset gather peaks at exactly 0.0000 s at every depth
while Lellouch's migrates 20.7 -> 11.5 ms.

**Corrected 2026-08-19.** This paragraph previously ended "Our picker reproduces
his Figure 9 from his own released traces at r = 0.948, so the method is not at
fault." That is **withdrawn as worded**. r = 0.948 is `corr(depth, velocity)` for
picks our picker made on his released 7d traces — a measure of how monotonic our
own picks are, **not** agreement with his published Figure 9 curve, which was
never compared value-for-value. The supportable version: on known-good input our
picker returns a monotonic velocity-depth trend of the expected form (2,416 m/s
at 50 m rising to 4,357 m/s at 700 m), so the picking stage is not what fails on
the 2024-25 data.

## 4. Two retractions, in order — read this before quoting any k = 0 number

This question produced a wrong claim and then a wrong retraction of it. Both are
recorded because the numbers appear in commit messages that are already pushed.

**4.1 Commit `8fa2266` overreached and is withdrawn in part.** It concluded that
"the contaminant and the target occupy the same unresolved region of wavenumber
space" and that reproducing Figure 7c "is not achievable on this recording". The
aperture arithmetic in it is correct and stands:

| freq | k at 3200 m/s | cells from k=0 (dk = 1/700) | inside Hann main lobe (+-2)? |
|---|---|---|---|
| 5 Hz | 0.00156 | 1.09 | yes, unresolved |
| 10 Hz | 0.00313 | 2.19 | marginal |
| 15 Hz | 0.00469 | 3.28 | no, resolved |
| 20 Hz | 0.00625 | 4.37 | no, resolved |

But the sweeping conclusion does not follow from it, for two reasons. First, by
its own table the target **is** resolved above ~12 Hz, so "unresolvable" was
never true across the band. Second, the prediction it implied was tested and
**failed**: at 12-20 Hz, where the target is resolved, 2024-25 shows *less* fan
energy (0.42 %) than 2017 (1.15 %) and k = 0 stays at 54.74 %. Better resolution
did not expose an arrival. The correct mechanism is section 1, which does not
depend on resolution at all.

**4.2 The withdrawal of the cross-epoch k = 0 comparison was itself wrong, and
that comparison is REINSTATED.** The 2017 release files `M1p33` and `M2p46` are
earthquake records, so comparing one of their windows against 2024-25 ambient
noise looked like the two-arms-processed-differently error behind five of the six
earlier withdrawals in this tree. `lellouch2017_window_audit.py` then **measured**
the arrival times rather than assuming them:

| record | arrival | census window | verdict |
|---|---|---|---|
| M1p33 | 4.100 s of 5.00 s | first 2.5 s | arrival is **outside** |
| M2p46 | 4.916 s of 5.00 s | first 2.5 s | arrival is **outside** |

Both arrivals fall after the 2.5 s cut, so the census's 2017 arm is genuine
pre-event noise, both arms were noise, and the comparison is not confounded by
signal presence. Its originally stated limits still apply: only ~5 s of 2017
data, and gauge length 10 m against 16.335 m (a 0.1-1.1 % effect at 3200 m/s).
Measured k = 0 share of the 2017 pre-event window: **0.13 %** (M1p33) and
**0.11 %** (M2p46), from `lellouch2017_window_audit.txt`.

> **Corrected 2026-08-20 (audit).** This sentence previously ended "against
> **~66 %** for 2024-25 ambient." **No product in this tree contains 66 %.** The
> matched measurement — both arms through
> `ambient_apparent_velocity_census.py`, same band, same k = 0 handling — gives
> **52.12 %** for 2024-25 against **0.04 %** for 2017
> (`ambient_apparent_velocity_census.txt:26`). Quoting 66 % against the 0.13 %
> from a *different* product with different processing was exactly the
> unmatched-arms error that caused five of this project's ten withdrawals. The
> direction of the contrast is unchanged and is, if anything, starker on the
> matched numbers; only the figure quoted was wrong.

## 5. Is it the interrogator?

**Partly answered 2026-08-19, and the answer is no — not inside the aperture that
matters.** `interrogator_and_illumination_v2.py` tested whether the dominant
spatial pattern is *stable across days*, which a fixed instrumental fingerprint
must be, on six days spanning 2024-05-21 to 2025-05-06:

| channels measured | cross-day median \|corr(u1,u1)\| | control \|corr(u1,u2+)\| | verdict | product |
|---|---:|---:|---|---|
| **0–699, lead-in INCLUDED** | **0.8426** | 0.3889 (95th) | stable | `interrogator_blame_test.txt` (v1) |
| 23–708 (the analysed aperture) | **0.0188** | 0.1282 (95th) | **not** stable — not an instrumental fingerprint | `interrogator_and_illumination_v2.txt` (v2) |

> **Corrected 2026-08-20 (audit).** The first row previously read
> "0–22 (surface lead-in) | 0.8426 | 0.3889 | stable — instrumental", implying two
> matched measurements over disjoint channel ranges. **There is no
> channels-0–22-only correlation measurement.** 0.8426 is v1's run over
> **channels 0–699** with the lead-in included; the row is relabelled above to
> what was actually computed. The inference that the stability *comes from* the
> lead-in rests on v1's separate **power-ratio** test (T2), not on a correlation
> over 0–22 — and that test's own output warns "a high ratio alone is NOT proof
> of an instrumental origin … does not decide on its own", with one of its six
> days (2025-03-13) giving an *inverted* ratio of 0.01 against the 2e4–1.3e5 band
> quoted elsewhere in this file.

So the interrogator imprints a stable spatial pattern on the array **as a whole
when the uncemented surface lead-in is included**, and the most likely seat of it
is that lead-in — **which the Figure 7c pipeline already excludes.** What is
directly measured, and what the conclusion rests on, is the second row: inside
channels 23–708 the dominant pattern is not distinguishable
from unrelated patterns, so an instrumental fingerprint is not supported there.
This also supersedes v1 of that script, whose "instrumental" verdict was
unrestricted and whose |A| = 1e-4 "no illumination" number was the algebraic
k-symmetry of a separable static pattern rather than a property of the wavefield.

Note this does **not** retract section 1: the static fixed-k pattern is still
there and still holds ~39 % of in-band energy. It says only that the pattern is
not a *time-invariant instrument signature*; a day-to-day-varying static spatial
pattern is fully consistent with section 1.

The settings that differ between the two acquisitions, retained for the record:

| | Lellouch 2017 | 2024-25 |
|---|---|---|
| channel spacing | 1.0 m | 1.0209 m |
| sample rate | 250 Hz | 500 Hz (decimated to 250 for all comparisons) |
| gauge length | 10 m | 16.335 m |
| output unit | — | `rad * 2PI/2^16`, i.e. absolute optical phase |

Channel spacing is effectively identical and **gauge length is not the
explanation** — the effect is measured at 0.1-1.1 % at 3200 m/s, three orders
below what would be needed. The last row remains the most plausible instrumental
candidate: absolute optical phase with no common-mode rejection carries a
laser/temperature/PSU term shared across channels, which is low-wavenumber energy
by construction, and it is consistent with the fixed-k finding in section 1.

**But it is now a weaker candidate than when this was written**, because a
laser/PSU term should be *stable across days* and the table above shows the
dominant pattern over channels 23-708 is not (0.0188 against a control 95th of
0.1282). Either the term drifts substantially on month timescales, or the
low-wavenumber energy inside the aperture is environmental rather than
instrumental. This is not resolved here.

**And it is not the load-bearing question any more.** Even a complete answer to
"which box made the pedestal" would not change the outcome, because the binding
constraint is upstream of the instrument: the 2024-25 wavefield carries no
downgoing/upgoing asymmetry in the fan, in any of 240 windows scanned across the
archive. See the banner at the top of this file.

**The test that would settle it cannot be run.** `matched_earthquake_census.py`
was written to put an earthquake in *both* arms, which controls for signal
presence. It refuses to run, correctly: the 2017 arrivals are at 4.100 s and
4.916 s of 5.00 s records, leaving **0.90 s** of post-arrival data in total,
which is less than one 2.0 s analysis window. Concatenating the full records
instead would give a 2017 arm that is ~90 % pre-event noise against a 2024-25 arm
that is mostly earthquake — the same confound inverted. Shortening the window to
fit 0.9 s is not a fix: at 5 Hz that is under 5 cycles from a single window.

Note also that `faultzone/repeaters/cache_all/` **cannot** be used for any k = 0
measurement. `extract_all.py` reads through `DASutils.readFile_HDF` without
`median=False`, and that default subtracts the per-sample median across channels,
removing the k = 0 component before the cache is written. This is the confound
that forced the withdrawal of `cross_epoch_noise_floor.py`. Read raw HDF5 with
`h5py` instead, as the census and `matched_earthquake_census.py` do.

## 6. Scripts

| script | role |
|---|---|
| `ambient_fixed_k_test.py` | **section 1** — the low-k centroid comparison, the primary result |
| `lellouch2017_window_audit.py` | **section 4.2** — measures where the 2017 arrivals are |
| `ambient_apparent_velocity_census.py` | the energy budget; `CENSUS_FMIN`/`CENSUS_FMAX`/`K0_REMOVE` env switches, band-tagged outputs |
| `matched_earthquake_census.py` | **section 5** — refuses to run, with the measured reason |

`ambient_fixed_k_test.py` reports a **secondary** peak-per-frequency diagnostic
that returns no verdict, because excluding +-2 cells around k = 0 pushes its peak
out to ~0.019 cyc/m — the genuine surface-wave energy at 267-640 m/s, a different
feature. Its own permutation control catches this: "static" wins on 1578/2000
frequency-shuffled sets, i.e. what random labels give. That section not
concluding does **not** retract section 1; it uses a different quantity in a
different region of wavenumber.
