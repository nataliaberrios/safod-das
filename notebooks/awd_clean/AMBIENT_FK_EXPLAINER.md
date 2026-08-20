# Why the ambient F–K fan produced a ridge that was not there

**Explainer, written 2026-08-20. NOT a status document and NOT authoritative.**

This file explains *why* the 2.5–4.5 km/s F–K fan produced a convincing but false
ridge. It introduces no new results and **overrules nothing**. Every number in it
is quoted from a product or a status doc, cited inline. Where this file and any
status doc disagree, **the status doc wins**:

- `Ambient_FK_QC_workflow.ipynb` (v11) — authoritative on the F–K fan verdict
- `AMBIENT_LOWK_MECHANISM.md` — authoritative on the fixed-*k* mechanism
- `AMBIENT_CC_LITERATURE_REVIEW.md` — authoritative on illumination and adaptive F–K
- `FIG7C_MULTIDAY_RESULT.md` — authoritative on the multi-day numbers
- `AMBIENT_FIG7C_STATUS.md` (+ `_ADDENDUM`) — authoritative on the single-day operator

Scope: SAFOD 2024–25 ambient DAS, 5–20 Hz, attempting Lellouch et al. (2019)
Figure 7c. All archive dates below are **UTC**; Parkfield local is UTC−8 (PST) in
winter and UTC−7 (PDT) in summer. Stack sizes are given as duration, since the
1-minute segment length is a property of this archive and does not transfer.

---

## 1. What the filter does

A 2-D FFT maps the record to frequency × wavenumber. Apparent velocity is
`v = f/k`, so a constant velocity is a **straight line through the origin** and a
velocity band is a **wedge**. From `ambient_fk_transfer_test.py:10`:

```python
K,F = np.meshgrid(k,f,indexing="ij"); v = |F|/max(|K|,1e-12)
mask = (|F|>=5)&(|F|<=20)&(v>=2500)&(v<=4500)&(|K|>0)
if mode=="positive": mask &= (F*K)>0     # one propagation direction
if mode=="negative": mask &= (F*K)<0     # the other
y = ifft2(fft2(xd)*mask).real
```

Three properties drive everything below:

1. **Hard binary mask** — no taper, brick-wall passband.
2. **`|K|>0` excludes exactly the k = 0 column** and nothing else near it.
3. **The wedge is defined by velocity**, so everything surviving it has
   `f/k ∈ [2500, 4500]` *by construction*.

`velocity_scores` then sweeps trial velocities 1,200–6,000 m/s, sampling the
correlation panel along `lag = distance/v`, and reports the maximum.

---

## 2. The QC programme

The advisor's 8-point spec (2026-08-13 23:44 UTC) became Sections 2–8 of the QC
notebook. `Ambient_FK_QC_workflow.ipynb` is the record; this is an index.

| § | Test | What it randomises / probes | Verdict |
|---|---|---|---|
| 4 | Duration-matched fan comparison | only the F–K mask changes | ridge appears **only** with the fan |
| 4a | Decimation / anti-alias | sampling artefacts | not the cause |
| 5 | Broadband white-noise input | order fabricated from IID noise | **passed** — no false ridge |
| 6 | **Pre-filter channel scramble** | destroys spatial order *before* F–K | **FAILED** — the deciding gate |
| 6 | Independent circular time shifts | destroys phase, preserves each spectrum | **FAILED** |
| 7 | "Is F–K required?" | visual vs scientific claim | visually yes, scientifically unresolved → rejected |
| 8 | Stack convergence, 5 h → 8 days | √N growth expected | **N^+0.042** vs theory **N^+0.50** |
| 13.1 | Adaptive F–K (Isken et al. 2022) | the one filter family never tried | 5 configs, none passed predeclared gates |

**The white-noise null passed.** That is why the programme did not look broken
early: the pipeline genuinely does not manufacture a ridge from unstructured
input. It manufactures one from *this* input, which is structured but is not a
wavefield.

---

## 3. Why it was convincing

Across **1,170 independent 10-minute windows** (`fk_controls_*.json`):

| branch | n | modal peak velocity | median p | windows p < 0.05 |
|---|---:|---:|---:|---:|
| negative | 1,170 | **3,075 m/s** (723 windows) | **0.002** | **100.0 %** |
| both | 1,155 | 2,850 m/s | 0.002 | 98.6 % |
| positive | 1,155 | 3,775 m/s | 0.966 | 0.0 % |

Four independent-looking reasons to believe it:

1. **Reproducible** — every one of 1,170 windows significant.
2. **Physically sensible velocity** — 3,075 m/s is a plausible formation Vp.
3. **Agrees with the literature** — within 4 % of Lellouch's published 3.2 km/s.
4. **An apparent directional control** — the opposite branch stayed silent at
   p = 0.97, matching his downgoing-only claim.

---

## 4. The null was in the wrong place

`ambient_fk_transfer_test.py:52`:

```python
null = [max(velocity_scores(top, lags, rng.permutation(dist))[1]) for _ in range(500)]
```

This permutes **receiver-distance labels on the finished, already-filtered
correlation panel**. It asks *"does this panel contain ordered moveout?"* — not
*"did the moveout come from the data or from the filter?"*

The filter runs upstream of the randomisation, so **the filter sits outside its
own null.** Per `Ambient_FK_QC_workflow.ipynb` §13.1:

> Any pre-correlation operator in this tree owes an input-level null, not merely
> a gather-level one.

### Move the randomisation upstream and the result inverts

Production pre-filter null, **5.0 hours** of 2024-12-20 UTC (local start
2024-12-19 16:00 PST), 20 realizations
(`fk_full_pipeline_null_v2_n300_r20/…_aggregate.json`):

| Null placement | Observed | Null 95th | p |
|---|---:|---:|---:|
| **After** the filter (distance permutation) | 0.222 | 0.098 | **0.002** |
| **Before** the filter (channel permutation) | 0.283 | **0.365** | **1.00** |
| **Before** the filter (circular time shift) | 0.283 | 0.333 | **1.00** |

The scrambled data scored **higher than the real data**, on both null families.

The diagnostic quantity is the null's own 95th percentile: it rises from ~0.10 to
~0.36, a factor of **3.7**, purely by letting the surrogate pass through the
filter. That increase *is* the structure the fan manufactures.

**Caveat, from the notebook itself (§7):** these two products use different
stacking pipelines and window counts and "cannot be combined into one acceptance
decision." The pre-filter result at p = 1.00 stands on its own regardless.

Note also that in the 5.0-hour run **both branches peak at 3,075 m/s with
near-identical scores** (0.279 / 0.283). The directional asymmetry of §3 does not
survive the larger stack, so reason 4 above never actually held.

---

## 5. The mechanism — the input has no velocity

`AMBIENT_LOWK_MECHANISM.md` §1. Since `v = f/k`, doubling frequency forces a
**wave** to double its wavenumber; a **static spatial pattern** holds `k` fixed:

| | 5–12 Hz | 12–20 Hz | ratio |
|---|---|---|---|
| band centre frequency | 8.50 Hz | 16.00 Hz | **1.882** |
| power-weighted \|k\| centroid | 0.00173 | 0.00175 cyc/m | **1.011** |

Frequency nearly doubled; wavenumber moved **1.1 %**. The dominant feature is at
fixed wavenumber — a **static spatial pattern, not a propagating wave**. It holds
~39 % of in-band energy; k = 0 holds **~66 %**, against **0.13 %** in Lellouch's
2017 records.

The circularity is then complete:

1. A static pattern is a near-vertical line at small `k`, broadband in `f`.
2. The wedge `2500 ≤ f/k ≤ 4500` slices that line and keeps part of it.
3. Whatever emerges has apparent velocity **inside the passband by construction**.
4. Scanning 1,200–6,000 m/s then finds its peak inside 2,500–4,500. Inevitably.

The filter did not *detect* a 3.1 km/s arrival. It *assigned* one.

This also explains why **eight** methods failed identically — F–K brick-wall and
tapered, tau-p/Radon, rank-*k* SVD, PCC/PWS, flat-event removal, common-mode
(table in `AMBIENT_LOWK_MECHANISM.md` §2). Every one separates energy *by
velocity*, and the contaminant has none. Tapering the wedge cut ringing 229× and
changed nothing (`8146dee`), because ringing was never the mechanism.

**One synthesis of mine, flagged as such:** the aperture table in
`AMBIENT_LOWK_MECHANISM.md` §4.1 gives `dk = 1/700 m⁻¹`, so a 3,200 m/s target at
5 Hz sits **1.09 cells from k = 0** — inside the leakage skirt of a pedestal
carrying ~66 % of the energy. The mask excludes the k = 0 column and passes cells
1, 2, 3…, where that leakage lives. Both facts are in the docs; joining them is
my reading, not an established result.

---

## 6. What actually explains the non-reproduction

Strict precedence — **do not cite (2) alone**:

1. **Illumination is the binding constraint** (QC §12). Lellouch's downgoing P is
   an *inference from a measured downgoing/upgoing asymmetry*. Same measurement,
   matched spatial rank: 2017 gives \|A\| = 0.348, p = 0.0050; 2024–25 gives
   \|A\| = 0.040, p = 0.7307. An archive-wide scan of **240 windows** spanning
   2024-05-21 → 2025-05-06 UTC (local 2024-05-20 17:00 PDT → 2025-05-05 17:00
   PDT), decision rule fixed before the run, found **11 significant against 12.0
   expected by chance** (Binomial 95th pct 18) — below chance. No filter can
   create a propagation direction absent from the recording.
2. **A static fixed-wavenumber pattern** (`AMBIENT_LOWK_MECHANISM.md`) explains
   the *processing* failures, not the absence of the arrival.

Corroborating, none requiring a multiplicity correction: six complete days give
min p = 0.1345, Fisher's p = 0.524; the coherent **96.0-hour**, 23,036-window
stack gives p = 0.9184; **0 of 181 velocities** clear the per-velocity null
(3,200 m/s: 2.752 against threshold 2.811); detectability grows as N^+0.042 where
a real arrival needs N^+0.50. Per Behm (2016), 30 s suffices under adequate
illumination — **a stack that fails to converge is the signature of an absent
signal, not a weak one**.

Ruled out rather than left open: gauge length (0.1–1.1 % at 3,200 m/s); the
interrogator inside the analysed aperture (stable only in channels 0–22, which
the pipeline excludes); day selection (the richest day scored *worst*).

---

## 7. What this does not say

From `Ambient_FK_QC_workflow.ipynb` §14:

- Not a criticism of Lellouch et al. (2019).
- Not a claim that ambient interferometry fails in general.
- **Not a claim that F–K filtering is invalid.** It is a standard, correct
  directional tool being asked a velocity question about an object that has no
  velocity.
- Not a claim that no arrival exists — a weak arrival can hide beneath a feature
  holding ~39 % of in-band energy.
- A statement about **this archive, this band (5–20 Hz), and these days**.

The velocity model remains reachable by another route: 206 cached earthquakes,
and G0 already matched the 2005 check shot to 0.2 %.

---

## 8. The transferable lesson

**Randomise at the input, not at the output.** A null that permutes the finished
gather cannot see any operator applied before the gather was formed. Section 13.1
is a second, independent instance of the same failure: adaptive F–K with no prior
static removal produced a spurious p = 0.0060 by amplifying the pedestal 6.6-fold,
while the cleanest configuration reached p = 0.9392.

Both times the filter passed the gather-level null and failed the input-level one.
