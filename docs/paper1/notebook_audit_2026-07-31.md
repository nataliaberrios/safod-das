# Notebook audit — `JULY24_fml_coupling.ipynb`, 2026-07-31

Mechanical audit of the live analysis notebook against
[`STATUS.md`](STATUS.md). Every claim below was checked against cell sources and
stored outputs in the `.ipynb` on disk, not from memory or from prior notes.

Notebook state at audit time: 180 cells (128 code, 52 markdown), 98 executed,
30 empty, **1 cell with an error**.

---

## 1. Roadmap item #3 is already complete

The roadmap asks to bump `TW_N_NULL` from 20 to 199 because "the 20-null
p=0.048 is the 1/21 floor, not a measured p."

That bump has already happened:

| Cell | Variable | Value |
|---|---|---|
| 158 (exec 86) | `TW_N_NULL` | **199** |
| 162 (exec 90) | `N_CONFIRM_NULLS` | **199** |

The summary in cell 173 reports the result consistent with 199 nulls: observed
semblance 0.326 against a shuffled-depth 99th percentile of 0.084, with none of
199 realizations reaching the observed value, giving `p_perm = 0.005`. That is
`1/200`, the correct floor for 199 nulls — not the old `1/21`.

**Action: strike item #3 from the roadmap.** The permutation screen stands as
published-quality.

## 2. Blocker #1 is half-complete

The roadmap names two cells to delete before the clean re-run.

| Target | Status |
|---|---|
| Cell starting `# Earth tide detection sensitivity using channel-specific timing noise` (quarantined ch-374 / 0.0000 ms zero-scatter) | **Already gone.** Not present in `JULY24`, and not in `JULY23` either — removed at or before the July 23 revision. |
| Cell starting `# Per-epoch OLS: dt(z) = alpha - eps * T(z)` with `SEP_ZMAX_M = 130.0` | **Still present** — cell 132, exec 68. |

Cell 132 already carries the header comment `# DEMO / DO NOT QUOTE:
shallow-aperture (0-130 m) slope fit.`, so the short-lever σ_α = 4.57 ms is
already fenced off in-place rather than circulating as a competing result.

**Action:** delete cell 132 (or convert to markdown) and the two-σ_α ambiguity is
closed. This is a one-cell edit, not a re-run.

## 3. The 48-vs-49 epoch drift is a labelling artifact, plus one real bug

The true epoch count is **49** (indices 0–48). Evidence:

```
cell  21 [exec 11]: reloaded 49 epochs, 49 begtimes (match: True)
cell 159 [exec 87]: Tracking candidate through 49 epochs...
cell 159 [exec 87]: epoch 00/48 ... epoch 48/48        <- 0-indexed display of 49
cell 163 [exec 91]: Filtering 49 epochs x 89 channels...
cell 175 [exec 97]: clean epochs: 47 of 49
```

The `48`s fall into three categories:

- **Superseded estimate.** Cell 15 prints `48 epochs estimated from file count
  (6 files each); true count finalized in Step 4`. Self-labelled as provisional;
  harmless.
- **Cosmetic 0-indexing.** `epoch 00/48 … epoch 48/48` is 49 epochs displayed
  from zero. Correct as written.
- **Stale prose.** The summary markdown (cell 173) says "48 half-hour epochs"
  and "46 of 48 clean epochs". The analysis that produced the headline numbers
  actually ran on 49, with **47 of 49** clean. Fix the text.

There is one genuine defect:

```
cell 155 [exec 84]: epoch_stacks has 49 epochs; contaminated has 49; t_hours has 48
```

`t_hours` is one element short of the data arrays. Any plot or fit pairing
`t_hours` against a 49-length series is either truncating an epoch or
misaligning the time axis by one half-hour slot. **This is the only finding in
this audit that can move a number**, and it should be fixed before the locked
table is rebuilt.

## 4. The single notebook error is dead code

Cell 175 (exec 97) — the trend-sensitivity table and null-vs-data figure for the
tidal non-detection — raises:

```
ValueError: matmul: Input operand 1 has a mismatch in its core dimension 0 ...
(size 1 is different from 2)
```

at this line:

```python
dfz = c_tr[0] + (Glin[:,1:]@c_tr[1:] - c_tr[0]) * 0  # baseline
```

`Glin` carries the trend columns *plus* the tide column; `c_tr` holds
trend-only coefficients, so the slice widths disagree. But the entire expression
is multiplied by `0`, and `dfz` is never read again anywhere in the cell — it is
assigned once and dropped.

**The line is dead code. Deleting it fixes the cell.** No result depends on it;
the De Fazio reference curve on the next line is computed independently as
`Xtr@c_tr + DEFAZIO*tau`.

This matters more than its size suggests: cell 175 is the robustness check
behind the paper's headline null (95% UL = 7.98×10⁻⁴). It is the last executed
cell in the notebook, so the notebook currently ends on a failure at exactly the
result the paper leads with.

## 5. Internal inconsistency in the depth claims

Cell 115 (markdown) states:

> CC of 0.027 at ch 500 (633m) confirms the AWD signal does not penetrate below
> **~100m**. All δt measurements at deep channels are measuring noise.

This contradicts the rest of the notebook, which repeatedly establishes a usable
interval extending to **373–374 m**: cell 39 (median SNR 22.5 dB, CC 0.987,
σ_t 0.434 ms over 127–374 m), cell 58 (within-burst CC 0.938 over the same
range), and the cell 173 summary (principal continuous usable interval
32–373 m).

The evidence cited in cell 115 — CC 0.027 at 633 m — supports "no usable signal
at 633 m", not "no penetration below 100 m". The `~100m` figure appears to be
stale text from an earlier aperture. A reviewer comparing §1 and §6c will catch
this immediately.

**Action:** restate cell 115 as "signal is incoherent by 633 m, consistent with
the 373 m usable limit established in §1."

## 6. Aperture reconciliation (roadmap item #2) — unchanged

Still genuinely open. Three apertures coexist:

| Aperture | Where it comes from |
|---|---|
| 127–374 m | the dv/v chain and the within-burst repeatability analysis |
| 32–373 m | the four-criteria gate as reported in the cell 173 summary |
| 76–373 m | the four-criteria gate with the `SIGMA_T_MAX` gate, per the roadmap |

The roadmap's recommendation — adopt the four-criteria aperture and make the
dv/v chain agree with it — still stands. Nothing in this audit changes that.

---

## Revised critical path

The roadmap's blocker is substantially smaller than assumed. Ordered by
effort-to-value:

1. Delete one dead line in cell 175 → restores the headline-null robustness
   figure. **Minutes.**
2. Fix `t_hours` length 48 → 49 in cell 155's upstream construction. **The only
   change that can move a number.**
3. Delete cell 132 → closes Blocker #1.
4. Correct stale prose: 48→49 epochs and 46→47 clean in cell 173; the `~100m`
   claim in cell 115.
5. Reconcile the aperture (roadmap #2) — still the largest genuine analysis task.
6. Re-run clean, rebuild the locked table.

Items 1–4 are edits, not re-runs. Only item 2 can change a reported value.
