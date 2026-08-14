# C2 follow-up — does tube-wave amplitude image permeable structure?

> **Phase 0 ran on 2026-08-13. Results in `C2_PHASE0_RESULTS.md`.**
>
> The seven candidates are **retired** — the count was a *deficit* against chance
> (7 observed, 36 expected), the list flips to 42 candidates under a robust σ, and
> nothing survives an autocorrelation-preserving surrogate null.
>
> But Phase 0 also showed the measurement had **no power at all** — a planted 95 %
> amplitude loss is detected at most 34 % of the time — so C2 was never tested
> rather than tested and failed. And 96 % of the fatal scatter is *static*
> per-channel response (ρ = 0.922 between burst halves), hence removable, worth
> **5.0×**: with it divided out a step would be detectable at an **18 % amplitude
> loss**.
>
> **So Phase 1 gained one mandatory first item — calibrate the static per-channel
> amplitude response — and that item alone decides whether C2 can exist.** The rest
> of the plan below stands as written.

**Status: unclaimed result, not yet a finding.** Written 2026-08-13.

The gate test `tube_wave_gate.py` returned C2 PASS: seven wireline channels more
than 2σ below a log-amplitude trend, at 129, 131, 137, 696, 1280, 1295 and
1609 m. Tube waves lose energy at permeable fractures, so localised amplitude
loss is the classic hydrophone-VSP permeable-fracture indicator (method
reference: Banerjee & Chatterjee 2021, in the paper collection).

**Read the pass criterion before believing the result.** It was

```python
c2 = any(F['drops'].size >= 3 for F in fibers.values())
```

— "at least three channels below −2σ on either fiber." That is not a statistical
test. For Gaussian residuals the expected count below −2σ is ~2.3% of channels,
so with hundreds of channels three is guaranteed by chance. C2 as it stands
establishes that *some* channels are low, not that any of them mean anything.

Everything below assumes the answer might be no, and is ordered so the cheapest
thing that could kill it comes first.

---

## Phase 0 — is there anything there? (half a day)

Do this before any interpretation. Any one of these could retire C2.

- [ ] **Count against chance.** How many channels are in the fit, and how many
      would fall below −2σ if residuals were Gaussian? With ~1372 wireline
      channels the expectation is ~31. **Seven is fewer than chance.** If that
      holds, the residual distribution is not Gaussian and −2σ is the wrong
      threshold — establish the actual null distribution before calling anything
      a candidate.
- [ ] **Collapse correlated channels.** Gauge length is 10 m and channel spacing
      2.04 m, so ~5 adjacent channels are not independent. 129/131/137 is one
      feature, 1280/1295 probably another. The effective count is ~4 features,
      not 7, and the effective N for the null is channels ÷ ~5.
- [ ] **Check σ is not self-inflating.** σ was estimated from residuals that
      contain the features being detected, so a strong feature raises the
      threshold and hides itself. Re-estimate σ robustly (MAD, or from a
      feature-excluded subset) and see whether the candidate list changes.
- [ ] **Test for the right shape.** The physics is a *step*: below a permeable
      fracture the tube wave has lost energy and stays weaker. The gate test
      detected *dips* — isolated excursions — which is a different thing. Fit a
      step-detection statistic (e.g. cumulative-sum or a piecewise-constant fit
      on log-amplitude) and ask whether any candidate is a step rather than a
      notch.

**Decision point.** If the candidates do not survive a proper null and do not
look like steps, write two paragraphs retiring C2 and stop. That is a legitimate
outcome and cheaper than the alternative.

---

## Phase 1 — a defensible amplitude measurement (2–3 days)

Only if Phase 0 survives. The gate measurement was quick and has known
weaknesses; redo it properly.

- [ ] **FIRST, AND DECISIVE — divide out the static per-channel amplitude
      response.** Phase 0 measured 96 % of the wireline log-amplitude scatter as
      static (ρ = 0.922 between burst halves), σ_static = 1.177 against
      σ_noise = 0.241. Until this is removed no threshold on the profile can
      work; once it is, a step is detectable at an 18 % amplitude loss. Estimate
      the gain from an **independent** window — a different frequency band, a
      pre-arrival time window, or ambient RMS — so the calibration does not come
      from the samples being tested. Then re-run `c2_phase0_significance.py` and
      check σ falls from 1.20 toward 0.24. **If it does not, the static term is
      not a simple gain and C2 should be retired for good.** Everything else in
      this phase is wasted effort before this works.
- [ ] **Use the frozen trajectory.** The gate test fitted its own tube velocity
      per fiber (1330 / 1440 m/s). The analysis now has a properly frozen
      outbound trajectory — 1544.6 m/s at t₀ = +0.100 s, selected on the
      discovery half and validated against 499 permutations. Measure amplitude
      along that, so the amplitude profile inherits the trajectory's provenance.
- [ ] **Freeze the amplitude definition before looking.** Window length relative
      to the wavelet, RMS versus envelope peak, band (15–30 Hz to match the
      validated mode). Write it down first.
- [ ] **Separate instrument from medium.** Amplitude varies with channel
      sensitivity, coupling and gauge-length response, none of which is
      permeability. Estimate the channel-response envelope from a band or a
      time window where no guided mode is present, and divide it out — or show
      the candidates survive without that correction.
- [ ] **Replace the linear log-trend.** A single exponential is a crude model of
      tube-wave decay. Fit something defensible (piecewise, or a physically
      motivated decay) and state why.

---

## Phase 2 — controls, in this project's usual style (3–4 days)

- [ ] **Pre-register.** Candidate definition, threshold, null construction and
      the decision rule, dated, before evaluating. This project already does
      this well — `DEEP_DVV_PREREGISTRATION.md` is the template.
- [ ] **Split-sample.** Measure candidates on the odd bursts, confirm on the
      even. A permeable fracture is a property of the rock and must appear in
      both halves; a processing artefact need not.
- [ ] **Surrogate null.** Channel-order permutation, as used for the mode
      validation. A real localised feature is destroyed by scrambling channels;
      a heavy-tailed residual distribution is not.
- [ ] **Both legs.** The outbound and return limbs pass the *same* depths in
      opposite order. A genuine fracture should appear at the same depth on
      both. This is the strongest internal check available and costs nothing
      extra.
- [ ] **The cemented fiber as a contrast.** It found zero channels below −2σ.
      If the wireline candidates are fluid-path features, that asymmetry is
      expected and is evidence. If the cemented fiber shows the same depths,
      they are more likely structural or instrumental.

---

## Phase 3 — independent corroboration (depends on Ettore)

Amplitude loss alone will not convince anyone. Tie candidates to something
measured by a different instrument.

- [ ] **Get the borehole logs.** Any of caliper, sonic, resistivity, temperature
      or flow would do. A permeable fracture at 696 m should be visible in at
      least one. This needs the completion and logging data — the same
      conversation as the outstanding interrogator-model and turnaround-depth
      questions.
- [ ] **Compare against known structure.** Li et al. (2004) map a 150 m wide
      low-velocity waveguide on the Parkfield SAF; Ellsworth & Malin place a
      damage channel at 2.7 km. Do any candidates coincide? Coincidence is not
      proof but non-coincidence is informative.
- [ ] **Resolve the depth mapping first.** Candidates are quoted in *distance
      along fiber*. Comparing them to logged depths requires the
      coordinate-to-depth transform, which this project still records as
      provisional. Without it, no candidate can be placed against a log.

---

## Phase 4 — what could be claimed, and what could not

- [ ] **Draft the claim and its ceiling together.** The most that this method can
      support is *localised tube-wave amplitude loss at specific along-fiber
      positions, reproducible across bursts and legs*. Calling that
      "permeability" requires the fluid-mechanical link, which needs forward
      modelling not done here — the same discipline the paper applies to
      apparent velocity.
- [ ] **Decide where it goes.** It is not part of the sensitivity paper's story.
      Either a short separate note, or a supplement section explicitly labelled
      exploratory.

---

## Effort and honest expectation

| Phase | Effort | Could it end here? |
|---|---|---|
| 0 — is there anything there | half a day | **yes, likely** |
| 1 — defensible measurement | 2–3 days | yes |
| 2 — controls | 3–4 days | yes |
| 3 — corroboration | blocked on logs | yes |
| 4 — write-up | 1–2 days | — |

**My expectation was that Phase 0 would retire it**, because seven candidates out
of 1603 channels is fewer than Gaussian chance gives, the features cluster into
four once channel correlation is accounted for, and the test looked for dips where
the physics predicts steps.

**Half right.** All three held, and the candidate list is gone. But the same run
found the measurement never had the power to detect a permeable fracture at all,
and that 96 % of what destroyed it is static and removable — so the honest verdict
is *untestable as run*, with 5.0× of headroom and an 18 % detection threshold
waiting behind one calibration step. Phase 1's new first item is now the whole
decision, and it is a day of work rather than three.
