## 6. Limitations, and claims we withdrew

We report this section in full because a negative result is only as credible as its
error bars, and because several of the intermediate claims made during this work
turned out to be wrong. Recording which ones, and why, is part of the result.

### 6.1 Limitations of the illumination conclusion

- **Window length.** The scan uses 30 s windows, one per sampled record. A burst of
  illumination shorter than 30 s, or one falling between sampled records, could be
  missed. The sampling covers 240 windows spread across 377,680 records, so this is
  a sparse sample of the archive in time even though it spans its full extent.
- **Rank removal is conservative against detection.** Projecting out the leading
  spatial patterns is necessary (section 3.3) but would also remove a genuinely
  low-rank plane-wave arrival. We report the full rank sweep for this reason, and
  note that the 2017 arm remains significant through rank 2 under the identical
  operation, so the comparison is fair at matched ranks.
- **The fan is a frozen selection.** 2,500-4,000 m/s was fixed in advance. An
  arrival outside it would not be counted.
- **|A| is along-fibre.** It measures a preferred direction along the fibre, which
  in this near-vertical borehole is close to, but not exactly, vertical.
- **Band.** All conclusions are for 5-20 Hz, the band of the original paper.

### 6.2 The cross-epoch comparison is not perfectly matched

Our 2017 arm uses the pre-event portions of two earthquake records, about 5 s in
total, because that is the only raw 2017 noise in the public release. The
correlograms he released are, per the release documentation, "a stack of 7 different
one-day correlations", and the raw ambient data behind them is not available. **The arm we measure is
therefore not the same acquisition that produced his figure.**

This does not affect the measured contrast, which uses identical processing on both
arms and finds p = 0.0050 versus p = 0.7307. It does mean we cannot attribute the
difference to any specific cause. In particular:

- We have **no evidence** about what generated the 2017 illumination. An earlier
  draft of this work speculated that a manned field deployment supplied it, on the
  grounds that an interrogator was installed at SAFOD in June 2017. There is no site
  log, no operational record, and no statement in Lellouch et al. (2019)
  attributing the noise to any source -- that paper infers surface origin from the
  same asymmetry we measure, which is not independent corroboration of a cause.
  **That speculation is withdrawn and should not be cited.**
- Settling it would require SAFOD/USGS operational records for both windows, or the
  raw ambient data behind Figure 7c. Neither is in hand.

We also cannot separate "illumination absent in 2024-2025" from "illumination
present but below our detection threshold". The positive control bounds the
threshold usefully -- the statistic detects the effect in ~5 s of 2017 data -- but a
weaker illumination than 2017's would not necessarily be found.

### 6.3 Claims made during this work and subsequently withdrawn

Four automated verdicts were printed from broken inputs during this study, and
several substantive claims were retracted. All are listed because they shaped the
analysis path.

| claim | status | reason |
|---|---|---|
| "81-99 % of energy below 1,500 m/s implies the arrival is absent" | withdrawn | no geometric baseline; 98.4 % of in-band cells lie below 1,500 m/s by construction, so white noise scores the same |
| "downgoing never exceeds 49.9 %" | withdrawn | factually false; the same outputs contain 51.6 % |
| "2024-25 ambient decorrelates in 4 m versus 26 m in 2017" | withdrawn | the two arms were processed differently: the 2024-25 arm passed through a reader whose default subtracts the per-sample median across channels, removing the common mode from one arm only. The re-run product gives 3.0-3.6 m for 2024-25 and 400 m for 2017 in band, so the "26 m" of the withdrawn claim is reproduced by nothing now on disk |
| "a detectable coherent component at 5,850 m/s" | withdrawn | that peak is the scan ceiling; a moveout-free control reproduces the curve to 3.5 % |
| "tapering the f-k mask buys no discriminating power" | withdrawn | the real and synthetic paths were not matched, so the real/white ratio of 0.63 that the claim rested on is not a valid comparison. The corrected ratios quoted when the claim was withdrawn (1.51-2.01) are carried by no product now on disk and should not be requoted; the surviving statement is the one in the section 4.5 table, that the taper suppresses the ringing floor and the fan still fails the pre-filter channel scramble |
| "the contaminant and target are unresolved at this aperture, so Figure 7c is unachievable" | withdrawn in part | the aperture arithmetic stands, but by its own table the target is resolved above ~12 Hz, and the 12-20 Hz test showed *less* fan energy than 2017 rather than a recovered arrival |
| withdrawal of the cross-epoch k = 0 comparison | **reinstated** | the withdrawal assumed the 2017 window contained an earthquake; measurement placed the arrivals at 4.100 s and 4.916 s of 5.00 s, outside the 2.5 s window, so both arms are noise |
| "no surface illumination", first version | superseded | the reported \|A\| = 1e-4 was the algebraic k-symmetry of a separable static pattern, not a property of the wavefield; the corrected measurement gives 0.040 and the same conclusion |
| "the interrogator is to blame" | **restricted** | true in the surface lead-in, where the spatial pattern is stable across a year; not supported inside the analysed aperture, which is where it would matter |
| "a manned field deployment supplied the 2017 illumination" | withdrawn | speculation with no supporting record (section 6.2) |
| "our picker recovers Lellouch's published Figure 9 model at r = 0.948" | **restated** | r = 0.948 is corr(depth, velocity) on picks from his released traces, not agreement with his published curve, which is not digitised here. The validation stands in its corrected form (section 4.1); this project's own planning document carried the overstatement too |

Five of the first six share one cause: **two arms of a comparison received
different processing.** Every comparison in the final analysis therefore asserts at
runtime that both arms passed through an identical operation sequence and refuses
to report a comparison otherwise. Three of the four false verdicts share a second
cause: a script concluding from an input that had failed a precondition. Every
script now gates its verdict on the pedestal diagnostic, the peak location, causal
dominance, and its own control, and prints "uninformative" rather than a
conclusion when a gate fails.

We report this history because the withdrawn claims were, in several cases,
superficially more interesting than the surviving ones -- and because the surviving
result is a null, which is exactly the situation in which the temptation to keep an
attractive intermediate finding is strongest.

## 7. Conclusions

1. Lellouch et al. (2019) Figure 7c does not reproduce on a 2024-2025 continuous
   DAS archive from the same fibre in the same borehole. Six days give a minimum
   p-value of 0.1345, a coherent 96-hour stack gives p = 0.9184, and the observed
   moveout score clears the per-velocity null at 0 of 181 trial velocities.
2. The implementation is validated: the same picker returns the published
   2,416-4,357 m/s velocity range from the original authors' released correlograms
   (depth-velocity r = 0.948) while yielding velocities spanning 146 to 1.6e7 m/s
   on ours. Note that r = 0.948 is a depth-velocity correlation, not a point-by-point
   match to the published Figure 9 curve, which is not digitised here.
3. The 2024-2025 field is dominated by a static spatial pattern at fixed
   wavenumber, holding about 39 % of the in-band energy. Because it has no apparent
   velocity, no velocity-domain filter separates it, which accounts for the
   identical failure of eight such methods.
4. **More data does not help, and that is measured rather than asserted.** Stacking
   one contiguous day in 1 to 24 hourly chunks grows detectability as N^+0.042
   without common-mode removal and N^+0.019 with it, against N^+0.50 for a coherent
   arrival accumulating against incoherent noise. The common-mode-removed exponent
   is measured with the pedestal diagnostic at -0.219, so the coherent contaminant
   is suppressed and the curve is still flat, which points to an absent arrival
   rather than a masked one.
5. The binding constraint is **illumination**. The downgoing/upgoing asymmetry that
   the original authors used as evidence for surface sources is significant in their
   records (p = 0.0050) and absent from ours (p = 0.7307), and a pre-registered scan
   of 240 windows across one year finds 11 significant windows against 12.0 expected
   by chance.
6. Neither the interrogator change nor the gauge-length change accounts for the
   null: the stable instrumental pattern is confined to the surface lead-in, which
   the analysis already excludes, and the gauge-length response attenuates the
   target by 0.1-1.1 %.
7. Borehole ambient-noise body-wave interferometry should be treated as
   illumination-limited. We recommend measuring the asymmetry as a feasibility gate
   before committing acquisition or compute, and reporting stack convergence, since
   a result that degrades with added data indicates an absent arrival rather than an
   under-sampled one.
