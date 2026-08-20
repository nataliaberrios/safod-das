### 4.8 Illumination is the binding constraint (Figures 6 and 7)

{{FIGURE:fig6_illumination}}

Lellouch's own evidence for surface sources is the downgoing/upgoing asymmetry. We
measure it under identical processing in both epochs, sweeping the number of spatial
patterns projected out:

| rank removed | 2017 pre-event | 2024-2025 |
|---|---|---|
| 0 | 0.348, **p = 0.0050** | 0.040, p = 0.7307 |
| 1 | 0.329, **p = 0.0050** | 0.040, p = 0.4314 |
| 2 | 0.250, **p = 0.0299** | 0.036, p = 0.4638 |
| 4 | 0.123, p = 0.2718 | 0.001, p = 0.9875 |
| 8 | 0.081, p = 0.4489 | 0.011, p = 0.8279 |

The measurement **independently recovers the asymmetry Lellouch reported, from his
own data**, and finds none in ours at any rank. This is a working positive control,
which is what makes the null interpretable.

The data-quantity argument runs the informative way: the 2017 arm succeeds on
approximately 5 s of noise, ours fails on 30 s. Behm (2016) found 30 s sufficient
under adequate illumination. Non-detection here is therefore not an argument for
more data.

Because illumination by surface activity is intermittent, a single window could
miss it, and any illuminated minority of windows would be diluted by indiscriminate
stacking -- which is what all our earlier stacks did, and which would explain their
degradation with more data. We therefore scanned **240 windows spanning
2024-05-21 to 2025-05-06**, with the decision rule fixed before the run: illumination
is claimed only if the count of windows significant at alpha = 0.05 exceeds the 95th
percentile of Binomial(N, 0.05).

| | |
|---|---|
| windows measured | **240 of 240** sampled |
| \|A\| median / 90th / max | 0.0296 / 0.0753 / 0.2385 |
| windows with p < 0.05 | **11** |
| expected by chance | **12.0** (95th percentile = 18) |

{{FIGURE:fig7_archive_scan}}

Eleven hits against twelve expected. The archive contains no illuminated windows,
and no diurnal concentration appears that would indicate weak cultural sources.

### 4.9 The interrogator is not responsible inside the analysed aperture

Because the two acquisitions used different interrogators, an instrumental origin
for the static pattern is a natural hypothesis, and an instrumental term makes a
testable prediction: a fixed optical or electronic response has a fixed spatial
fingerprint, whereas an earth or site pattern varies with conditions. We correlated
the leading spatial pattern between six days spanning a year (2024-05-21, 08-29,
11-24, 2025-01-17, 03-13, 05-06), against a control of correlations between the
leading pattern of one day and the sub-leading patterns of another. The two channel
ranges come from two separate runs, and their controls use four and five sub-leading
patterns respectively, so the control percentiles are not interchangeable.

| channel range | cross-day median \|corr\| | control 95th | reading |
|---|---|---|---|
| 0-699 (lead-in included) | **0.8426** | 0.3889 (n = 120) | stable: instrumental |
| 23-708 (as analysed) | **0.0188** | 0.1282 (n = 150) | not stable |

The distinction is the whole result. There *is* a highly stable instrumental
pattern, but it lives in the surface lead-in. Its lead-in to deep power ratio
reaches 2 x 10^4 to 1.3 x 10^5 on four of the six days, with a median across days of
2.3 x 10^4; the remaining two days give 323 and 0.01, so the concentration in the
lead-in is strong but not uniform, and we quote the range with that qualification
rather than as a property of every record. Channel 23 is the wellhead, and the
analysis begins there, so that pattern is already excluded. Inside the analysed
aperture the dominant pattern is indistinguishable from unrelated patterns, so a
fixed instrumental fingerprint is not supported where it would matter.

Two limits on this section. The lead-in records real surface ground motion, so a
high lead-in to deep ratio is by itself not proof of an instrumental origin, and the
product that measures it says so. And the test is one file per day, so a pattern
that is stable within a day but drifts over weeks could read either way.

Gauge length is likewise ruled out: the `sinc(pi k L)` response attenuates a
3,200 m/s arrival by 0.1-1.1 % across the band for a 16.335 m gauge.
