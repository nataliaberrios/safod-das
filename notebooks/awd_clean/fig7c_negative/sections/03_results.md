## 4. Results

### 4.1 The implementation is correct (Figure 2)

{{FIGURE:fig2_method_validation}}

Before reporting a null we establish that the processing works. Lellouch's release
includes stacked correlograms for the constant-offset geometry. Applying our
picker -- a three-sample maximum with parabolic interpolation, following Nakata &
Snieder (2012) as his paper specifies -- to *his* correlograms yields a
monotonically increasing profile rising from **2,416 m/s at 50 m to 4,357 m/s at
700 m**, with a depth-velocity correlation of **r = 0.948**. That range and that
monotonicity are what his Figure 9 reports.

We are deliberately precise about what this does and does not show, because an
earlier draft of this work -- and this project's own planning document --
described the same number as our picker "recovering his published Figure 9
velocity model". It is not that. The 0.948 is the correlation between depth and
the velocities our picker returns from his traces; his published Figure 9 curve is
not digitised in our repository, so no point-by-point agreement with it is
established.

What matters for the null is the contrast, and it is stark. The *same picker,
unchanged*, returns a physically sensible profile from his data and nothing usable
from ours:

| | picks returned | velocity range | median |
|---|---|---|---|
| his released traces | 14 of 14 | 2,416 - 4,357 m/s | 3,268 m/s |
| our archive | 114 of 171 fibre positions | 146 to 1.6e7 m/s | 26,380 m/s |

Our failure is not that the picker returns nothing; it is that what it returns is
physically impossible -- five orders of magnitude of scatter about a median an
order of magnitude above any crustal P velocity. That is what a picker does when
handed a correlation with no arrival in it: it reports the position of whatever
noise maximum happens to fall in the search window. The picker, the geometry, the band and the moveout logic are
therefore not at fault. Anything that follows is a statement about the data.

### 4.2 Figure 7c does not reproduce (Figure 1)

{{FIGURE:fig1_no_reproduction}}

Applying the same processing to the 2024-2025 archive yields no arrival that
survives its controls.

| test | result |
|---|---|
| best single day, familywise p | **0.1345** (2024-11-30, 21.3 h block) |
| six independent days | none reach p < 0.05 |
| Fisher combination over five independent days | chi-sq = 9.08 on 10 df, **p = 0.524** |
| coherent four-day stack, 23,036 windows | **p = 0.9184** |
| per-velocity null, 181 trial velocities | **0 clear the 95th percentile** |
| at 3,200 m/s specifically | score 2.752 against a threshold of 2.811 |
| detectability of the four-day stack, peak / null 95th | **0.72** (needs 1.00) |

The Fisher combination is worth its own line: chi-sq approximately equal to its
degrees of freedom is exactly what independent pure-noise p-values produce, so
there is no residual signal sitting below per-day significance.

Every day's peak falls at 5,850-5,925 m/s, at the top of the scan, rather than at
the published 3,200 m/s, and the causal-to-acausal ratio at 3,200 m/s is below 1
on every day. The peak velocity is itself an artefact of the statistic rather than
a measurement: the score samples each trace in a gate at `t = offset/v`, so as `v`
rises the gates migrate toward zero lag where the pedestal lives.

**One configuration does reach p < 0.05, and it is instructive rather than
positive.** Stacking the two days that individually came closest, 2024-11-30 and
2024-12-20, returns p = 0.039. Those two days were selected *because* they had the
lowest p-values. Running all ten pairs of the five usable 500 Hz days shows the
selected pair to be the extreme of a distribution otherwise consistent with noise,
with pair p-values ranging up to 0.83. We report it because a reader who selected
the same way would find the same number, and because it is the precise shape of
the error this study exists to avoid.

The per-velocity result is the strongest single statement and it needs no
multiplicity correction: at no trial velocity, including the published 3,200 m/s,
does the observed score exceed what receiver-order permutation produces.

Two features of this table deserve emphasis. First, the coherent 96-hour stack is
*worse* than the best single day. Second, the baseline moveout statistic is almost
entirely pedestal: `corr(trial velocity, score) = +0.976`, rising monotonically
across the scan, which means the statistic is measuring proximity to the zero-lag
lobe rather than moveout. Common-mode removal fixes the diagnostic (`corr = -0.381`)
but does not produce a detection.

A statistic that improves with data is a weak signal. A statistic that degrades
with data is an absent one plus an accumulating artefact, and section 4.4 identifies
the artefact.

### 4.3 The constant-offset geometry isolates the failure (Figure 3)

{{FIGURE:fig3_fig7d_isolation}}

Figure 7c is a moveout gather, but the velocity model comes from Figure 7d, the
constant-offset geometry in which source and receiver are held 50 m apart and slid
down the array together. That geometry needs coherence over 50 m rather than the
full 700 m aperture, so it is the easier target and the sharper diagnostic.

In Lellouch's released 7d traces the correlation peak **migrates with depth**, from
20.7 ms at shallow depth to 11.5 ms at depth, exactly as `50 m / v(z)` requires for
a rising velocity profile.

In ours the peak sits at **exactly 0.0000 s at every depth**, and the parabolic
picker returns NaN because there is no off-zero maximum to interpolate. Removing the
common mode reduces the amplitude by a factor of about 120 and the peak stays at
zero lag.

A correlation that peaks at zero lag for every receiver pair, independent of
separation, is the signature of a component common to the channels rather than a
wave propagating between them. This observation motivated everything that follows.

### 4.4 The contaminant is a static pattern at fixed wavenumber (Figure 5)

{{FIGURE:fig5_static_pattern}}

Apparent velocity is `v = f / k`, so how a feature behaves as the analysis band
moves identifies what it is. A wave at fixed velocity must have its wavenumber
track frequency; a static spatial pattern sits at fixed wavenumber and its
*apparent* velocity therefore rises with frequency.

Measured on the first 16 non-zero wavenumber cells of a 700 m aperture, with k = 0
itself excluded so the answer is not fixed by construction:

| | 5-12 Hz | 12-20 Hz | ratio |
|---|---|---|---|
| band centre frequency | 8.50 Hz | 16.00 Hz | **1.882** |
| power-weighted \|k\| centroid | 0.00173 cyc/m | 0.00175 cyc/m | **1.011** |
| cells from k = 0 | 1.21 | 1.23 | |
| share of in-band energy | 39.33 % | 39.08 % | |

A wave predicts 1.882. A static pattern predicts 1.000. We observe **1.011**: the
centre frequency nearly doubled and the wavenumber moved by 1.1 %.

The measurement had the range to see a wave and did not. A 3,200 m/s arrival sits
at 1.86 cells at 8.5 Hz and 3.50 cells at 16.0 Hz, both inside the window measured,
so a dominant arrival would have shifted the centroid.

The same behaviour is visible in the coarser energy budget: raising the band moves
about 28 % of the in-band energy from the 3,000-4,000 m/s bin to the
6,000-10,000 m/s bin. Energy that changes apparent velocity when the band changes
is at fixed wavenumber, not fixed velocity.

### 4.5 Why eight velocity-domain methods fail identically (Figure 4)

{{FIGURE:fig4_eight_methods}}

A static pattern at fixed wavenumber has **no velocity**. Every method we applied
separates energy by velocity, so all of them are asking a question the contaminant
does not answer:

| method | outcome |
|---|---|
| fixed-fan f-k, brick-wall mask | apparent ridge; fails the pre-filter channel-scramble gate |
| fixed-fan f-k, raised-cosine taper | no discriminating power gained |
| tau-p slant stack, slowness mute swept to 8 km/s | pedestal survives every mute |
| rank-k subspace projection, global and windowed | pedestal not suppressed |
| phase cross-correlation | amplitude-blind; does not remove a phase-coherent pedestal |
| phase-weighted stacking | no recovery |
| offset-axis median flat-event removal | no recovery |
| median / mean common-mode removal | diagnostic improves, no detection |

This is one failure, not eight. Figure 4 shows the pedestal diagnostic for each.

The mean is the exact projection that annihilates k = 0, and it performs *worse*
than the robust median (`corr = +0.888` versus `-0.381`) because it is not robust to
the glitched channels that the 2017 release documentation itself warns require
trace editing: subtracting a mean contaminated by glitches injects them into every
channel. Order matters -- outliers must be interpolated before any exact projection.

We also note a resolution limit that constrains fixed-fan filtering specifically.
At a 700 m aperture, `dk = 1/700` and a Hann taper's main lobe spans +-2 cells, so a
3,200 m/s target lies 1.09 cells from k = 0 at 5 Hz and 4.37 cells at 20 Hz. Below
about 12 Hz the target is inside the main lobe of the zero-wavenumber peak and is
not resolvable from it at this aperture. Above 12 Hz it is resolved -- and we tested
that: in 12-20 Hz our fan energy is *lower* (0.42 %) than 2017's (1.15 %) and the
zero-wavenumber share remains 54.74 %. Better resolution exposed no arrival, which
is why resolution is a contributing limitation and not the explanation.

### 4.7 Illumination is the binding constraint (Figures 6 and 7)

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

### 4.8 The interrogator is not responsible inside the analysed aperture

Because the two acquisitions used different interrogators, an instrumental origin
for the static pattern is a natural hypothesis, and an instrumental term makes a
testable prediction: a fixed optical or electronic response has a fixed spatial
fingerprint, whereas an earth or site pattern varies with conditions. We correlated
the leading spatial pattern between six days spanning a year, against a control of
correlations between the leading pattern of one day and the second-to-fifth patterns
of another.

| channel range | cross-day median \|corr\| | control 95th | reading |
|---|---|---|---|
| 0-699 (lead-in included) | **0.8426** | 0.3889 | stable: instrumental |
| 23-708 (as analysed) | **0.0188** | 0.1282 | not stable |

The distinction is the whole result. There *is* a highly stable instrumental
pattern, but it lives in the surface lead-in: its lead-in to deep power ratio is
2 x 10^4 to 1.3 x 10^5. Channel 23 is the wellhead, and the analysis begins there,
so that pattern is already excluded. Inside the analysed aperture the dominant
pattern is indistinguishable from unrelated patterns, so a fixed instrumental
fingerprint is not supported where it would matter.

Gauge length is likewise ruled out: the `sinc(pi k L)` response attenuates a
3,200 m/s arrival by 0.1-1.1 % across the band for a 16.335 m gauge.
