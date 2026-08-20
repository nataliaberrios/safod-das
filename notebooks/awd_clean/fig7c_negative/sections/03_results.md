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
| our archive | 114 of 171 fibre positions | 146 to 1.6e7 m/s | 26,382 m/s |

Our failure is not that the picker returns nothing; it is that what it returns is
physically impossible -- five orders of magnitude of scatter about a median an
order of magnitude above any crustal P velocity. That is what a picker does when
handed a correlation with no arrival in it: it reports the position of whatever
noise maximum happens to fall in the search window. The picker, the geometry, the band and the moveout logic are
therefore not at fault. Anything that follows is a statement about the data.

Two points of precision about that range, so it is not read as stronger than it
is. The extremes come from the parabolic sub-sample correction, which is unbounded:
the argmax is constrained to the 5-45 ms search window, but the interpolated offset
added to it is not, so 146 m/s corresponds to an extrapolated 0.343 s and 1.6e7 m/s
to 3.15e-6 s. The distribution-level statement is the robust one, and it says the
same thing: the interquartile range of our finite picks is 7,614 to 85,352 m/s, and
only 10 of the 114 fall inside the 1,500-6,000 m/s scan the rest of this study uses.

### 4.2 Figure 7c does not reproduce (Figure 1)

{{FIGURE:fig1_no_reproduction}}

Applying the same processing to the 2024-2025 archive yields no arrival that
survives its controls.

| test | result |
|---|---|
| best single day, familywise p | **0.1345** (2024-11-30, 21.3 h block) |
| six days measured | none reach p < 0.05 |
| Fisher combination over the five independent days | chi-sq = 9.08 on 10 df, **p = 0.524** |
| coherent four-day stack, 23,036 windows | **p = 0.9184** |
| per-velocity null, 181 trial velocities | **0 clear the 95th percentile** |
| at 3,200 m/s specifically | score 2.752 against a threshold of 2.811 |
| detectability of the four-day stack, peak / null 95th | **0.72** (needs 1.00) |

The Fisher combination is worth its own line: chi-sq approximately equal to its
degrees of freedom is exactly what independent pure-noise p-values produce, so
there is no residual signal sitting below per-day significance.

On each of the five 500 Hz days the peak falls at 5,850-5,925 m/s, at the top of
the scan, rather than at the published 3,200 m/s. (The sixth day, 2024-05-11, was
recorded at 5,000 Hz and is excluded from the stack; its peak is at 1,675 m/s, at
the other end of the scan.) The peak velocity is itself an artefact of the statistic
rather than a measurement: the score samples each trace in a gate at
`t = offset/v`, so as `v` rises the gates migrate toward zero lag where the pedestal
lives. The causal-to-acausal ratio at 3,200 m/s never exceeds 1.13 on any day and
is below 1 on the three days closest to significance, so there is no day on which
the causal side carries the excess a downgoing arrival would require.

**Two pairs of days do reach p < 0.05, and they are instructive rather than
positive.** Stacking the two days that individually came closest, 2024-11-30 and
2024-12-20, returns p = 0.0390; the pair 2024-12-20 + 2025-02-24 returns p = 0.0165.
Over all ten pairs of the five usable 500 Hz days, 2 reach alpha = 0.05 where 0.5
are expected by chance, and pair p-values run up to 0.8251. **Multiplicity is
therefore not the objection -- velocity is.** Both pairs that reach p < 0.05 peak
at 5,850 m/s, with causal/acausal ratios at 3,200 m/s of 0.97 and 1.01: the scan
ceiling and no causal preference, which is the pedestal signature of section 3.2
rather than an arrival. We report these pairs because a reader who selected the same way
would find the same numbers, and because a p-value that survives only while its
peak sits at the edge of the scan is the precise shape of the error this study
exists to avoid.

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
the artefact. Section 4.7 measures the convergence directly, on a single day at
fixed geometry, and finds it flat.

### 4.3 The constant-offset geometry isolates the failure (Figure 3)

{{FIGURE:fig3_fig7d_isolation}}

Figure 7c is a moveout gather, but the velocity model comes from Figure 7d, the
constant-offset geometry in which source and receiver are held 50 m apart and slid
down the array together. That geometry needs coherence over 50 m rather than the
full 700 m aperture, so it is the easier target and the sharper diagnostic.

In Lellouch's released 7d traces the correlation peak **migrates with depth**, from
20.7 ms at shallow depth to 11.5 ms at depth, exactly as `50 m / v(z)` requires for
a rising velocity profile.

In ours -- with the common mode already removed, which is the branch this product
was computed in and the branch Figure 3(a) shows -- the peak sits at **exactly
0.0000 s at 13 of the 14 depths**, and the parabolic picker returns NaN at those 13
because there is no off-zero maximum to interpolate. The single exception is the
shallowest, 50 m pair: its global maximum lies at +0.114 s, and inside the picker's
5-45 ms window it returns 1,010 m/s, an implied lag of 49.5 ms that has been pushed
past the window edge by the unbounded parabolic correction and is in any case more
than twice the 20.7 ms his 50 m pair requires. Neither lag is a candidate travel
time. Removing the common mode therefore does not move the peak off zero lag; it
only lowers its amplitude.

A correlation that peaks at zero lag for essentially every receiver pair,
independent of
separation, is the signature of a component common to the channels rather than a
wave propagating between them. This observation motivated everything that follows.

### 4.4 The contaminant is a static pattern at fixed wavenumber (Figure 5)

{{FIGURE:fig5_static_pattern}}

Apparent velocity is `v = f / k`, so how a feature behaves as the analysis band
moves identifies what it is. A wave at fixed velocity must have its wavenumber
track frequency; a static spatial pattern sits at fixed wavenumber and its
*apparent* velocity therefore rises with frequency.

Measured on the first 16 non-zero wavenumber cells of a 700 m aperture -- `+-1` to
`+-8` on the two-sided wavenumber axis -- with k = 0 itself excluded so the answer
is not fixed by construction:

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

**What this measurement does not show**, stated because the product that produced it
states it. It identifies the *dominant* low-wavenumber feature as a static pattern.
It does not show that no arrival exists: an arrival weak relative to a feature
holding about 39 % of the in-band energy would not move the centroid. The
implication is about method rather than absence -- because the contaminant is a
fixed *spatial* pattern, removing it has to be a spatial operation, an estimate of
the static per-channel response divided out, rather than a velocity filter. The
evidence for absence is sections 4.7 and 4.8, not this one.

### 4.5 Why eight velocity-domain methods fail identically (Figure 4)

{{FIGURE:fig4_eight_methods}}

A static pattern at fixed wavenumber has **no velocity**. Every method we applied
separates energy by velocity, so all of them are asking a question the contaminant
does not answer:

| method | outcome |
|---|---|
| fixed-fan f-k, brick-wall mask | apparent ridge; fails the pre-filter channel-scramble gate |
| fixed-fan f-k, raised-cosine taper | ringing floor suppressed; the fan still fails the pre-filter channel scramble |
| tau-p slant stack, 6 km/s slowness mute | pedestal worsens, +0.740 before the mute to +0.940 after |
| rank-k subspace projection, global and windowed | pedestal not suppressed |
| phase cross-correlation | amplitude-blind; does not remove a phase-coherent pedestal |
| phase-weighted stacking | no recovery |
| offset-axis median flat-event removal | no recovery |
| median / mean common-mode removal | diagnostic improves, no detection |

This is one failure, not eight. Figure 4 shows the pedestal diagnostic for the five
of these that have one in a saved product, alongside the untreated baseline; the two
fixed-fan mask families and rank-k are characterised in their own products instead.

Two rows carry a provenance caveat and are marked here rather than left to be
found. The rank-k outcome is recorded in a source comment
(`ambient_lellouch2019_exact_stack.py:284`, `corr` stays at +0.96 to +0.98 for ranks
1, 2, 4 and 8) and not in a standalone product. The tau-p mute is a single declared
value, not a sweep: an earlier version of this table claimed a sweep to 8 km/s,
which was the script's default rather than the value the run used.

The mean is the exact projection that annihilates k = 0, and it performs *worse*
than the robust median (`corr = +0.888` versus `-0.381`) because it is not robust to
the glitched channels that the 2017 release documentation itself warns require
trace editing: subtracting a mean contaminated by glitches injects them into every
channel. Order matters -- outliers must be interpolated before any exact projection.
The `+0.888` is a third provenance caveat: it survives only as a measurement
recorded in a source comment (`ambient_apparent_velocity_census.py:136`) and is not
carried by any saved product, so it cannot be recomputed from what is on disk.

We also note a resolution limit that constrains fixed-fan filtering specifically.
At a 700 m aperture, `dk = 1/700` and a Hann taper's main lobe spans +-2 cells, so a
3,200 m/s target lies 1.09 cells from k = 0 at 5 Hz and 4.38 cells at 20 Hz. Below
about 12 Hz the target is inside the main lobe of the zero-wavenumber peak and is
not resolvable from it at this aperture. Above 12 Hz it is resolved -- and we tested
that: in 12-20 Hz our fan energy is *lower* (0.42 %) than 2017's (1.15 %) and the
zero-wavenumber share remains 54.74 %. Better resolution exposed no arrival, which
is why resolution is a contributing limitation and not the explanation.
