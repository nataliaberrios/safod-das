## 2. Data

### 2.1 The two acquisitions

Both datasets come from the same fibre in the SAFOD main hole. Lellouch et al.
(2019) report a fibre under approximately 1 N tension in a 0.9 mm steel tube
cemented between casing strings, reaching 864 m, with analysis limited to 800 m by
a loop failure at the fibre end.

| | Lellouch et al. (2019) | This study |
|---|---|---|
| recorded | June 2017 | 2024-05-21 to 2025-05-06 |
| interrogator | OptaSense ODH3.1 | `SAFOD-Zumberge-16mGL-500HzFs` |
| channel spacing | 1.0 m | 1.0209 m |
| sample rate | 250 Hz | 500 Hz |
| gauge length | 10 m | 16.335 m |
| recorded quantity | — | optical phase, `rad * 2PI/2^16` |
| available to us | 2 earthquake records (800 x 1250, 250 Hz) plus stacked correlograms | 377,680 continuous 60 s records |

The acquisitions differ in interrogator and gauge length. Channel spacing is
effectively identical. We quantify the gauge-length difference in section 3.4 and
find it far too small to matter: the `sinc(pi k L)` response attenuates a 3,200 m/s
arrival by 0.1-1.1 % across 5-20 Hz.

An important asymmetry in what is available: Lellouch's public release contains raw
data only for two earthquake records (`M1p33`, `M2p46`) plus **stacked**
correlograms, which the release documentation describes as "a stack of 7 different
one-day correlations". The raw ambient data behind his Figure 7c is not public.
Our 2017 measurements therefore use the pre-event portions of the two earthquake
records, and we state explicitly where that limits the comparison (section 6.2).

### 2.2 Verifying that the 2017 comparison windows are noise

Because `M1p33` and `M2p46` are earthquake records, using them as an ambient proxy
requires knowing that the analysed portion precedes the arrival; otherwise a
comparison against our ambient data would contrast a loud coherent arrival with
noise. We picked the arrivals rather than assuming: a short-term-average trigger at
3x the quiet-quarter level places them at **4.100 s** and **4.916 s** of 5.00 s
records. Our comparison window is the first 2.5 s, so both arrivals fall outside
it and both arms of every cross-epoch comparison below are noise.

This check was not incidental. We initially withdrew our central cross-epoch
result on the suspicion that this window contained an earthquake, and the
measurement showed the withdrawal was unnecessary. Section 6.3 records that.

## 3. Methods

### 3.1 Faithful reproduction

We implement only the operations reported in Lellouch et al. (2019) section 4.1:
time differentiation of recorded phase to a strain-rate proxy; running-absolute-mean
temporal normalisation; 30 s correlation windows with 15 s overlap over a
continuous day; a fixed shallow virtual source with receiver centres every 50 m
from 50 to 700 m; the published neighbour sum `C_{S,R} = sum_{Z=R-10}^{R+10} C_{S,Z}`;
simple unshifted stacking; and a 5-20 Hz bandpass applied to the stacked
correlations. The bandpass is applied to the full correlation and only then cropped
in lag, since the reverse order aliases energy into the retained window.

Cross-spectra are accumulated in the frequency domain with
`n_fft = 2^ceil(log2(2W-1))`, so the correlation is linear rather than circular.
The paper does not report a running-absolute-mean duration; we use 0.1 s following
Bensen et al. (2007), who recommend approximately half the maximum period of the
analysis band, and we treat it as a declared sensitivity rather than a recovered
constant. Common-mode removal and source-power stabilisation are available only as
labelled sensitivity branches, and f-k filtering is absent from the baseline
because it is not reported in the paper's ambient-interferometry section.

Day-level products are stored as summed cross-spectra in resumable chunks, which
allows exact full-day aggregation and convergence analysis without recomputing
correlations. We verified that an independent re-sum of one day's 24 chunks
reproduces the stored aggregate to a maximum absolute difference of 0.000e+00, and
that window counts match expectation (5,759 for a contiguous day).

### 3.2 Statistics and controls

Our decision statistic is the maximum envelope moveout score over a declared
1,500-6,000 m/s scan. The null is a **receiver-order permutation**: offsets are
held fixed and receiver identities permuted, so each null realisation carries the
same geometry as the observation, and the velocity selection is repeated
independently in every realisation to account for choosing a velocity after
inspecting the published figure.

Two additional controls matter. First, we report a **per-velocity** null in
addition to the familywise one, since a familywise maximum can hide the fact that
no individual velocity is anomalous. Second, we report a **pedestal diagnostic**,
`corr(trial velocity, moveout score)`. A moveout score that rises monotonically
with trial velocity is measuring proximity to the zero-lag lobe rather than
moveout, and its p-value is then not interpretable as a detection. We use
|corr| < 0.5 as the threshold for a statistic that measures what it claims to.

Every script in this study refuses to declare a recovery unless the pedestal is
suppressed, the peak lies inside 2,500-4,000 m/s, the peak is not at a scan edge,
the causal side dominates at 3,200 m/s, and p < 0.05. These gates exist because
four automated verdicts were printed from broken inputs during this work; all four
are recorded in section 6.3.

### 3.3 The illumination measurement

Lellouch's evidence for surface sources is the amplitude asymmetry between
downgoing and upgoing energy. We measure it as

    A = ( E(+k) - E(-k) ) / ( E(+k) + E(-k) )

over the frozen 2,500-4,000 m/s fan at **positive frequencies only**. The
restriction is essential: for real input `P(-k,-f) = P(k,f)`, so a budget over all
frequencies would make the two signs identical by construction and force A = 0. We
report |A| only, since which sign corresponds to downgoing depends on a
channel-ordering convention, whereas the existence of a preferred direction does
not.

The null shuffles fan-cell powers between the +k and -k halves, preserving total
fan energy and the marginal distribution of cell powers while destroying any
directional preference.

One correction is central to interpreting our numbers. A separable field
`x(ch,t) = a(ch) s(t)` has `P = |A(k)|^2 |S(f)|^2`, and for real `a(ch)` this gives
`|A(-k)| = |A(k)|` exactly -- perfect k-symmetry, so |A| = 0 regardless of
illumination. Because our field is dominated by exactly such a static pattern
(section 4.4), |A| must be measured **after** projecting out the leading spatial
subspace. We sweep the rank so the reader can see the answer emerge, and we note
that high ranks are conservative against detection, since a genuinely low-rank
plane-wave arrival would also be removed.

### 3.4 Adaptive f-k

We implement the adaptive f-k filter of Isken et al. (2022), which differs from
every other f-k filter used here in assuming no velocity. The mask is the data's own
amplitude spectrum raised to an exponent, computed in sliding time-space windows
with 50 % overlap and Bartlett-tapered recombination, and the amplitude spectrum is
deliberately not smoothed, following their finding that smoothing "distorts the
signal amplitudes and narrows the filter band". A normalised variant retains the
amplitude of the most dominant f-k component.

Before use on real data we verified four properties: an exponent of zero is a
bit-identical no-op; Bartlett 50 %-overlap recombination reconstructs the input to
a maximum relative error of 0.000e+00; the filter raises the coherence of a
synthetic 3,200 m/s plane wave in noise from 0.0523 to 0.1980 at the exponent
alpha = 1 used in production, and to between 0.1231 and 0.1968 across the other five
exponent and normalisation settings checked; and, critically, it does **not**
manufacture moveout from pure noise (0.0000 to 0.0000 at both alpha = 1 and
alpha = 2). The last check targets the failure mode that recurred throughout this
study.

Because an adaptive filter enhances the *dominant* coherent component, and ours is
the static pattern, we predicted in advance that applying it without prior spatial
removal would amplify the contaminant. That prediction was recorded before the run
and is tested in section 4.6.
