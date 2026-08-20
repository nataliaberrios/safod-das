# Ambient-noise body-wave interferometry is illumination-limited, not processing-limited: a failed reproduction of Lellouch et al. (2019) Figure 7c on a later SAFOD DAS deployment

## Abstract

Lellouch et al. (2019) recovered a downgoing ~3.2 km/s P-wave packet from one day
of ambient noise on a distributed acoustic sensing (DAS) array in the San Andreas
Fault Observatory at Depth (SAFOD) main hole, and used the resulting virtual-source
gathers to estimate a P-velocity profile. We attempted to reproduce that result on
a 2024-2025 continuous DAS archive from the same borehole and the same fibre, and
we could not. This paper reports the failure and, more usefully, identifies its
cause.

We first establish that our implementation is correct: applying our picker to
Lellouch's own released correlograms returns a monotonic profile from 2,416 to
4,357 m/s with a depth-velocity correlation of r = 0.948, the range and behaviour
his Figure 9 reports, while the same picker on our data returns velocities
spanning 146 to 1.6e7 m/s with a median of 26,380 m/s -- five orders of magnitude
of scatter about a value an order of magnitude above any crustal velocity. Applied to our archive, the same processing yields no arrival. Six
independent days give a minimum p-value of 0.1345 against a receiver-order
permutation null, with a Fisher combination of p = 0.524; a coherent four-day
stack gives p = 0.9184; and the observed moveout score exceeds the per-velocity
null at 0 of 181 trial velocities. In the
constant-offset geometry that produced the published velocity model, his
correlation peaks migrate from 20.7 to 11.5 ms with depth while ours sit at
exactly zero lag at every depth.

Eight velocity-domain methods fail identically: fixed-fan f-k filtering with
brick-wall and raised-cosine masks, tau-p slant stacking with a slowness-mute
sweep, rank-k subspace projection, phase cross-correlation, phase-weighted
stacking, offset-axis median flat-event removal, and median and mean common-mode
removal. We show why. The field is dominated by a **static spatial pattern at fixed
wavenumber**: raising the analysis band from 5-12 Hz to 12-20 Hz doubles the centre
frequency (ratio 1.882) but moves the low-wavenumber power centroid by 1.1 %
(ratio 1.011), where a wave at fixed velocity requires the two ratios to be equal.
A pattern at fixed wavenumber has no velocity, so no velocity-domain filter can
separate it. We additionally implement the adaptive f-k filter of Isken et al.
(2022), the one f-k family not previously tried, which makes no velocity
assumption.

The decisive result is upstream of all of this. Lellouch's evidence for
surface-origin sources is the amplitude asymmetry between downgoing and upgoing
energy. We measure that asymmetry directly. It is significant in his 2017 records
(|A| = 0.348, p = 0.005) and absent from ours (|A| = 0.040, p = 0.73) under
identical processing at matched spatial ranks; a pre-registered scan of 240
windows spanning one year of archive finds 11 significant windows against 12.0
expected by chance. The wavefield we recorded contains no net downgoing component
in this band, and no filter can create a propagation direction that is not in the
data.

We conclude that borehole ambient-noise body-wave interferometry is limited first
by **illumination** and only then by processing. For two receivers on a vertical
line the stationary-phase sources lie in a narrow cone directly above the
wellhead, so the method depends on local, often cultural, near-vertical-incidence
noise that may simply be absent. Because the same fibre in the same borehole
yielded a positive result in 2017 and a null in 2024-2025, medium properties
cannot account for the difference, and the reproducibility of published
ambient-noise body-wave results should not be assumed to transfer across
deployments at a single site. We recommend that the asymmetry measurement reported
here be run as an inexpensive feasibility gate before any such campaign, and we
report the negative alongside every control, including several of our own claims
that we withdrew during the work.
