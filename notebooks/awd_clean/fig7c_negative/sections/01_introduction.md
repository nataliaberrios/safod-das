## 1. Introduction

Seismic interferometry retrieves the Green's function between two receivers by
cross-correlating the ambient wavefield they both record. In borehole geometry the
appeal is direct: a vertical array of receivers, correlated against one another,
should yield the vertical travel times and hence a velocity profile, with no
active source and no additional acquisition cost. Distributed acoustic sensing
(DAS) makes the geometry unusually favourable, since a single fibre provides
hundreds of channels at metre spacing down the hole.

Lellouch et al. (2019) demonstrated this at the San Andreas Fault Observatory at
Depth (SAFOD). Using an OptaSense ODH3.1 interrogator on the main-hole fibre in
June 2017, they cross-correlated ambient noise against a fixed shallow virtual
source and, after a 5-20 Hz bandpass, recovered an impulsive arrival on the causal
side of the correlation functions with an apparent velocity of about 3,200 m/s,
which they interpreted as a downgoing P wave. Their Figure 7c shows the
fixed-source gather; their Figure 7d shows a constant-offset variant; and their
Figure 9 presents the P-velocity profile derived from it. They attribute the
dominant ambient sources to the surface, inferring this from the amplitude
asymmetry between the downgoing and upgoing energy, and they note that the
approach "can be especially helpful in seismically quiet areas."

That result is a natural reproduction target. The SAFOD fibre remained in place,
and a continuous DAS archive was recorded on it from May 2024 to May 2025 -- the
same borehole, the same fibre, the same rock, roughly seven years later. If
ambient-noise body-wave interferometry is a robust monitoring tool, as the
time-lapse literature assumes when it proposes using it for reservoir or
fault-zone surveillance, it should reproduce.

It does not. This paper documents that failure in full, because the reason turns
out to be more informative than a successful reproduction would have been.

Our contribution is threefold.

**First, we separate implementation error from data.** Reproduction failures are
uninformative unless the implementation is independently validated. We validate
ours against Lellouch's own released correlograms: the picker described in his
paper, applied to his data, returns the monotonic 2,416-4,357 m/s profile his
Figure 9 reports (depth-velocity r = 0.948), and the identical picker applied to
our data returns physically impossible velocities spanning five orders of
magnitude. The processing is therefore not
at fault.

**Second, we explain why an entire class of remedies cannot work.** Eight
velocity-domain methods, including two f-k mask families, tau-p slant stacking and
rank-k subspace projection, fail in the same way. We show that this is not eight
independent failures but one: the field is dominated by a static spatial pattern at
fixed wavenumber, which has no apparent velocity and therefore cannot be separated
by any operator that discriminates on velocity. We also implement and test the
adaptive f-k filter of Isken et al. (2022), which makes no velocity assumption and
which was the one f-k family we had not tried.

**Third, and most usefully for practitioners, we identify the binding
constraint.** Body-wave retrieval on a vertical array requires noise sources in
the stationary-phase zone for a vertical path, which is a narrow cone directly
above the wellhead. That illumination is local and frequently cultural. We
measure the downgoing/upgoing asymmetry that Lellouch used as his own evidence for
surface sources, find it in his data and not in ours, and confirm the absence with
a pre-registered scan over one year of archive. The reproduction fails because the
wavefield lacks the required component, not because the processing lacks skill.

This matters beyond one figure. Behm (2016) retrieved borehole P and S waves by
the virtual-source method in a producing oil field and found that "ambient noise
from time periods as short as 30 seconds is sufficient to obtain robust
interferograms" -- illumination there came from surface industrial activity. If 30 s
suffices when illumination is adequate, then a campaign that fails after days of
stacking is not under-sampled; it is unilluminated. Our own stacks behave exactly
this way: they get *worse* with more data. Treating non-detection as a signal-to-noise
problem, and answering it with longer records or more aggressive filtering, is
then a category error, and one that consumes substantial acquisition and compute
budget before it is discovered.

We therefore recommend the asymmetry measurement as an inexpensive feasibility
gate: it requires seconds of data, no velocity assumption, and no correlation
stack, and it would have told us at the outset that this particular reproduction
was not available.
