#!/usr/bin/env python3
"""The velocities the figures draw, in one place, with their provenance.

WHY THIS MODULE EXISTS. `1675.0` was hard-coded in five scripts as the Deep
arrival velocity. It is wrong, and because it was copied rather than imported,
every figure in the repository drew a wrong line and each had to be found by
grep. Anything a figure annotates belongs here.

--------------------------------------------------------------------------
DEEP FIBRE ARRIVAL -- 1675 m/s is RETRACTED (2026-08-20)
--------------------------------------------------------------------------
The `+-0.35 s` lag window cannot hold a 700 m offset below 1934 m/s. On a
300-6000 m/s grid that left 71 of 229 trial velocities with their far-offset
gates outside the window; the median then ran over the surviving NEAR offsets,
which sit closest to the zero-lag lobe and score high. The score therefore rose
as velocity fell for a purely geometric reason -- the "pedestal" -- and the peak
was dragged down to 1675 m/s.

Re-aggregated at +-2.5 s, which holds 700 m down to 280 m/s, three geometrically
valid source channels give:

    src 211 (TVD    2 m)  1350 m/s   score 18.82   max-null p = 0.0002
    src 400 (TVD  389 m)  1525 m/s   score 27.08   max-null p = 0.0002
    src 800 (TVD 1208 m)  1550 m/s   score  7.66   max-null p = 0.0002

THIS IS A DEPTH TREND, NOT SCATTER. The velocity rises monotonically with source
depth, which is what a mode whose speed tracks formation stiffness should do, and
it means NO SINGLE NUMBER is right for every figure -- a figure must draw the
velocity measured for ITS OWN source. `deep_velocity_for_source()` exists so that
happens by construction rather than by remembering.

An independent envelope-peak regression (lag on separation, 100-450 m) over four
gathers spanning BOTH fibre limbs gives 1443-1537 m/s, inside this range.

Treat the trend as provisional: three sources is three points, the three gathers
have very different scores (7.7 to 27.1), and no per-velocity null has yet been
run on the corrected aggregates. It is reported because it is what the numbers
say, not because it is established.

WHAT IT IS NOT. This is not Lellouch's 3200 m/s P wave, and it is not a body
wave. It is ~2x slower, it appears on a WIRELINE (fluid-coupled) fibre, and the
per-velocity p at 3200 m/s is 0.20-0.36. It is consistent with a fluid-guided
mode, but tube-wave IDENTITY is unverified: no borehole diameter, fluid, casing
or cement data exist in this workspace. Do not upgrade "consistent with" to "is".

The active-source comparison at 1547 m/s is NOT independent corroboration --
`deep_dvv_frozen_trajectory.json` shows that semblance search was confined to
1300-1800 m/s on tube-wave grounds, so it could not have returned anything else.

--------------------------------------------------------------------------
NANO FIBRE -- which number is the a-priori target depends on the BAND
--------------------------------------------------------------------------
`nano_mode_identification.txt` measures 3300 m/s at 15-25 Hz and 2950 m/s at
35-80 Hz, a slowness trend of 0.483 us m^-1 Hz^-1. Extrapolated into the 5-20 Hz
analysis band that gives 3102-3238 m/s.

So 2950 m/s is a 30-60 Hz number and was the WRONG a-priori target for 5-20 Hz
ambient work; the in-band expectation is ~3200 m/s, which coincides with
Lellouch's value. Testing 3200 m/s at 5-20 Hz on the cemented fibre was
band-appropriate, so the Figure 7c negative is not an artefact of band choice.

See AUDIT_2026-08-20.md.
"""
from __future__ import annotations

# Deep fibre, 5-20 Hz, per source channel: scan peak at +-2.5 s over 300-6000 m/s.
# Keyed by source so a figure cannot draw another source's velocity.
V_DEEP_BY_SOURCE = {211: 1350.0, 400: 1525.0, 800: 1550.0}

# Default for contexts with no single source (whole-archive stacks, schematics).
# The wellhead is the Lellouch Figure 7c geometry, so it is the canonical one.
V_DEEP_ARRIVAL = V_DEEP_BY_SOURCE[211]
V_DEEP_ARRIVAL_RANGE = (1350.0, 1550.0)

# Nano / cemented fibre. Use the in-band value for 5-20 Hz work.
V_NANO_INBAND = 3200.0            # 5-20 Hz, extrapolated 3102-3238
V_NANO_HIGHBAND = 2950.0          # 30-60 Hz only -- not an ambient a-priori target

# Lellouch et al. (2019) Figure 7c, for reference lines.
V_LELLOUCH = 3200.0

# Deep active-source comparison. Reported for context only: see the docstring,
# its search was confined to 1300-1800 m/s so it cannot corroborate.
V_DEEP_ACTIVE = 1547.0
V_DEEP_ACTIVE_BANDMATCHED = 1592.0


def deep_label() -> str:
    """Figure annotation for the Deep arrival, carrying its uncertainty."""
    return "%.0f m/s (%.0f-%.0f)" % (V_DEEP_ARRIVAL, *V_DEEP_ARRIVAL_RANGE)


def deep_velocity_for_source(source_channel):
    """Measured arrival velocity for THIS source channel, or None if unmeasured.

    The velocity rises with source depth (1350 -> 1550 m/s over TVD 2 -> 1208 m),
    so a figure that draws a global constant draws the wrong line for two of the
    three sources. Returning None rather than a fallback is deliberate: a caller
    with no measurement should say so on the figure, not silently annotate it
    with a number measured somewhere else.
    """
    return V_DEEP_BY_SOURCE.get(int(source_channel))


def plot_lag_limit(lags, max_offset_m, v_m_s, margin=1.3):
    """Lag half-range a figure should SHOW, so the drawn moveout is on the page.

    Several figures were hard-limited to +-0.35 s while annotating a moveout
    curve that leaves it: 700 m at 1350 m/s arrives at 0.52 s, so the curve ran
    off the axis and the arrival it marks was simply not in the picture. Derive
    the limit from the curve being drawn, and never show more than was retained.
    """
    import numpy as _np
    have = float(_np.abs(_np.asarray(lags)).max())
    want = margin * float(max_offset_m) / float(v_m_s)
    return min(have, want) if want > 0 else have


def required_lag_s(max_offset_m, v_min_m_s, gate_half_s=0.012, safety=1.5):
    """Smallest lag half-window that can hold every offset at every velocity.

    Pick the window from the geometry, not by hand. The whole 1675 m/s error came
    from a hand-picked 0.35 s that could not hold a 700 m offset below 1934 m/s;
    the far gates then fell outside the window and the score was taken over the
    near offsets only, which sit closest to the zero-lag lobe and read high. The
    bias grows as velocity falls, which manufactures a low-velocity peak.

    There is no meaningful cost to headroom -- a 30 s correlation window already
    contains lags to +-30 s, so a wider crop only keeps more of what was computed.
    The one real consequence is that `normalized_envelope` divides each trace by
    its median envelope OVER THE RETAINED WINDOW, so widening the window changes
    that denominator and the score scale. Scores are therefore comparable only
    between runs using the same window, which is why this is derived and recorded
    rather than maximised.

    `safety` is deliberately modest: the point is to clear the requirement, not to
    inflate the window until the normalisation drifts.
    """
    return float(safety * (float(max_offset_m) / float(v_min_m_s) + gate_half_s))


if __name__ == "__main__":
    print("Deep arrival   : %s" % deep_label())
    print("Nano, 5-20 Hz  : %.0f m/s" % V_NANO_INBAND)
    print("Lellouch 7c    : %.0f m/s" % V_LELLOUCH)
