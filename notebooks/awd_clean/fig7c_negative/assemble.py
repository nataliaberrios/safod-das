#!/usr/bin/env python3
"""Concatenate sections/*.md into DOCUMENT.md. Never edit DOCUMENT.md by hand.

Follows the convention of awd_clean/manuscript/assemble_manuscript.py: the
sections are the source, the assembled file is a build product. Figure
placeholders of the form {{FIGURE:name}} are replaced with a markdown image plus
caption drawn from CAPTIONS below, so a figure cannot appear in the text without
a caption or vice versa -- the script fails if either is missing.
"""
from pathlib import Path
import re, sys

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"

CAPTIONS = {
 "fig1_no_reproduction":
  "**Figure 1. Figure 7c does not reproduce.** (a) Our fixed-source R+-10 gather "
  "from 24 h of 2024-12-20, with the 3,200 m/s trajectory Lellouch reports "
  "overlaid; no arrival follows it. (b) The moveout scan for the baseline and "
  "common-mode-removed branches. The baseline rises monotonically with trial "
  "velocity, the signature of a statistic measuring proximity to the zero-lag "
  "lobe rather than moveout. (c) The observed maximum against 5,000 "
  "receiver-order permutations; the observation sits inside the null.",
 "fig2_method_validation":
  "**Figure 2. The implementation is correct; the data is not.** (a) Our picker, "
  "unchanged, applied to Lellouch's own released constant-offset correlograms: 14 "
  "of 14 depths return a pick, rising monotonically from 2,416 to 4,357 m/s, the "
  "range his Figure 9 reports, with a depth-velocity correlation of r = 0.948. "
  "Note that r = 0.948 is the correlation between depth and velocity on these "
  "picks, not a point-by-point match to his published curve, which is not "
  "digitised here. (b) The identical picker on our archive, log scale: 114 of 171 "
  "fibre positions return a pick, but they span 146 to 1.6e7 m/s with a median of "
  "26,380 m/s, against the physical range from (a) shaded for comparison. The "
  "failure is not that the picker returns nothing; it is that it returns whatever "
  "noise maximum falls in its search window, which is what a picker does when the "
  "correlation contains no arrival.",
 "fig3_fig7d_isolation":
  "**Figure 3. The constant-offset geometry isolates the failure.** (a) Our "
  "constant-offset gather, source and receiver 50 m apart and slid down the "
  "array together. Every trace peaks at zero lag. (b) Peak lag against depth: "
  "Lellouch's peaks migrate from 20.7 to 11.5 ms as 50 m / v(z) requires for a "
  "rising velocity profile, while ours remain at exactly 0.0000 s at every "
  "depth. A correlation peaking at zero lag independent of receiver separation "
  "indicates a component common to the channels, not a wave propagating between "
  "them.",
 "fig4_eight_methods":
  "**Figure 4. One failure, not eight.** The pedestal diagnostic corr(trial "
  "velocity, moveout score) for every velocity-domain method applied. Values "
  "beyond +-0.5 indicate a statistic dominated by proximity to the zero-lag lobe "
  "rather than by moveout, so its p-value is not interpretable as a detection. "
  "All of these methods discriminate on velocity, and the contaminant has none.",
 "fig5_static_pattern":
  "**Figure 5. The contaminant is at fixed wavenumber, not fixed velocity.** "
  "(a) The low-wavenumber power marginal over the first 16 non-zero cells, for "
  "two disjoint bands, with the power-weighted centroids marked. (b) Raising the "
  "band centre by a factor of 1.882 moves the centroid by 1.011. A wave at fixed "
  "velocity requires the wavenumber ratio to equal the frequency ratio; a static "
  "spatial pattern requires 1.000. The measurement had the range to see a wave: "
  "a 3,200 m/s arrival would sit at 1.86 and 3.50 cells respectively.",
 "fig6_illumination":
  "**Figure 6. The illumination test, with a working positive control.** "
  "(a) Downgoing/upgoing fan asymmetry |A| against the number of spatial "
  "patterns projected out, for both epochs under identical processing, with each "
  "arm's own permutation null. The statistic recovers the asymmetry Lellouch "
  "reported from his own records and finds none in ours. (b) The same at ranks "
  "0-2 with p-values. Rank removal is required because a separable static "
  "pattern forces |A| to zero by algebra, and is conservative against detection.",
 "fig7_archive_scan":
  "**Figure 7. No illuminated windows in one year of archive.** (a) |A| for 240 "
  "windows spanning 2024-05-21 to 2025-05-06. (b) The distribution against the "
  "median null 95th percentile: 11 windows reach p < 0.05 where 12.0 are expected "
  "by chance, with a pre-registered ceiling of 18. (c) Median |A| by hour of day; "
  "cultural surface sources would concentrate in working hours and do not.",
}


def main():
    order = sorted(p for p in (HERE / "sections").glob("*.md"))
    if not order:
        sys.exit("no sections found")
    parts, used = [], set()
    for p in order:
        text = p.read_text()
        for m in re.finditer(r"\{\{FIGURE:([a-z0-9_]+)\}\}", text):
            name = m.group(1)
            if name not in CAPTIONS:
                sys.exit("figure %s has no caption in assemble.py" % name)
            png = FIG / (name + ".png")
            if not png.is_file():
                sys.exit("figure %s referenced but figures/%s.png is missing; "
                         "run make_figures.py" % (name, name))
            used.add(name)
            text = text.replace(
                m.group(0),
                "![%s](figures/%s.png)\n\n%s" % (name, name, CAPTIONS[name]))
        parts.append(text.rstrip() + "\n")
    missing = set(CAPTIONS) - used
    if missing:
        sys.exit("captioned but never referenced: %s" % sorted(missing))
    out = HERE / "DOCUMENT.md"
    out.write_text(
        "<!-- GENERATED by assemble.py from sections/*.md -- do not edit -->\n\n"
        + "\n".join(parts))
    words = len(out.read_text().split())
    print("wrote DOCUMENT.md (%d words, %d sections, %d figures)"
          % (words, len(order), len(used)))


if __name__ == "__main__":
    main()
