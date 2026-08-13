"""Build a one-document overview of the processing chain, raw data to result.

Intended for someone who has not followed the analysis: each page is one stage,
saying what goes in, what happens, what comes out, and showing the figure that
stage produced.

    python make_overview_pdf.py   ->  PROCESSING_OVERVIEW.pdf

Landscape, two columns: narrow text on the left, large figure on the right — so
the figure gets real space rather than being squeezed under full-width prose.

Uses matplotlib's PDF backend; no LaTeX required (none is installed on the
cluster). Figures come from figures/, so rerun make_paper_figures.py and
make_poster_figures.py first if anything upstream changed.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
OUT = HERE / "PROCESSING_OVERVIEW.pdf"

PAGE = (11.0, 8.5)
INK, MUTED, RULE = "#0b0b0b", "#52514e", "#c9c8c3"
NANO, OUTB, RETN = "#2a78d6", "#eb6834", "#1baf7a"

TEXT_L, TEXT_R = 0.055, 0.40          # text column
FIG_L, FIG_R = 0.435, 0.965           # figure column
FIG_B, FIG_T = 0.075, 0.80


def page(pdf: PdfPages, *, kicker: str, title: str, body: list[str],
         image: str | None = None, numbers: list[tuple[str, str]] | None = None):
    fig = plt.figure(figsize=PAGE)
    fig.patch.set_facecolor("white")

    fig.text(TEXT_L, 0.945, "  ".join(kicker), fontsize=10.5, color=OUTB, weight="bold")
    fig.text(TEXT_L, 0.905, title, fontsize=21, color=INK, weight="bold", va="top")
    fig.add_artist(plt.Line2D([TEXT_L, 0.965], [0.862, 0.862], color=RULE, lw=1.2))

    y = 0.815
    for line in body:
        bold = line.startswith("**")
        fig.text(TEXT_L, y, line.replace("**", ""), fontsize=11.5,
                 color=INK if bold else MUTED, va="top",
                 weight="bold" if bold else "normal")
        y -= 0.034 if line else 0.018

    if numbers:
        ny = 0.235
        for label, value in numbers:
            fig.text(TEXT_L, ny, value, fontsize=17, color=INK, weight="bold")
            fig.text(TEXT_L + 0.115, ny + 0.004, label, fontsize=10.5, color=MUTED)
            ny -= 0.048

    if image and (FIG / image).exists():
        data = mpimg.imread(FIG / image)
        avail_w = (FIG_R - FIG_L) * PAGE[0]
        avail_h = (FIG_T - FIG_B) * PAGE[1]
        ratio = data.shape[1] / data.shape[0]
        if avail_w / avail_h > ratio:
            h_in, w_in = avail_h, avail_h * ratio
        else:
            w_in, h_in = avail_w, avail_w / ratio
        w, h = w_in / PAGE[0], h_in / PAGE[1]
        left = FIG_L + ((FIG_R - FIG_L) - w) / 2
        bottom = FIG_B + ((FIG_T - FIG_B) - h) / 2
        ax = fig.add_axes([left, bottom, w, h])
        ax.imshow(data)
        ax.axis("off")

    pdf.savefig(fig)
    plt.close(fig)


def flow_page(pdf: PdfPages):
    fig = plt.figure(figsize=PAGE)
    fig.patch.set_facecolor("white")
    fig.text(0.5, 0.94, "SAFOD repeated-source borehole DAS", fontsize=12.5,
             color=OUTB, weight="bold", ha="center")
    fig.text(0.5, 0.885, "From raw fibre data to a sensitivity number",
             fontsize=27, color=INK, weight="bold", ha="center", va="top")
    fig.text(0.5, 0.818, "Two fibres in one borehole, one repeated source, 24 hours",
             fontsize=15, color=MUTED, ha="center", va="top")

    ax = fig.add_axes([0.03, 0.30, 0.94, 0.44])
    ax.set_xlim(0, 10); ax.set_ylim(0, 4); ax.axis("off")
    steps = [
        ("1  Raw data", "989 GPS drop times\n2 fibres, 24 h", "#e8e8e4"),
        ("2  Stack", "859 drops on both\n46 bursts", "#e8e8e4"),
        ("3  Find arrival", "semblance scan,\nhalf the bursts", "#ffe2d3"),
        ("4  Measure", "delay vs travel time,\nslope = change", "#ffe2d3"),
        ("5  Inject", "known change in,\nblind recovery out", "#d6f0e4"),
        ("6  Sensitivity", "smallest change\nreliably detected", "#d6f0e4"),
    ]
    for i, (head, sub, colour) in enumerate(steps):
        x = 0.12 + i * 1.635
        ax.add_patch(FancyBboxPatch((x, 1.5), 1.40, 1.65,
                                    boxstyle="round,pad=0.06,rounding_size=0.12",
                                    facecolor=colour, edgecolor="#b8b7b2", lw=1.1))
        ax.text(x + 0.70, 2.83, head, fontsize=11, weight="bold", ha="center", color=INK)
        ax.text(x + 0.70, 2.18, sub, fontsize=9.5, ha="center", va="center",
                color=MUTED, linespacing=1.5)
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + 1.47, 2.32), (x + 1.60, 2.32),
                                         arrowstyle="-|>", mutation_scale=16,
                                         color="#8f8e89", lw=1.6))

    fig.text(0.28, 0.235, "Nano   cemented, shallow", fontsize=13.5, color=NANO,
             weight="bold", ha="center")
    fig.text(0.28, 0.196, "927 m of fibre · 30–60 Hz · ~2950 m/s", fontsize=11.5,
             color=MUTED, ha="center")
    fig.text(0.72, 0.235, "Deep   wireline, reversing", fontsize=13.5, color=OUTB,
             weight="bold", ha="center")
    fig.text(0.72, 0.196, "6532 m of fibre, two legs · 15–30 Hz · ~1547 m/s",
             fontsize=11.5, color=MUTED, ha="center")

    fig.add_artist(plt.Line2D([0.15, 0.85], [0.145, 0.145], color=RULE, lw=1.2))
    fig.text(0.5, 0.095, "The cemented fibre resolves a 1% velocity change; "
             "the wireline outbound branch resolves 0.5%",
             fontsize=14, color=INK, ha="center", weight="bold")
    fig.text(0.5, 0.052, "The wireline arrival travels 11.8× longer, but its timing "
             "is 7.6× noisier — so the gain is 2×, not 12×",
             fontsize=12, color=MUTED, ha="center")
    pdf.savefig(fig)
    plt.close(fig)


def main() -> None:
    with PdfPages(OUT) as pdf:
        flow_page(pdf)

        page(pdf, kicker="STEP 1–2", title="Raw data to burst stacks",
             body=[
                 "**In: raw fibre recordings, plus 989",
                 "**GPS times of the weight drops.**",
                 "",
                 "For each drop a 3.5 s window is cut",
                 "from every channel — 0.5 s before,",
                 "3.0 s after. A drop is kept only if",
                 "that window fits inside its file.",
                 "",
                 "**Out: one averaged waveform per",
                 "**burst, per channel.**",
                 "",
                 "Single drops vary by ~20× in signal",
                 "power, so they are averaged in groups.",
                 "Drops >60 s apart are separate bursts.",
                 "",
                 "The Deep fibre stopped after burst 45,",
                 "so the last three bursts are Nano-only.",
                 "That is why 49 bursts becomes 46 — not",
                 "any quality cut.",
             ],
             numbers=[("bursts analysed", "46"), ("drops, both fibres", "859"),
                      ("survey duration", "23.96 h"), ("burst cadence", "30 min")],
             image="fig05_fig12_repeatability_publication.png")

        page(pdf, kicker="STEP 3", title="Finding the arrival",
             body=[
                 "**An arrival moving along the fibre is a",
                 "**straight line in distance vs time.**",
                 "",
                 "Describing it takes two numbers: a start",
                 "time and a speed.",
                 "",
                 "Shift every channel by its predicted",
                 "arrival time and sum. A correct guess",
                 "adds up; a wrong one cancels. The peak",
                 "of that ratio picks the trajectory.",
                 "",
                 "**The search uses half the bursts only,**",
                 "so sensitivity is not measured on the",
                 "same data that chose what to look at.",
                 "",
                 "Scrambling channel order 499 times",
                 "never reproduces the real arrival.",
             ],
             numbers=[("Nano", "2975 m/s"), ("Deep outbound", "1544.6 m/s"),
                      ("Deep return", "1549.7 m/s"), ("legs agree to", "0.3%")],
             image="fig04_deep_tube_null_tests.png")

        page(pdf, kicker="STEP 4", title="Measuring a velocity change",
             body=[
                 "**Compare one burst against the average",
                 "**of all the others.**",
                 "",
                 "Measure the delay at several positions",
                 "along the fibre — per channel for Nano,",
                 "per 400 m aperture beam for Deep, which",
                 "needs the extra stacking at 3 km.",
                 "",
                 "**If the source fired early, every delay",
                 "**shifts by the same amount.**",
                 "**If the medium sped up, the delay grows",
                 "**with distance travelled.**",
                 "",
                 "So fit a line through delay vs travel",
                 "time. The intercept absorbs source",
                 "timing; the slope is the velocity change.",
                 "",
                 "Precision scales with the span in travel",
                 "time — the lever arm.",
             ],
             numbers=[("Nano lever arm", "0.121 s"), ("Deep lever arm", "1.428 s"),
                      ("ratio", "11.8×")],
             image="fig01b_lever_arm.png")

        page(pdf, kicker="STEP 5", title="How small a change is detectable",
             body=[
                 "**Do not estimate this from theory —",
                 "**inject known changes into real data.**",
                 "",
                 "Each channel is shifted exactly as a real",
                 "velocity change would shift it, at 15",
                 "levels from ±0.01% to ±1%, plus zero.",
                 "",
                 "**The recovery code never sees the",
                 "**injected value.** A sealed answer key is",
                 "joined only at the end.",
                 "",
                 "The noise floor is what comes back when",
                 "nothing was injected. The reliable level",
                 "is the smallest injection recovered with",
                 "the right sign, above that floor, 95% of",
                 "the time.",
             ],
             numbers=[("injection levels", "15"), ("trials per leg", "345"),
                      ("held-out bursts", "23"), ("outbound floor", "0.18%")],
             image="poster_p1_recovery.png")

        page(pdf, kicker="RESULT", title="What each installation resolves",
             body=[
                 "**Cemented fibre: 1%.**",
                 "**Wireline outbound branch: 0.5%.**",
                 "**Wireline return branch: 1%.**",
                 "",
                 "The improvement therefore belongs to the",
                 "outbound branch, not to the wireline",
                 "installation as a whole.",
                 "",
                 "**Why only 2×, when the lever arm is",
                 "**11.8× longer?**",
                 "",
                 "Because the deep arrival's timing is 7.6×",
                 "noisier burst to burst. 11.8 ÷ 7.6 = 1.55,",
                 "exactly the observed ratio of scatter.",
                 "",
                 "Geometry helps only in proportion to",
                 "timing repeatability. That is the result",
                 "that transfers to other experiments.",
             ],
             numbers=[("Nano", "1.0%"), ("Deep outbound", "0.5%"),
                      ("Deep return", "1.0%"), ("timing penalty", "7.6×")],
             image="poster_p2_sensitivity.png")

        page(pdf, kicker="CHECKS", title="Why we believe it",
             body=[
                 "**Synthetic test.** A noiseless wave",
                 "perturbed two independent ways recovers",
                 "the same answer to five significant",
                 "figures — the injection is faithful.",
                 "",
                 "**Timing control.** Shift everything by",
                 "5 ms: the intercept absorbs 4.87 ms and",
                 "the velocity estimate moves 7×10⁻⁷.",
                 "Source timing does not leak in.",
                 "",
                 "**Permutation test.** Scramble channel",
                 "order 499 times: zero scrambles reach the",
                 "observed beam power, in all four tests.",
                 "",
                 "**Jackknife.** Drop each aperture, then",
                 "each burst. Nothing moves the result by",
                 "as much as the noise floor.",
                 "",
                 "**One control failed, and is reported.**",
                 "The perturbed-trajectory test recovers",
                 "injections even pointed at the wrong",
                 "arrival, because the injection shifts",
                 "whole traces. It fires the same on both",
                 "legs, so it cannot discriminate — the",
                 "return branch is left unclassified rather",
                 "than explained after the fact.",
             ],
             image="s1_deep_tube_repeatability.png")

        page(pdf, kicker="LIMITS", title="What this does not show",
             body=[
                 "**Not a formation velocity.** The measured",
                 "quantity is the apparent speed of a guided",
                 "mode along the fibre. Both observables are",
                 "dispersive; neither is a clean body wave.",
                 "",
                 "**No depth resolution.** Position is",
                 "distance along fibre. The depth mapping is",
                 "provisional and no result depends on it.",
                 "",
                 "**No tidal detection.** Fitting a solid-",
                 "Earth-tide model to the per-burst",
                 "velocities gives a clean null, with a 95%",
                 "upper limit ~18× above the tidal signal",
                 "expected at this site.",
                 "",
                 "The limit is the 24-hour survey length,",
                 "not the instrument: over one cycle a daily",
                 "signal is hard to separate from drift.",
             ],
             image="s3_deep_dvv_tidal_fit.png")

    size = OUT.stat().st_size / 1e6
    print(f"wrote {OUT.name}  ({size:.1f} MB, 8 pages)")


if __name__ == "__main__":
    main()
