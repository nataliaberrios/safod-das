#!/usr/bin/env python3
"""Promote the standard directional F-K audit and correct Figure 45."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "AWD_results_dashboard.ipynb"
TEX = ROOT / "AWD_advisor_figure_guide.tex"
TEX_V49 = ROOT / "AWD_advisor_figure_guide_v49.tex"
AUDIT_DIR = ROOT / "ambient_transfer" / "fk_directional_audit_v1"
REPORT = AUDIT_DIR / "ambient_fk_directional_audit_v1.json"
FIGURE = AUDIT_DIR / "ambient_fk_directional_audit_v1.png"


def lines(text: str) -> list[str]:
    return text.splitlines(keepends=True)


def unique_cell(notebook: dict, needle: str) -> dict:
    matches = [
        cell
        for cell in notebook["cells"]
        if needle in "".join(cell.get("source", []))
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one notebook cell containing {needle!r}, found {len(matches)}")
    return matches[0]


def validate(report: dict) -> None:
    if report.get("workflow_version") != "ambient_fk_directional_audit_v1":
        raise ValueError("unexpected audit workflow version")
    decision = report["decision"]
    if decision["signed_fk_implementation_validated_by_one_way_synthetics"] is not True:
        raise ValueError("one-way synthetic sign validation did not pass")
    if decision["standard_direction_only_filter_shows_independent_3p2_km_s_result"]:
        raise ValueError("audit unexpectedly claims direction-only 3.2 km/s recovery")
    if decision["production_velocity_fan_result_is_operator_conditioned"] is not True:
        raise ValueError("audit did not flag velocity-fan conditioning")


def main() -> None:
    for path in (NOTEBOOK, TEX, REPORT, FIGURE):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    report = json.loads(REPORT.read_text())
    validate(report)

    held = report["mask_comparisons"]["held_out"]
    production_corr = held["production_2p5_4p5"]["mirror_correlation"]
    direction_corr = held["direction_only"]["mirror_correlation"]
    direction_negative = held["direction_only"]["negative_score_3200"]
    direction_positive = held["direction_only"]["positive_score_3200"]
    power_ratio = report["real_pre_filter_direction_balance"][
        "positive_to_negative_total_power_ratio"
    ]
    legacy = report["legacy_figure45_mirror_check"]

    corrected_notebook = f"""## Velocity-fan-conditioned signed-lag symmetry diagnostic

The aggregate contains 11,523 one-minute files from eight preselected dates.
The two panels were produced by applying the two complementary signed masks
independently and evaluating their nominal physical lag sides. They are not
copies of one another, and the one-way synthetic test verifies the code's sign
mapping. However, a direct numerical audit shows that the two real-data panels
become almost exact lag reflections after the 2.5–4.5 km/s fan:
mirror correlation {legacy['mirror_correlation']:.4f}, relative RMS difference
{legacy['mirror_relative_rms']:.3f}, and norm ratio {legacy['mirrored_norm_ratio_positive_to_negative']:.3f}.

The complementary signed fans do not select Hermitian-conjugate pairs; the
conjugate coefficients of a real gather remain within the same F×K sign. The
near-mirroring is therefore an empirical property of this processed dataset,
not a mathematical requirement. Both fans restrict the magnitude of slowness
and use the same channel-0 correlation geometry, so balanced opposite-direction
energy or shared coherent artifacts can yield similar reflected kernels. The
direction-only audit is required to decide whether the 3.2 km/s structure exists
outside that conditioning. Receiver permutation
after selection tests order within this conditioned output; it does not make
the selected apparent velocity independent of the fan.

**Corrected evidence flag — operator-conditioned diagnostic.** Figure 45 is
retained for provenance and for verifying implementation symmetry, but it is
not evidence for two physical propagation directions or for an independently
measured 3.1–3.2 km/s wave. The standard direction-only held-out result and the
pre-filter spectrum are evaluated separately in the v49 audit below.
"""
    corrected_caption = f"""**Publication-style caption. Velocity-fan-conditioned
signed-lag symmetry diagnostic.** The aggregate uses 11,523 one-minute records
from eight preselected dates, Nano channel 0 as the virtual source, increasing
receiver offset along the cemented fiber, 5–20 Hz preprocessing, 5 s
running-absolute-mean temporal normalization, and complementary signed
2.5–4.5 km/s F–K fans. The upper-left and upper-right panels were reconstructed
and correlated independently; they are not duplicated arrays. Black dashed
trajectories mark the conditional 3.075 km/s score peak, and gray dotted
trajectories mark 3.2 km/s. After reversing the lag axis of the right panel,
the two sections have correlation {legacy['mirror_correlation']:.4f}, relative
RMS difference {legacy['mirror_relative_rms']:.3f}, and norm ratio
{legacy['mirrored_norm_ratio_positive_to_negative']:.3f}. This near-mirror
near-mirror symmetry is an empirical property of the fan-filtered output, not
a requirement of Hermitian Fourier symmetry: conjugate coefficients remain
inside the same F×K sign. Balanced opposite-direction energy or shared coherent
artifacts can produce similar reflected kernels when both fans impose the same
slowness magnitude and channel-0 correlation geometry. The receiver-order
probabilities shown below the
sections are conditional on that selection. The figure therefore validates
consistent implementation of the two complementary fan operators but does not
demonstrate two physical ambient Green's-function branches, upgoing energy,
downgoing energy, wave type, or formation velocity.
"""
    notebook_summary = f"""## Standard directional F–K audit

**Question.** Is the F–K implementation valid, and does the real held-out
wavefield retain the apparent 3.2 km/s arrival when the velocity restriction is
removed?

F–K wavefield separation is standard VSP processing. The crucial distinction
is between (1) selecting a signed quadrant to separate coordinate directions
and (2) selecting a narrow velocity fan and then using the filtered output to
argue for that same velocity. The second calculation is allowed as a
conditioned display, but it is not an independent velocity measurement.

The code was tested with one-way plane waves in both coordinate directions.
Each synthetic is retained by its expected signed mask, rejected by the
opposite mask by more than the frozen factor-of-eight requirement, and produces
the predicted cross-correlation lag. Thus the FFT sign convention and mask
implementation are validated.

For real data, the held-out mask ladder gives:

- production 2.5–4.5 km/s fan: reflected-branch correlation
  {production_corr:.3f};
- direction-only signed quadrant, with no velocity restriction:
  reflected-branch correlation {direction_corr:.3f};
- direction-only scores at 3.2 km/s: {direction_negative:.4f} and
  {direction_positive:.4f};
- pre-filter positive/negative total F–K power ratio: {power_ratio:.4f};
- independent target-corridor enrichment: fails both frozen familywise nulls.

The direction-only held-out sections are not mirror copies, showing that the
implementation can distinguish the signed branches in real data. But neither
direction-only branch contains a stable 3.2 km/s correlation score. The narrow
fan is what makes the two branches almost exact lag mirrors and makes the
3.1–3.2 km/s trajectory visually strong.

**Decision.** F–K filtering is valid and correctly signed. The existing
ambient 3.1–3.2 km/s result is operator-conditioned and is not independently
validated as a physical arrival. It must not be called an upgoing/downgoing
pair or a formation velocity. A physical branch would require corroboration
outside the selected fan—for example a known event/source, an independently
picked ridge in the unmasked F–K spectrum, or a slant-stack/Radon result with a
predeclared decision statistic.

**Relation to Lellouch et al. (2019).** Their reported F–K step was applied to
the S phase of a known earthquake after P-moveout correction, retaining
upgoing events below 11 km/s and rejecting flattened P energy and downgoing
free-surface reflections. Their ambient interferometry section instead used a
top-channel virtual source, 50 m receiver increments, one-day correlations,
and running-absolute-mean normalization. Our signed F–K filtering of ambient
data is therefore an added experiment, not a literal reproduction of their
ambient workflow.
"""
    notebook_caption = f"""**Publication-style caption. Standard directional
F–K audit of the ambient SAFOD workflow.** Panel (a) shows mean 5–20 Hz
pre-filter F–K power from 300 one-minute Nano records before any signed or
velocity mask; total positive- and negative-wavenumber powers are balanced to
a ratio of {power_ratio:.4f}. Panels (b) and (c) show the two held-out
coordinate directions after signed-quadrant selection only, with no
apparent-velocity restriction; dashed trajectories mark 3.2 km/s solely as a
fixed diagnostic. Panel (d) compares correlation between complementary
branches after reflecting one lag axis for the development date and seven
held-out dates. The held-out mirror correlation is {production_corr:.3f} for
the frozen 2.5–4.5 km/s fan but only {direction_corr:.3f} for direction-only
selection. Panel (e) shows held-out moveout scores for the direction-only
branches; their fixed 3.2 km/s values are {direction_negative:.4f} and
{direction_positive:.4f}, providing no independent recovery of that
trajectory. Panel (f) summarizes the balanced pre-filter directional power,
the near-mirror fan response, and the non-mirrored direction-only response.
One-way synthetic tests independently verify mask sign and lag conventions.
The audit therefore supports the correctness of the F–K implementation but
classifies the earlier narrow-fan ambient ridge as operator-conditioned rather
than as an independently demonstrated physical wave.
"""

    notebook = json.loads(NOTEBOOK.read_text())
    status = notebook.setdefault("metadata", {}).setdefault("awd_dashboard", {})
    if status.get("version") != "v48":
        raise ValueError(f"expected v48 notebook, found {status.get('version')}")
    unique_cell(notebook, "## Corrected eight-day physical-lag branch comparison")[
        "source"
    ] = lines(corrected_notebook)
    unique_cell(notebook, "**Publication-style caption. Corrected eight-day signed-lag")[
        "source"
    ] = lines(corrected_caption)
    code = """from pathlib import Path
import json
from IPython.display import Image, display

cwd = Path.cwd()
dashboard_root = next(
    (item for item in (cwd, cwd / "awd_clean")
     if (item / "AWD_results_dashboard.ipynb").exists()),
    None,
)
if dashboard_root is None:
    raise FileNotFoundError("Run from awd_clean or its parent")
product = dashboard_root / "ambient_transfer" / "fk_directional_audit_v1"
report = product / "ambient_fk_directional_audit_v1.json"
figure = product / "ambient_fk_directional_audit_v1.png"
for required in (report, figure):
    if not required.is_file():
        raise FileNotFoundError(required)
print(json.dumps(json.loads(report.read_text()), indent=2))
display(Image(filename=str(figure), width=1200))
"""
    notebook["cells"].extend(
        [
            {"cell_type": "markdown", "metadata": {}, "source": ["# v49\n"]},
            {"cell_type": "markdown", "metadata": {}, "source": lines(notebook_summary)},
            {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": lines(code)},
            {"cell_type": "markdown", "metadata": {}, "source": lines(notebook_caption)},
        ]
    )
    status["version"] = "v49"
    status["last_status_sync"] = "2026-08-13"
    NOTEBOOK.write_text(json.dumps(notebook, indent=1) + "\n")

    tex = TEX.read_text()
    if tex.count("version v48") != 1:
        raise ValueError("authoritative TeX does not have exactly one v48 marker")
    tex = tex.replace("version v48", "version v49", 1)
    tex = tex.replace(r"\date{11 August 2026}", r"\date{13 August 2026}", 1)
    start = r"\subsection{Corrected eight-day signed-lag comparison (v47)}"
    stop = r"\section{Ambient F--K injection--recovery calibration (v48)}"
    if tex.count(start) != 1 or tex.count(stop) != 1:
        raise ValueError("could not isolate Figure 45 section")
    before, remainder = tex.split(start, 1)
    _, after = remainder.split(stop, 1)
    corrected_tex = rf"""\subsection{{Velocity-fan-conditioned signed-lag symmetry diagnostic (corrected v49)}}
\figquestion{{What does the near-equality of the two signed, velocity-fan-filtered
correlation panels establish?}}

The two branches were filtered and correlated independently, and one-way
synthetics validate their sign mapping.  Nevertheless, reversing the lag axis
of the right panel makes the two real-data sections nearly identical: their
correlation is {legacy['mirror_correlation']:.4f}, relative RMS difference is
{legacy['mirror_relative_rms']:.3f}, and norm ratio is
{legacy['mirrored_norm_ratio_positive_to_negative']:.3f}.  This is a property
of the selected operator, not evidence for two independently observed waves.
The complementary signed fans do not select Hermitian-conjugate pairs; the
conjugate coefficients of a real gather remain within the same $F K$ sign.
The near-mirroring is therefore an empirical property of this processed data
set, not a mathematical requirement.  Balanced opposite-direction energy or
shared coherent artifacts can produce similar reflected kernels when both
fans impose the same slowness magnitude and channel-0 correlation geometry.
The direction-only audit is needed to test whether the 3.2-km-s\( ^{{-1}} \)
structure persists outside that conditioning.  Receiver
permutation after filtering tests spatial order in that conditioned output but
does not make the selected velocity independent of the fan.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.99\linewidth]{{ambient_transfer/signed_lag_v2/seasonal_signed_fk_v2_aggregate.png}}
\caption{{\textbf{{Velocity-fan-conditioned signed-lag symmetry diagnostic.}} The aggregate uses 11,523 one-minute records from eight preselected dates, Nano channel~0 as the virtual source, increasing receiver offset along the cemented fiber, 5--20~Hz preprocessing, 5-s running-absolute-mean temporal normalization, and complementary signed 2.5--4.5~km~s\( ^{{-1}} \) F--K fans. The two upper panels were reconstructed and correlated independently; they are not duplicated arrays. Black dashed trajectories mark the conditional 3.075~km~s\( ^{{-1}} \) score peak, and gray dotted trajectories mark 3.2~km~s\( ^{{-1}} \). After reversing the lag axis of the right panel, the sections have correlation {legacy['mirror_correlation']:.4f}, relative RMS difference {legacy['mirror_relative_rms']:.3f}, and norm ratio {legacy['mirrored_norm_ratio_positive_to_negative']:.3f}. The near-mirror symmetry is an empirical property of the fan-filtered output, not a requirement of Hermitian Fourier symmetry: conjugate coefficients remain inside the same $F K$ sign. Balanced opposite-direction energy or shared coherent artifacts can produce similar reflected kernels when both fans impose the same slowness magnitude and channel-0 correlation geometry. The receiver-order probabilities below the sections are conditional on that selection. The figure validates consistent implementation of the complementary fan operators but does not demonstrate two physical ambient Green's-function branches, upgoing energy, downgoing energy, wave type, or formation velocity.}}
\end{{figure}}

\takeaway{{Figure 45 is retained as an operator-symmetry and provenance
diagnostic.  It is not evidence that both physical propagation directions are
present.  The standard direction-only test is reported in Figure 47.}}


{stop}"""
    tex = before + corrected_tex + after

    end_marker = r"\end{document}"
    if tex.count(end_marker) != 1:
        raise ValueError("could not isolate TeX document end")
    audit_tex = rf"""
\section{{Standard directional F--K audit (v49)}}
\figquestion{{Is the signed F--K implementation correct, and does the held-out
ambient wavefield retain a 3.2-km-s\( ^{{-1}} \) arrival without a velocity
restriction?}}

F--K wavefield separation is standard VSP processing.  The evidentiary issue
is narrower: signed-quadrant selection separates coordinate directions, while
a velocity fan additionally imposes an allowed range of slopes.  A fan-filtered
display is valid, but the same fan cannot supply independent evidence for a
velocity inside its own passband.

The present code passes one-way plane-wave tests in both coordinate directions:
each synthetic is retained by the expected signed mask, rejected by the
opposite mask by more than a factor of eight, and produces the predicted
correlation-lag sign.  The FFT convention, channel direction, mask sign, and
lag extraction are therefore validated.

In the held-out real data, the 2.5--4.5~km~s\( ^{{-1}} \) production fan gives
a reflected-branch correlation of {production_corr:.3f}, whereas standard
direction-only selection gives {direction_corr:.3f}.  Direction-only scores
at 3.2~km~s\( ^{{-1}} \) are {direction_negative:.4f} and
{direction_positive:.4f}.  Before masking, the ratio of total positive- to
negative-wavenumber power is {power_ratio:.4f}, and the frozen target-corridor
energy statistic fails both familywise null tests.  Thus the implementation
distinguishes directions, but the current held-out data do not retain an
independent 3.2-km-s\( ^{{-1}} \) correlation after the velocity fan is removed.

This distinction matches standard practice.  Lellouch et al. (2019,
\href{{https://doi.org/10.1029/2019JB017533}}{{doi:10.1029/2019JB017533}})
used F--K filtering on a known earthquake S-phase gather after P-moveout
correction, passing upgoing events below 11~km~s\( ^{{-1}} \) and suppressing
flattened P energy and downgoing free-surface reflections.  Their ambient
interferometry used a top-channel virtual source, receivers at 50-m increments,
one-day correlations, and running-absolute-mean normalization; signed F--K
filtering of the ambient correlations was not the reported step.  Rao and Wang
(2016, \href{{https://doi.org/10.1088/1742-2132/13/3/412}}{{doi:10.1088/1742-2132/13/3/412}})
likewise describe F--K separation as standard while documenting overlap of
opposite wavefields near zero wavenumber and limitations of a simple fan mask.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.99\linewidth]{{ambient_transfer/fk_directional_audit_v1/ambient_fk_directional_audit_v1.png}}
\caption{{\textbf{{Standard directional F--K audit of the ambient SAFOD workflow.}} (a) Mean 5--20~Hz pre-filter F--K power from 300 one-minute Nano records before any signed or velocity mask; total positive- and negative-wavenumber powers have ratio {power_ratio:.4f}. (b,c) Held-out coordinate-direction sections after signed-quadrant selection only, with no apparent-velocity restriction; dashed trajectories mark 3.2~km~s\( ^{{-1}} \) solely as a fixed diagnostic. (d) Correlation between complementary branches after reflecting one lag axis, shown for the development date and seven held-out dates. The held-out mirror correlation is {production_corr:.3f} for the frozen 2.5--4.5~km~s\( ^{{-1}} \) fan but {direction_corr:.3f} for direction-only selection. (e) Held-out moveout scores for the direction-only branches; fixed 3.2-km-s\( ^{{-1}} \) values are {direction_negative:.4f} and {direction_positive:.4f}, providing no independent recovery of that trajectory. (f) Summary of balanced pre-filter directional power, near-mirror fan response, and non-mirrored direction-only response. One-way synthetic tests independently verify mask sign and lag conventions. The audit supports the correctness of the F--K implementation but classifies the earlier narrow-fan ambient ridge as operator-conditioned rather than as an independently demonstrated physical wave.}}
\end{{figure}}

\takeaway{{The filtering code is correct.  The current ambient 3.1--3.2-km-s
\( ^{{-1}} \) interpretation is not.  It remains a selected visualization until
an unmasked or independently constrained method detects the same event.}}

"""
    tex = tex.replace(end_marker, audit_tex + end_marker)
    TEX.write_text(tex)
    TEX_V49.write_text(tex)
    print(
        json.dumps(
            {
                "version": "v49",
                "figure45_mirror_correlation": legacy["mirror_correlation"],
                "held_out_production_mirror_correlation": production_corr,
                "held_out_direction_only_mirror_correlation": direction_corr,
                "held_out_direction_only_scores_3200": [
                    direction_negative,
                    direction_positive,
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
