#!/usr/bin/env python3
"""Promote the consolidated F-K validation and faithful Lellouch transfer test."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "AWD_results_dashboard.ipynb"
TEX = ROOT / "AWD_advisor_figure_guide.tex"
TEX_V50 = ROOT / "AWD_advisor_figure_guide_v50.tex"
VALIDATION_DIR = ROOT / "ambient_transfer" / "fk_validation_evidence_v1"
VALIDATION_REPORT = VALIDATION_DIR / "fk_validation_evidence_v1.json"
VALIDATION_FIGURE = VALIDATION_DIR / "fk_validation_evidence_v1.png"
LELLOUCH_DIR = ROOT / "ambient_transfer" / "lellouch2019_reproduction_v1"


def lines(text: str) -> list[str]:
    return text.splitlines(keepends=True)


def unique_cell(notebook: dict, needle: str) -> dict:
    matches = [
        cell for cell in notebook["cells"]
        if needle in "".join(cell.get("source", []))
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one notebook cell containing {needle!r}, found {len(matches)}"
        )
    return matches[0]


def latest_complete_day() -> tuple[Path, Path, dict]:
    reports = sorted(LELLOUCH_DIR.glob("lellouch2019_2024-12-20_start0_requestedall_used*.json"))
    if len(reports) != 1:
        raise ValueError(
            "expected exactly one completed full-day Lellouch report; "
            f"found {len(reports)}"
        )
    report_path = reports[0]
    report = json.loads(report_path.read_text())
    if report.get("workflow_version") != "ambient_lellouch2019_reproduction_v1":
        raise ValueError("unexpected Lellouch workflow version")
    if report.get("used_files") != 1440 or report.get("used_30_s_windows") != 5759:
        raise ValueError(
            "full-day completion requirement failed: expected 1440 files and 5759 windows"
        )
    if len(report.get("target_channels", [])) != 16:
        raise ValueError("full-day product does not contain 50--800 m targets")
    figure_path = report_path.with_suffix(".png")
    if not figure_path.is_file() or figure_path.stat().st_size == 0:
        raise FileNotFoundError(figure_path)
    return report_path, figure_path, report


def validate(validation: dict) -> None:
    if validation.get("workflow_version") != "fk_validation_evidence_v1":
        raise ValueError("unexpected consolidated F-K report version")
    decision = validation["decision"]
    required = (
        "operator_implementation_validated",
        "both_coordinate_directions_recovered_on_real_background",
        "real_ambient_input_is_nearly_direction_balanced",
    )
    for key in required:
        if decision.get(key) is not True:
            raise ValueError(f"consolidated F-K validation failed: {key}")


def main() -> None:
    for path in (NOTEBOOK, TEX, VALIDATION_REPORT, VALIDATION_FIGURE):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    validation = json.loads(VALIDATION_REPORT.read_text())
    validate(validation)
    lellouch_report_path, lellouch_figure_path, lellouch = latest_complete_day()

    synthetic = validation["exact_5_20_hz_synthetic_correct_to_wrong_rms_ratio"]
    active = validation["awd_empirical_direction"][
        "expected_to_opposite_fixed_tube_energy_ratio"
    ]
    ambient = validation["real_ambient_direction_balance"]
    metrics = lellouch["metrics"]
    isolated = metrics["isolated_receiver"]
    simple = metrics["local_21_channel_simple"]
    aligned = metrics["local_21_channel_aligned_3200"]
    alignment_gain = abs(aligned["causal_score_at_3200_m_s"]) / max(
        abs(isolated["causal_score_at_3200_m_s"]), 1e-30
    )

    notebook = json.loads(NOTEBOOK.read_text())
    status = notebook.setdefault("metadata", {}).setdefault("awd_dashboard", {})
    if status.get("version") != "v49":
        raise ValueError(f"expected v49 notebook, found {status.get('version')}")

    symmetry = unique_cell(
        notebook, "## Velocity-fan-conditioned signed-lag symmetry diagnostic"
    )
    symmetry_text = "".join(symmetry["source"])
    symmetry_text = symmetry_text.replace(
        "**Corrected evidence flag — operator-conditioned diagnostic.** Figure 45 is\n"
        "retained for provenance and for verifying implementation symmetry, but it is\n"
        "not evidence for two physical propagation directions or for an independently\n"
        "measured 3.1–3.2 km/s wave. The standard direction-only held-out result and the\n"
        "pre-filter spectrum are evaluated separately in the v49 audit below.",
        "**Status update v50 — valid selected-operator diagnostic.** Figure 45 is\n"
        "retained for provenance. The v50 validation shows that both signed operators\n"
        "recover known waves in the correct coordinate direction and that the real\n"
        "ambient input is nearly balanced between those directions. The near-mirror\n"
        "correlations still do not identify physical upgoing/downgoing body waves or\n"
        "provide an independent velocity estimate.",
    )
    symmetry["source"] = lines(symmetry_text)

    v49 = unique_cell(notebook, "## Standard directional F–K audit")
    v49_text = "".join(v49["source"])
    old_decision = """**Decision.** F–K filtering is valid and correctly signed. The existing
ambient 3.1–3.2 km/s result is operator-conditioned and is not independently
validated as a physical arrival. It must not be called an upgoing/downgoing
pair or a formation velocity. A physical branch would require corroboration
outside the selected fan—for example a known event/source, an independently
picked ridge in the unmasked F–K spectrum, or a slant-stack/Radon result with a
predeclared decision statistic.
"""
    new_decision = """**Status update v50 — interpretation superseded.** F–K filtering is valid,
correctly signed, and recovers injected waves in both coordinate directions on
real background. Failure of an intentionally broader direction-only filter to
show the same ridge is a sensitivity result, not a rejection of the standard
velocity fan. The fan-assisted 3.1–3.2 km/s ridge is a reproducible selected
observable in an independently motivated P-wave interval; the fan itself
cannot supply an unbiased velocity estimate or establish wave type. Physical
labels remain “increasing-coordinate” and “decreasing-coordinate” until source
distribution or reflection geometry resolves upgoing versus downgoing.
"""
    if old_decision not in v49_text:
        raise ValueError("could not update stale v49 notebook decision")
    v49["source"] = lines(v49_text.replace(old_decision, new_decision, 1))

    validation_summary = f"""## Signed F–K validation ladder

**Question.** Can the signed 2.5–4.5 km/s F–K operator be trusted, and what
physical language is supported for the two ambient branches?

F–K filtering is a standard VSP wavefield-separation operation. A defensible
application must separate operator validation from geological interpretation.
Four independent checks now provide that separation:

1. Exact 5–20 Hz one-way synthetics pass through the production sampling and
   masks. The correct/wrong output-RMS ratios are
   {synthetic['increasing_z']:.2f}:1 for increasing-coordinate propagation and
   {synthetic['decreasing_z']:.2f}:1 for decreasing-coordinate propagation.
2. Broadband known-direction waves injected before nonlinear normalization
   into 300 real one-minute records are independently recovered in the correct
   branch for both directions and for all in-wedge test speeds (2.75, 3.2, and
   4.0 km/s). The off-wedge 1.8 km/s specificity control is not recovered.
3. The known surface AWD empirically fixes the physical coordinate sign where
   its arrival is directional. Expected/opposite energy in a fixed 2.975 km/s
   moveout tube is {active['5_20_hz']:.2f} at 5–20 Hz but
   {active['25_60_hz']:.1f} at 25–60 Hz. Thus the active source strongly
   calibrates the sign at higher frequency; its low-frequency gather is too
   balanced to serve as a decisive 5–20 Hz direction standard.
4. The unmasked real ambient input is itself nearly balanced: the
   F×K<0/F×K>0 target-fan power ratio is
   {ambient['pre_filter_fan_power']:.3f}, and the corresponding filtered
   correlation-score ratio is {ambient['post_filter_correlation_score']:.3f}.

**Decision.** The filtering implementation is validated. The two outputs are
opposite propagation directions along the fiber coordinate and are not
duplicated arrays. The 5–20 Hz ambient wavefield contains comparable energy in
both coordinate directions. That is compatible with illumination from both
sides, reflections, scattering, and F–K leakage where events overlap near
zero wavenumber. It does not by itself identify two body-wave Green's-function
branches. Use the labels **increasing-coordinate (F×K<0)** and
**decreasing-coordinate (F×K>0)**; reserve “downgoing” and “upgoing” for a
source- or reflection-calibrated interpretation.

**What the fan result answers.** A 2.5–4.5 km/s fan is scientifically
reasonable because its passband was motivated independently by Lellouch et
al.'s approximately 3.2 km/s SAFOD result and the AWD's approximately 2.98
km/s mode. The filtered correlation demonstrates coherent energy compatible
with that predeclared band. It does not independently discover or precisely
estimate velocity inside the band.
"""
    validation_code = """from pathlib import Path
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
product = dashboard_root / "ambient_transfer" / "fk_validation_evidence_v1"
report = product / "fk_validation_evidence_v1.json"
figure = product / "fk_validation_evidence_v1.png"
for required in (report, figure):
    if not required.is_file():
        raise FileNotFoundError(required)
print(json.dumps(json.loads(report.read_text()), indent=2))
display(Image(filename=str(figure), width=1200))
"""
    validation_caption = f"""**Publication-style caption. Validation ladder for
the signed SAFOD ambient F–K operator.** (a) Exact 10 Hz, 3.2 km/s plane waves
are propagated separately toward increasing and decreasing fiber coordinate,
decimated identically to the production workflow, and passed through the
complementary 5–20 Hz, 2.5–4.5 km/s signed masks. The ratio of correct-branch
to wrong-branch output RMS is {synthetic['increasing_z']:.2f}:1 for the
increasing-coordinate wave and {synthetic['decreasing_z']:.2f}:1 for the
decreasing-coordinate wave, validating the Fourier-sign implementation.
(b) Minimum injected broadband RMS, expressed relative to median real-channel
5–20 Hz RMS, required for independent pre-filter recovery after injection into
the same 300 real one-minute records. Both coordinate directions are recovered
at 2.75, 3.2, and 4.0 km/s; the independently audited calculation requires
correct-branch target enrichment above the stricter frozen null threshold and
a positive paired-bootstrap lower bound. The 1.8 km/s off-wedge control is
excluded from this panel because it is never recovered, demonstrating tested
velocity specificity. (c) The exact signed fan is applied at full resolution
to the known surface-AWD Nano stack over 80–440 m. Energy is measured within a
fixed ±12 ms tube around the independently estimated 2.975 km/s
increasing-coordinate moveout. Expected/opposite energy ratios are
{active['5_20_hz']:.2f} at 5–20 Hz and {active['25_60_hz']:.1f} at 25–60 Hz,
showing that empirical direction calibration is strong at higher frequency but
not decisive in the ambient band. (d) In the unmasked ambient input, the
F×K<0/F×K>0 target-fan power ratio is
{ambient['pre_filter_fan_power']:.3f}; after independent signed filtering and
physical-lag correlation, the 3.2 km/s score ratio is
{ambient['post_filter_correlation_score']:.3f}. The real 5–20 Hz wavefield is
therefore nearly balanced between coordinate directions. Together the panels
validate the signed F–K implementation and show that the paired ambient panels
are not a software mirror. They do not determine whether decreasing-coordinate
energy originates below the array, is reflected or scattered, or leaks where
wavefields overlap near zero wavenumber; nor do they provide an independent
formation-velocity estimate.
"""

    lellouch_summary = f"""## Published ambient-workflow transfer test

**Question.** Does a close implementation of the ambient workflow reported by
Lellouch et al. (2019) recover the approximately 3.2 km/s causal arrival in the
2024–2025 acquisition without F–K filtering, and did the previously omitted
receiver-neighborhood stack matter?

One complete UTC day ({lellouch['used_files']:,} one-minute files;
{lellouch['used_30_s_windows']:,} overlapping 30 s windows) was processed with
Nano channel 0 as the fixed top virtual source and target receivers every 50 m
from 50 to 800 m. Raw phase/strain was differentiated to a strain-rate proxy,
divided by a 5 s running absolute mean, correlated in 30 s windows with 15 s
overlap, stacked, and filtered at 5–20 Hz. The 5 s normalization duration is an
explicit project choice because the paper does not report that duration.

The panels isolate the effect of the published spatial enhancement:

- individual receivers at 50 m increments;
- the reported local stack over receiver channels R−10 through R+10, with the
  same channel-0 virtual source;
- the same local correlations shifted at the published average 3.2 km/s before
  stacking.

At the fixed 3.2 km/s trajectory, the isolated, simple-local, and aligned-local
causal scores are {isolated['causal_score_at_3200_m_s']:.4f},
{simple['causal_score_at_3200_m_s']:.4f}, and
{aligned['causal_score_at_3200_m_s']:.4f}, respectively. The aligned-to-isolated
absolute-score ratio is {alignment_gain:.2f}. The aligned section's descriptive
causal peak is {aligned['descriptive_causal_peak_velocity_m_s']/1000:.3f} km/s;
its causal/anti-causal absolute-score ratio at 3.2 km/s is
{aligned['abs_causal_to_anti_causal_ratio_at_3200']:.2f}.

**Interpretive boundary.** This is the closest documented transfer of the
published ambient subsection to the available 2024–2025 data, not a
bit-for-bit reproduction of the June–July 2017 acquisition. No F–K filter is
used here because none is reported in that ambient subsection. The comparison
therefore separates failure of the unfiltered transfer from validity of the
added F–K-assisted observable.
"""
    lellouch_code = f"""from pathlib import Path
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
product = dashboard_root / "ambient_transfer" / "lellouch2019_reproduction_v1"
report = product / "{lellouch_report_path.name}"
figure = product / "{lellouch_figure_path.name}"
for required in (report, figure):
    if not required.is_file():
        raise FileNotFoundError(required)
print(json.dumps(json.loads(report.read_text()), indent=2))
display(Image(filename=str(figure), width=1300))
"""
    lellouch_caption = f"""**Publication-style caption. Transfer of the Lellouch
et al. (2019) ambient-interferometry workflow to the 2024–2025 SAFOD Nano
acquisition.** A complete UTC day ({lellouch['used_files']:,} contiguous
one-minute files) is divided into {lellouch['used_30_s_windows']:,} 30 s
correlation windows with 15 s overlap. Nano channel 0 is the fixed top virtual
source; receiver coordinates increase down the cemented main-hole fiber and
are displayed every 50 m from 50 to 800 m using the measured
{lellouch['channel_spacing_m']:.6f} m channel spacing. Recorded phase,
proportional to strain, is differentiated in time to form a strain-rate proxy;
earthquake transients and other amplitude excursions are suppressed by 5 s
running-absolute-mean normalization. Correlations are averaged over the day
and bandpass filtered at 5–20 Hz. Black variable-area wiggles show positive
correlation lag, and the red line is the fixed 3.2 km/s published reference,
not a fitted velocity. (a) Individual channel-0-to-receiver correlations at
50 m increments. (b) The exact spatial enhancement reported in the paper,
formed by stacking correlations from the same virtual source to receiver
channels R−10 through R+10. (c) The same 21 local correlations shifted at the
fixed 3.2 km/s travel-time difference before stacking, applied only after the
simple-stack reference as prescribed in the paper. Fixed-trajectory causal
scores for (a–c) are {isolated['causal_score_at_3200_m_s']:.4f},
{simple['causal_score_at_3200_m_s']:.4f}, and
{aligned['causal_score_at_3200_m_s']:.4f}; the aligned value is
{alignment_gain:.2f} times the isolated absolute score. The aligned section's
descriptive causal maximum occurs at
{aligned['descriptive_causal_peak_velocity_m_s']/1000:.3f} km/s and its
causal/anti-causal absolute-score ratio at 3.2 km/s is
{aligned['abs_causal_to_anti_causal_ratio_at_3200']:.2f}. These metrics state
what the transferred workflow produces; they are not a formation-velocity
inversion. Differences from the published figure can reflect the distinct
2024–2025 acquisition, interrogator/gauge settings, source distribution, and
the assumed inheritance of the 2017 channel-0 coordinate. No F–K filter is
used in this figure, because Lellouch et al. did not report F–K filtering in
their ambient-interferometry subsection.
"""

    notebook["cells"].extend(
        [
            {"cell_type": "markdown", "metadata": {}, "source": ["# v50\n"]},
            {"cell_type": "markdown", "metadata": {}, "source": lines(validation_summary)},
            {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": lines(validation_code)},
            {"cell_type": "markdown", "metadata": {}, "source": lines(validation_caption)},
            {"cell_type": "markdown", "metadata": {}, "source": lines(lellouch_summary)},
            {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": lines(lellouch_code)},
            {"cell_type": "markdown", "metadata": {}, "source": lines(lellouch_caption)},
        ]
    )
    status["version"] = "v50"
    status["last_status_sync"] = "2026-08-13"
    NOTEBOOK.write_text(json.dumps(notebook, indent=1) + "\n")

    tex = TEX.read_text()
    if tex.count("version v49") != 1:
        raise ValueError("authoritative TeX does not contain exactly one v49 marker")
    tex = tex.replace("version v49", "version v50", 1)
    tex = tex.replace(
        r"\subsection{Seasonal robustness of the ambient F--K result (running)}",
        r"\subsection{Seasonal robustness of the ambient F--K result (completed)}",
        1,
    )
    stale_seasonal = r"""At the time of this guide revision, the processing is still running. The status is therefore \planned, and no seasonal persistence claim is made. The live notebook section records completion fraction and will display the day-level and across-day products once the aggregation file exists. A consistent peak velocity and branch asymmetry across seasons would support a persistent propagating wavefield; systematic seasonal changes would instead identify environmental or coupling dependence that must be included in the interpretation."""
    completed_seasonal = r"""\textbf{Status update v50: completed.} The corrected physical-lag aggregate contains 11,523 one-minute files from all eight preselected dates. The \(F K<0\) and \(F K>0\) branches both peak at 3.075~km~s\(^{-1}\), with scores 0.279 and 0.288 and receiver-order probabilities \(p=0.0002\) for each branch. These probabilities demonstrate seasonal repeatability of the fan-selected ordered observable; they remain conditional on the fixed 2.5--4.5~km~s\(^{-1}\) passband and do not independently determine velocity or physical upgoing/downgoing identity."""
    if tex.count(stale_seasonal) != 1:
        raise ValueError("could not replace stale seasonal-running paragraph")
    tex = tex.replace(stale_seasonal, completed_seasonal, 1)
    tex = tex.replace(
        r"\subsection{Independent unfiltered verification (running)}",
        r"\subsection{Independent unfiltered verification (completed)}",
        1,
    )
    stale_unfiltered = r"""The verification is intentionally resumable and produces day-level products before any across-day aggregate. Agreement between this independent full-resolution path and the F--K routine's unfiltered control would argue against an implementation bug. A strong F--K branch alongside a weak independently verified unfiltered result would instead indicate that F--K isolation materially changes observability. Until the unfiltered day-level products and permutation nulls exist, no conclusion is assigned to this branch."""
    completed_unfiltered = r"""\textbf{Status update v50: completed.} The exact-resolution aggregate contains the same 11,523 weighted one-minute records. At the fixed 3.2~km~s\(^{-1}\) trajectory, the causal score is 0.002711 with receiver-order \(p=0.341\), and the anti-causal score is 0.0000929 with \(p=0.202\). This independently confirms that the broad unfiltered correlation remains weak for the selected statistic. It does not diagnose an F--K implementation bug and does not invalidate a standard passband used to isolate low-SNR coherent energy."""
    if tex.count(stale_unfiltered) != 1:
        raise ValueError("could not replace stale unfiltered-running paragraph")
    tex = tex.replace(stale_unfiltered, completed_unfiltered, 1)
    stale_takeaway = r"""\takeaway{The filtering code is correct.  The current ambient 3.1--3.2-km-s
\( ^{-1} \) interpretation is not.  It remains a selected visualization until
an unmasked or independently constrained method detects the same event.}"""
    corrected_takeaway = r"""\takeaway{The filtering implementation is correct.  The fan-assisted
3.1--3.2-km-s\( ^{-1} \) ridge is a reproducible observable selected in an
independently motivated apparent-velocity interval. Direction-only
non-recovery constrains robustness but does not reject standard fan filtering.
The fan does not independently estimate velocity or identify physical wave
type; v50 supplies the completed validation ladder.}"""
    if tex.count(stale_takeaway) != 1:
        raise ValueError("could not replace the stale v49 TeX takeaway")
    tex = tex.replace(stale_takeaway, corrected_takeaway, 1)

    end_marker = r"\end{document}"
    if tex.count(end_marker) != 1:
        raise ValueError("could not isolate TeX document end")
    validation_tex = rf"""
\section{{Signed F--K validation ladder (v50)}}
\figquestion{{Can the signed 2.5--4.5~km~s\( ^{{-1}} \) F--K operator be
trusted, and what physical language is supported for the ambient branches?}}

F--K filtering is standard VSP wavefield separation.  Here its implementation
and physical interpretation are evaluated separately. Exact one-way
5--20~Hz synthetics give correct/wrong branch RMS ratios of
{synthetic['increasing_z']:.2f} and {synthetic['decreasing_z']:.2f} for
increasing- and decreasing-coordinate propagation. Known-direction broadband
waves injected into 300 real records are independently recovered in the
correct branch for both directions at 2.75, 3.2, and 4.0~km~s\( ^{{-1}} \),
while the 1.8~km~s\( ^{{-1}} \) off-wedge control is rejected. The known
surface AWD gives expected/opposite fixed-tube energy ratios of
{active['5_20_hz']:.2f} at 5--20~Hz and {active['25_60_hz']:.1f} at
25--60~Hz, so empirical physical-sign calibration is strong only in the higher
band. The real ambient F--K target-fan power ratio is
{ambient['pre_filter_fan_power']:.3f}, and the filtered score ratio is
{ambient['post_filter_correlation_score']:.3f}; the recorded low-frequency
field is therefore nearly balanced between coordinate directions.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.99\linewidth]{{ambient_transfer/fk_validation_evidence_v1/fk_validation_evidence_v1.png}}
\caption{{\textbf{{Validation ladder for the signed SAFOD ambient F--K operator.}} (a) Exact 10-Hz, 3.2-km-s\( ^{{-1}} \) plane waves are propagated separately toward increasing and decreasing fiber coordinate, decimated identically to the production workflow, and passed through complementary 5--20-Hz, 2.5--4.5-km-s\( ^{{-1}} \) signed masks. Correct/wrong output-RMS ratios are {synthetic['increasing_z']:.2f}:1 and {synthetic['decreasing_z']:.2f}:1, validating Fourier sign. (b) Minimum injected broadband RMS relative to median real-channel 5--20-Hz RMS required for independent pre-filter recovery after injection into 300 real one-minute records. Both directions are recovered for all in-wedge test speeds; the omitted 1.8-km-s\( ^{{-1}} \) off-wedge control is never recovered. (c) Full-resolution application to the known surface-AWD Nano stack over 80--440~m. Energy within a fixed \(\pm12\)-ms tube around the independently estimated 2.975-km-s\( ^{{-1}} \) increasing-coordinate moveout gives expected/opposite ratios {active['5_20_hz']:.2f} at 5--20~Hz and {active['25_60_hz']:.1f} at 25--60~Hz, demonstrating frequency-dependent empirical directionality. (d) The unmasked ambient target-fan power ratio is {ambient['pre_filter_fan_power']:.3f}, and the independently filtered 3.2-km-s\( ^{{-1}} \) score ratio is {ambient['post_filter_correlation_score']:.3f}. Thus the two branches are valid coordinate-direction outputs rather than a software mirror, while the physical origin of the decreasing-coordinate energy remains unresolved. The figure does not establish two body-wave Green's-function branches, wave type, or an independent formation velocity.}}
\end{{figure}}

\takeaway{{Use increasing-coordinate (\(F K<0\)) and decreasing-coordinate
(\(F K>0\)) labels for the 5--20-Hz ambient branches.  Reserve upgoing and
downgoing for a source- or reflection-calibrated physical interpretation.  The
standard fan is valid as an independently motivated passband, but cannot be
used alone to estimate velocity within that passband.}}

\section{{Published ambient-workflow transfer test (v50)}}
\figquestion{{Does the unfiltered Lellouch et al. (2019) ambient workflow
recover the causal 3.2-km-s\( ^{{-1}} \) arrival in the 2024--2025 data, and
what is gained by the previously omitted local receiver stack?}}

One complete UTC day ({lellouch['used_files']:,} one-minute files;
{lellouch['used_30_s_windows']:,} overlapping 30-s windows) is processed with
channel~0 as the fixed top virtual source, receivers every 50~m from 50 to
800~m, time differentiation to a strain-rate proxy, 5-s
running-absolute-mean normalization, 15-s window overlap, and a final
5--20-Hz correlation filter. The 5-s normalization duration is declared as a
project choice because the paper does not report it. The spatial operator
stacks correlations from the same source to receivers \(R-10\) through
\(R+10\), with a second version aligned at the published 3.2~km~s\( ^{{-1}} \).
No F--K filter is applied because none is reported in the paper's ambient
subsection.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.99\linewidth]{{ambient_transfer/lellouch2019_reproduction_v1/{lellouch_figure_path.name}}}
\caption{{\textbf{{Transfer of the Lellouch et al. (2019) ambient-interferometry workflow to the 2024--2025 SAFOD Nano acquisition.}} A complete UTC day ({lellouch['used_files']:,} contiguous one-minute files) is divided into {lellouch['used_30_s_windows']:,} 30-s windows with 15-s overlap. Nano channel~0 is the fixed top virtual source; target receiver positions increase down the cemented main-hole fiber at 50-m intervals from 50 to 800~m using measured {lellouch['channel_spacing_m']:.6f}-m spacing. Phase/strain is differentiated to a strain-rate proxy and divided by a 5-s running absolute mean before correlation; the daily stack is filtered at 5--20~Hz. Black variable-area wiggles show positive correlation lag, and the red line is a fixed 3.2-km-s\( ^{{-1}} \) reference rather than a fit. (a) Individual receiver correlations. (b) The published local stack over receiver channels \(R-10\) through \(R+10\), retaining the same channel-0 source. (c) The same 21 correlations shifted at fixed 3.2-km-s\( ^{{-1}} \) travel-time differences before stacking. Fixed-trajectory causal scores are {isolated['causal_score_at_3200_m_s']:.4f}, {simple['causal_score_at_3200_m_s']:.4f}, and {aligned['causal_score_at_3200_m_s']:.4f}; the aligned absolute score is {alignment_gain:.2f} times the isolated value. The aligned descriptive causal maximum is {aligned['descriptive_causal_peak_velocity_m_s']/1000:.3f}~km~s\( ^{{-1}} \), with causal/anti-causal absolute-score ratio {aligned['abs_causal_to_anti_causal_ratio_at_3200']:.2f} at 3.2~km~s\( ^{{-1}} \). This is a documented transfer, not a bit-for-bit reproduction of the 2017 records; acquisition epoch, interrogator/gauge settings, ambient source distribution, and channel-origin inheritance can affect the comparison.}}
\end{{figure}}

\takeaway{{The unfiltered transfer and the added F--K-assisted workflow answer
different questions.  The former tests reproducibility of the published 2017
processing on the 2024--2025 acquisition; the latter tests whether a
predeclared apparent-velocity/direction passband can recover a coherent
observable.  Neither alone provides a formation-\(V_P\) inversion.}}

"""
    tex = tex.replace(end_marker, validation_tex + end_marker)
    TEX.write_text(tex)
    TEX_V50.write_text(tex)
    print(
        json.dumps(
            {
                "version": "v50",
                "lellouch_report": str(lellouch_report_path.relative_to(ROOT)),
                "lellouch_figure": str(lellouch_figure_path.relative_to(ROOT)),
                "used_files": lellouch["used_files"],
                "used_windows": lellouch["used_30_s_windows"],
                "alignment_gain_at_3200": alignment_gain,
                "ambient_fan_power_ratio": ambient["pre_filter_fan_power"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
