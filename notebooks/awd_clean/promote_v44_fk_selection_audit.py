#!/usr/bin/env python3
"""Promote completed ambient F-K selection audits into the v44 dashboard.

This script is intentionally fail-closed.  It validates the held-out mask
sensitivity product and both n=30 and n=300 pre-F-K full-pipeline null
products before preparing either document.  The corrected seasonal signed-lag
aggregate is optional because its retry jobs may still be running.

Running this script modifies ``AWD_results_dashboard.ipynb`` and
``AWD_advisor_figure_guide.tex`` once, then writes the versioned TeX copy
``AWD_advisor_figure_guide_v44.tex``.  A second invocation refuses to promote
the same version again.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
NOTEBOOK = HERE / "AWD_results_dashboard.ipynb"
TEX = HERE / "AWD_advisor_figure_guide.tex"
TEX_VERSIONED = HERE / "AWD_advisor_figure_guide_v44.tex"
MARKER = "# v44"
PROMOTION_DATE = "2026-08-10"

TRANSFER = HERE / "ambient_transfer"
MASK_DIR = TRANSFER / "fk_mask_sensitivity_v2"
MASK_JSON = MASK_DIR / "ambient_fk_mask_sensitivity_v2.json"
MASK_PNG = MASK_DIR / "ambient_fk_mask_sensitivity_v2.png"

N30_DIR = TRANSFER / "fk_full_pipeline_null_v2_pilot_n30_r20"
N30_JSON = N30_DIR / "fk_full_pipeline_null_v2_aggregate.json"
N30_PNG = N30_DIR / "fk_full_pipeline_null_v2_aggregate_figure.png"

N300_DIR = TRANSFER / "fk_full_pipeline_null_v2_n300_r20"
N300_JSON = N300_DIR / "fk_full_pipeline_null_v2_aggregate.json"
N300_PNG = N300_DIR / "fk_full_pipeline_null_v2_aggregate_figure.png"

SEASONAL_DIR = TRANSFER / "signed_lag_v2"
SEASONAL_JSON = SEASONAL_DIR / "seasonal_signed_fk_v2_aggregate.json"
SEASONAL_PNG = SEASONAL_DIR / "seasonal_signed_fk_v2_aggregate.png"
SEASONAL_ARRAY_JOB = "38562360"
SEASONAL_POSTPROCESS_JOB = "38562403"


def source_text(value: Any) -> str:
    """Flatten a notebook source field for marker checks."""
    if isinstance(value, str):
        return value
    if isinstance(value, list):
        return "".join(source_text(item) for item in value)
    return str(value)


def markdown(source: str) -> dict[str, Any]:
    """Return a clean Markdown notebook cell."""
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": (source.rstrip() + "\n").splitlines(True),
    }


def code(source: str) -> dict[str, Any]:
    """Return an unexecuted code notebook cell."""
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": (source.rstrip() + "\n").splitlines(True),
    }


def read_json(path: Path) -> dict[str, Any]:
    """Load one JSON object and reject missing, empty, or malformed products."""
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"required JSON product is missing or empty: {path}")
    try:
        value = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"malformed JSON product: {path}") from exc
    if not isinstance(value, dict):
        raise TypeError(f"expected a JSON object in {path}")
    return value


def require_png(path: Path) -> None:
    """Require a nontrivial PNG before any document is modified."""
    if not path.is_file() or path.stat().st_size < 1024:
        raise FileNotFoundError(f"required PNG product is missing or too small: {path}")
    if path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"required figure is not a valid PNG: {path}")


def require_keys(mapping: dict[str, Any], keys: tuple[str, ...], context: str) -> None:
    """Raise a useful error for an incomplete analysis report."""
    missing = [key for key in keys if key not in mapping]
    if missing:
        raise KeyError(f"{context} is missing keys: {', '.join(missing)}")


def validate_mask_report(report: dict[str, Any]) -> None:
    """Validate the frozen development/held-out mask comparison."""
    if report.get("workflow_version") != "ambient_fk_mask_sensitivity_v2_aggregate":
        raise ValueError("unexpected held-out mask-sensitivity workflow version")
    if report.get("partial") is not False or report.get("split_is_frozen") is not True:
        raise ValueError("mask sensitivity must be complete and use the frozen split")
    heldout_dates = report.get("heldout_dates")
    available_dates = report.get("available_heldout_dates")
    if not isinstance(heldout_dates, list) or len(heldout_dates) != 7:
        raise ValueError("mask sensitivity must declare exactly seven held-out dates")
    if set(available_dates or []) != set(heldout_dates):
        raise ValueError("not all seven held-out dates are available")
    heldout = report.get("group_results", {}).get("held_out", {})
    for mask in (
        "production_2p5_4p5",
        "narrow_2p8_3p8",
        "broad_2p0_5p5",
        "direction_only",
    ):
        for branch in ("negative", "positive"):
            metrics = heldout.get(mask, {}).get(branch, {})
            require_keys(
                metrics,
                (
                    "conditional_peak_velocity_m_s",
                    "conditional_peak_absolute_score",
                    "absolute_score_3200",
                    "conditional_p_fixed_3200",
                    "conditional_p_peak_scan",
                ),
                f"held-out {mask}/{branch}",
            )


def validate_null_report(report: dict[str, Any], expected_files: int, label: str) -> None:
    """Validate one complete pre-F-K full-pipeline null aggregate."""
    if report.get("workflow_version") != "ambient_fk_full_pipeline_null_v2_aggregate":
        raise ValueError(f"unexpected workflow version for {label}")
    if int(report.get("requested_files", -1)) != expected_files:
        raise ValueError(f"{label} requested_files is not {expected_files}")
    if int(report.get("used_files", -1)) != expected_files:
        raise ValueError(f"{label} used_files is not {expected_files}")
    n_null = int(report.get("null_realizations", 0))
    ids = report.get("null_realization_ids")
    if n_null < 1 or not isinstance(ids, list) or len(ids) != n_null:
        raise ValueError(f"{label} has an incomplete null-realization inventory")
    if len(set(ids)) != n_null:
        raise ValueError(f"{label} contains duplicate null realization IDs")
    methods = report.get("null_methods")
    if set(methods or []) != {"channel_permutation", "circular_time_shift"}:
        raise ValueError(f"{label} must contain both pre-F-K null constructions")
    for method in methods:
        result = report.get("null_results", {}).get(method, {})
        require_keys(result.get("familywise", {}), ("null95", "p"), f"{label}/{method}")
        for branch in ("negative", "positive"):
            require_keys(
                result.get("branches", {}).get(branch, {}),
                (
                    "observed_peak_absolute_score_in_wedge",
                    "observed_peak_velocity_m_s_in_wedge",
                    "null95_peak_score_in_wedge",
                    "p_peak_in_wedge",
                ),
                f"{label}/{method}/{branch}",
            )


def validate_seasonal_report(report: dict[str, Any]) -> None:
    """Validate the optional corrected causal/anti-causal seasonal aggregate."""
    if report.get("workflow_version") != "signed_lag_v2":
        raise ValueError("unexpected corrected seasonal workflow version")
    if len(report.get("days", [])) != 8:
        raise ValueError("corrected seasonal aggregate must contain eight dates")
    if int(report.get("weighted_files", 0)) < 1:
        raise ValueError("corrected seasonal aggregate has no weighted files")
    for branch in ("negative", "positive"):
        require_keys(
            report.get("branches", {}).get(branch, {}),
            (
                "physical_lag_sign",
                "peak_velocity_m_s",
                "peak_absolute_score",
                "absolute_score_3200",
                "absolute_opposite_lag_leakage_3200",
                "null95",
                "p_peak",
            ),
            f"corrected seasonal/{branch}",
        )


def p_string(value: float) -> str:
    """Format a probability without implying unsupported precision."""
    value = float(value)
    if value < 0.001:
        return f"{value:.4f}"
    if value < 0.01:
        return f"{value:.3f}"
    return f"{value:.3g}"


def null_summary(report: dict[str, Any]) -> str:
    """Return a compact advisor-facing summary of one null ensemble."""
    pieces = []
    for method in ("channel_permutation", "circular_time_shift"):
        result = report["null_results"][method]
        neg = result["branches"]["negative"]["p_peak_in_wedge"]
        pos = result["branches"]["positive"]["p_peak_in_wedge"]
        familywise = result["familywise"]["p"]
        label = "channel permutation" if method == "channel_permutation" else "circular time shift"
        pieces.append(
            f"{label}: negative p = {p_string(neg)}, positive p = {p_string(pos)}, "
            f"familywise p = {p_string(familywise)}"
        )
    return "; ".join(pieces)


def passes_familywise_nulls(report: dict[str, Any], alpha: float = 0.05) -> bool:
    """Return whether both predeclared full-pipeline familywise nulls are rejected."""
    return all(
        float(report["null_results"][method]["familywise"]["p"]) <= alpha
        for method in ("channel_permutation", "circular_time_shift")
    )


def notebook_root_code(relative_path: str, variable: str) -> str:
    """Build a self-contained, working-directory-tolerant display cell."""
    return f'''from pathlib import Path
from IPython.display import Image, display

cwd = Path.cwd()
dashboard_root = next(
    (candidate for candidate in (cwd, cwd / "awd_clean")
     if (candidate / "AWD_results_dashboard.ipynb").exists()),
    None,
)
if dashboard_root is None:
    raise FileNotFoundError("Run this notebook from awd_clean/ or its parent directory")
{variable} = dashboard_root / "{relative_path}"
if not {variable}.is_file():
    raise FileNotFoundError({variable})
display(Image(filename=str({variable}), width=1200))'''


def seasonal_notebook_code() -> str:
    """Build a runtime-aware seasonal status/display cell."""
    return f'''from pathlib import Path
import json
from IPython.display import Image, display

cwd = Path.cwd()
dashboard_root = next(
    (candidate for candidate in (cwd, cwd / "awd_clean")
     if (candidate / "AWD_results_dashboard.ipynb").exists()),
    None,
)
if dashboard_root is None:
    raise FileNotFoundError("Run this notebook from awd_clean/ or its parent directory")
seasonal_dir = dashboard_root / "ambient_transfer" / "signed_lag_v2"
seasonal_json = seasonal_dir / "seasonal_signed_fk_v2_aggregate.json"
seasonal_png = seasonal_dir / "seasonal_signed_fk_v2_aggregate.png"
if seasonal_json.is_file() and seasonal_png.is_file():
    seasonal = json.loads(seasonal_json.read_text())
    print("Corrected eight-day signed-lag aggregate")
    print(json.dumps(seasonal, indent=2))
    display(Image(filename=str(seasonal_png), width=1200))
else:
    print(
        "Corrected seasonal retry still pending: SLURM array {SEASONAL_ARRAY_JOB}; "
        "after-success postprocess {SEASONAL_POSTPROCESS_JOB}."
    )'''


def make_notebook_cells(
    mask_report: dict[str, Any],
    n30_report: dict[str, Any],
    n300_report: dict[str, Any],
    seasonal_report: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    """Construct the additive v44 notebook section."""
    heldout = mask_report["group_results"]["held_out"]
    prod_neg = heldout["production_2p5_4p5"]["negative"]
    prod_pos = heldout["production_2p5_4p5"]["positive"]
    narrow_neg = heldout["narrow_2p8_3p8"]["negative"]
    broad_neg = heldout["broad_2p0_5p5"]["negative"]
    broad_pos = heldout["broad_2p0_5p5"]["positive"]
    direction_neg = heldout["direction_only"]["negative"]
    direction_pos = heldout["direction_only"]["positive"]
    n300_passes = passes_familywise_nulls(n300_report)

    if n300_passes:
        evidence_flag = (
            "**Current evidence flag — provisionally retained.** The production-wedge "
            "observable repeats on seven frozen held-out dates and the five-hour observed "
            "familywise statistic exceeds both stated pre-F–K surrogate ensembles. This "
            "supports preserved directional coherence under the tested nulls, but the "
            "approximately 3.1 km/s ridge remains mask-selected and is not an independent "
            "formation-velocity estimate."
        )
        null_decision = (
            "The five-hour test rejects both predeclared familywise nulls at the 0.05 "
            "threshold. The selected ambient observable is therefore retained for a "
            "held-out replication test; this is not yet mode identification or a velocity "
            "measurement."
        )
    else:
        evidence_flag = (
            "**Current evidence flag — withdrawn as physical-wave evidence.** The "
            "production-wedge observable repeats on seven frozen held-out dates, but the "
            "five-hour full-pipeline test does not reject both coherence-destroying nulls. "
            "The approximately 3.1 km/s ridge is retained only as a processing-sensitive "
            "diagnostic and must not be presented as ambient Green's-function convergence, "
            "a propagating mode, or formation velocity."
        )
        null_decision = (
            "The five-hour test fails the predeclared requirement that both familywise "
            "nulls be rejected at the 0.05 threshold. Ambient F–K escalation therefore "
            "stops here unless a new, independently justified observable is defined; longer "
            "stacking of the same selected statistic would not repair this failed test."
        )

    if seasonal_report is None:
        seasonal_status = (
            f"The corrected eight-day signed-lag aggregate is still pending under retry "
            f"array `{SEASONAL_ARRAY_JOB}` and after-success postprocess "
            f"`{SEASONAL_POSTPROCESS_JOB}`. The following code cell checks dynamically "
            "for the completed product each time the notebook is run."
        )
        seasonal_caption = (
            "**Status caption.** The corrected seasonal comparison is pending. When available, "
            "the figure will place the `F×K < 0` branch at positive lag and the `F×K > 0` "
            "branch at negative lag, use a common display scale, and report branch-specific "
            "permutation nulls. No causal/anti-causal amplitude asymmetry is inferred while "
            f"retry jobs `{SEASONAL_ARRAY_JOB}` and `{SEASONAL_POSTPROCESS_JOB}` are incomplete."
        )
    else:
        neg = seasonal_report["branches"]["negative"]
        pos = seasonal_report["branches"]["positive"]
        seasonal_status = (
            f"The corrected eight-day signed-lag aggregate contains "
            f"{int(seasonal_report['weighted_files']):,} weighted one-minute files. The "
            f"`F×K < 0`/positive-lag branch peaks at {neg['peak_velocity_m_s']/1000:.3f} "
            f"km/s with p = {p_string(neg['p_peak'])}; the `F×K > 0`/negative-lag branch "
            f"peaks at {pos['peak_velocity_m_s']/1000:.3f} km/s with "
            f"p = {p_string(pos['p_peak'])}."
        )
        seasonal_caption = (
            "**Publication-style caption. Corrected seasonal signed-lag comparison.** "
            f"Eight selected dates ({int(seasonal_report['weighted_files']):,} weighted "
            "one-minute files) are processed with the channel-0 virtual source, 5–20 Hz "
            "band, running-absolute-mean temporal normalization, and the frozen 2.5–4.5 "
            "km/s signed F–K wedge. The `F×K < 0` branch is evaluated on positive correlation "
            "lag and the `F×K > 0` branch on negative lag, following the synthetic sign and "
            "lag-reversal tests. Dashed moveout lines are references, not independent fits. "
            "The figure measures branch-specific correlation coherence and opposite-lag "
            "leakage; it does not by itself establish Green's-function convergence, wave "
            "type, causality, or formation velocity."
        )

    return [
        markdown(
            f"""# v44

## Ambient transfer validation

{evidence_flag}

All products below use Nano channel 0 as the virtual source, receiver offset increasing along the fiber, a 5–20 Hz band, and the same signed F–K convention established by the deterministic synthetic tests. Under the project convention, `F×K < 0` is scored on positive lag for propagation toward increasing fiber coordinate; `F×K > 0` is scored on negative lag for the opposite coordinate direction. The seven held-out dates were frozen before this comparison. The older same-positive-lag statement that the positive branch was absent is superseded."""
        ),
        markdown(
            f"""### Does the result survive held-out dates and different F–K masks?

For the seven held-out dates, the production wedge peaks at {prod_neg['conditional_peak_velocity_m_s']/1000:.3f} km/s on the negative branch and {prod_pos['conditional_peak_velocity_m_s']/1000:.3f} km/s on the positive branch. At 3.2 km/s the equal-day scores are {prod_neg['absolute_score_3200']:.3f} and {prod_pos['absolute_score_3200']:.3f}, respectively, with conditional receiver-order permutation probabilities p = {p_string(prod_neg['conditional_p_fixed_3200'])} and p = {p_string(prod_pos['conditional_p_fixed_3200'])}.

The dependence on the mask is scientifically decisive. The narrow 2.8–3.8 km/s wedge forces the negative-branch conditional peak to {narrow_neg['conditional_peak_velocity_m_s']/1000:.3f} km/s and raises its score to {narrow_neg['conditional_peak_absolute_score']:.3f}. The broad 2.0–5.5 km/s wedge instead peaks at {broad_neg['conditional_peak_velocity_m_s']/1000:.3f} km/s and {broad_pos['conditional_peak_velocity_m_s']/1000:.3f} km/s; its 3.2 km/s fixed-velocity probabilities are p = {p_string(broad_neg['conditional_p_fixed_3200'])} and p = {p_string(broad_pos['conditional_p_fixed_3200'])}. With direction-only masks, the 3.2 km/s scores fall to {direction_neg['absolute_score_3200']:.4f} and {direction_pos['absolute_score_3200']:.4f}. Thus the production observable generalizes across time, but the approximately 3.1 km/s estimate is conditional on velocity selection."""
        ),
        code(
            notebook_root_code(
                "ambient_transfer/fk_mask_sensitivity_v2/ambient_fk_mask_sensitivity_v2.png",
                "mask_figure",
            )
        ),
        markdown(
            """**Publication-style caption. Held-out sensitivity of the ambient transfer observable to the F–K mask.** The development date (2024-12-20) is separated from seven preselected held-out dates spanning the available seasonal archive. Each date is processed independently with Nano channel 0 as the virtual source, 5–20 Hz preprocessing and temporal normalization; dates receive equal weight in the held-out aggregate. Columns compare the frozen production wedge (2.5–4.5 km/s), a narrow wedge (2.8–3.8 km/s), a broad wedge (2.0–5.5 km/s), and signed-direction-only filtering with no velocity restriction. Negative and positive F–K branches are evaluated at positive and negative physical lags, respectively. Receiver-order permutations quantify ordered moveout only after F–K selection. Persistence of the production-wedge score on held-out dates establishes temporal repeatability of that selected observable. The large changes in peak position and score across masks show that the approximately 3.1 km/s ridge is not an independent or unbiased estimate of formation velocity; the receiver-order null also is not a full-pipeline test of mask selection."""
        ),
        markdown(
            f"""### Does the fixed wedge create a ridge from coherence-destroyed data?

The full-pipeline null modifies the wavefield *before* F–K filtering and then repeats the entire production sequence: signed 2.5–4.5 km/s filtering, channel-0 correlation, stacking, and velocity scanning. Channel permutation destroys spatial order while retaining individual trace spectra; independent circular time shifts disrupt interchannel phase alignment while approximately retaining each trace's spectrum and amplitude distribution. These are more stringent tests of F–K selection than permuting receivers after a filtered correlation section.

For the 30-file pilot, {null_summary(n30_report)}. For the 300-file test, {null_summary(n300_report)}. A large p value means the observed maximum is common among these coherence-destroyed surrogates and therefore does **not** reject the relevant null; it is not evidence that the physical signal is absent.

**Decision:** {null_decision}"""
        ),
        code(
            notebook_root_code(
                "ambient_transfer/fk_full_pipeline_null_v2_pilot_n30_r20/fk_full_pipeline_null_v2_aggregate_figure.png",
                "null_n30_figure",
            )
        ),
        markdown(
            f"""**Publication-style caption. Thirty-minute pre-F–K full-pipeline null pilot.** The observed branch-specific peak median-correlation statistics and the maximum across both branches are compared with {int(n30_report['null_realizations'])} channel-permutation and {int(n30_report['null_realizations'])} independent-circular-time-shift realizations. Each surrogate is generated before F–K filtering and is passed through the same channel-0 virtual-source workflow, 5–20 Hz preprocessing, frozen 2.5–4.5 km/s signed wedge, physical-lag correlation, stack, and velocity scan as the unmodified 30-file record. Vertical black lines show observed statistics, dotted lines show the surrogate 95th percentiles, and empirical probabilities include a plus-one correction. The familywise test includes both signed branches and the full velocity scan. The pilot yields {null_summary(n30_report)}; it therefore does not reject either coherence-destroying null and cannot by itself establish that the selected ridge originates from preserved wave propagation."""
        ),
        code(
            notebook_root_code(
                "ambient_transfer/fk_full_pipeline_null_v2_n300_r20/fk_full_pipeline_null_v2_aggregate_figure.png",
                "null_n300_figure",
            )
        ),
        markdown(
            f"""**Publication-style caption. Five-hour pre-F–K full-pipeline null test.** This analysis repeats the complete null construction for 300 one-minute files, increasing the stack duration tenfold relative to the pilot while preserving the same Nano channel-0 geometry, 5–20 Hz band, temporal normalization, 2.5–4.5 km/s signed wedge, physical-lag definitions, and branch-plus-velocity familywise statistic. Channel permutations test dependence on spatial ordering, whereas independent circular shifts test dependence on interchannel phase alignment; neither surrogate represents every possible spatially correlated ambient source field. The recorded results are {null_summary(n300_report)}. These probabilities test whether the observed selected coherence exceeds these explicitly constructed surrogates. They do not identify the wave type, prove Green's-function convergence, or turn the conditional ridge peak into formation velocity."""
        ),
        markdown(
            f"""### Corrected seasonal signed-lag status

{seasonal_status}"""
        ),
        code(seasonal_notebook_code()),
        markdown(seasonal_caption),
    ]


def replace_once(text: str, old: str, new: str, label: str) -> str:
    """Replace one expected block and fail if the source has drifted."""
    count = text.count(old)
    if count != 1:
        raise ValueError(f"expected exactly one {label}; found {count}")
    return text.replace(old, new, 1)


def tex_section(
    mask_report: dict[str, Any],
    n30_report: dict[str, Any],
    n300_report: dict[str, Any],
    seasonal_report: dict[str, Any] | None,
) -> str:
    """Construct the discreet advisor-facing ambient validation section."""
    heldout = mask_report["group_results"]["held_out"]
    prod_neg = heldout["production_2p5_4p5"]["negative"]
    prod_pos = heldout["production_2p5_4p5"]["positive"]
    narrow_neg = heldout["narrow_2p8_3p8"]["negative"]
    broad_neg = heldout["broad_2p0_5p5"]["negative"]
    broad_pos = heldout["broad_2p0_5p5"]["positive"]
    direction_neg = heldout["direction_only"]["negative"]
    direction_pos = heldout["direction_only"]["positive"]
    n300_passes = passes_familywise_nulls(n300_report)

    if n300_passes:
        decision_tex = (
            "The five-hour result rejects both predeclared familywise nulls at "
            "$\\alpha=0.05$.  The selected ambient observable is therefore retained for "
            "held-out replication, but not yet as mode identification or a formation-velocity "
            "measurement."
        )
        takeaway_tex = (
            "The production-wedge observable repeats on seven held-out dates and exceeds both "
            "stated five-hour pre-F--K surrogate ensembles.  It is provisionally retained as "
            "directional coherence under those tests.  Its approximately "
            "3.1~km~s$^{-1}$ peak remains mask selected and is not formation $V_P$."
        )
    else:
        decision_tex = (
            "The five-hour result fails the predeclared requirement that both familywise "
            "nulls be rejected at $\\alpha=0.05$.  The ambient ridge is therefore withdrawn "
            "as physical-wave evidence and retained only as a processing-sensitive diagnostic. "
            "Longer stacking of the same selected statistic is not treated as a remedy."
        )
        takeaway_tex = (
            "The production-wedge observable repeats on seven held-out dates, but it does not "
            "exceed both stated five-hour pre-F--K surrogate ensembles.  It is therefore "
            "withdrawn as physical-wave evidence.  The approximately 3.1~km~s$^{-1}$ feature "
            "must not be presented as Green's-function convergence, a propagating mode, or "
            "formation $V_P$."
        )

    if seasonal_report is None:
        seasonal_tex = rf"""
\subsection{{Corrected seasonal signed-lag status}}
The corrected eight-day causal/anti-causal aggregation remains pending under
retry array {SEASONAL_ARRAY_JOB} and after-success postprocess
{SEASONAL_POSTPROCESS_JOB}.  Consequently, no seasonal branch-amplitude
asymmetry is promoted in v44.  This status is deliberately separate from the
completed held-out mask-sensitivity test, which used its own versioned
products.
"""
    else:
        neg = seasonal_report["branches"]["negative"]
        pos = seasonal_report["branches"]["positive"]
        seasonal_tex = rf"""
\subsection{{Corrected seasonal signed-lag comparison}}
The corrected aggregate contains {int(seasonal_report['weighted_files']):,}
weighted one-minute files.  The $F K<0$ branch evaluated at positive lag peaks
at {neg['peak_velocity_m_s']/1000:.3f}~km~s$^{{-1}} with
$p={p_string(neg['p_peak'])}$; the $F K>0$ branch evaluated at negative lag
peaks at {pos['peak_velocity_m_s']/1000:.3f}~km~s$^{{-1}} with
$p={p_string(pos['p_peak'])}$.  These are branch-specific selected coherence
diagnostics, not independent estimates of formation velocity.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.98\linewidth]{{ambient_transfer/signed_lag_v2/seasonal_signed_fk_v2_aggregate.png}}
\caption{{\textbf{{Corrected eight-day causal and anti-causal signed-lag comparison.}} Eight selected dates ({int(seasonal_report['weighted_files']):,} weighted one-minute files) are processed with Nano channel~0 as the virtual source, a 5--20~Hz band, running-absolute-mean temporal normalization, and the frozen 2.5--4.5~km~s$^{{-1}}$ signed F--K wedge. The $F K<0$ branch is scored at positive correlation lag and the $F K>0$ branch at negative lag, as established by the synthetic direction and impulse-lag tests. A common display scale permits visual comparison, while the lower panels report branch-specific velocity scores, receiver-order null levels, and opposite-lag leakage. Dashed moveout lines are references rather than independent fits. The figure measures corrected physical-lag coherence; it does not establish Green's-function convergence, wave type, causality, or formation $V_P$.}}
\end{{figure}}
"""

    return rf"""
\section{{Ambient transfer validation}}

\subsection{{Held-out mask sensitivity}}
\figquestion{{Does the ambient transfer observable persist on dates excluded
from development, and is its apparent velocity independent of the selected
F--K mask?}}

Nano channel~0 is the virtual source.  All eight dates use a 5--20~Hz band and
the same temporal normalization.  The development date is 2024-12-20; the
remaining seven dates were frozen as held out before aggregation.  Under the
documented Fourier and coordinate convention, $F K<0$ is evaluated at positive
lag for propagation toward increasing fiber coordinate, whereas $F K>0$ is
evaluated at negative lag for the opposite coordinate direction.

On the held-out dates, the production 2.5--4.5~km~s$^{{-1}}$ wedge peaks at
{prod_neg['conditional_peak_velocity_m_s']/1000:.3f} and
{prod_pos['conditional_peak_velocity_m_s']/1000:.3f}~km~s$^{{-1}}$ on the
negative and positive branches.  Their equal-day scores at
3.2~km~s$^{{-1}}$ are {prod_neg['absolute_score_3200']:.3f} and
{prod_pos['absolute_score_3200']:.3f}, with conditional receiver-order
probabilities $p={p_string(prod_neg['conditional_p_fixed_3200'])}$ and
$p={p_string(prod_pos['conditional_p_fixed_3200'])}$.  The observable therefore
repeats beyond the development date under the frozen production processing.

The inferred velocity is not independent of the selection.  The narrow
2.8--3.8~km~s$^{{-1}}$ wedge moves the negative-branch conditional peak to
{narrow_neg['conditional_peak_velocity_m_s']/1000:.3f}~km~s$^{{-1}}$ and raises
its score to {narrow_neg['conditional_peak_absolute_score']:.3f}.  The broad
2.0--5.5~km~s$^{{-1}}$ wedge instead peaks at
{broad_neg['conditional_peak_velocity_m_s']/1000:.3f} and
{broad_pos['conditional_peak_velocity_m_s']/1000:.3f}~km~s$^{{-1}}$, and its
fixed 3.2~km~s$^{{-1}}$ probabilities are
$p={p_string(broad_neg['conditional_p_fixed_3200'])}$ and
$p={p_string(broad_pos['conditional_p_fixed_3200'])}$.  Direction-only scores
at 3.2~km~s$^{{-1}}$ fall to {direction_neg['absolute_score_3200']:.4f} and
{direction_pos['absolute_score_3200']:.4f}.  Thus the measured result is a
repeatable, directionally and kinematically selected correlation observable;
the approximately 3.1~km~s$^{{-1}}$ ridge is not an unbiased formation-velocity
measurement.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.99\linewidth]{{ambient_transfer/fk_mask_sensitivity_v2/ambient_fk_mask_sensitivity_v2.png}}
\caption{{\textbf{{Held-out sensitivity of the ambient transfer observable to the F--K mask.}} The 2024-12-20 development date is separated from seven preselected held-out dates spanning the available seasonal archive. Each date is processed independently with Nano channel~0 as the virtual source, 5--20~Hz preprocessing, and running-absolute-mean temporal normalization; dates receive equal weight in the held-out aggregate. Columns compare the frozen production wedge (2.5--4.5~km~s$^{{-1}}$), a narrow wedge (2.8--3.8~km~s$^{{-1}}$), a broad wedge (2.0--5.5~km~s$^{{-1}}$), and signed-direction-only filtering without a velocity restriction. The $F K<0$ and $F K>0$ branches are evaluated on positive and negative physical lags, respectively. Receiver-order permutations quantify ordered moveout after F--K selection; they are not full-pipeline selection nulls. Persistence of the production-wedge score establishes temporal repeatability of the selected observable, whereas the large mask-dependent changes in ridge position and strength demonstrate that the approximately 3.1~km~s$^{{-1}}$ value is not an independent or unbiased estimate of formation velocity.}}
\end{{figure}}

\subsection{{Pre-F--K full-pipeline null tests}}
\figquestion{{Does the frozen velocity wedge produce an apparently coherent
ridge after the input wavefield's spatial or phase coherence has been
destroyed?}}

The nulls are applied before F--K selection.  Channel permutation preserves
individual trace spectra but destroys the receiver ordering.  Independent
circular time shifts approximately preserve each trace's spectrum and
amplitude distribution but disrupt interchannel phase alignment.  Every
surrogate is then passed through the same 5--20~Hz processing, frozen
2.5--4.5~km~s$^{{-1}}$ signed wedge, channel-0 correlation, stacking, and
velocity scan.  The familywise statistic is the maximum across both physical
lag branches and the scanned velocity range.

For the 30-file pilot, {null_summary(n30_report)}.  A large probability means
that the observed maximum is common among the stated surrogates and therefore
does not reject that null; it does not show that a physical signal is absent.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.99\linewidth]{{ambient_transfer/fk_full_pipeline_null_v2_pilot_n30_r20/fk_full_pipeline_null_v2_aggregate_figure.png}}
\caption{{\textbf{{Thirty-minute pre-F--K full-pipeline null pilot.}} The observed negative-branch, positive-branch, and branch-plus-velocity familywise peak median-correlation statistics are compared with {int(n30_report['null_realizations'])} channel-permutation and {int(n30_report['null_realizations'])} independent-circular-time-shift realizations. Each surrogate is constructed before F--K filtering and passes through the same Nano channel-0 virtual-source workflow, 5--20~Hz preprocessing, frozen 2.5--4.5~km~s$^{{-1}}$ signed wedge, physical-lag correlation, stack, and velocity scan as the unmodified 30-file record. Black vertical lines show observed statistics, dotted lines show surrogate 95th percentiles, and one-sided empirical probabilities use a plus-one correction. The pilot does not reject either coherence-destroying null. It therefore cannot establish that the selected ridge depends on preserved wave propagation, and it does not identify wave type or formation velocity.}}
\end{{figure}}

The 300-file test gives {null_summary(n300_report)}.  It increases stack duration
tenfold while retaining the identical mask and test statistic, so differences
from the pilot describe stack-length sensitivity rather than a changed
selection rule.

\textbf{{Decision.}} {decision_tex}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.99\linewidth]{{ambient_transfer/fk_full_pipeline_null_v2_n300_r20/fk_full_pipeline_null_v2_aggregate_figure.png}}
\caption{{\textbf{{Five-hour pre-F--K full-pipeline null test.}} The full null construction is repeated for 300 one-minute files with the same Nano channel-0 geometry, 5--20~Hz band, temporal normalization, 2.5--4.5~km~s$^{{-1}}$ signed wedge, branch-specific physical lag conventions, and maximum-across-branches-and-velocity familywise statistic used in the pilot. Channel permutations test reliance on spatial ordering; independent circular shifts test reliance on interchannel phase alignment. Both surrogate classes precede F--K filtering, so this is a direct test of whether the fixed selection can recover similar statistics from these coherence-destroyed inputs. The nulls do not span every possible spatially correlated ambient source field. The figure constrains processing-selection risk; it does not prove Green's-function convergence, identify the mode, establish causality, or convert a conditional ridge peak into formation $V_P$.}}
\end{{figure}}

\takeaway{{{takeaway_tex}}}

{seasonal_tex.rstrip()}
"""


def update_tex(
    original: str,
    section: str,
    seasonal_report: dict[str, Any] | None,
) -> str:
    """Apply narrow corrections and append the v44 section before the document end."""
    if "version v44" in original or "\\section{Ambient transfer validation}" in original:
        raise SystemExit("v44 TeX marker already present; refusing a second promotion")
    if original.count("\\end{document}") != 1:
        raise ValueError("expected exactly one \\end{document}")
    text = replace_once(original, "version v43", "version v44", "document version")
    text, date_count = re.subn(r"\\date\{[^}]*\}", r"\\date{10 August 2026}", text, count=1)
    if date_count != 1:
        raise ValueError("could not update the TeX date")

    text = replace_once(
        text,
        "the positive branch is not significant, while the two-sided branch remains enhanced.",
        "the legacy positive-mask curve was evaluated at the same positive lag and is superseded as a physical opposite-branch test, while the two-sided branch remains enhanced.",
        "stale positive-branch statement",
    )
    text = replace_once(
        text,
        "\\subsection{Completed eight-day seasonal F--K validation and unfiltered control}",
        "\\subsection{Superseded legacy eight-day F--K summary and unfiltered control}",
        "legacy seasonal subsection title",
    )
    text = replace_once(
        text,
        "\\caption{\\textbf{Eight-day seasonal repeatability of the negative signed F--K correlation section.}",
        "\\caption{\\textbf{Superseded legacy view of eight-day negative-mask positive-lag repeatability.}",
        "legacy seasonal figure caption label",
    )

    old_running = (
        "The corrected routine now passes three deterministic regression tests: a delayed impulse returns "
        "$+0.20$~s, the reversed pair returns $-0.20$~s, and broadband waves traveling toward increasing "
        "and decreasing fiber coordinate are retained by opposite signed F--K masks and peak at opposite "
        "correlation lags. A one-file real-data smoke test shows comparable candidate branches near "
        "3.1~km~s$^{-1}$ on opposite lag signs. That smoke result is diagnostic only. The versioned "
        "eight-day recomputation is running as SLURM array 37717572, with after-success aggregation job "
        "37717613. No physical causal/anti-causal asymmetry is claimed until those jobs complete."
    )
    if seasonal_report is None:
        new_status = (
            "The corrected routine passes three deterministic regression tests: delayed and reversed "
            "impulses return opposite lags, and broadband waves traveling toward increasing and "
            "decreasing fiber coordinate are retained by opposite signed F--K masks and peak at opposite "
            f"correlation lags. The corrected seasonal retry is currently array {SEASONAL_ARRAY_JOB}, "
            f"with after-success postprocess {SEASONAL_POSTPROCESS_JOB}. No seasonal physical "
            "causal/anti-causal asymmetry is claimed while those products remain pending."
        )
    else:
        neg = seasonal_report["branches"]["negative"]
        pos = seasonal_report["branches"]["positive"]
        new_status = (
            "The corrected routine passes three deterministic regression tests: delayed and reversed "
            "impulses return opposite lags, and broadband waves traveling toward increasing and "
            "decreasing fiber coordinate are retained by opposite signed F--K masks and peak at opposite "
            "correlation lags. The completed corrected aggregate evaluates the $F K<0$ branch at positive "
            f"lag (peak {neg['peak_velocity_m_s']/1000:.3f}~km~s$^{{-1}}$, "
            f"$p={p_string(neg['p_peak'])}$) and the $F K>0$ branch at negative lag "
            f"(peak {pos['peak_velocity_m_s']/1000:.3f}~km~s$^{{-1}}$, "
            f"$p={p_string(pos['p_peak'])}$)."
        )
    text = replace_once(text, old_running, new_status, "stale v43 running-job paragraph")

    old_takeaway = (
        "\\takeaway{v43 preserves the repeatable negative-mask positive-lag observable near "
        "3.075~km~s$^{-1}$, withdraws any claim that the physical opposite branch is absent, and records "
        "the versioned corrected seasonal rerun. The result remains conditional on the selected F--K "
        "wedge and is not an independent formation-velocity estimate.}"
    )
    new_takeaway = (
        "\\takeaway{This legacy summary is retained for provenance but is superseded for branch physics "
        "and F--K-selection inference by the v44 ambient transfer validation below. The old positive-mask "
        "same-positive-lag result must not be interpreted as absence of the physical opposite branch.}"
    )
    text = replace_once(text, old_takeaway, new_takeaway, "legacy v43 takeaway")

    return text.replace("\\end{document}", section.rstrip() + "\n\n\\end{document}", 1)


def atomic_write(path: Path, content: str) -> None:
    """Write one prepared text product through a sibling temporary file."""
    temporary = path.with_name(f".{path.name}.v44.tmp")
    temporary.write_text(content)
    os.replace(temporary, path)


def main() -> None:
    """Validate all required products and promote the authoritative documents once."""
    if not NOTEBOOK.is_file() or not TEX.is_file():
        raise FileNotFoundError("authoritative notebook or TeX guide is missing")
    if TEX_VERSIONED.exists():
        raise FileExistsError(f"versioned TeX already exists: {TEX_VERSIONED}")

    notebook = json.loads(NOTEBOOK.read_text())
    if not isinstance(notebook.get("cells"), list):
        raise ValueError("authoritative notebook has no cell list")
    if any(MARKER in source_text(cell.get("source", [])) for cell in notebook["cells"]):
        raise SystemExit("v44 notebook marker already present; refusing a second promotion")
    current_version = notebook.get("metadata", {}).get("awd_dashboard", {}).get("version")
    if current_version != "v43":
        raise ValueError(f"expected notebook metadata version v43, found {current_version!r}")

    # Fail closed before constructing or writing either document.
    for figure in (MASK_PNG, N30_PNG, N300_PNG):
        require_png(figure)
    mask_report = read_json(MASK_JSON)
    n30_report = read_json(N30_JSON)
    n300_report = read_json(N300_JSON)
    validate_mask_report(mask_report)
    validate_null_report(n30_report, 30, "n30 full-pipeline null")
    validate_null_report(n300_report, 300, "n300 full-pipeline null")

    seasonal_json_exists = SEASONAL_JSON.is_file()
    seasonal_png_exists = SEASONAL_PNG.is_file()
    if seasonal_json_exists != seasonal_png_exists:
        raise FileNotFoundError(
            "corrected seasonal aggregate is only partially present; require both JSON and PNG or neither"
        )
    seasonal_report = read_json(SEASONAL_JSON) if seasonal_json_exists else None
    if seasonal_report is not None:
        require_png(SEASONAL_PNG)
        validate_seasonal_report(seasonal_report)

    notebook["cells"].extend(
        make_notebook_cells(mask_report, n30_report, n300_report, seasonal_report)
    )
    dashboard_metadata = notebook.setdefault("metadata", {}).setdefault("awd_dashboard", {})
    dashboard_metadata["version"] = "v44"
    dashboard_metadata["last_status_sync"] = PROMOTION_DATE
    notebook_text = json.dumps(notebook, indent=1, ensure_ascii=False) + "\n"

    tex_original = TEX.read_text()
    tex_updated = update_tex(
        tex_original,
        tex_section(mask_report, n30_report, n300_report, seasonal_report),
        seasonal_report,
    )

    # All validation and text transformations have succeeded before these writes.
    atomic_write(NOTEBOOK, notebook_text)
    atomic_write(TEX, tex_updated)
    atomic_write(TEX_VERSIONED, tex_updated)
    print(f"promoted {NOTEBOOK.name} and {TEX.name} to v44")
    if seasonal_report is None:
        print(
            f"corrected seasonal aggregate pending: array {SEASONAL_ARRAY_JOB}; "
            f"postprocess {SEASONAL_POSTPROCESS_JOB}"
        )
    else:
        print(f"included corrected seasonal figure: {SEASONAL_PNG}")


if __name__ == "__main__":
    main()
