#!/usr/bin/env python3
"""Promote the independently audited signed-lag aggregate to dashboard v46."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "AWD_results_dashboard.ipynb"
TEX = ROOT / "AWD_advisor_figure_guide.tex"
TEX_V46 = ROOT / "AWD_advisor_figure_guide_v46.tex"
OUT = ROOT / "ambient_transfer" / "signed_lag_v2"
AGGREGATE = OUT / "seasonal_signed_fk_v2_aggregate.json"
AUDIT = OUT / "seasonal_signed_fk_v2_completion_audit.json"
FIGURE = OUT / "seasonal_signed_fk_v2_aggregate.png"


def lines(text: str) -> list[str]:
    return text.splitlines(keepends=True)


def conditional_result(p_value: float) -> str:
    return (
        "rejects the conditional receiver-order null"
        if p_value <= 0.05
        else "does not reject the conditional receiver-order null"
    )


def main() -> None:
    for path in (NOTEBOOK, TEX, AGGREGATE, AUDIT, FIGURE):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    aggregate = json.loads(AGGREGATE.read_text())
    audit = json.loads(AUDIT.read_text())
    required = {
        "frozen_manifest_complete",
        "branch_lag_signs_correct",
        "saved_metrics_reproduced_from_sections",
        "aggregate_figure_present",
    }
    if audit.get("workflow_version") != "signed_lag_v2_completion_audit":
        raise ValueError("unexpected audit workflow version")
    if not all(audit.get("checks", {}).get(key) is True for key in required):
        raise ValueError(f"audit is not complete: {audit.get('checks')}")
    if aggregate.get("branch_lag_sign") != {"negative": 1.0, "positive": -1.0}:
        raise ValueError("corrected branch-lag convention is absent")
    if aggregate.get("weighted_files") != audit.get("expected_and_used_files"):
        raise ValueError("aggregate and audit file counts disagree")

    neg = audit["branches"]["negative"]
    pos = audit["branches"]["positive"]
    ratio = float(audit["positive_to_negative_absolute_score_ratio_3200"])
    asymmetry = float(audit["signed_branch_asymmetry_3200"])
    relation = (
        "comparable in magnitude"
        if 0.67 <= ratio <= 1.5
        else ("positive-branch dominant" if ratio > 1.5 else "negative-branch dominant")
    )
    used = int(audit["expected_and_used_files"])
    ndays = len(audit["days"])

    notebook = json.loads(NOTEBOOK.read_text())
    status = notebook.setdefault("metadata", {}).setdefault("awd_dashboard", {})
    if status.get("version") == "v46":
        raise RuntimeError("notebook is already v46")
    if status.get("version") != "v45":
        raise ValueError(f"expected v45 notebook, found {status.get('version')}")
    notebook["cells"][170]["source"] = lines(
        """# v43 (completed in v46)

## Critical signed-lag audit

The v43 audit found that the legacy negative-lag extraction was invalid and
prohibited interpreting the positive-mask same-positive-lag control as the
physical opposite branch. Corrected synthetic tests and the complete
eight-day comparison are reported in v46.
"""
    )
    notebook["cells"][182]["source"] = lines(
        """### Corrected seasonal signed-lag result

The corrected eight-day aggregation and its independent completion audit are
complete. The dynamic code cell below loads the versioned product; quantitative
interpretation and the detailed caption are provided in v46.
"""
    )
    notebook["cells"][184]["source"] = lines(
        """**Status caption.** Superseded by the independently audited v46 result.
The corrected convention evaluates F×K < 0 at positive lag and F×K > 0 at
negative lag. The legacy same-positive-lag comparison remains only an
opposite-mask leakage control.
"""
    )
    summary = f"""## Corrected eight-day physical-lag branch comparison

The aggregate contains {used:,} one-minute files from {ndays} preselected
dates. F×K < 0 is evaluated at positive lag for propagation toward increasing
fiber coordinate, whereas F×K > 0 is evaluated at negative lag for the
opposite coordinate direction. These are coordinate-direction labels, not
automatic upgoing or downgoing assignments.

- **F×K < 0, positive lag:** peak {neg['peak_velocity_m_s']/1000:.3f} km/s;
  absolute score at 3.2 km/s {neg['absolute_score_3200']:.4f};
  p = {neg['p_peak']:.4g}; physical/leakage ratio
  {neg['physical_to_leakage_ratio_3200']:.2f}. It
  {conditional_result(neg['p_peak'])}.
- **F×K > 0, negative lag:** peak {pos['peak_velocity_m_s']/1000:.3f} km/s;
  absolute score at 3.2 km/s {pos['absolute_score_3200']:.4f};
  p = {pos['p_peak']:.4g}; physical/leakage ratio
  {pos['physical_to_leakage_ratio_3200']:.2f}. It
  {conditional_result(pos['p_peak'])}.
- The positive/negative absolute-score ratio at 3.2 km/s is {ratio:.3f}, and
  normalized asymmetry is {asymmetry:.3f}. The branches are **{relation}**
  after the frozen signed F–K selection.

**Interpretive boundary.** Both physical lag sides are now computed with the
correct FFT extraction and sign convention. Receiver-order permutation tests
ordered moveout in each selected correlation section, but remains conditional
on the F–K operator. The v45 pre-filter test found no target-power enrichment
and the strict unfiltered 3.2 km/s trajectories were weak. The corrected
branches are therefore repeatable filtered observables, not independently
validated Green's functions, upgoing or downgoing waves, identified modes, or
formation velocity.
"""
    caption = f"""**Publication-style caption. Corrected eight-day signed-lag
comparison after independent completion audit.** The aggregate uses {used:,}
one-minute records from {ndays} preselected dates, Nano channel 0 as the
virtual source, increasing receiver offset along the cemented fiber, 5–20 Hz
filtering, 5 s running-absolute-mean temporal normalization, and a frozen
2.5–4.5 km/s signed F–K wedge. The upper-left panel shows F×K < 0 on its
physical positive-lag side; the upper-right shows F×K > 0 on its physical
negative-lag side, using a common correlation scale. Black dashed trajectories
mark 3.2 km/s and the yellow star marks the virtual source. The lower-left
shows absolute median moveout scores and branch-specific receiver-permutation
95% thresholds. The lower-right compares physical-side and opposite-lag
leakage at 3.2 km/s. Negative and corrected positive branches peak at
{neg['peak_velocity_m_s']/1000:.3f} and {pos['peak_velocity_m_s']/1000:.3f}
km/s, with p = {neg['p_peak']:.4g} and p = {pos['p_peak']:.4g}; their
3.2-km/s absolute-score ratio is {ratio:.3f}. The audit independently
reproduces the scores, verifies the lag signs, and confirms every selected
file. Permutation significance quantifies ordered moveout after frozen F–K
selection; it does not establish Green's-function convergence, causality,
upgoing or downgoing identity, wave type, or formation velocity.
"""
    notebook["cells"].extend(
        [
            {"cell_type": "markdown", "metadata": {}, "source": ["# v46\n"]},
            {"cell_type": "markdown", "metadata": {}, "source": lines(summary)},
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": lines(
                    """from pathlib import Path
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
product = dashboard_root / "ambient_transfer" / "signed_lag_v2"
aggregate_json = product / "seasonal_signed_fk_v2_aggregate.json"
audit_json = product / "seasonal_signed_fk_v2_completion_audit.json"
figure = product / "seasonal_signed_fk_v2_aggregate.png"
for required in (aggregate_json, audit_json, figure):
    if not required.is_file():
        raise FileNotFoundError(required)
print(json.dumps(json.loads(aggregate_json.read_text()), indent=2))
print(json.dumps(json.loads(audit_json.read_text()), indent=2))
display(Image(filename=str(figure), width=1200))
"""
                ),
            },
            {"cell_type": "markdown", "metadata": {}, "source": lines(caption)},
        ]
    )
    status["version"] = "v46"
    status["last_status_sync"] = "2026-08-10"
    NOTEBOOK.write_text(json.dumps(notebook, indent=1) + "\n")

    tex = TEX.read_text()
    if tex.count("version v45") != 1:
        raise ValueError("authoritative TeX does not have one v45 marker")
    tex = tex.replace("version v45", "version v46", 1)
    start = r"\subsection{Corrected seasonal signed-lag status}"
    end = r"\end{document}"
    if tex.count(start) != 1 or tex.count(end) != 1:
        raise ValueError("could not isolate final TeX status subsection")
    section = rf"""\subsection{{Corrected eight-day signed-lag comparison (v46)}}
\figquestion{{Does corrected \(F K>0\) contain ordered moveout on its physical
negative-lag side, and how does it compare with \(F K<0\) at positive lag?}}

The independently audited aggregate contains {used:,} one-minute files from
{ndays} frozen dates. \(F K<0\) is evaluated at positive lag for propagation
toward increasing fiber coordinate; \(F K>0\) is evaluated at negative lag
for the opposite coordinate direction. These are coordinate directions, not
automatic upgoing or downgoing labels.

The negative branch peaks at {neg['peak_velocity_m_s']/1000:.3f}~km~s\( ^{{-1}} \),
has absolute score {neg['absolute_score_3200']:.4f} at 3.2~km~s\( ^{{-1}} \),
and receiver-permutation \(p={neg['p_peak']:.4g}\). The corrected positive
branch peaks at {pos['peak_velocity_m_s']/1000:.3f}~km~s\( ^{{-1}} \), has
absolute score {pos['absolute_score_3200']:.4f}, and
\(p={pos['p_peak']:.4g}\). Their positive/negative absolute-score ratio is
{ratio:.3f}, with normalized asymmetry {asymmetry:.3f}; the branches are
{relation} after signed F--K selection.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.99\linewidth]{{ambient_transfer/signed_lag_v2/seasonal_signed_fk_v2_aggregate.png}}
\caption{{\textbf{{Corrected eight-day signed-lag comparison after independent completion audit.}} The aggregate uses {used:,} one-minute records from {ndays} preselected dates, Nano channel~0 as the virtual source, increasing offset along the cemented fiber, 5--20~Hz filtering, 5-s running-absolute-mean normalization, and a frozen 2.5--4.5~km~s\( ^{{-1}} \) signed F--K wedge. The upper-left panel shows \(F K<0\) on positive lag; the upper-right shows \(F K>0\) on negative lag, with a common correlation scale. Dashed trajectories mark 3.2~km~s\( ^{{-1}} \). The lower-left shows moveout scores and receiver-permutation 95\% thresholds; the lower-right compares physical and opposite-lag leakage at 3.2~km~s\( ^{{-1}} \). Branches peak at {neg['peak_velocity_m_s']/1000:.3f} and {pos['peak_velocity_m_s']/1000:.3f}~km~s\( ^{{-1}} \), with \(p={neg['p_peak']:.4g}\) and \(p={pos['p_peak']:.4g}\); their score ratio is {ratio:.3f}. The audit reproduces scores, verifies signs, and confirms all files. Significance is conditional on F--K selection and does not independently establish Green's-function convergence, causality, propagation identity, wave type, or formation velocity.}}
\end{{figure}}

\takeaway{{Both physical lag sides have now been computed correctly and are
{relation}. Their significance remains conditional on F--K selection. The v45
pre-filter test and strict unfiltered controls do not independently corroborate
the branches as ambient Green's functions, identified modes, or formation
\(V_P\).}}

"""
    tex = tex[: tex.index(start)] + section + tex[tex.index(end) :]
    TEX.write_text(tex)
    TEX_V46.write_text(tex)
    print(json.dumps({"version": "v46", "ratio_3200": ratio, "asymmetry": asymmetry}, indent=2))


if __name__ == "__main__":
    main()
