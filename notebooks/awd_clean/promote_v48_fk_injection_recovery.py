#!/usr/bin/env python3
"""Promote the audited ambient F-K injection-recovery result to v48."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "AWD_results_dashboard.ipynb"
TEX = ROOT / "AWD_advisor_figure_guide.tex"
TEX_V48 = ROOT / "AWD_advisor_figure_guide_v48.tex"
OUT = ROOT / "ambient_transfer" / "fk_injection_recovery_v1_n300"
AGGREGATE = OUT / "ambient_fk_injection_recovery_v1_aggregate.json"
FIGURE = OUT / "ambient_fk_injection_recovery_v1_aggregate.png"
EXPECTED_AMPLITUDES = [0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]
EXPECTED_SCENARIOS = {
    (velocity, direction)
    for velocity in (1800.0, 2750.0, 3200.0, 4000.0)
    for direction in (1, -1)
}


def lines(text: str) -> list[str]:
    return text.splitlines(keepends=True)


def format_ratio(value: float | None) -> str:
    return "not reached" if value is None else f"{value:g}"


def validate(summary: dict) -> None:
    if summary.get("workflow_version") != "ambient_fk_injection_recovery_v1_aggregate":
        raise ValueError("unexpected aggregate workflow version")
    if summary.get("used_files") != 300:
        raise ValueError(f"expected 300 files, found {summary.get('used_files')}")
    scenarios = summary.get("scenarios", [])
    keys = {(float(item["velocity_m_s"]), int(item["direction"])) for item in scenarios}
    if keys != EXPECTED_SCENARIOS or len(scenarios) != len(EXPECTED_SCENARIOS):
        raise ValueError(f"scenario grid is incomplete: {keys}")
    for scenario in scenarios:
        amplitudes = [float(item["amplitude_ratio"]) for item in scenario["amplitudes"]]
        if amplitudes != EXPECTED_AMPLITUDES:
            raise ValueError(f"amplitude grid changed for {scenario['velocity_m_s']}")
        zero = scenario["amplitudes"][0]
        if zero["independent_recovery"] or zero["postfilter_recovery"]:
            raise ValueError("zero-injection baseline was classified as recovered")
        for point in scenario["amplitudes"]:
            if len(point["prefilter_uplift_ci95"]) != 2:
                raise ValueError("missing pre-filter confidence interval")
            if len(point["postfilter_score_uplift_ci95"]) != 2:
                raise ValueError("missing post-filter confidence interval")


def scenario_rows(summary: dict) -> str:
    rows = []
    for scenario in sorted(
        summary["scenarios"], key=lambda item: (item["velocity_m_s"], -item["direction"])
    ):
        direction = "increasing" if scenario["direction"] == 1 else "decreasing"
        rows.append(
            f"{scenario['velocity_m_s']/1000:g} km/s, {direction}: "
            f"independent {format_ratio(scenario['minimum_independent_recovery_ratio'])}; "
            f"post-filter {format_ratio(scenario['minimum_postfilter_recovery_ratio'])}"
        )
    return "\n".join(f"- {row}" for row in rows)


def latex_rows(summary: dict) -> str:
    rows = []
    for scenario in sorted(
        summary["scenarios"], key=lambda item: (item["velocity_m_s"], -item["direction"])
    ):
        direction = "increasing" if scenario["direction"] == 1 else "decreasing"
        rows.append(
            f"{scenario['velocity_m_s']/1000:g} & {direction} & "
            f"{format_ratio(scenario['minimum_independent_recovery_ratio'])} & "
            f"{format_ratio(scenario['minimum_postfilter_recovery_ratio'])} \\\\"
        )
    return "\n".join(rows)


def main() -> None:
    for path in (NOTEBOOK, TEX, AGGREGATE, FIGURE):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    summary = json.loads(AGGREGATE.read_text())
    validate(summary)
    notebook = json.loads(NOTEBOOK.read_text())
    status = notebook.setdefault("metadata", {}).setdefault("awd_dashboard", {})
    if status.get("version") == "v48":
        raise RuntimeError("notebook is already v48")
    if status.get("version") != "v47":
        raise ValueError(f"expected v47 notebook, found {status.get('version')}")

    rows = scenario_rows(summary)
    in_wedge = [
        item for item in summary["scenarios"] if item["inside_production_wedge"]
    ]
    independent = [item["minimum_independent_recovery_ratio"] for item in in_wedge]
    postfilter = [item["minimum_postfilter_recovery_ratio"] for item in in_wedge]
    independent_text = ", ".join(format_ratio(value) for value in independent)
    postfilter_text = ", ".join(format_ratio(value) for value in postfilter)
    off_target = [
        item for item in summary["scenarios"] if not item["inside_production_wedge"]
    ]
    off_target_pass = any(
        item["minimum_independent_recovery_ratio"] is not None for item in off_target
    )
    specificity = (
        "The 1.8 km/s control entered the target-energy decision and therefore reveals "
        "limited velocity specificity under this statistic."
        if off_target_pass
        else "The 1.8 km/s off-wedge control never passed the target-energy decision, "
        "supporting velocity specificity over the tested amplitude ladder."
    )

    summary_markdown = f"""## Real-noise plane-wave injection–recovery

This calibration injects deterministic broadband 5–20 Hz plane waves into the
original one-minute samples before detrending, bandpass filtering, 5 s temporal
normalization, decimation, F–K selection, or correlation. Each amplitude is
the injected RMS divided by the median real-channel 5–20 Hz RMS for that file.
The same 300 files are paired at every amplitude.

The table below reports the smallest tested positive amplitude ratio satisfying
each frozen decision. “Independent” requires the pre-filter target-energy
statistic to exceed the stricter original null threshold and have a positive
paired-bootstrap lower bound. “Post-filter” additionally requires a positive
signed-correlation uplift and a conditional velocity peak within 0.25 km/s of
the injected in-wedge velocity.

{rows}

Across the six in-wedge direction–velocity scenarios, independent thresholds
are {independent_text}; post-filter thresholds are {postfilter_text}.
{specificity}

**Interpretive boundary.** This is a sensitivity calibration in real noise,
not evidence that the uninjected ridge is a Green's function. If post-filter
recovery occurs below independent pre-filter recovery, that difference
quantifies the leverage of the frozen F–K operator; it does not by itself show
that the operator either recovered or manufactured the uninjected feature.
"""
    caption = f"""**Publication-style caption. Ambient F–K injection–recovery in
real SAFOD noise.** Synthetic random broadband plane waves spanning 5–20 Hz
were added to 300 original one-minute Nano records before all production
preprocessing. Injected RMS was scaled independently in each record to the
median real-channel 5–20 Hz RMS, with ratios 0, 0.003, 0.01, 0.03, 0.1, 0.3,
and 1.0. Colors distinguish apparent velocities of 1.8, 2.75, 3.2, and
4.0 km/s; solid and dashed curves distinguish propagation toward increasing
and decreasing fiber coordinate. Panel (a) gives the correct-branch mean
pre-filter log target/reference power statistic. Panel (b) subtracts the
identical zero-injection files and shows paired 95% bootstrap intervals.
Panel (c) gives the corresponding uplift in the physical-lag correlation score
after the frozen signed 2.5–4.5 km/s F–K wedge. Panel (d) shows the conditional
peak velocity after selection; horizontal guides mark the injected velocities.
The independent decision uses the stricter branch-specific 95th percentile
from the previously frozen channel-permutation and circular-time-shift nulls.
The off-wedge 1.8 km/s case tests velocity specificity. {specificity} Because
the production wedge restricts output slowness by construction, post-filter
recovery is interpreted only relative to the independent pre-filter threshold;
neither result proves Green's-function convergence or identifies the
uninjected ambient feature as a physical mode.
"""
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
product = dashboard_root / "ambient_transfer" / "fk_injection_recovery_v1_n300"
aggregate_json = product / "ambient_fk_injection_recovery_v1_aggregate.json"
figure = product / "ambient_fk_injection_recovery_v1_aggregate.png"
for required in (aggregate_json, figure):
    if not required.is_file():
        raise FileNotFoundError(required)
print(json.dumps(json.loads(aggregate_json.read_text()), indent=2))
display(Image(filename=str(figure), width=1200))
"""
    notebook["cells"].extend([
        {"cell_type": "markdown", "metadata": {}, "source": ["# v48\n"]},
        {"cell_type": "markdown", "metadata": {}, "source": lines(summary_markdown)},
        {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": lines(code)},
        {"cell_type": "markdown", "metadata": {}, "source": lines(caption)},
    ])
    status["version"] = "v48"
    status["last_status_sync"] = "2026-08-11"
    NOTEBOOK.write_text(json.dumps(notebook, indent=1) + "\n")

    tex = TEX.read_text()
    if tex.count("version v47") != 1:
        raise ValueError("authoritative TeX does not have exactly one v47 marker")
    tex = tex.replace("version v47", "version v48", 1)
    marker = r"\end{document}"
    if tex.count(marker) != 1:
        raise ValueError("could not isolate TeX document end")
    section = rf"""
\section{{Ambient F--K injection--recovery calibration (v48)}}
\figquestion{{At what input amplitude do the independent pre-filter statistic
and the production F--K correlation workflow recover known plane waves in the
actual SAFOD noise field?}}

Synthetic random broadband 5--20~Hz plane waves were added to the original
samples before detrending, bandpass filtering, 5-s running-absolute-mean
normalization, decimation, F--K selection, or correlation. The amplitude is
injected RMS divided by median real-channel 5--20~Hz RMS within each file. All
amplitudes use the same 300 one-minute files, so confidence intervals are
paired against an exactly matched zero-injection baseline.

\begin{{table}}[H]
\centering
\begin{{tabular}}{{llll}}
\toprule
Injected speed (km/s) & Coordinate direction & Independent ratio & Post-filter ratio \\
\midrule
{latex_rows(summary)}
\bottomrule
\end{{tabular}}
\caption{{Smallest tested positive injection ratio satisfying each frozen
decision; ``not reached'' means no amplitude through unity passed. Independent
recovery requires both the stricter original pre-filter null threshold and a
positive paired-bootstrap lower bound. Post-filter recovery additionally
requires positive score uplift and a conditional peak within
0.25~km~s\( ^{{-1}} \) of an in-wedge injection.}}
\end{{table}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.99\linewidth]{{ambient_transfer/fk_injection_recovery_v1_n300/ambient_fk_injection_recovery_v1_aggregate.png}}
\caption{{\textbf{{Ambient F--K injection--recovery in real SAFOD noise.}} Synthetic random broadband plane waves spanning 5--20~Hz were added to 300 original one-minute Nano records before all production preprocessing. Injected RMS was scaled in each record to median real-channel 5--20~Hz RMS at ratios 0, 0.003, 0.01, 0.03, 0.1, 0.3, and 1.0. Colors distinguish 1.8, 2.75, 3.2, and 4.0~km~s\( ^{{-1}} \); line style distinguishes increasing- and decreasing-coordinate propagation. (a) Correct-branch mean pre-filter log target/reference power. (b) Paired pre-filter uplift with 95\% bootstrap intervals. (c) Physical-lag score uplift after the frozen signed 2.5--4.5~km~s\( ^{{-1}} \) F--K wedge. (d) Conditional peak velocity, with horizontal injected-velocity guides. The independent decision uses the stricter branch-specific 95th percentile from the frozen channel-permutation and circular-time-shift nulls; 1.8~km~s\( ^{{-1}} \) is an off-wedge specificity control. {specificity} Post-filter recovery is meaningful only relative to the independent threshold because the selected corridor constrains output slowness. This calibration does not prove Green's-function convergence or identify the uninjected feature as a physical mode.}}
\end{{figure}}

\takeaway{{The injection experiment quantifies the gap between independent
input-wavefield observability and conditional visibility after the frozen F--K
operator. It calibrates processing sensitivity, not the physical identity of
the uninjected ambient ridge.}}

"""
    tex = tex.replace(marker, section + marker)
    TEX.write_text(tex)
    TEX_V48.write_text(tex)
    print(json.dumps({
        "version": "v48",
        "used_files": summary["used_files"],
        "in_wedge_independent_thresholds": independent,
        "in_wedge_postfilter_thresholds": postfilter,
        "off_target_passed": off_target_pass,
    }, indent=2))


if __name__ == "__main__":
    main()
