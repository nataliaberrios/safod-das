#!/usr/bin/env python3
"""Promote completed ambient-noise validation products into v41 records."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTEBOOK = HERE / "AWD_results_dashboard.ipynb"
GUIDE_V40 = HERE / "AWD_advisor_figure_guide_v40.tex"
GUIDE_V41 = HERE / "AWD_advisor_figure_guide_v41.tex"
GUIDE_CURRENT = HERE / "AWD_advisor_figure_guide.tex"


def md(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(True)}


def code(text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": text.splitlines(True),
    }


def promote_notebook() -> None:
    notebook = json.loads(NOTEBOOK.read_text())
    metadata = notebook.setdefault("metadata", {}).setdefault("awd_dashboard", {})
    metadata["version"] = "v41"
    metadata["last_status_sync"] = "2026-08-05"

    replacements = {
        "The next validation is to hold the negative wedge fixed on additional days, run the same signed controls on those days, and test nearby velocity wedges that do not move the peak artificially. Only after that should this be promoted as a reproduction of Lellouch’s ambient result.": "The completed eight-day aggregate now holds the negative wedge fixed across independent dates and signed controls. The result remains a conditional transfer test because the wedge contains the target velocity interval; nearby-wedge and held-out-mask tests remain the appropriate next robustness checks.",
        "but complete-day and across-season products are still required.": "and the completed seasonal products are now reported in v41.",
        "This is a preliminary one-day partial stack; the complete seasonal jobs remain in progress.": "This remains an interim one-day partial stack; the completed seasonal aggregate is reported in the v41 section below.",
    }
    for cell in notebook["cells"]:
        source = "".join(cell.get("source", []))
        for old, new in replacements.items():
            source = source.replace(old, new)
        if source:
            cell["source"] = source.splitlines(True)

    notebook["cells"].extend([
        md("""# v41

## Completed seasonal F–K validation, frequency-band diagnostics, and anti-alias pilot

The seasonal and preprocessing-validation jobs have completed. This section promotes their saved products into the authoritative dashboard and separates robust conclusions from diagnostics that remain conditional.
"""),
        code("""from pathlib import Path
import json
import pandas as pd
from IPython.display import Image, display

here = Path.cwd()
out = here / "ambient_transfer" if (here / "ambient_transfer").exists() else here / "awd_clean" / "ambient_transfer"

fk = json.loads((out / "fk_seasonal_aggregate.json").read_text())
unfiltered = json.loads((out / "seasonal_unfiltered_aggregate.json").read_text())
alias = json.loads((out / "alias_sensitivity_2024-12-20_start0_n30.json").read_text())

summary = pd.DataFrame([
    {"path": "8-day negative signed F–K", **fk["modes"]["negative"]},
    {"path": "8-day positive signed F–K", **fk["modes"]["positive"]},
    {"path": "8-day two-sided F–K", **fk["modes"]["both"]},
    {"path": "8-day exact-resolution unfiltered", "peak_velocity_m_s": unfiltered["peak_v_mps"], "peak_score": unfiltered["peak_score"], "score_3200": unfiltered["score_3200"], "null95": unfiltered["null95"], "p_peak": unfiltered["p_peak"]},
])
display(summary)
"""),
        code("""display(Image(filename=str(out / "fk_seasonal_day_comparison.png"), width=1100))
"""),
        md("""**Publication-style caption.** Eight-day seasonal validation of the Lellouch-style ambient transfer observable. Each day was processed independently from complete one-minute records using the channel-0 top-source geometry, measured 1.020952 m spacing, 5–20 Hz temporal normalization, and the fixed signed F–K masks. The negative branch reaches a peak apparent velocity of 3.075 km/s with score 0.279 and a 500-permutation receiver null probability of 0.0002 across 11,523 weighted files. The positive branch is not significant (p = 0.978), while the two-sided branch remains enhanced (p = 0.0026). The exact-resolution unfiltered aggregate is weak (p = 0.343). Because the negative wedge was defined to include 2.5–4.5 km/s, this is a statistically repeatable, directionally selected transfer observable rather than an independent velocity estimate or an exact reproduction of the 2017 acquisition."""),
        md("""## Frequency-band anomaly diagnostic

The completed frequency-band scan processed four ten-file blocks per selected day in 3–8, 5–12, 8–20, and 15–30 Hz bands. The low-frequency 3–8 Hz control gives stable early and late apparent velocities near 3.15 and 3.12 km/s, respectively, with median breakpoint near 425 m. Higher-frequency row-maximum fits are less stable and can produce nonphysical late-branch values; they are therefore treated as diagnostics of bandwidth and ridge-picking sensitivity, not as independent geological constraints. The full anomaly postprocess contains 1,167 products across nine dates and consistently places the descriptive transition near 350–425 m, not uniquely at 500–520 m."""),
        code("""band = json.loads((out / "frequency_band_anomaly_test.json").read_text())
rows = []
for b in band["bands_hz"]:
    subset = [r for r in band["results"] if r["band_hz"] == b]
    rows.append({
        "band_hz": f"{b[0]}–{b[1]}",
        "n_blocks": len(subset),
        "median_early_m_s": pd.Series([r["early_velocity_m_s"] for r in subset]).median(),
        "median_late_m_s": pd.Series([r["late_velocity_m_s"] for r in subset]).median(),
        "median_breakpoint_m": pd.Series([r["best_breakpoint_m"] for r in subset]).median(),
    })
display(pd.DataFrame(rows))
"""),
        md("""**Diagnostic caption.** Frequency-band sensitivity of the conditional negative-branch moveout anomaly. Each row summarizes 32 ten-file products from the selected seasonal days. The 3–8 Hz band is the most stable control; the higher bands show broader or unstable row maxima and are retained to expose processing sensitivity. Breakpoints are descriptive fiber coordinates under the inherited same-fiber convention, not measured lithologic boundaries or a V_P inversion."""),
        md("""## Anti-alias preprocessing sensitivity

The anti-alias pilot compares the current direct factor-of-two decimation, explicit anti-aliased polyphase resampling in time and space, and full-resolution processing. All three paths use identical normalized correlations, geometry, 5–20 Hz bandpass, signed 2.5–4.5 km/s wedge, and receiver-permutation nulls. On the first 30 one-minute records from 2024-12-20, all three paths are weak and nonsignificant: direct peak 1.875 km/s (p = 1.00), anti-aliased peak 1.850 km/s (p = 0.990), and full-resolution peak 1.925 km/s (p = 0.910). This pilot does not test the long-stack seasonal signal; it shows that the short-stack control is not materially changed by the preprocessing choice. A five-hour or complete-day anti-alias comparison remains the appropriate final sensitivity check."""),
        code("""display(Image(filename=str(out / "alias_sensitivity_2024-12-20_start0_n30.png"), width=1100))
"""),
        md("""**Publication-style caption.** Anti-alias preprocessing sensitivity for a 30-minute ambient transfer pilot. Columns compare direct factor-of-two slicing, explicit anti-aliased polyphase resampling along both time and channel axes, and full-resolution processing. The upper row shows the corresponding top-source correlation sections; the lower row shows trial-velocity scores with receiver-permutation 95th-percentile thresholds. All paths are weak and nonsignificant over this short interval, and their peak velocities differ by only 75 m/s. The result is a preprocessing control, not a rejection of the long-stack F–K observable; the longer-stack comparison is still required before declaring the seasonal result insensitive to decimation."""),
        md("""## v41 conclusion and remaining qualification

The completed seasonal products establish that the negative signed F–K observable is repeatable across eight dates and is absent from the positive branch and exact-resolution unfiltered control. The result is conditional on the selected velocity wedge and should not yet be called an independent formation-velocity measurement. The completed frequency-band products show that the apparent spatial transition is broad and processing-sensitive. The anti-alias pilot is reassuring but short. The next robustness run is a longer anti-aliased-versus-direct comparison using the same five-hour or complete-day stack."""),
    ])
    NOTEBOOK.write_text(json.dumps(notebook, indent=1, ensure_ascii=True) + "\n")


def promote_guide() -> None:
    text = GUIDE_V40.read_text()
    text = text.replace("\\date{2 August 2026}", "\\date{5 August 2026}")
    text = text.replace("version v30", "version v41")
    text = text.replace("the complete seasonal jobs remain running.", "the completed seasonal aggregate is reported below.")
    text = text.replace("the complete seasonal jobs remain in progress.", "the completed seasonal aggregate is reported in the v41 section below.")
    text = text.replace(
        "The next tests are a complete-day stack and held-out-day/wedge-sensitivity validation.",
        "The complete-day and seasonal products are now available below; held-out-day and nearby-wedge tests remain the appropriate robustness checks.",
    )
    text = text.replace(
        "The next validation is to repeat the scan in independent frequency bands and in the completed seasonal day-level aggregates, then compare the registered coordinate with the Lellouch $V_P/V_S$ and mode-conversion evidence.",
        "The completed frequency-band and seasonal products now provide that comparison. They support a broad transition-like feature, but the registered coordinate remains conditional and does not uniquely identify the Lellouch anomaly or a lithologic boundary.",
    )
    text = text.replace(
        "The ongoing seasonal and frequency-band controls test whether the F--K enhancement is persistent and whether it survives changes in passband.",
        "The completed seasonal and frequency-band controls show that the F--K enhancement persists across dates, while the apparent transition is sensitive to passband and row-maximum stability.",
    )
    addition = r"""
\subsection{Completed eight-day seasonal F--K validation and unfiltered control}
\figquestion{Does the signed F--K observable persist across independent dates, and does it remain absent from the exact-resolution unfiltered control?}

The seasonal jobs processed eight selected dates independently and then formed an aggregate from 11,523 weighted one-minute files. All records were detrended, filtered from 5--20~Hz, divided by a 5-s running mean of absolute amplitude, and correlated with the channel-0 top-source geometry using the measured 1.020952-m spacing. The fixed negative signed wedge has a peak at 3.075~km~s$^{-1}$, score 0.279, and a 500-permutation receiver null probability $p=0.0002$. The positive branch peaks at 3.825~km~s$^{-1}$ but is not significant ($p=0.978$). The two-sided branch is enhanced ($p=0.0026$), while the exact-resolution unfiltered aggregate is weak ($p=0.343$).

This establishes a repeatable directionally selected transfer observable across seasons, not an independent velocity estimate: the negative wedge was defined to include 2.5--4.5~km~s$^{-1}$. Earthquake-containing intervals were retained; temporal normalization suppresses their amplitude leverage but does not remove their waveforms.

\begin{figure}[H]
\centering
\includegraphics[width=0.98\linewidth]{ambient_transfer/fk_seasonal_day_comparison.png}
\caption{\textbf{Eight-day seasonal validation of the ambient signed F--K observable.} Each selected date was processed independently with the channel-0 top-source geometry, measured 1.020952-m spacing, 5--20-Hz preprocessing, and fixed signed velocity masks before aggregation. The negative branch reaches 3.075~km~s$^{-1}$ with score 0.279 and receiver-permutation $p=0.0002$ across 11,523 weighted files. The positive branch is not significant ($p=0.978$), the two-sided branch remains enhanced ($p=0.0026$), and the exact-resolution unfiltered control is weak ($p=0.343$). The figure demonstrates seasonal repeatability of a directionally selected transfer observable; because the negative wedge contains the target velocity interval by construction, it is not an independent formation-velocity estimate or an exact reproduction of the 2017 acquisition.}\end{figure}

\subsection{Completed frequency-band and anomaly diagnostics}
The frequency-band scan processed four ten-file blocks per selected day in 3--8, 5--12, 8--20, and 15--30~Hz bands. The 3--8-Hz control is the most stable, with median early and late apparent velocities of approximately 3.15 and 3.12~km~s$^{-1}$ and a median descriptive breakpoint near 425~m. Higher-frequency row-maximum fits are less stable and sometimes yield nonphysical late-branch values, so they are retained as sensitivity diagnostics rather than independent geological constraints. The full anomaly postprocess contains 1,167 products across nine dates and places the descriptive transition broadly near 350--425~m, not uniquely at 500--520~m.

\takeaway{The seasonal result is robust as a selected directional observable. The apparent spatial transition is repeatable but broad, frequency-dependent, and conditional on fiber-coordinate registration; it is not yet a lithologic boundary or $V_P$ inversion.}

\subsection{Anti-alias preprocessing sensitivity pilot}
\figquestion{Does the short-stack F--K result depend materially on direct factor-of-two decimation?}

The pilot compares the current direct slicing path, explicit anti-aliased polyphase resampling along both time and channel axes, and full-resolution processing. The same normalized correlations, 5--20-Hz bandpass, 2.5--4.5~km~s$^{-1}$ signed wedge, and receiver-permutation nulls are used in all paths. For the first 30 one-minute records on 2024-12-20, the direct path peaks at 1.875~km~s$^{-1}$ ($p=1.00$), the anti-aliased path at 1.850~km~s$^{-1}$ ($p=0.990$), and the full-resolution path at 1.925~km~s$^{-1}$ ($p=0.910$). These are all weak short-stack controls and differ by only 75~m~s$^{-1}$.

The pilot is reassuring about gross preprocessing sensitivity but does not validate the long-stack seasonal signal: 30 minutes is much shorter than the five-hour and complete-day products. A longer anti-aliased-versus-direct comparison remains the next robustness run.

\begin{figure}[H]
\centering
\includegraphics[width=0.98\linewidth]{ambient_transfer/alias_sensitivity_2024-12-20_start0_n30.png}
\caption{\textbf{Anti-alias preprocessing sensitivity for a 30-minute ambient pilot.} Columns compare direct factor-of-two slicing, explicit anti-aliased polyphase resampling along the time and channel axes, and full-resolution processing. The upper row shows top-source correlation sections; the lower row shows trial-velocity scores with receiver-permutation 95th-percentile thresholds. All three paths are weak and nonsignificant over this short interval, with peak velocities from 1.850 to 1.925~km~s$^{-1}$. The result is a preprocessing control rather than a test of the full seasonal F--K signal; the longer-stack comparison is required before declaring the seasonal result insensitive to decimation.}\end{figure}

\takeaway{v41 promotes the completed seasonal, unfiltered, frequency-band, and anti-alias products. The defensible claim is a repeatable, conditional negative signed F--K observable near 3.075~km~s$^{-1}$, with a broad processing-sensitive spatial transition and an anti-alias pilot that is reassuring but not yet long-stack validation.}

"""
    text = text.replace("\\end{document}", addition + "\\end{document}")
    GUIDE_V41.write_text(text)
    shutil.copyfile(GUIDE_V41, GUIDE_CURRENT)


if __name__ == "__main__":
    promote_notebook()
    promote_guide()
    print(f"Promoted dashboard and guide to v41: {NOTEBOOK}, {GUIDE_V41}")

