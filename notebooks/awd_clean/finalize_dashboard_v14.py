"""Finalize v14 wording with the completed target-scan result."""
from pathlib import Path
import json

P = Path(__file__).resolve().parent / "AWD_results_dashboard.ipynb"
d = json.loads(P.read_text())

def src(c):
    return "".join(c.get("source", []))

for c in d["cells"]:
    s = src(c)
    if s.startswith("## Advisor-ready core results — current dashboard v13"):
        c["source"] = s.replace("current dashboard v13", "current dashboard v14").splitlines(True)
    if s.startswith("## Current status registry"):
        old = "| Deep target-depth observability | **Targeted test added; depth conditional** | 3000–3450 m provisional scan; explicit SDZ/CDZ windows; 15–30 Hz target nulls and 3–15 Hz controls |"
        new = "| Deep target-depth observability | **Accepted conditionally; depth mapping unverified** | At both assumed SDZ/CDZ depths, positive 15–30 Hz branch validates on outbound and return legs (all four p=0.002); 3–15 Hz control is weak |"
        c["source"] = s.replace(old, new).splitlines(True)
    if s.startswith("### Targeted Deep observability at provisional SDZ and CDZ depths"):
        s = s.replace(
            "**Answer.** This branch answers whether the existing observable is present in the target windows; it does not claim that the experiment measured creep. The 15–30 Hz band is the primary tube-wave candidate test, while 3–15 Hz is a control. The numerical conclusion is taken from the target table and null report, not from visual continuity alone.",
            "**Answer.** Under the provisional mapping, the positive 15–30 Hz branch validates at both assumed strand depths on both legs: outbound and return tests at 3192 m and 3302 m each give permutation p=0.002. The 3–15 Hz control is generally weak, and isolated negative-direction p-values do not form a coherent opposing branch. This establishes target-window observability of the slow mode under the mapping hypothesis; it does not claim that the experiment measured creep or a casing reflection."
        )
        c["source"] = s.splitlines(True)
    if s.startswith("**Figure 21. Targeted Deep tube-wave observability"):
        c["source"] = s.replace(
            "Accordingly, the figure establishes target-window observability under a stated registration hypothesis, not a measured tube-wave reflection, casing deformation, permeability anomaly, or creep rate.",
            "The positive 15–30 Hz branch is validated at both provisional strand depths on both legs (each explicit channel-order test has p=0.002), whereas the 3–15 Hz control is weak and the isolated negative-direction results are not a coherent opposing branch. Accordingly, the figure establishes target-window observability under a stated registration hypothesis, not a measured tube-wave reflection, casing deformation, permeability anomaly, or creep rate."
        ).splitlines(True)

P.write_text(json.dumps(d, indent=1) + "\n")
print("Finalized dashboard v14 wording")
