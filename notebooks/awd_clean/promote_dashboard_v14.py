"""Promote the authoritative AWD dashboard with the targeted Deep scan."""
from pathlib import Path
import json

P = Path(__file__).resolve().parent / "AWD_results_dashboard.ipynb"
d = json.loads(P.read_text())

def src(cell):
    return "".join(cell.get("source", []))

def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(True)}

def code(text):
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": text.splitlines(True)}

# Increment the single visible version marker exactly once.
for cell in d["cells"]:
    if src(cell).strip() == "v13":
        cell["source"] = ["v14\n"]
        break

# Update the opening narrative and the status registry without changing prior
# numerical claims.  The new branch is explicitly conditional on registration.
for cell in d["cells"]:
    s = src(cell)
    if s.startswith("# SAFOD June 2026 AWD-DAS Results Dashboard"):
        text = s.replace(
            "A separate 1.4–1.56 km s⁻¹ low-frequency Deep mode is accepted as a repeatable signed-slowness observation after independent burst-split and spatial-permutation tests. Its physical interpretation as a tube wave remains preliminary and is not a permeability detection.",
            "A separate 1.4–1.56 km s⁻¹ low-frequency Deep mode is accepted as a repeatable signed-slowness observation after independent burst-split and spatial-permutation tests. A new targeted scan tests whether that observable persists in windows covering the provisional SDZ/CDZ depths; the depth mapping remains conditional and its physical interpretation as a tube wave remains preliminary."
        )
        text = text.replace(
            "narrative I am trying to follow:",
            "working project direction:\n\n> Establish a downhole active-source baseline for testing whether repeatable Deep tube-wave observables are present at SAFOD's creeping strands. Quantify the AWD detection floor and the registration assumptions required before making a future deformation-monitoring claim.\n\nThe original Nano mode-identification narrative remains an active branch; the Deep target scan does not relabel the Nano mode or claim that this 24-hour experiment detected creep.\n\nnarrative I am trying to follow:"
        )
        cell["source"] = text.splitlines(True)
    elif s.startswith("## Current status registry"):
        cell["source"] = s.replace(
            "*Status flags last synchronized: 2 August 2026; dashboard promotion v13; repeatability, Deep registration/forward products, source-history regression, and physically phased CCA synchronized.*",
            "*Status flags last synchronized: 2 August 2026; dashboard promotion v14; targeted Deep strand-depth scan added with provisional registration explicitly flagged.*"
        )
        marker = "| Deep tube-wave interpretation | **Preliminary** | Consistent slow directed mode; physical tube-wave label remains unproven |"
        addition = marker + "\n| Deep target-depth observability | **Targeted test added; depth conditional** | 3000–3450 m provisional scan; explicit SDZ/CDZ windows; 15–30 Hz target nulls and 3–15 Hz controls |"
        cell["source"] = cell["source"].replace(marker, addition)
    elif s.startswith("## 13. Deep fiber-coordinate registration"):
        cell["source"] = s.replace(
            "### Limitation and next step",
            "### Target-depth observability branch — conditional baseline test\n\n**Question.** If the plausible surface-lead-in interpretation is used provisionally, is the validated Deep mode observable in windows that cover the SDZ (3192 m) and CDZ (3302 m)?\n\n**Method.** `deep_target_scan.py` scans 200 m windows at 50 m target spacing across provisional 3000–3450 m depth. Outbound coordinate is treated as depth; the reversed return coordinate is mapped as `turnaround-depth`, with turnaround 3475.31 m. Odd/even non-empty epochs provide discovery and validation. Both signed slowness directions and 3–15 and 15–30 Hz bands are tested. At the two strand depths, 499 channel-order permutations provide target-specific validation nulls.\n\n**Answer.** The targeted branch is an observability test, not a creep detection. The 15–30 Hz branch should be interpreted only through the generated target table and null report; 3–15 Hz is retained as a lower-frequency control. A positive result would justify a depth-registered baseline study; a negative result would limit the proposed creep-monitoring observable.\n\n**Figure.** `deep_target_scan.png` shows target-depth validation semblance, strand windows, and permutation nulls. Every depth label is explicitly provisional.\n\n### Limitation and next step"
        )
    elif s.startswith("## 19. Current conclusion"):
        cell["source"] = s.replace(
            "No geological `V_P` inversion or permeability inference is claimed.",
            "No geological `V_P` inversion or permeability inference is claimed. The project direction is now a conditional downhole baseline study: the targeted Deep strand-depth scan tests whether the repeatable slow observable is present where the SDZ/CDZ would lie if the surface-lead-in mapping is correct. The present 24-hour data establish observability and sensitivity requirements, not a creep rate; repeated surveys would be required for deformation monitoring."
        )

# Insert the target-scan figure and caption after the existing Deep registration
# figure caption, but do not create a second notebook.
if not any("Targeted Deep observability at provisional SDZ and CDZ depths" in src(c) for c in d["cells"]):
    idx = next(i for i,c in enumerate(d["cells"]) if src(c).startswith("**Figure 17. Deep fiber-coordinate registration")) + 1
    cells = [
        md("""### Targeted Deep observability at provisional SDZ and CDZ depths — conditional baseline result\n\n**Question.** Does the repeatable Deep slow mode remain observable in the provisional measured-depth interval containing the Southwest and Central Deforming Zones?\n\n**Method.** The scan uses 200 m local windows sampled every 50 m across provisional 3000–3450 m. Outbound and reversed-return legs are mapped separately under the explicit conditional transform `s_out = depth` and `s_return = 3475.31 m − depth`. Odd non-empty epochs select signed slowness and arrival time; even epochs validate the fixed trajectory. The two strand depths receive 499 channel-order permutations in each direction and frequency band.\n\n**Answer.** This branch answers whether the existing observable is present in the target windows; it does not claim that the experiment measured creep. The 15–30 Hz band is the primary tube-wave candidate test, while 3–15 Hz is a control. The numerical conclusion is taken from the target table and null report, not from visual continuity alone.\n\n**Limitation.** The depth transform is unverified, the current record is one 24-hour survey, and the statistic is a phase-coherence measure rather than a calibrated casing-deformation amplitude. A future repeated survey and deployment/depth documentation are required before making a deformation-monitoring claim.\n"""),
        code("""show_image('deep_target_scan.png')\n+print((HERE/'deep_target_scan.txt').read_text())\n+"""),
        md("""**Figure 21. Targeted Deep tube-wave observability at the provisional SAFOD deforming-strand depths.** The upper panel plots even-epoch validation semblance for 200-m Deep windows sampled every 50 m across the provisional 3000–3450 m interval. Solid curves are the outbound leg and dashed curves are the reversed-return leg; blue and red distinguish positive and negative signed slowness. The orange band marks the independently reported 3150–3414 m damage-zone interval, while vertical markers identify the provisional SDZ (3192 m) and CDZ (3302 m). The middle panels show display-scaled validation record sections for windows centered near those two depths; scaling is applied only for visualization and does not enter the coherence statistics. The lower panel overlays the 499 channel-order permutation distributions for the explicit strand-depth tests; red lines mark the observed fixed-trajectory validation statistics and labels report empirical tail probabilities. The mapping `s_out = depth` and `s_return = 3475.31 m − depth` assumes that the interrogator lead-in is surface cable and that the hairpin coordinate is the turnaround depth; neither assumption is independently surveyed in the present files. Accordingly, the figure establishes target-window observability under a stated registration hypothesis, not a measured tube-wave reflection, casing deformation, permeability anomaly, or creep rate.\n""")
    ]
    d["cells"][idx:idx] = cells

P.write_text(json.dumps(d, indent=1) + "\n")
print(f"Promoted {P} to v14 with {len(d['cells'])} cells")
