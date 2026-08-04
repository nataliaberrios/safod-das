"""Audit Deep tube-validation products and stage an unreviewed dashboard candidate."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import csv
import json
import os
import re
import subprocess

import numpy as np


HERE = Path(__file__).resolve().parent
PARENT_JOB = os.environ.get("AWD_PARENT_JOB", "37217153")
MAIN_NOTEBOOK = HERE / "AWD_results_dashboard.ipynb"
PROPOSAL_FILE = HERE / "PROPOSED_NOTEBOOK_UPDATE.md"
STATUS_FILE = HERE / "NIGHTLY_STATUS.md"
AUDIT_FILE = HERE / "deep_tube_validation_audit.json"

EXPECTED = {
    "record_sections": HERE / "deep_tube_record_sections.png",
    "repeatability_figure": HERE / "deep_tube_repeatability.png",
    "null_figure": HERE / "deep_tube_null_tests.png",
    "candidate_table": HERE / "deep_tube_candidates.csv",
    "burst_table": HERE / "deep_tube_burst_repeatability.csv",
    "numerical_product": HERE / "deep_tube_validation.npz",
    "text_report": HERE / "deep_tube_validation.txt",
}


def job_state(job_id: str) -> str:
    try:
        result = subprocess.run(
            ["sacct", "-X", "-j", job_id, "--noheader", "--format=State,ExitCode"],
            check=False, capture_output=True, text=True,
        )
        return " ".join(result.stdout.split()) or "unknown"
    except Exception as error:
        return f"unavailable ({error})"


def parse_report(path: Path):
    rows = []
    pattern = re.compile(
        r"^(outbound|return) (\d+)-(\d+) Hz: validation fixed-p power=([0-9.eE+-]+), "
        r"permutation p=([0-9.eE+-]+), bursts \+p>-p=([0-9.eE+-]+)$"
    )
    for line in path.read_text().splitlines():
        match = pattern.match(line.strip())
        if match:
            rows.append({
                "leg": match.group(1), "band_low_hz": float(match.group(2)),
                "band_high_hz": float(match.group(3)),
                "validation_power": float(match.group(4)),
                "permutation_p": float(match.group(5)),
                "direction_fraction": float(match.group(6)),
            })
    return rows


def markdown_cell(text: str):
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(True)}


def code_cell(text: str):
    return {
        "cell_type": "code", "execution_count": None, "metadata": {},
        "outputs": [], "source": text.splitlines(True),
    }


def result_sentence(metrics):
    pieces = []
    for row in metrics:
        pieces.append(
            f"{row['leg']} {row['band_low_hz']:.0f}–{row['band_high_hz']:.0f} Hz: "
            f"permutation p={row['permutation_p']:.3f}, "
            f"{100*row['direction_fraction']:.0f}% of bursts have +p power > −p power"
        )
    return "; ".join(pieces)


def main():
    timestamp = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
    files = {
        name: {
            "path": str(path), "exists": path.exists(),
            "bytes": path.stat().st_size if path.exists() else 0,
        }
        for name, path in EXPECTED.items()
    }
    problems = [name for name, item in files.items() if not item["exists"] or item["bytes"] == 0]
    metrics = []
    candidate_rows = burst_rows = 0
    numerical_valid = False
    if not problems:
        try:
            metrics = parse_report(EXPECTED["text_report"])
            with EXPECTED["candidate_table"].open() as stream:
                candidate_rows = sum(1 for _ in csv.DictReader(stream))
            with EXPECTED["burst_table"].open() as stream:
                burst_rows = sum(1 for _ in csv.DictReader(stream))
            with np.load(EXPECTED["numerical_product"]) as product:
                numerical_valid = all(np.size(product[key]) > 0 for key in product.files)
        except Exception as error:
            problems.append(f"parse_error: {error}")
    if len(metrics) != 4:
        problems.append(f"expected 4 report metrics, found {len(metrics)}")
    if candidate_rows != 24:
        problems.append(f"expected 24 candidate rows, found {candidate_rows}")
    if burst_rows != 184:
        problems.append(f"expected 184 burst rows, found {burst_rows}")
    if not numerical_valid:
        problems.append("NPZ integrity check failed")

    computational_status = "PASS" if not problems else "FAIL"
    screens = []
    for row in metrics:
        if row["permutation_p"] <= 0.05 and row["direction_fraction"] >= 0.70:
            screen = "PASS"
        elif row["permutation_p"] <= 0.10 or row["direction_fraction"] >= 0.60:
            screen = "REVIEW"
        else:
            screen = "FAIL"
        screens.append({**row, "screen": screen})
    audit = {
        "timestamp": timestamp, "parent_job": PARENT_JOB,
        "parent_job_state": job_state(PARENT_JOB),
        "computational_status": computational_status,
        "problems": problems, "files": files,
        "candidate_rows": candidate_rows, "burst_rows": burst_rows,
        "metrics": screens,
        "interpretation_policy": "Automated screens are not scientific phase labels.",
    }
    AUDIT_FILE.write_text(json.dumps(audit, indent=2) + "\n")

    lines = [
        "# SAFOD AWD overnight status", "",
        f"- Generated: `{timestamp}`",
        f"- Parent validation job: `{PARENT_JOB}` ({audit['parent_job_state']})",
        f"- Computational audit: **{computational_status}**",
        f"- Candidate table rows: `{candidate_rows}` (expected 24)",
        f"- Burst table rows: `{burst_rows}` (expected 184)", "",
        "## Statistical screens", "",
        "| Leg | Band | Fixed-p validation power | Permutation p | Bursts +p > −p | Automated screen |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in screens:
        lines.append(
            f"| {row['leg']} | {row['band_low_hz']:.0f}–{row['band_high_hz']:.0f} Hz | "
            f"{row['validation_power']:.5f} | {row['permutation_p']:.4f} | "
            f"{100*row['direction_fraction']:.1f}% | **{row['screen']}** |"
        )
    lines += ["", "## Integrity problems", ""]
    lines += [f"- {problem}" for problem in problems] if problems else ["- None detected."]
    lines += [
        "", "## Interpretation boundary", "",
        "The automated screen tests computational integrity, split-sample spatial repeatability, fixed-p channel ordering, and burst directionality. It does **not** establish a tube-wave phase, locate a permeable fracture, or promote any notebook status flag.",
        "", "## Morning review checklist", "",
        "1. Inspect all three figures for artifacts and branch continuity.",
        "2. Verify candidate trajectories against the time-domain wavefield.",
        "3. Review the provisional hairpin orientation and depth registration.",
        "4. Decide whether the result remains preliminary, becomes supported, or is rejected.",
        "5. If scientifically accepted, update the single authoritative notebook in place, increment its version once, and synchronize every status flag.",
    ]
    if computational_status == "PASS":
        proposal = [
            "# Proposed update to the single AWD dashboard", "",
            "**Automated report only. Do not treat this as a promoted scientific interpretation.**", "",
            "Authoritative notebook: `AWD_results_dashboard.ipynb`", "",
            "## Machine-extracted result", "", result_sentence(metrics), "",
            "## Required review before updating the notebook", "",
            "1. Inspect all generated figures.",
            "2. Assign accepted, preliminary, inconclusive, or rejected status.",
            "3. Update the authoritative notebook in place.",
            "4. Increment its internal version cell by exactly one.",
            "5. Synchronize the status registry, captions, evidence table, and conclusion.",
            "6. Clear or regenerate affected rendered outputs.",
        ]
        PROPOSAL_FILE.write_text("\n".join(proposal) + "\n")
        lines += ["", f"Proposed update report: `{PROPOSAL_FILE.name}`. No notebook was created or modified automatically."]
    else:
        lines += ["", "No notebook update was proposed because the computational audit failed."]
    STATUS_FILE.write_text("\n".join(lines) + "\n")
    print(f"Computational audit: {computational_status}")
    print(f"Wrote {STATUS_FILE}, {AUDIT_FILE}")
    if computational_status == "PASS":
        print(f"Wrote {PROPOSAL_FILE}; authoritative notebook unchanged")


if __name__ == "__main__":
    main()
