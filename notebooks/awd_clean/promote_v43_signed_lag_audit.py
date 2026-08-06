#!/usr/bin/env python3
"""Append the signed-lag audit status to the authoritative dashboard."""
from __future__ import annotations

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
NOTEBOOK = HERE / "AWD_results_dashboard.ipynb"
MARKER = "# v43"


def source_text(value):
    if isinstance(value, str):
        return value
    if isinstance(value, list):
        return "".join(source_text(item) for item in value)
    return str(value)


def markdown(source):
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": (source.rstrip() + "\n").splitlines(True),
    }


def code(source):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": (source.rstrip() + "\n").splitlines(True),
    }


def main():
    notebook = json.loads(NOTEBOOK.read_text())
    if any(MARKER in source_text(cell.get("source", [])) for cell in notebook["cells"]):
        raise SystemExit("v43 marker already present")
    notebook["cells"].extend([
        markdown(
            """# v43

## Critical signed-lag audit: physical positive branch pending corrected seasonal rerun

The legacy FFT correlation assembled the negative-lag half from a truncated array. The zero/positive-lag half was correct, so the reported negative-mask target-moveout scores and their receiver-permutation nulls remain reproducible. The negative-lag half was not valid. Consequently, the legacy positive-mask value `p = 0.978` is a same-positive-lag leakage control; it must not be described as evidence that the physical opposite-propagating branch is absent.

The corrected convention has now passed three deterministic regressions: a delayed impulse returns `+0.20 s`, the reversed pair returns `-0.20 s`, and broadband waves propagating toward increasing/decreasing fiber coordinate are retained by opposite signed F–K masks and peak at opposite correlation lags. A one-file real-data smoke test gives comparable physical branches near 3.1 km/s on opposite lag signs. This is diagnostic only; the versioned eight-day recomputation is running as SLURM array `37717572`, with dependent aggregation job `37717613`.

Until that aggregation completes, Figure 39 supports only the repeatability of the negative-mask positive-lag observable. Figure 34's positive-mask score is an opposite-mask leakage control, not a physical anti-causal branch test."""
        ),
        code(
            """from pathlib import Path
import json
from IPython.display import Image, display

root = Path.cwd()
v2 = (root / "ambient_transfer" / "signed_lag_v2"
      if (root / "ambient_transfer").exists()
      else root / "awd_clean" / "ambient_transfer" / "signed_lag_v2")
aggregate = v2 / "seasonal_signed_fk_v2_aggregate.json"
figure = v2 / "seasonal_signed_fk_v2_aggregate.png"
if aggregate.exists():
    print(json.dumps(json.loads(aggregate.read_text()), indent=2))
    if figure.exists():
        display(Image(filename=str(figure), width=1200))
else:
    print("Corrected eight-day signed-lag aggregation pending: array 37717572; postprocess 37717613.")"""
        ),
        markdown(
            """**Status caption.** The corrected seasonal figure will compare the negative F–K mask at positive lag with the positive F–K mask at negative lag, using the same channel-0 virtual source, common display scale, signed 3.2 km/s references, branch-specific velocity scores, and receiver-permutation nulls. No conclusion about causal/anti-causal asymmetry should be drawn until this cell finds the completed v2 aggregate."""
        ),
    ])
    notebook.setdefault("metadata", {}).setdefault("awd_dashboard", {})["version"] = "v43"
    notebook["metadata"]["awd_dashboard"]["last_status_sync"] = "2026-08-05"
    NOTEBOOK.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n")
    print("promoted", NOTEBOOK, "to v43")


if __name__ == "__main__":
    main()
