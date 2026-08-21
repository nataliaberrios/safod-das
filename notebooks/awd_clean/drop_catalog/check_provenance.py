"""Fail if any section of AWD_drop_catalog.ipynb does not name what produced it.

Written after an audit found five of seven figure sections naming no source, and
a second pass found the virtual-source section still missing one because the
first check only looked at headings containing "Figure" or "numbers". So this
checks EVERY numbered section that has output, rather than a hand-picked list.

    python3 awd_clean/drop_catalog/check_provenance.py

Exits non-zero if a section that produces a figure or printed output has no
`**Source.**` marker in the markdown immediately preceding it.
"""
import json
import re
import sys
from pathlib import Path

NB = Path(__file__).resolve().parent / "AWD_drop_catalog.ipynb"
# Sections that legitimately have no source: they configure or explain, and
# produce nothing a reader would need to reproduce.
EXEMPT = ("1. Configuration", "2. Load the products",
          "11. Provenance", "12. What to be careful about")


def main():
    nb = json.loads(NB.read_text())
    cells = nb["cells"]
    section, has_source, produces = None, False, False
    problems, checked = [], 0

    def close():
        nonlocal checked
        if section is None:
            return
        if any(e in section for e in EXEMPT):
            return
        checked += 1
        if produces and not has_source:
            problems.append(section)

    for c in cells:
        src = "".join(c["source"])
        if c["cell_type"] == "markdown":
            heads = [l for l in src.split("\n") if re.match(r"^##\s+\d+\.", l)]
            if heads:
                close()
                section = heads[0].lstrip("# ").strip()
                has_source, produces = False, False
            if "**Source.**" in src or "**Source. " in src:
                has_source = True
        else:
            outs = c.get("outputs", [])
            if outs:
                produces = True
    close()

    print(f"checked {checked} sections that produce output")
    for p in problems:
        print(f"  MISSING **Source.**  {p}")
    if problems:
        print(f"\nFAIL: {len(problems)} section(s) do not say what produced them")
        return 1
    print("PASS: every section that produces output names its source")
    return 0


if __name__ == "__main__":
    sys.exit(main())
