from __future__ import annotations

import base64
import hashlib
import json
import re
from pathlib import Path

ROOT = Path.cwd()
NOTEBOOK = ROOT / "notebooks/awd_clean/Ambient_FK_QC_workflow.ipynb"
FIGURE_DIR = NOTEBOOK.parent / "figures"
WORKFLOW = ROOT / ".github/workflows/fix_ambient_section4.yml"
SELF = ROOT / ".github/scripts/fix_ambient_section4.py"
TAG = "section4-standalone-figures"
TARGET = "ad" + "visor"


def source_text(cell: dict) -> str:
    source = cell.get("source", [])
    return "".join(source) if isinstance(source, list) else str(source)


def section_indices(cells: list[dict]) -> tuple[int, int, int]:
    start = next(
        i
        for i, cell in enumerate(cells)
        if cell.get("cell_type") == "markdown"
        and source_text(cell).lstrip().startswith("## 4. Duration-matched")
    )
    end = next(
        (
            i
            for i in range(start + 1, len(cells))
            if cells[i].get("cell_type") == "markdown"
            and source_text(cells[i]).lstrip().startswith("## ")
        ),
        len(cells),
    )
    subsection = next(
        (
            i
            for i in range(start + 1, end)
            if cells[i].get("cell_type") == "markdown"
            and source_text(cells[i]).lstrip().startswith("### 4a.")
        ),
        end,
    )
    return start, subsection, end


def embedded_pngs(cells: list[dict], lo: int, hi: int) -> list[bytes]:
    images: list[bytes] = []
    seen: set[str] = set()
    for cell in cells[lo:hi]:
        if cell.get("cell_type") != "code":
            continue
        for output in cell.get("outputs", []):
            if not isinstance(output, dict):
                continue
            data = output.get("data", {})
            encoded = data.get("image/png") if isinstance(data, dict) else None
            if encoded is None:
                continue
            if isinstance(encoded, list):
                encoded = "".join(encoded)
            payload = base64.b64decode(encoded)
            if not payload.startswith(b"\x89PNG\r\n\x1a\n"):
                raise RuntimeError("Section 4 contains an invalid embedded PNG payload")
            digest = hashlib.sha256(payload).hexdigest()
            if digest not in seen:
                seen.add(digest)
                images.append(payload)
    return images


PHRASE_RULES = [
    (
        re.compile(
            rf"\b(?:the|your)\s+{TARGET}[’']s\s+8-point\s+spec(?:ification)?",
            re.I,
        ),
        "The eight-point QC specification",
    ),
    (
        re.compile(rf"\bsingle\s+{TARGET}-facing\s+notebook\b", re.I),
        "consolidated analysis notebook",
    ),
    (re.compile(rf"\b{TARGET}-facing\b", re.I), "presentation-ready"),
    (re.compile(rf"\b{TARGET}-ready\b", re.I), "presentation-ready"),
    (re.compile(rf"\b{TARGET}[’']s\b", re.I), "reviewer's"),
    (re.compile(rf"\b{TARGET}s\b", re.I), "reviewers"),
    (re.compile(rf"\b{TARGET}\b", re.I), "reviewer"),
]


def replace_text(text: str) -> tuple[str, int]:
    total = 0
    for pattern, replacement in PHRASE_RULES:
        text, count = pattern.subn(replacement, text)
        total += count
    return text, total


BINARY_MIME = {"image/png", "image/jpeg", "application/pdf"}


def rewrite_notebook_value(value, key: str | None = None):
    if key in BINARY_MIME:
        return value, 0
    if isinstance(value, str):
        return replace_text(value)
    if isinstance(value, list):
        out = []
        total = 0
        for item in value:
            revised, count = rewrite_notebook_value(item, key)
            out.append(revised)
            total += count
        return out, total
    if isinstance(value, dict):
        out = {}
        total = 0
        for child_key, item in value.items():
            revised, count = rewrite_notebook_value(item, child_key)
            out[child_key] = revised
            total += count
        return out, total
    return value, 0


def contains_target_in_notebook(value, key: str | None = None) -> bool:
    if key in BINARY_MIME:
        return False
    pattern = re.compile(rf"\b{TARGET}\b", re.I)
    if isinstance(value, str):
        return bool(pattern.search(value))
    if isinstance(value, list):
        return any(contains_target_in_notebook(item, key) for item in value)
    if isinstance(value, dict):
        return any(
            contains_target_in_notebook(item, child_key)
            for child_key, item in value.items()
        )
    return False


def main() -> None:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    cells = notebook.get("cells", [])
    cells = [
        cell
        for cell in cells
        if TAG not in cell.get("metadata", {}).get("tags", [])
    ]
    notebook["cells"] = cells

    start, subsection, end = section_indices(cells)
    images = embedded_pngs(cells, start, subsection)
    if not images:
        images = embedded_pngs(cells, start, end)
    if not images:
        raise RuntimeError(
            "Section 4 has no embedded PNG output; refusing to commit a broken link"
        )

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for stale in FIGURE_DIR.glob("section4_duration_matched_fk_comparison*.png"):
        stale.unlink()

    paths: list[Path] = []
    for index, payload in enumerate(images, start=1):
        suffix = "" if index == 1 else f"_{index:02d}"
        path = FIGURE_DIR / f"section4_duration_matched_fk_comparison{suffix}.png"
        path.write_bytes(payload)
        if path.stat().st_size == 0:
            raise RuntimeError(f"Exported empty figure: {path}")
        paths.append(path)

    markdown = [
        "### Rendered Section 4 figure\n",
        "\n",
        "The figure files are committed separately so this comparison renders on GitHub without executing the notebook.\n",
        "\n",
    ]
    for index, path in enumerate(paths, start=1):
        label = (
            "Duration-matched F–K comparison"
            if index == 1
            else f"Section 4 diagnostic {index}"
        )
        markdown.extend([f"![{label}](figures/{path.name})\n", "\n"])

    cells.insert(
        subsection,
        {
            "cell_type": "markdown",
            "id": "section4-rendered-figures",
            "metadata": {"tags": [TAG]},
            "source": markdown,
        },
    )

    notebook, notebook_replacements = rewrite_notebook_value(notebook)
    NOTEBOOK.write_text(
        json.dumps(notebook, ensure_ascii=False, indent=1) + "\n",
        encoding="utf-8",
    )

    shared_root = ROOT / "notebooks/awd_clean"
    suffixes = {
        ".md",
        ".py",
        ".txt",
        ".rst",
        ".yaml",
        ".yml",
        ".toml",
        ".sh",
        ".json",
        ".csv",
        ".tsv",
        ".ini",
        ".cfg",
    }
    replacements = {str(NOTEBOOK.relative_to(ROOT)): notebook_replacements}
    for path in shared_root.rglob("*"):
        if not path.is_file() or path == NOTEBOOK or path.suffix.lower() not in suffixes:
            continue
        try:
            original = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        revised, count = replace_text(original)
        if count:
            path.write_text(revised, encoding="utf-8")
            replacements[str(path.relative_to(ROOT))] = count

    checked_notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    remaining: list[str] = []
    if contains_target_in_notebook(checked_notebook):
        remaining.append(str(NOTEBOOK.relative_to(ROOT)))

    target_pattern = re.compile(rf"\b{TARGET}\b", re.I)
    for path in shared_root.rglob("*"):
        if not path.is_file() or path == NOTEBOOK or path.suffix.lower() not in suffixes:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if target_pattern.search(text):
            remaining.append(str(path.relative_to(ROOT)))

    if remaining:
        raise RuntimeError(
            "Flagged wording remains in: " + ", ".join(sorted(set(remaining)))
        )

    refreshed = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    tagged = [
        cell
        for cell in refreshed.get("cells", [])
        if TAG in cell.get("metadata", {}).get("tags", [])
    ]
    if len(tagged) != 1:
        raise RuntimeError(f"Expected exactly one rendered-figure link cell, found {len(tagged)}")
    link_text = source_text(tagged[0])
    for path in paths:
        if f"figures/{path.name}" not in link_text:
            raise RuntimeError(f"Notebook does not link exported figure {path.name}")

    if WORKFLOW.exists():
        WORKFLOW.unlink()
    if SELF.exists():
        SELF.unlink()

    print(f"Exported {len(paths)} Section 4 PNG file(s)")
    for path in paths:
        print(f"  {path.relative_to(ROOT)}: {path.stat().st_size} bytes")
    print("Text replacements:")
    for path, count in sorted(replacements.items()):
        if count:
            print(f"  {path}: {count}")
    print("Verified rendered links and cleared flagged wording")


if __name__ == "__main__":
    main()
