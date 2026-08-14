# Manuscript — start here

Everything paper-facing for *Contrasting Borehole DAS Installations at SAFOD*.
The analysis itself lives one level up in `awd_clean/`; this folder holds only
what goes into, or explains, the paper.

## The three files that matter

| File | What it is |
|---|---|
| **[`AWD_reproduce_analysis.ipynb`](AWD_reproduce_analysis.ipynb)** | **The reproducibility notebook.** Every analysis step, start to finish, runnable and inspectable. Start here to follow or check the analysis. |
| **[`MANUSCRIPT.md`](MANUSCRIPT.md)** | The whole paper, abstract → conclusions, ~9,400 words. **Read this one.** Generated — do not edit it directly. |
| **[`REPRODUCE.md`](REPRODUCE.md)** | What code produced every number, in the order it must run. Includes a two-minute plain-language explanation of the experiment. |
| **[`figures/FIGURES.md`](figures/FIGURES.md)** | The figure set with manuscript numbering, and which script made each one. |

## Two commands

```bash
cd notebooks/awd_clean/manuscript
python assemble_manuscript.py     # sections/ -> MANUSCRIPT.md
python make_paper_figures.py      # ../*.png  -> figures/
python build_notebook.py          # rebuild the notebook from its builder
```

The notebook needs ~8 GB (the Deep stack array is 2.2 GB), so run it under
`sh_dev --mem=16G` or in a job, not on a login node. Smoke-test it headless with
`python run_nb_cells.py AWD_reproduce_analysis.ipynb`.

Run the first after editing any section. It fails loudly if a heading marker has
moved, rather than quietly emitting a truncated paper.

## Layout

```
manuscript/
├── README.md                 you are here
├── MANUSCRIPT.md             generated — the paper
├── REPRODUCE.md              what was run
├── assemble_manuscript.py    sections/ -> MANUSCRIPT.md
├── make_paper_figures.py     analysis PNGs -> figures/
├── sections/                 edit these, not MANUSCRIPT.md
│   ├── MANUSCRIPT_ABSTRACT_CONCLUSIONS.md
│   ├── MANUSCRIPT_INTRODUCTION.md
│   ├── MANUSCRIPT_METHODS.md          §2, §3.1–3.2, §3.5–3.6
│   ├── DEEP_DVV_METHODS_DRAFT.md      §3.3–3.4, §3.7–3.8, §4.6
│   ├── MANUSCRIPT_RESULTS.md          §4.1–4.5
│   └── MANUSCRIPT_DISCUSSION.md       §5
└── figures/                  generated
```

Sections are split by where they were drafted, not by manuscript order. The
assembler interleaves them; `MANUSCRIPT.md` reads in the right order.

## What deliberately stayed in `awd_clean/`

These are the analysis record, not the paper:

- `DEEP_DVV_STATUS.md` — the authoritative claim/evidence table. **Where a number
  here disagrees with the manuscript, this file wins.**
- `DEEP_DVV_PREREGISTRATION.md` — the frozen analysis design and decision rules,
  timestamped before results existed.
- All analysis scripts and their outputs.

## Still open

| Item | Who |
|---|---|
| Figure 1 — experiment geometry schematic | never made; see below |
| Interrogator model numbers; physical depth at channel 1702; which cable is which | Ettore |
| Model B to the supplement as a model-dependence bound | not yet written in; Model A **is** integrated (§3.9, §5.7) |
| Section numbering gaps if sections are added | re-run the assembler |
