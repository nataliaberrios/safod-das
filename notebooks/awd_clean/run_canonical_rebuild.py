"""Execute the legacy builder with corrected AWD bounds, without overwriting it."""

from pathlib import Path

root = Path(__file__).resolve().parents[1]
legacy = root / "paired_stack_job_deep_all.py"
source = legacy.read_text()
replacements = {
    'OUT_NPZ = os.path.join(FIG_DIR, "epoch_stacks_paired_deep_all.npz")':
        'OUT_NPZ = str(Path(__file__).resolve().parent / '
        '"canonical_epoch_stacks_paired_deep_all.npz")',
    'awd_files = [f for f in all_nano if in_nano_survey_window(f)]':
        'awd_files = all_nano',
    'print(f"setup: {len(awd_files)} AWD files -> {n_epochs} bursts")':
        'print(f"setup: indexing all {len(awd_files)} Nano files")',
}
for old, new in replacements.items():
    if old not in source:
        raise RuntimeError(f"Expected legacy line not found: {old}")
    source = source.replace(old, new, 1)

exec(compile(source, str(legacy), "exec"), {"__name__": "__main__", "__file__": __file__})
