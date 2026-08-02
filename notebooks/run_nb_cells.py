"""Execute a notebook's code cells in order to verify it runs clean.

The das env has no nbconvert, and a notebook that raises on cell 4 is worse than
no notebook, so this reads the .ipynb as JSON and execs the code cells in one
shared namespace -- the same thing Run All does, minus the kernel.
"""
import json, sys, traceback

path = sys.argv[1]
nb = json.load(open(path))
cells = [c for c in nb['cells'] if c['cell_type'] == 'code']
ns = {'__name__': '__main__'}
for i, c in enumerate(cells, 1):
    src = ''.join(c['source'])
    if not src.strip():
        continue
    print(f'\n----- cell {i}/{len(cells)} -----', flush=True)
    try:
        exec(compile(src, f'<cell {i}>', 'exec'), ns)
    except Exception:
        print(f'FAILED on cell {i}:', flush=True)
        traceback.print_exc()
        sys.exit(1)
print('\nALL CELLS OK')
