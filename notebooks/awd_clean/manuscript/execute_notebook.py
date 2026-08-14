"""Execute a notebook and save its outputs, without nbconvert.

The `das` env has no nbconvert, nbclient, or nbformat -- which is why
`run_nb_cells.py` exists.  But `run_nb_cells.py` only proves the cells run: it
`exec`s them in a shared namespace and writes nothing back, so the committed
notebook has no figures in it.  Anyone opening it sees code and blank space until
they run it themselves.  For an advisor reading the analysis, that is the same as
having no figures.

This drives the `das` kernel over jupyter_client (which IS installed) and captures
the real outputs -- stdout, results, and inline PNGs -- into the notebook JSON,
in nbformat v4 shape.  No new dependencies.

Deliberate choices:

  * `%matplotlib inline` and the figure dpi are set in a setup cell that is
    executed but NOT written into the notebook, so the committed source stays
    exactly what a reader runs, and only the outputs are added.
  * A cell that raises stops the run and exits non-zero, after writing what it
    got, so the failure is diagnosable rather than silently half-executed.
  * dpi is capped because a notebook with a dozen 16x6-inch figures at
    matplotlib's default dpi is tens of MB, and this repo already has 44 MB
    notebooks that are painful to review.

Usage:
    python execute_notebook.py IN.ipynb [-o OUT.ipynb] [--kernel das]
                               [--timeout 1800] [--dpi 110]
"""
import argparse
import json
import queue
import sys
import time
from pathlib import Path

from jupyter_client.manager import KernelManager

SETUP = """\
%matplotlib inline
%config InlineBackend.figure_formats = ['png']
import matplotlib as _mpl
_mpl.rcParams['figure.dpi'] = {dpi}
_mpl.rcParams['savefig.dpi'] = {dpi}
"""


def collect(kc, msg_id, timeout):
    """Drain iopub for one execution, returning nbformat-v4 outputs."""
    outputs, streams = [], {}
    idle = False
    got_reply = False
    exec_count = None
    error = None
    deadline = time.time() + timeout

    while not (idle and got_reply):
        if time.time() > deadline:
            raise TimeoutError(f'cell exceeded {timeout}s')
        # shell reply carries execution_count and the error status
        try:
            r = kc.get_shell_msg(timeout=0.2)
            if r['parent_header'].get('msg_id') == msg_id:
                got_reply = True
                exec_count = r['content'].get('execution_count')
                if r['content'].get('status') == 'error':
                    error = r['content']
        except queue.Empty:
            pass
        try:
            m = kc.get_iopub_msg(timeout=0.2)
        except queue.Empty:
            continue
        if m['parent_header'].get('msg_id') != msg_id:
            continue
        t, c = m['msg_type'], m['content']
        if t == 'status':
            idle = idle or c['execution_state'] == 'idle'
        elif t == 'stream':
            # coalesce consecutive writes to the same stream, as Jupyter does
            if c['name'] in streams:
                streams[c['name']]['text'].append(c['text'])
            else:
                o = {'output_type': 'stream', 'name': c['name'],
                     'text': [c['text']]}
                streams[c['name']] = o
                outputs.append(o)
        elif t in ('display_data', 'execute_result'):
            o = {'output_type': t, 'data': c['data'],
                 'metadata': c.get('metadata', {})}
            if t == 'execute_result':
                o['execution_count'] = c.get('execution_count')
            outputs.append(o)
            streams.clear()          # a figure breaks the stream run
        elif t == 'error':
            outputs.append({'output_type': 'error', 'ename': c['ename'],
                            'evalue': c['evalue'], 'traceback': c['traceback']})
    return outputs, exec_count, error


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('notebook')
    ap.add_argument('-o', '--output')
    ap.add_argument('--kernel', default='das')
    ap.add_argument('--timeout', type=int, default=1800)
    ap.add_argument('--dpi', type=int, default=110)
    a = ap.parse_args()

    src = Path(a.notebook)
    dst = Path(a.output) if a.output else src
    nb = json.loads(src.read_text())
    code = [c for c in nb['cells'] if c['cell_type'] == 'code']
    print(f'{src.name}: {len(nb["cells"])} cells, {len(code)} code, '
          f'kernel {a.kernel}, dpi {a.dpi}\n')

    km = KernelManager(kernel_name=a.kernel)
    km.start_kernel(cwd=str(src.parent.resolve()))
    kc = km.client()
    kc.start_channels()
    failed = None
    try:
        kc.wait_for_ready(timeout=120)
        # setup is executed but not recorded in the notebook
        mid = kc.execute(SETUP.format(dpi=a.dpi), store_history=False)
        _, _, err = collect(kc, mid, 120)
        if err:
            raise RuntimeError(f'setup cell failed: {err["ename"]}: '
                               f'{err["evalue"]}')

        n = 0
        for i, cell in enumerate(nb['cells']):
            if cell['cell_type'] != 'code':
                continue
            n += 1
            code_str = ''.join(cell['source'])
            if not code_str.strip():
                cell['outputs'], cell['execution_count'] = [], None
                continue
            t0 = time.time()
            mid = kc.execute(code_str)
            outs, cnt, err = collect(kc, mid, a.timeout)
            cell['outputs'] = outs
            cell['execution_count'] = cnt
            imgs = sum(1 for o in outs
                       if 'image/png' in o.get('data', {}))
            print(f'  [{n:2d}/{len(code)}] cell {i:2d}  {time.time()-t0:6.1f}s  '
                  f'{len(outs)} outputs, {imgs} figure(s)'
                  f'{"   <-- ERROR" if err else ""}')
            if err:
                print('\n'.join(err.get('traceback', [])))
                failed = (i, err)
                break
    finally:
        kc.stop_channels()
        km.shutdown_kernel(now=True)

    dst.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + '\n')
    mb = dst.stat().st_size / 1e6
    total_imgs = sum(1 for c in nb['cells'] for o in c.get('outputs', [])
                     if 'image/png' in o.get('data', {}))
    print(f'\nWrote {dst}  ({mb:.1f} MB, {total_imgs} embedded figures)')
    if failed:
        print(f'FAILED at cell {failed[0]}: {failed[1]["ename"]}: '
              f'{failed[1]["evalue"]}')
        sys.exit(1)

    # A run where every cell "succeeds" and no figure is captured is the exact
    # failure this script was written to fix: an imported module had forced the
    # Agg backend, so plt.show() was a silent no-op.  Cells execute cleanly and
    # the notebook looks fine until someone opens it.  Check for it explicitly.
    drew = [i for i, c in enumerate(nb['cells'])
            if c['cell_type'] == 'code'
            and any(k in ''.join(c['source'])
                    for k in ('plt.subplots', 'plt.figure'))]
    silent = [i for i in drew
              if not any('image/png' in o.get('data', {})
                         for o in nb['cells'][i].get('outputs', []))]
    print(f'{len(drew)} cells draw figures, {len(drew)-len(silent)} produced one')
    if silent:
        print(f'\nERROR: cells {silent} create figures but captured none.')
        print('Almost certainly a non-interactive backend: check whether an')
        print(f'imported module calls matplotlib.use("Agg").  Backend in use is')
        print('reported by the first code cell.')
        sys.exit(2)
    print('All cells executed cleanly, and every drawing cell produced a figure.')


if __name__ == '__main__':
    main()
