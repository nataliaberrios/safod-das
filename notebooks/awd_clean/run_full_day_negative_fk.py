
from pathlib import Path
import subprocess,sys,json,time
date='2025-02-20'; n=1443; chunk=10
root=Path(__file__).resolve().parent; out=root/'ambient_transfer'; log=out/'full_day_negative_progress.log'
das=sys.executable
with log.open('w') as lf:
    for start in range(0,n,chunk):
        count=min(chunk,n-start)
        cmd=[das,str(root/'ambient_fk_transfer_test.py'),'--date',date,'--start',str(start),'--nfiles',str(count),'--modes','negative']
        t=time.time(); lf.write(f'START {start} {count}\n'); lf.flush()
        r=subprocess.run(cmd,cwd=root.parent,capture_output=True,text=True)
        lf.write(r.stdout); lf.write(r.stderr); lf.write(f'END {start} rc={r.returncode} elapsed={time.time()-t:.1f}\n'); lf.flush()
        if r.returncode!=0: raise SystemExit(f'chunk {start} failed rc={r.returncode}')
print('all chunks complete')
