"""New solve mathematics against NumPy; never original-DLL differential evidence."""
import argparse
import hashlib
import json
import os
import pathlib
import subprocess
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
p=argparse.ArgumentParser()
p.add_argument('--probe',required=True)
p.add_argument('--output')
args=p.parse_args()
rng=np.random.default_rng(20261026)
cases=[]
for n in range(1,33):
    for r in (1,3,5):
        a=rng.uniform(-1,1,(n,n))+np.eye(n)*(2*n)
        b=rng.uniform(-1,1,(n,r))
        cases.append((a,b))
cases.append((np.array([[0.,1.],[2.,3.]]),np.array([[5.,5.],[11.,9.]])))
text=''.join(f'{len(a)} {b.shape[1]}\n'+' '.join(format(float(x),'.17g') for x in a.ravel(order='F'))+'\n'+' '.join(format(float(x),'.17g') for x in b.ravel(order='F'))+'\n' for a,b in cases)
text+='2 1\n1 2 2 4\n1 2\n'
run=subprocess.run([args.probe],input=text,text=True,capture_output=True,check=True,timeout=120)
rows=[json.loads(line) for line in run.stdout.splitlines()]
if len(rows)!=len(cases)+1 or rows[-1]!={'singular':True}:raise RuntimeError('Missing results or singularity mismatch')
max_error=0.;max_backward=0.
for (a,b),row in zip(cases,rows):
    x=np.array(row['x']).reshape(b.shape,order='F');reference=np.linalg.solve(a,b)
    error=float(np.max(np.abs(x-reference)))
    backward=float(np.linalg.norm(a@x-b,np.inf)/(np.linalg.norm(a,np.inf)*np.linalg.norm(x,np.inf)+np.linalg.norm(b,np.inf)))
    if error>2e-13 or backward>64*len(a)*np.finfo(float).eps:raise RuntimeError(f'Accuracy mismatch {len(a)} {error} {backward}')
    max_error=max(max_error,error);max_backward=max(max_backward,backward)
report={'seed':20261026,'cases':len(rows),'numpy_version':np.__version__,'max_absolute_solution_difference':max_error,'max_scaled_backward_error':max_backward,'probe_sha256':hashlib.sha256(pathlib.Path(args.probe).read_bytes()).hexdigest(),'scope':'Independent NumPy reference, not original DLL'}
if args.output:pathlib.Path(args.output).write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
