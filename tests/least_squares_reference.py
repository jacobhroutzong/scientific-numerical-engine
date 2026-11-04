"""Seeded mathematical references and optional API benchmarks; no original DLL."""
import argparse, hashlib, json, os, pathlib, platform, subprocess, time
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--probe',required=True)
p.add_argument('--output')
p.add_argument('--benchmark',action='store_true')
p.add_argument('--compiler',default='not recorded (validation only)')
p.add_argument('--flags',default='not recorded (validation only)')
p.add_argument('--cpu',default=platform.processor())
args=p.parse_args()
if args.benchmark and args.compiler.startswith('not recorded'):raise SystemExit('Benchmark requires --compiler and --flags metadata')
rng=np.random.default_rng(20261104)
cases=[]
for n in range(1,21):
    for nrhs in (1,3):
        a=rng.uniform(-1,1,(2*n+3,n));a[:n]+=np.eye(n)*n
        b=rng.uniform(-1,1,(len(a),nrhs))
        cases.append((f'tall-{n}-rhs{nrhs}',a,b,'ok'))
for n in (4,8,16,32,64):
    a=rng.uniform(-1,1,(n,n))+np.eye(n)*(2*n)
    cases.append((f'square-{n}',a,rng.uniform(-1,1,(n,3)),'ok'))
for n in (4,8):
    h=1./(np.arange(n)[:,None]+np.arange(n)[None,:]+1.)
    cases.append((f'hilbert-{n}',h,h@np.ones((n,1)),'ok'))
cases += [('near-singular-accepted',np.diag([1.,1e-12]),np.array([[1.],[1e-12]]),'ok'),
          ('near-singular-rejected',np.diag([1.,1e-18]),np.ones((2,1)),'rank_rejected'),
          ('singular',np.array([[1.,2.],[2.,4.],[3.,6.]]),np.ones((3,1)),'rank_rejected')]

def wire(a,b,estimate=1,repetitions=1):
    return f'{a.shape[0]} {a.shape[1]} {b.shape[1]} -1 {estimate} {repetitions}\n'+' '.join(format(x,'.17g') for x in a.ravel(order='F'))+'\n'+' '.join(format(x,'.17g') for x in b.ravel(order='F'))+'\n'
def run(text):
    result=subprocess.run([args.probe],input=text,text=True,capture_output=True,timeout=180,check=True)
    return [json.loads(line) for line in result.stdout.splitlines()]
rows=run(''.join(wire(a,b) for _,a,b,_ in cases));assert len(rows)==len(cases)
results=[]
for (name,a,b,status),row in zip(cases,rows):
    assert row['status']==status,(name,row)
    reference,_,rank,_=np.linalg.lstsq(a,b,rcond=np.finfo(float).eps*max(a.shape))
    record={'case':name,'m':a.shape[0],'n':a.shape[1],'rhs':b.shape[1],'status':status,'numpy_svd_rank':int(rank),'input_sha256':hashlib.sha256(wire(a,b).encode()).hexdigest()}
    if status=='ok':
        x=np.array(row['x']).reshape((a.shape[1],b.shape[1]),order='F')
        cond=float(np.linalg.cond(a))
        relative=float(np.linalg.norm(x-reference)/max(np.linalg.norm(reference),1e-300))
        tolerance=min(.05,256*np.finfo(float).eps*max(a.shape)*max(cond,1))
        assert relative<=tolerance,(name,relative,tolerance)
        residual=a@x-b; norm_a=np.linalg.norm(a,'fro')
        for j,d in enumerate(row['diagnostics']):
            residual_norm=float(np.linalg.norm(residual[:,j]))
            assert abs(d[0]-residual_norm)<=1e-11*max(1,np.linalg.norm(b[:,j])),(name,d,residual_norm)
            expected=residual_norm/(norm_a*np.linalg.norm(x[:,j])+np.linalg.norm(b[:,j]))
            assert abs(d[1]-expected)<1e-12,(name,d,expected)
            if residual_norm>1e-10*np.linalg.norm(b[:,j]):assert d[2]<1e-10,(name,d)
        _,r=np.linalg.qr(a,mode='reduced');r_condition=float(np.linalg.cond(r,np.inf))
        assert abs(row['condition_inf_r']-r_condition)/r_condition<max(1e-11,tolerance),(name,row['condition_inf_r'],r_condition)
        record.update(relative_solution_difference=relative,comparison_tolerance=tolerance,condition2_a_reference=cond,condition_inf_r=row['condition_inf_r'],max_scaled_residual=max(d[1] for d in row['diagnostics']))
    results.append(record)

timings=[]
if args.benchmark:
    selected=[c for c in cases if c[0] in {'tall-4-rhs1','tall-16-rhs3','square-4','square-8','square-16','square-32','square-64','hilbert-8','near-singular-accepted','near-singular-rejected','singular'}]
    for name,a,b,status in selected:
        if status!='ok':
            timings.append({'case':name,'status':status,'timing':'omitted: rejection and NumPy minimum-norm solve are different contracts'});continue
        repetitions=32 if a.shape[1]<=16 else 4
        for _ in range(3):np.linalg.lstsq(a,b,rcond=np.finfo(float).eps*max(a.shape))
        run(wire(a,b,estimate=0,repetitions=repetitions))
        native=run(wire(a,b,estimate=0,repetitions=repetitions)*7)
        samples=[]
        for _ in range(7):
            start=time.perf_counter_ns()
            for __ in range(repetitions):np.linalg.lstsq(a,b,rcond=np.finfo(float).eps*max(a.shape))
            samples.append((time.perf_counter_ns()-start)/repetitions)
        ns=[x['nanoseconds_per_call'] for x in native]
        timings.append({'case':name,'m':len(a),'n':a.shape[1],'rhs':b.shape[1],'calls_per_sample':repetitions,'cpp_samples_ns':ns,'numpy_samples_ns':samples,'cpp_median_ns':float(np.median(ns)),'numpy_median_ns':float(np.median(samples))})
config=getattr(np.__config__,'CONFIG',{}).get('Build Dependencies',{}).get('blas',{})
report={'seed':20261104,'scope':'Independent NumPy/SVD and QR mathematics; not original-binary differential evidence','numpy':np.__version__,'python':platform.python_version(),'os':platform.platform(),'cpu':args.cpu,'compiler':args.compiler,'flags':args.flags,'blas':{k:config.get(k) for k in ['name','version','openblas configuration']},'requested_blas_threads':1,'probe_sha256':hashlib.sha256(pathlib.Path(args.probe).read_bytes()).hexdigest(),'cases':results,'benchmarks':timings,'timing_boundary':'Native timer excludes process/serialization, includes validation, QR, explicit thin Q, projection, solve and residual diagnostics; optional R condition disabled for timing. NumPy times lstsq API using SVD and its API overhead; does not compute our diagnostics. No speed-superiority or kernel-equivalence claim.'}
if args.output:pathlib.Path(args.output).write_text(json.dumps(report,indent=2),encoding='utf-8',newline='\n')
print(json.dumps({'cases':len(results),'accepted':sum(r['status']=='ok' for r in results),'rank_rejected':sum(r['status']=='rank_rejected' for r in results),'benchmarks':len(timings),'max_relative_solution_difference':max(r.get('relative_solution_difference',0) for r in results)},indent=2))
