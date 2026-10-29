"""Build the actual engine sources. No emulation, networking or dependency download."""
import argparse
import hashlib
import json
import os
import pathlib
import subprocess
import sys

p=argparse.ArgumentParser()
p.add_argument('--emsdk',required=True,type=pathlib.Path)
p.add_argument('--engine',required=True,choices=['numerical','symbolic'])
p.add_argument('--source',required=True,type=pathlib.Path)
p.add_argument('--out',required=True,type=pathlib.Path)
args=p.parse_args()
args.source=args.source.resolve();args.emsdk=args.emsdk.resolve();args.out=args.out.resolve()
compiler=args.emsdk/'upstream/emscripten/em++.py'
env=os.environ.copy();env['EM_CONFIG']=str((args.emsdk/'.emscripten').resolve())
temporary=args.source/'.wasm-tmp';temporary.mkdir(parents=True,exist_ok=True)
env['TEMP']=env['TMP']=env['EMCC_TEMP_DIR']=str(temporary)
version=subprocess.check_output([sys.executable,str(compiler),'--version'],env=env,text=True).splitlines()[0]
if '6.0.13' not in version:raise SystemExit('Expected pinned Emscripten 6.0.13, got '+version)
sources=['src/lu.cpp','src/solve.cpp'] if args.engine=='numerical' else ['src/integer.cpp','src/rational.cpp','src/expression.cpp']
sources+=['browser/adapter.cpp']
args.out.mkdir(parents=True,exist_ok=True)
command=[sys.executable,str(compiler),'-std=c++20','-O2','-ffp-contract=off','-fexceptions','--bind','--no-entry',
         '-sMODULARIZE=1','-sEXPORT_ES6=1','-sENVIRONMENT=worker,node','-sALLOW_MEMORY_GROWTH=1',
         '-sINITIAL_MEMORY=16777216','-sMAXIMUM_MEMORY=134217728','-sSTACK_SIZE=1048576','-sFILESYSTEM=0',
         '-I',str(args.source/'include'),'-ffile-prefix-map='+str(args.source)+'=/engine',
         '-ffile-prefix-map='+str(args.emsdk)+'=/emsdk',*[str(args.source/f) for f in sources],'-o',str(args.out/(args.engine+'.mjs'))]
subprocess.run(command,env=env,check=True)
manifest={'compiler':version,'engine':args.engine,'sources':{f:hashlib.sha256((args.source/f).read_bytes()).hexdigest() for f in sources},
          'headers':{str(f.relative_to(args.source)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((args.source/'include').rglob('*.hpp'))},
          'buildScriptSha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
          'outputs':{f.name:{'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size} for f in [args.out/(args.engine+'.mjs'),args.out/(args.engine+'.wasm')]},
          'flags':command[2:command.index('-I')],'scope':'Private development build; publication licensing approval pending'}
def canonical_hash(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
manifest['sources']={f:canonical_hash(args.source/f) for f in sources}
manifest['headers']={str(f.relative_to(args.source)).replace('\\','/'):canonical_hash(f) for f in sorted((args.source/'include').rglob('*.hpp'))}
manifest['buildScriptSha256']=canonical_hash(pathlib.Path(__file__))
manifest['sourceHashing']='SHA-256 of UTF-8 source bytes with CRLF normalized to LF'
manifest['sourceCommit']=subprocess.check_output(['git','-C',str(args.source),'rev-parse','HEAD'],text=True).strip()
manifest['sourceTreeClean']=not subprocess.check_output(['git','-C',str(args.source),'status','--porcelain','--untracked-files=normal'],text=True).strip()
manifest['scope']='Source-only project license and third-party notices apply; retain runtime notices when distributing generated artifacts'
(args.out/(args.engine+'-manifest.json')).write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps(manifest,indent=2))
