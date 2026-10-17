"""Record commands, environment and logs for reproducible validation runs."""
import argparse,json,os,subprocess,sys,time,hashlib,zipfile
from pathlib import Path

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--log-dir',required=True)
    parser.add_argument('--timeout',type=float,default=300.)
    parser.add_argument('command',nargs=argparse.REMAINDER)
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    out=(root/args.log_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    local=root.parent/'native_libs/root/usr/lib/x86_64-linux-gnu'
    native=Path(env.get('QT_NATIVE_LIB_DIR',str(local)))
    if native.exists():
        env['LD_LIBRARY_PATH']=str(native)+':'+str(native/'openblas-pthread')+':'+env.get('LD_LIBRARY_PATH','')
        env.setdefault('DEVSIM_MATH_LIBS','libopenblas.so.0')
    env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg',PYTHONPATH=str(root))
    cmd=args.command
    if cmd and cmd[0]=='--':cmd=cmd[1:]
    if not cmd:parser.error('a command is required')
    if cmd[0]=='python':cmd[0]=sys.executable
    def manifest():
        paths=list(root.rglob('*.py'))+list((root/'configs').glob('*.json'))+[root/'requirements.txt']
        return {p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                for p in paths if 'results' not in p.relative_to(root).parts}
    before=manifest()
    (out/'source_manifest.json').write_text(json.dumps(before,indent=2)+'\n')
    with zipfile.ZipFile(out/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as snapshot:
        for path in before:snapshot.write(root/path,path)
    start=time.time()
    with (out/'execution.log').open('w') as log:
        try:
            process=subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout)
            code=process.returncode
        except subprocess.TimeoutExpired:
            code=124;log.write('\nVALIDATION_TIMEOUT: run stopped; no convergence claimed.\n')
    info={'command':cmd,'cwd':str(root),'python':sys.version,'exit_code':code,
          'runtime_s':time.time()-start,'native_library_dir':str(native) if native.exists() else None,
          'threads':1,'timeout_s':args.timeout,'source_unchanged_during_run':before==manifest()}
    (out/'execution.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info),flush=True)
    return code

if __name__=='__main__':raise SystemExit(main())
