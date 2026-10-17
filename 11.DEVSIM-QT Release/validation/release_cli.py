"""Release CLI: environment, existing device API, refinement, export and certificate."""
import argparse,os,sys,json,subprocess,importlib.metadata,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

def bootstrap():
    native=Path(os.environ.get('QT_NATIVE_LIB_DIR',str(ROOT.parent/'native_libs/root/usr/lib/x86_64-linux-gnu')))
    if native.exists() and not os.environ.get('QT_NATIVE_BOOTSTRAPPED'):
        env=os.environ.copy();env.update(QT_NATIVE_BOOTSTRAPPED='1',LD_LIBRARY_PATH=str(native)+':'+str(native/'openblas-pthread')+':'+env.get('LD_LIBRARY_PATH',''))
        env.setdefault('DEVSIM_MATH_LIBS','libopenblas.so.0')
        env.setdefault('OPENBLAS_NUM_THREADS','1');env.setdefault('OMP_NUM_THREADS','1');env['MPLBACKEND']='Agg'
        os.execvpe(sys.executable,[sys.executable,*sys.argv],env)

def doctor(out=None):
    info={'python':sys.version,'cpu_affinity':len(os.sched_getaffinity(0)),'packages':{},'native':{}}
    for name in ('numpy','scipy','pytest','gmsh','devsim','matplotlib'):
        try:info['packages'][name]=importlib.metadata.version(name)
        except Exception as exc:info['packages'][name]=repr(exc)
    for name in ('gmsh','devsim'):
        try:__import__(name);info['native'][name]='PASS'
        except Exception as exc:info['native'][name]={'status':'FAIL','reason':repr(exc)}
    for name in ('memory.max','cpu.max'):
        p=Path('/sys/fs/cgroup')/name;info[name]=p.read_text().strip() if p.exists() else None
    print(json.dumps(info,indent=2))
    if out:Path(out).write_text(json.dumps(info,indent=2))
    return 0 if all(v=='PASS' for v in info['native'].values()) else 2

def compare(paths):
    import numpy as np
    from config import Config
    from scipy.spatial import cKDTree
    pairs=[]
    for old,new in zip(paths[:-1],paths[1:]):
        a=np.load(old);b=np.load(new);ca=Config(**json.loads(str(a['config_json'])));cb=Config(**json.loads(str(b['config_json'])))
        # This comparison intentionally refuses physical S/D geometry changes.
        shape=('ls','lg','ld','spacer','width','thickness','sd_width','sd_thickness','tox','geometry_profile')
        if any(getattr(ca,k)!=getattr(cb,k) for k in shape):raise ValueError('Physical geometry changed: category C; not numerical buffer convergence')
        xyz=b['xyz_nm'];dist,ind=cKDTree(a['xyz_nm']).query(xyz)
        na=a['electron_density_cm3'][ind];nb=b['electron_density_cm3'];x=xyz[:,0]
        regions={'channel':(x>=cb.gate_start)&(x<=cb.gate_end),'junction_extension':((x>=cb.ls)&(x<cb.gate_start))|((x>cb.gate_end)&(x<=cb.drain_start)), 'reservoir':(x<cb.ls)|(x>cb.drain_start)}
        row={'old':str(old),'new':str(new),'density_transfer':'nearest cell value; spatial discretization diagnostic',
             'max_sample_distance_nm':float(max(dist)),'mu_difference_ev':abs(float(a['mu_l_ev']-b['mu_l_ev'])),
             'current_absolute_difference_a':abs(float(a['id_a']-b['id_a'])),
             'current_relative_difference':abs(float(a['id_a']-b['id_a']))/max(abs(float(b['id_a'])),1e-15),
             'barrier_phi_max_difference_v':float(max(abs(a['phi_v'][ind]-b['phi_v'])))}
        for name,mask in regions.items():row[name+'_density_l2']=float(np.linalg.norm((na-nb)[mask])/max(np.linalg.norm(nb[mask]),1e-30))
        row['density_target_pass']=all(row[k+'_density_l2']<.01 for k in regions)
        row['lead_mu_target_pass']=row['mu_difference_ev']<1e-4
        row['transport_target_pass']=row['current_relative_difference']<.01 and max(abs(float(a['vd_v'])),abs(float(b['vd_v'])))>0
        row['status']='PASS' if all(row[k] for k in ('density_target_pass','lead_mu_target_pass','transport_target_pass')) else 'FAIL'
        pairs.append(row)
    return {'comparisons':pairs,'largest_basis_is_not_an_exact_reference':True,'status':'PASS' if pairs and all(p['status']=='PASS' for p in pairs) else 'INCONCLUSIVE' if not pairs else 'FAIL'}

def refine(args):
    from config import Config
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True);base=json.loads(Path(args.config).read_text());states=[];runs=[]
    for pair in args.pairs.split(','):
        m,sd=map(int,pair.split('/'));cfg=dict(base,modes=m,sd_modes=sd);config=out/f'modes_{m}_{sd}.json';config.write_text(json.dumps(cfg,indent=2));runout=out/f'modes_{m}_{sd}'
        cmd=[sys.executable,str(ROOT/'run.py'),'point','--config',str(config.resolve()),'--out',str(runout.resolve()),'--vg',str(args.vg),'--vd',str(args.vd)]
        tic=time.time()
        with (out/f'modes_{m}_{sd}.log').open('w') as log:process=subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        runs.append({'requested_modes':m,'requested_sd_modes':sd,'exit_code':process.returncode,'runtime_s':time.time()-tic})
        candidates=list((runout/'states').glob('*.npz'))
        if process.returncode==0:states.extend(candidates)
        else:break
    report=compare(states);report['runs']=runs
    (out/'convergence.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));return 0 if report['status']=='PASS' else 2

def export(args):
    import numpy as np
    from export_fields import export_state
    state=np.load(args.state);out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    table=np.column_stack([state['xyz_nm'],state['phi_v'],state['electron_density_cm3']])
    np.savetxt(out/'density_potential.csv',table,delimiter=',',header='x_nm,y_nm,z_nm,phi_v,electron_density_cm3',comments='')
    export_state(Path(args.state),out/'fields.vtu')
    if args.transmission:
        from config import Config
        from modespace.hamiltonian import Grid,ModeBuilder
        from modespace.junction_basis import JunctionCompatibleBuilder
        from device.definition import DeviceDefinition
        from transport.registry import create_transport_backend
        from transport.ms_postprocess import transmission
        c=Config(**json.loads(str(state['config_json'])));g=Grid.build(c);cls=JunctionCompatibleBuilder if c.junction_basis=='nested' else ModeBuilder
        builders=[cls(c,g,v) for v in DeviceDefinition(transport_axis=c.transport_axis).frame.valleys()];backend=create_transport_backend(c.backend,c);backend.calibrate_contacts({'builders':builders,'grid':g,'config':c})
        mu=float(state['mu_l_ev']);vd=float(state['vd_v']);energy=np.linspace(mu-vd-.15,mu+.4,201);data={}
        for builder in builders:
            builder.block_only=True;s=builder.build(state['phi_v']);data[s.valley.name]=transmission(s,c,energy,vd)
        np.savez_compressed(out/'transmission.npz',energy_ev=energy,**data)
        np.savetxt(out/'transmission.csv',np.column_stack([energy,*data.values()]),delimiter=',',header=','.join(['energy_ev',*data]),comments='')
    return 0

def main():
    bootstrap();p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    d=sub.add_parser('doctor');d.add_argument('--out')
    for command in ('run','sweep','resume'):
        sub.add_parser(command,add_help=False).add_argument('args',nargs=argparse.REMAINDER)
    c=sub.add_parser('compare');c.add_argument('states',nargs='+');c.add_argument('--out',required=True)
    r=sub.add_parser('refine');r.add_argument('--config',required=True);r.add_argument('--out',required=True);r.add_argument('--pairs',default='8/32,12/48,16/64');r.add_argument('--vg',type=float,default=.4);r.add_argument('--vd',type=float,default=.01)
    e=sub.add_parser('export');e.add_argument('--state',required=True);e.add_argument('--out',required=True);e.add_argument('--transmission',action='store_true')
    m=sub.add_parser('select-modes');m.add_argument('--config',required=True);m.add_argument('--out',required=True)
    s=sub.add_parser('certify');s.add_argument('--root',default=str(ROOT))
    a=p.parse_args()
    if a.command=='doctor':return doctor(a.out)
    if a.command in ('run','sweep','resume'):
        extra=a.args[1:] if a.args and a.args[0]=='--' else a.args
        if a.command=='resume' and '--initial-state' not in extra and '--resume' not in extra:p.error('resume requires --initial-state or --resume')
        return subprocess.call([sys.executable,str(ROOT/'run.py'),*extra],cwd=ROOT)
    if a.command=='compare':report=compare(a.states);Path(a.out).write_text(json.dumps(report,indent=2));return 0 if report['status']=='PASS' else 2
    if a.command=='refine':return refine(a)
    if a.command=='export':return export(a)
    if a.command=='select-modes':
        from validation.mode_selection import select
        print(json.dumps(select(a.config,a.out),indent=2));return 0
    if a.command=='certify':
        from validation.release_certificate import certify
        result=certify(Path(a.root));print(json.dumps(result,indent=2));return 0 if result['RELEASE_READY'] else 2

if __name__=='__main__':raise SystemExit(main())
