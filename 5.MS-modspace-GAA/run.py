#!/usr/bin/env python3
"""STEP05 command line; refuse stale meshes, conflicting resumes and fake I-V points."""
import argparse,csv,hashlib,json,os,sys,traceback
import importlib.metadata
from pathlib import Path
from dataclasses import fields
import numpy as np
from config import Config,VALLEYS
from geometry.build import build_mesh
from modespace.hamiltonian import Grid,ModeBuilder
from transport.reservoir import neutral_mu
from solver import run_scf


def voltages(spec):
    if ':' not in spec: return [float(v) for v in spec.split(',')]
    lo,hi,step=map(float,spec.split(':'))
    if step<=0 or hi<lo: raise ValueError('Voltage range must be lo:hi:positive_step')
    return [float(v) for v in np.round(np.arange(lo,hi+step*1e-6,step),10)]


def parser():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--phase',choices=['mesh','poisson','point','idvgs','idvds','plot'],default='point')
    a.add_argument('--out',default='output/run');a.add_argument('--config',help='JSON Config overrides')
    a.add_argument('--vg',type=float,default=.4);a.add_argument('--vd',type=float,default=.1)
    a.add_argument('--vgs',default='0:1.2:0.05');a.add_argument('--vds',default='0:1.2:0.05')
    a.add_argument('--gates',default='0.4,0.8,1.2');a.add_argument('--max-points',type=int)
    a.add_argument('--resume',action='store_true');a.add_argument('--keep-going',action='store_true')
    a.add_argument('--spectra',action='store_true',help='Expensive finite-reservoir POST-SCF spectral functions')
    a.add_argument('--warm-start-model',help='Optional checked surrogate initial guess; never skips reference SCF')
    a.add_argument('--smoke',action='store_true',help='Explicitly coarse numerical smoke preset, NOT production convergence')
    a.add_argument('--backend',choices=['spectral','lindblad','ms'])
    for f in fields(Config):
        if f.name=='backend': continue
        a.add_argument('--'+f.name.replace('_','-'),type=type(f.default),default=None)
    return a


def configuration(args):
    d=Config().to_dict()
    if args.smoke: d.update(dx=3.,modes=2,sd_modes=4,mesh_h=3.5,channel_h=1.8,junction_h=1.3)
    if args.config:
        overrides=json.loads(Path(args.config).read_text())
        unknown=set(overrides)-set(d)
        if unknown: raise ValueError(f'Unknown Config keys: {unknown}')
        d.update(overrides)
    for k in d:
        v=getattr(args,k,None)
        if v is not None: d[k]=v
    return Config(**d).validate()


def fingerprint(c):
    return hashlib.sha256(json.dumps(c.to_dict(),sort_keys=True).encode()).hexdigest()


def source_digest():
    root=Path(__file__).resolve().parent;h=hashlib.sha256()
    for f in sorted(root.rglob('*.py')):
        if 'output' in f.relative_to(root).parts or 'tests' in f.relative_to(root).parts: continue
        h.update(str(f.relative_to(root)).encode());h.update(f.read_bytes())
    return h.hexdigest()


def save_settings(out,c,args):
    f=out/'settings.json';record={'config':c.to_dict(),'fingerprint':fingerprint(c),
                               'source_sha256':source_digest(),
                               'versions':{name:importlib.metadata.version(name) for name in ['numpy','scipy','matplotlib','gmsh','devsim']},
                               'source_commit':os.environ.get('GIT_COMMIT','see repository git log'),
                               'preset':'coarse_smoke' if args.smoke else 'configured',
                               'argv':sys.argv}
    if f.exists():
        old=json.loads(f.read_text())
        if old['fingerprint']!=record['fingerprint']: raise ValueError('Configuration changed; choose a NEW output folder')
        if args.phase!='plot' and old.get('source_sha256',record['source_sha256'])!=record['source_sha256']:
            raise ValueError('Simulation source changed; choose a NEW output folder to keep results reproducible')
        if not args.resume and args.phase not in ('plot','mesh'): raise ValueError('Output exists: use --resume or a new directory')
    else: f.write_text(json.dumps(record,indent=2))


def write_points(path,rows):
    if not rows: return
    temp=path.with_suffix('.tmp')
    with temp.open('w',newline='') as f:
        keys=list(dict.fromkeys(k for r in rows for k in r))
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    os.replace(temp,path)


def main():
    args=parser().parse_args();c=configuration(args);out=Path(args.out).resolve()
    out.mkdir(parents=True,exist_ok=True);save_settings(out,c,args)
    if args.phase=='plot':
        from plotting import plot_run
        plot_run(out)
        # A plotting/export interruption must not leave an old partial summary
        # after all bias snapshots have already been written successfully.
        phases={}
        for phase in ['point','idvgs','idvds']:
            table=out/(phase+'.csv')
            if table.exists():
                saved=list(csv.DictReader(table.open()))
                phases[phase]={'finished_points':len(saved),
                    'converged_points':sum(r['status']=='SCF_CONVERGED_PROTOTYPE' for r in saved),
                    'failed_points':sum(r['status']!='SCF_CONVERGED_PROTOTYPE' for r in saved)}
        report={'summary_kind':'aggregate_from_saved_bias_tables','phase_summaries':phases,
                'quantitative_device_validation':False}
        if len(phases)==1: report.update(next(iter(phases.values())))
        (out/'run_summary.json').write_text(json.dumps(report,indent=2))
        return 0
    meshdir=out/'mesh';mesh=meshdir/'gaa_cm.msh'
    if not mesh.exists(): mesh=build_mesh(c,meshdir)
    if args.phase=='mesh':
        from plotting import plot_geometry
        plot_geometry(c,meshdir,out/'figures');return 0
    from poisson.device import Poisson
    from coupling.mapping import Mapping
    p=Poisson(c,mesh)
    try:
        grid=Grid.build(c);mapping=Mapping.build(meshdir/'gaa_nm.msh',p.xyz['silicon'],p.volume,grid)
        report={**p.diagnostics(),**mapping.diagnostics()}
        (out/'mapping.json').write_text(json.dumps(report,indent=2))
        print('Mapping:',json.dumps(report),flush=True)
        if args.phase=='poisson':
            # Neutral doping cancellation + all contacts equal is an EXACT constant solution.
            phi=p.solve(c.flatband,0.)
            zero=float(max(abs(phi)));report['neutral_laplace_error_v']=zero
            phi_low=p.solve(0.,args.vd);phi_high=p.solve(1.2,args.vd)
            report['gate_response_v']=float(np.max(phi_high-phi_low))
            (out/'poisson_checks.json').write_text(json.dumps(report,indent=2))
            p.export(out/'poisson.vtk')
            if zero>1e-8: raise RuntimeError('Neutral zero-field Poisson check failed')
            print('POISSON_CHECKS_PASS',report);return 0
        builders=[ModeBuilder(c,grid,v) for v in VALLEYS]
        mu,calibration=neutral_mu(builders,c,grid)
        (out/'contact_calibration.json').write_text(json.dumps(calibration,indent=2))
        print('Contact neutrality:',calibration,flush=True)
        if args.phase=='point': biases=[(args.vg,args.vd)]
        elif args.phase=='idvgs': biases=[(v,args.vd) for v in voltages(args.vgs)]
        else: biases=[(g,d) for g in voltages(args.gates) for d in voltages(args.vds)]
        if args.max_points: biases=biases[:args.max_points]
        path=out/(args.phase+'.csv');rows=list(csv.DictReader(path.open())) if path.exists() else []
        initial=None;failures=0;previous_vg=None
        for vg,vd in biases:
            if args.phase=='idvds' and previous_vg is not None and vg!=previous_vg:
                initial=None
            previous_vg=vg
            key=f'vg{vg:.6f}_vd{vd:.6f}';prior=next((r for r in rows if r['key']==key),None)
            state=out/'states'/f'{key}.npz'
            if args.resume and prior and prior['status']=='SCF_CONVERGED_PROTOTYPE' and state.exists():
                initial=np.load(state)['poisson_density_cm3'];print('RESUME skip',key);continue
            try:
                if args.warm_start_model and initial is None:
                    from surrogate import predict
                    channel=(grid.xyz[:,0]>=c.ls)&(grid.xyz[:,0]<=c.drain_start)
                    try:
                        prediction=predict(args.warm_start_model,vg,vd,fingerprint(c),grid.xyz[channel])
                        seed=c.doping(grid.xyz);seed[channel]=prediction;initial=mapping.charge(seed)
                    except ValueError as exc:
                        print('Surrogate rejected; using full reference initialisation:',exc,flush=True)
                # A completed FAILED bias has no accepted checkpoint. Prefer
                # the newly restored neighbouring converged bias for continuation.
                # Interrupted biases without a completed row still resume in-place.
                resume_point=args.resume and not(prior and prior['status']!='SCF_CONVERGED_PROTOTYPE' and initial is not None)
                result,next_density,final=run_scf(c,p,mapping,builders,grid,mu,vg,vd,out,key,initial,resume_point)
                rows=[r for r in rows if r['key']!=key]+[result];write_points(path,rows)
                ok=result['status']=='SCF_CONVERGED_PROTOTYPE'
                if ok: initial=next_density
                else: failures+=1
                if args.spectra and ok:
                    from transport.spectral import spectrum
                    E=np.arange(min(mu,mu-vd)-.5,max(mu,mu-vd)+.8,.002)
                    totalT=np.zeros(len(E));totalJ=totalT.copy()
                    for system,rho,baths in final['matrices']:
                        t,j=spectrum(system.H,baths,E,2*system.valley.multiplicity,c.backend,mu,mu-vd,c.temperature)
                        totalT+=2*system.valley.multiplicity*t;totalJ+=j
                    np.savez_compressed(out/'states'/f'{key}_spectra.npz',energy_ev=E,transmission=totalT,
                                        spectral_current_a_per_ev=totalJ,integrated_current_a=np.trapezoid(totalJ,E))
                    integrated=float(np.trapezoid(totalJ,E))
                    (out/'states'/f'{key}_spectra.json').write_text(json.dumps({
                        'postprocessing_only':True,'energy_step_ev':.002,'integrated_current_a':integrated,
                        'bond_current_a':result['id_a'],'relative_mismatch':abs(integrated-result['id_a'])/max(abs(result['id_a']),1e-18),
                        'note':'sampled spectrum is not certified; refine energy spacing/window for narrow resonances'},indent=2))
                if not ok and not args.keep_going: break
            except Exception:
                failures+=1;(out/f'{key}_failure.txt').write_text(traceback.format_exc())
                print(traceback.format_exc(),file=sys.stderr)
                if not args.keep_going: break
        from plotting import plot_run,plot_geometry
        plot_geometry(c,meshdir,out/'figures');plot_run(out)
        p.export(out/'electrostatics.vtk')
        summary={'requested_points':len(biases),'finished_points':len(rows),
                 'converged_points':sum(r['status']=='SCF_CONVERGED_PROTOTYPE' for r in rows),
                 'failures_this_run':failures,'quantitative_device_validation':False,
                 'reason':'mesh/modes/reservoir-length/gamma convergence must be assessed independently'}
        (out/'run_summary.json').write_text(json.dumps(summary,indent=2))
        if args.phase=='idvgs':
            from analysis import transfer_metrics
            (out/'transfer_metrics.json').write_text(json.dumps(transfer_metrics(rows,c),indent=2))
        print(json.dumps(summary,indent=2));return 0 if failures==0 else 2
    finally: p.close()


if __name__=='__main__':
    try: raise SystemExit(main())
    except Exception: traceback.print_exc();raise SystemExit(2)
