#!/usr/bin/env python3
"""Standalone 3D GAA driver; arbitrary 0–6 V gate/drain ranges and failure logs."""
import argparse,csv,json,sys,time,traceback
from pathlib import Path
import numpy as np
from config import Config
from provenance import source_digest,fingerprint


def voltage_range(value):
    if ':' not in value:return [float(x) for x in value.split(',')]
    lo,hi,step=map(float,value.split(':'))
    if step<=0 or hi<lo:raise ValueError('Use start:stop:positive_step')
    count=int(np.floor((hi-lo)/step+1e-9));values=lo+np.arange(count+1)*step
    if abs(values[-1]-hi)<1e-8:values[-1]=hi
    return values.tolist()


def parser():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('phase',choices=['point','idvgs','idvds','stress','frozen'])
    a.add_argument('--config');a.add_argument('--out',required=True)
    a.add_argument('--backend',choices=['ms','negf','spectral'])
    a.add_argument('--vg',type=float,default=.8);a.add_argument('--vd',type=float,default=.1)
    a.add_argument('--gate-range',default='0:1.2:.05');a.add_argument('--drain-range',default='0:1.2:.05')
    a.add_argument('--gates',default='.4,.8,1.2');a.add_argument('--axis',choices=['gate','drain'],default='drain')
    a.add_argument('--minimum-bias-step',type=float,default=.00625)
    a.add_argument('--initial-state');a.add_argument('--state');a.add_argument('--resume',action='store_true')
    a.add_argument('--dry-run',action='store_true');a.add_argument('--max-iterations',type=int)
    return a


def configuration(args):
    values=json.loads(Path(args.config).read_text()) if args.config else {}
    if args.backend:values['backend']=args.backend
    if args.max_iterations:values['max_iterations']=args.max_iterations
    return Config(**values).validate()


def prepare(c,out):
    from geometry.build import build_mesh
    from geometry.validate import validate_mesh
    from poisson.device import Poisson
    from coupling.mapping import Mapping
    from modespace.hamiltonian import Grid,ModeBuilder
    from modespace.junction_basis import JunctionCompatibleBuilder
    from device.definition import DeviceDefinition
    from transport.negf import neutral_lead_mu
    from solver import write_json
    mesh=out/'mesh'/'gaa_cm.msh'
    if not mesh.exists():mesh=build_mesh(c,out/'mesh')
    write_json(out/'geometry_checks.json',validate_mesh(out/'mesh'/'gaa_nm.msh',c))
    p=Poisson(c,mesh)
    try:
        grid=Grid.build(c);mapping=Mapping.build(out/'mesh'/'gaa_nm.msh',p.xyz['silicon'],p.volume,grid,c.mapping)
        write_json(out/'mapping.json',{**p.diagnostics(),**mapping.diagnostics()})
        device=DeviceDefinition(transport_axis=c.transport_axis)
        cls=JunctionCompatibleBuilder if c.junction_basis=='nested' else ModeBuilder
        builders=[cls(c,grid,v) for v in device.frame.valleys()]
        if c.mu_source_ev is not None:mu=c.mu_source_ev;cal={'calibration':'explicit fixed source μ','mu_source_ev':mu}
        elif c.backend=='spectral':
            from transport.reservoir import neutral_mu
            mu,cal=neutral_mu(builders,c,grid)
        else:mu,cal=neutral_lead_mu(builders,c,grid)
        write_json(out/'contact_calibration.json',cal);write_json(out/'device.json',device.describe())
        return p,grid,mapping,builders,mu
    except Exception:p.close();raise


def write_rows(path,rows):
    if not rows:return
    columns=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,columns);w.writeheader();w.writerows(rows)


def main():
    args=parser().parse_args();c=configuration(args);out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    if args.phase=='frozen':
        from validation.negf_comparison import frozen
        if not args.state:raise ValueError('--state is required')
        frozen(Path(args.state),out);return 0
    gates=voltage_range(args.gate_range);drains=voltage_range(args.drain_range)
    if args.phase=='point':points=[(args.vg,args.vd)]
    elif args.phase=='idvgs':points=[(v,args.vd) for v in gates]
    elif args.phase=='idvds':points=[(g,v) for g in voltage_range(args.gates) for v in drains]
    else:points=[(v,args.vd) if args.axis=='gate' else (args.vg,v) for v in (gates if args.axis=='gate' else drains)]
    if any(min(g,d)<0 or max(g,d)>6 for g,d in points):raise ValueError('Supported gate/drain requested range is 0–6 V')
    print(json.dumps({'requested_points':len(points),'configuration':c.to_dict(),'source_sha256':source_digest()}),flush=True)
    if args.dry_run:return 0
    from solver import scf,write_json,atomic_npz
    settings=out/'settings.json';identity=fingerprint(c);sha=source_digest()
    if settings.exists():
        old=json.loads(settings.read_text())
        if old['configuration_fingerprint']!=identity or old['source_sha256']!=sha:
            raise ValueError('Output directory configuration/source mismatch; choose a fresh directory')
    write_json(settings,{'config':c.to_dict(),'configuration_fingerprint':identity,'source_sha256':sha,
                        'command':sys.argv,'python':sys.version})
    rows=[];initial=None;previous=None
    p,grid,mapping,builders,mu=prepare(c,out)
    if args.initial_state:
        saved=np.load(args.initial_state)
        if saved['node_xyz_nm'].shape!=p.xyz['silicon'].shape or not np.allclose(saved['node_xyz_nm'],p.xyz['silicon'],atol=1e-9,rtol=0):
            p.close();raise ValueError('Seed has different Poisson nodes')
        initial=saved['poisson_density_cm3']
    def attempt(g,d,seed,kind='requested'):
        key=f'vg{g:.6f}_vd{d:.6f}'
        meta=out/'states'/f'{key}.json'
        if args.resume and meta.exists():
            row=json.loads(meta.read_text())
            if row['numerical_status']=='SCF_CONVERGED':
                s=np.load(meta.with_suffix('.npz'));return {k:v for k,v in row.items() if not isinstance(v,(dict,list))},s['poisson_density_cm3'],None
        try:
            row,n,q=scf(c,p,mapping,builders,grid,mu,g,d,out,key,seed)
            row['point_kind']=kind;return row,n,q
        except Exception as exc:
            failure=traceback.format_exc();(out/f'{key}_failure.txt').write_text(failure)
            if hasattr(exc,'state'):
                state=exc.state
                atomic_npz(out/'failures'/f'{key}_occupation.npz',**{k:v for k,v in state.items() if isinstance(v,np.ndarray)})
            if hasattr(exc,'bound_diagnostics'):
                write_json(out/'failures'/f'{key}_bound_poles.json',exc.bound_diagnostics)
            row={'key':key,'backend':c.backend,'vg_v':g,'vd_v':d,'id_a':float('nan'),
                'numerical_status':'MS_OCCUPATION_UNRESOLVED' if hasattr(exc,'state') else 'BOUND_PREPARATION_REQUIRED' if hasattr(exc,'bound_diagnostics') else 'ERROR',
                'physical_status':'NOT_EVALUATED','transport_status':'NOT_RUN','error':repr(exc),
                'source_sha256':sha,'configuration_fingerprint':identity,'point_kind':kind}
            # Report actual fields of the failed iterate, clearly distinguished
            # from a converged physical prediction. This matters for stress runs.
            try:
                from analysis.validity import element_fields
                fields=element_fields(p.mesh_path,p)
                high=any(f['maximum_v_nm']>c.high_field_v_nm for f in fields.values())
                row['physical_status']='MODEL_OUT_OF_SCOPE' if high else 'NOT_EVALUATED'
                row['maximum_failed_iterate_field_v_nm']=max(f['maximum_v_nm'] for f in fields.values())
                write_json(out/'failures'/f'{key}_field.json',{'field_state':'LAST_AVAILABLE_POISSON_FIELD_NOT_SCF_SOLUTION_OR_UNEVALUATED_NEWTON_TRIAL',
                    'fields':fields,'high_field_warning':high,'threshold_v_nm':c.high_field_v_nm})
            except Exception:pass
            write_json(out/'failures'/f'{key}.json',row);print(failure,flush=True)
            return row,None,None
    try:
        for g,d in points:
            if previous is not None and abs(previous[0]-g)>1e-12 and args.phase=='idvds':initial=None;previous=None
            row,n,q=attempt(g,d,initial)
            # On numerical SCF failure, bisect continuation from last accepted point.
            # An unresolved occupation is missing model evidence; do not bypass it.
            if previous is not None and row['numerical_status'] in ('ERROR','SCF_NOT_CONVERGED') and 'integration' not in str(row.get('error','')).lower():
                old=np.asarray(previous);new=np.array([g,d]);distance=float(max(abs(new-old)))
                intervals=2
                while distance/intervals>=args.minimum_bias_step and row['numerical_status']!='SCF_CONVERGED':
                    seed=initial;ok=True
                    for f in np.linspace(0,1,intervals+1)[1:]:
                        bias=old+(new-old)*f;rr,nn,qq=attempt(*bias,seed,'adaptive')
                        rows.append(rr)
                        if rr['numerical_status']!='SCF_CONVERGED':ok=False;break
                        seed=nn
                    if ok:row,n=rr,nn
                    intervals*=2
            rows.append(row);write_rows(out/(args.phase+'.csv'),rows)
            if row['numerical_status']=='SCF_CONVERGED':initial=n;previous=(g,d)
        return 0 if all(r['numerical_status']=='SCF_CONVERGED' for r in rows) else 2
    finally:p.close()


if __name__=='__main__':raise SystemExit(main())
