#!/usr/bin/env python3
"""Executed research alternative: nested wide/narrow transverse subspaces.

The low wide modes are ENRICHED with the zero-extended adjacent narrow modes.
No hopping rescaling/fitting is used: every onsite, shoulder wall and overlap
is reprojected from the same full 3D FD operator. This keeps a common trace
space at an abrupt section change. Extra retained S/D ranks are reported.
The extension has its own SHA/config identity, separate from baseline runs.
"""
import hashlib,json,sys,time,traceback
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.linalg import orth,eigh
from config import KIN,VALLEYS
from modespace.hamiltonian import ModeBuilder,System,Grid,transverse_matrix


class JunctionCompatibleBuilder(ModeBuilder):
    def build(self,phi):
        original=super().build(phi);c,g,v=self.c,self.g,self.v
        n=len(g.x);channel=np.where((g.x>=c.ls)&(g.x<=c.drain_start))[0]
        targets={0:original.bases[channel[0]],2:original.bases[channel[-1]]}
        bases=[];onsite=[];subbands=[];regions=[];tx=KIN/(v.mx*c.dx**2)
        for i,x in enumerate(g.x):
            mask=g.mask[i];region=0 if x<c.ls else 2 if x>c.drain_start else 1
            U=original.bases[i].reshape(-1,original.bases[i].shape[-1])
            if region!=1:
                narrow=targets[region].reshape(-1,targets[region].shape[-1])
                U=orth(np.column_stack([U,narrow]),rcond=1e-11)
            u=U[mask.ravel()];m=u.shape[1]
            transverse=transverse_matrix(mask,KIN/(v.my*c.dy**2),KIN/(v.mz*c.dz**2))
            ht=u.T@transverse@u-u.T@(g.slice_values(phi,i)[:,None]*u)
            left=g.mask[i-1] if i else np.zeros_like(mask);right=g.mask[i+1] if i+1<n else np.zeros_like(mask)
            wall=(~left).astype(float)+(~right).astype(float)
            h=ht+2*tx*np.eye(m)+tx*(u.T@(wall[mask,None]*u))
            bases.append(U.reshape(*mask.shape,m));onsite.append(h)
            subbands.append(eigh(ht,eigvals_only=True,check_finite=False));regions.extend([region]*m)
        off=np.r_[0,np.cumsum([len(h) for h in onsite])];H=np.zeros((off[-1],off[-1]));hop=[]
        for i,h in enumerate(onsite):H[off[i]:off[i+1],off[i]:off[i+1]]=h
        for i in range(n-1):
            a=bases[i].reshape(-1,bases[i].shape[-1]);b=bases[i+1].reshape(-1,bases[i+1].shape[-1]);t=-tx*a.T@b;hop.append(t)
            H[off[i]:off[i+1],off[i+1]:off[i+2]]=t;H[off[i+1]:off[i+2],off[i]:off[i+1]]=t.T
        return System(H,bases,subbands,off,np.array(regions),g,hop,v)


def main():
    from run import parser,configuration,save_settings,write_points
    from provenance import fingerprint,source_digest
    from geometry.build import build_mesh
    from geometry.validate import validate_mesh
    from poisson.device import Poisson
    from coupling.mapping import Mapping
    from transport.reservoir import neutral_mu
    from solver import run_scf,atomic_npz
    from resources import estimate
    args=parser().parse_args()
    if args.phase!='point' or args.resume:raise ValueError('Research extension supports fresh single-point runs only')
    c=configuration(args);out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    extension=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    identity=hashlib.sha256((fingerprint(c)+extension).encode()).hexdigest()
    plan={**estimate(c),'extra_sd_rank_up_to':c.modes,'note':'base resource estimate excludes extra S/D enrichment; dimension and measured runtime are reported per valley'}
    print(json.dumps(plan),flush=True)
    if args.dry_run:return 0
    save_settings(out,c,args)
    settings=json.loads((out/'settings.json').read_text());settings.update(extension_source_sha256=extension,extension_configuration_fingerprint=identity,basis_algorithm='nested wide/narrow subspaces; SD low modes plus zero-extended narrow modes')
    (out/'settings.json').write_text(json.dumps(settings,indent=2))
    meshdir=out/'mesh';mesh=meshdir/'gaa_cm.msh'
    if not mesh.exists():mesh=build_mesh(c,meshdir)
    (out/'geometry_checks.json').write_text(json.dumps(validate_mesh(meshdir/'gaa_nm.msh',c),indent=2))
    p=Poisson(c,mesh);key=f'vg{args.vg:.6f}_vd{args.vd:.6f}';tic=time.perf_counter()
    try:
        grid=Grid.build(c);mapping=Mapping.build(meshdir/'gaa_nm.msh',p.xyz['silicon'],p.volume,grid,c.mapping)
        (out/'mapping.json').write_text(json.dumps({**p.diagnostics(),**mapping.diagnostics()},indent=2))
        builders=[JunctionCompatibleBuilder(c,grid,v) for v in VALLEYS]
        mu,calibration=neutral_mu(builders,c,grid);(out/'contact_calibration.json').write_text(json.dumps(calibration,indent=2))
        initial=None
        if args.initial_state:
            s=np.load(args.initial_state)
            if s['node_xyz_nm'].shape!=p.xyz['silicon'].shape or not np.allclose(s['node_xyz_nm'],p.xyz['silicon'],atol=1e-9,rtol=0):raise ValueError('Initial Poisson nodes mismatch')
            initial=s['poisson_density_cm3']
        result,_,final=run_scf(c,p,mapping,builders,grid,mu,args.vg,args.vd,out,key,initial)
        result.update(extension_source_sha256=extension,extension_configuration_fingerprint=identity,sd_basis_enriched=True)
        write_points(out/'point.csv',[result])
        file=out/'states'/(key+'.npz');s=np.load(file);data={k:s[k] for k in s.files}
        data.update(extension_source_sha256=extension,base_configuration_fingerprint=data['configuration_fingerprint'],configuration_fingerprint=identity,sd_basis_enriched=True)
        atomic_npz(file,**data)
        j=out/'states'/(key+'.json');meta=json.loads(j.read_text());meta.update(result,base_configuration_fingerprint=meta['configuration_fingerprint'],configuration_fingerprint=identity)
        j.write_text(json.dumps(meta,indent=2));p.export(out/'electrostatics.vtk')
        from plotting import plot_run,plot_geometry
        from export_fields import export_last
        plot_geometry(c,meshdir,out/'figures');plot_run(out);export_last(out)
        print(json.dumps(result,indent=2),flush=True)
        return 0 if result['status']=='SCF_CONVERGED_PROTOTYPE' else 2
    except Exception:
        (out/(key+'_failure.txt')).write_text(traceback.format_exc());write_points(out/'point.csv',[dict(key=key,vg_v=args.vg,vd_v=args.vd,id_a=float('nan'),status='ERROR',extension_source_sha256=extension)])
        print(traceback.format_exc(),flush=True);return 2
    finally:p.close()


if __name__=='__main__':raise SystemExit(main())
