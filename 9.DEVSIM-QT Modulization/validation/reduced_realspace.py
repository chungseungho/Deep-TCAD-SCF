"""Independent full 3D stencil vs mode projection on three reduced grids."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.linalg import eigh
from config import Config,VALLEYS,KIN
from modespace.hamiltonian import Grid,ModeBuilder
from modespace.junction_basis import JunctionCompatibleBuilder
from multispace.occupation import fermi
from solver import write_json,atomic_npz
from provenance import source_digest


def full_operator(grid,c,v,phi):
    mask=grid.mask;ids=np.full(mask.shape,-1,int);ids[mask]=np.arange(mask.sum())
    H=np.diag(-phi);hops=[KIN/(mass*step**2) for mass,step in zip((v.mx,v.my,v.mz),(c.dx,c.dy,c.dz))]
    for position in zip(*np.where(mask)):
        a=ids[position]
        for axis,t in enumerate(hops):
            for direction in (-1,1):
                q=list(position);q[axis]+=direction
                inside=all(0<=x<s for x,s in zip(q,mask.shape)) and mask[tuple(q)]
                H[a,a]+=t if inside else 2*t
                if inside:H[a,ids[tuple(q)]]=-t
    return H


def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[]
    for dx in (1.,.5,.25):
        base=Config(ls=2,ld=2,lg=2,spacer=1,width=2,thickness=1,sd_width=4,sd_thickness=3,
                    tox=.5,dx=dx,dy=1,dz=1,modes=2,sd_modes=12,ms_boundary_shift_slices=0)
        grid=Grid.build(base);phi=.03*np.sin(grid.xyz[:,0]);v=VALLEYS[0]
        H=full_operator(grid,base,v,phi);e,u=eigh(H);rho=(u*fermi(e,.15,300))@u.T;truth=np.diag(rho)
        for name,m,sd,basis,junction in [('full',2,12,'fixed','plain'),('fixed_low',1,2,'fixed','plain'),
                ('nested_low',1,2,'fixed','nested'),('local_nested_low',1,2,'local','nested')]:
            from dataclasses import replace
            c=replace(base,modes=m,sd_modes=sd,basis=basis,junction_basis=junction)
            cls=JunctionCompatibleBuilder if junction=='nested' else ModeBuilder
            s=cls(c,grid,v).build(phi);U=s.realspace_transform()
            er,ur=eigh(s.H);rr=(ur*fermi(er,.15,300))@ur.T
            projected=np.diag(U@rr@U.T)
            defect=float(np.max(abs(s.H-U.T@H@U)))
            row={'dx_nm':dx,'basis':name,'dimension':s.dimension,'full_dimension':len(H),
                'operator_projection_max_error_ev':defect,
                'density_relative_l2_to_full':float(np.linalg.norm(projected-truth)/np.linalg.norm(truth)),
                'electrons_per_degeneracy':float(projected.sum()),'config':c.to_dict()}
            rows.append(row)
            atomic_npz(out/f'dx{dx:g}_{name}.npz',xyz_nm=grid.xyz,density_per_cell=projected,
                full_density_per_cell=truth,source_sha256=source_digest(),config_json=json.dumps(c.to_dict()))
    result={'source_sha256':source_digest(),'case_scope':'reduced frozen 3D operator and equilibrium density; no GAA SCF plateau',
        'cases':rows,'pass':bool(all(x['operator_projection_max_error_ev']<1e-11 for x in rows) and
           all(x['density_relative_l2_to_full']<1e-11 for x in rows if x['basis']=='full'))}
    write_json(out/'summary.json',result);print(json.dumps(result,indent=2));return result


if __name__=='__main__':run(sys.argv[1] if len(sys.argv)>1 else 'results/reduced_realspace')
