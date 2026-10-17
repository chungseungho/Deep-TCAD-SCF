"""Full real-space operator, closed density and independent open inversion."""
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import replace
import numpy as np
from scipy.linalg import eigh
from config import Config,VALLEYS,KIN,Q,HBAR_EV_S
from modespace.hamiltonian import Grid,ModeBuilder
from modespace.junction_basis import JunctionCompatibleBuilder
from multispace.occupation import fermi
from transport.negf_reference import device_blocks,lead_matrices,solve_negf,blocks_density
from transport.ms_postprocess import transmission
from validation.reduced_realspace import full_operator
from validation.step08_modes import density_eigen
from solver import write_json,atomic_npz
from provenance import source_digest


def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[];mu=.3;v=VALLEYS[0]
    for dx in (1.,.5,.25):
        for length in (2.,4.,8.):
            c=Config(ls=length,ld=length,lg=2,spacer=1,width=2,thickness=1,sd_width=4,sd_thickness=3,
                tox=.5,dx=dx,dy=1,dz=1,modes=2,sd_modes=12,ms_boundary_shift_slices=0)
            grid=Grid.build(c);phi=np.zeros(len(grid.xyz));H=full_operator(grid,c,v,phi)
            er,ur=eigh(H);truth=(abs(ur)**2@fermi(er,mu,300))*4/grid.cell_volume*1e21
            full=ModeBuilder(c,grid,v).build(phi);U=full.realspace_transform()
            Es=np.linspace(.08,1.1,101);ons,hop,tx,leads=device_blocks(full,c,(0.,0.))
            # Independent real-space open operator from the stencil, including
            # explicit half-cell terminal correction; only analytic lead g reused.
            n0=grid.mask[0].sum();nn=grid.mask[-1].sum();Ho=H.copy()
            Ho[:n0,:n0]-=tx*np.eye(n0);Ho[-nn:,-nn:]-=tx*np.eye(nn)
            UL=full.bases[0][grid.mask[0]];UR=full.bases[-1][grid.mask[-1]];Td=[]
            for E in Es:
                sl,sr=lead_matrices(np.array([complex(E)]),tx,leads)
                SL=UL@sl[0]@UL.T;SR=UR@sr[0]@UR.T;A=E*np.eye(len(H),dtype=complex)-Ho
                A[:n0,:n0]-=SL;A[-nn:,-nn:]-=SR
                rhs=np.zeros((len(H),nn));rhs[-nn:]=np.eye(nn)
                G=np.linalg.solve(A,rhs)[:n0]
                Td.append(np.trace((1j*(SL-SL.conj().T))@G@(1j*(SR-SR.conj().T))@G.conj().T).real)
            Td=np.array(Td);fullblocks,_,_,_=solve_negf(full,c,mu,mu,0)
            openfull=blocks_density(full,fullblocks,4)
            # Occupation-only difference on the SAME equilibrium H, declared
            # frozen transport diagnostic, not electrostatically biased device.
            window=fermi(Es,mu,300)-fermi(Es,mu-.05,300)
            Id=np.trapezoid(Td*window,Es)*Q/(2*np.pi*HBAR_EV_S)*4
            for name,m,sd,junction in [('low_plain',1,2,'plain'),('low_nested',1,2,'nested'),
                    ('M2',2,4,'nested'),('M3',2,8,'nested'),('full',2,12,'plain')]:
                cc=replace(c,modes=m,sd_modes=sd,junction_basis=junction)
                builder=JunctionCompatibleBuilder if junction=='nested' else ModeBuilder
                s=builder(cc,grid,v).build(phi);V=s.realspace_transform();e,u=eigh(s.H)
                density=density_eigen(s,u,fermi(e,mu,300));t=transmission(s,cc,Es,0)
                I=np.trapezoid(t*window,Es)*Q/(2*np.pi*HBAR_EV_S)*4
                blocks,_,_,_=solve_negf(s,cc,mu,mu,0);op=blocks_density(s,blocks,4)
                ch=(grid.xyz[:,0]>=c.gate_start)&(grid.xyz[:,0]<=c.gate_end)
                row={'dx_nm':dx,'electrode_length_nm':length,'basis':name,'dimension':s.dimension,'full_dimension':len(H),
                    'projection_error_ev':float(max(abs(s.H-V.T@H@V).ravel())),
                    'orthogonality_error':float(max(abs(V.T@V-np.eye(s.dimension)).ravel())),
                    'hermiticity_error':float(max(abs(s.H-s.H.T).ravel())),
                    'closed_density_relative_to_full':float(np.linalg.norm(density-truth)/np.linalg.norm(truth)),
                    'open_density_relative_to_full':float(np.linalg.norm(op-openfull)/np.linalg.norm(openfull)),
                    'closed_open_channel_relative_l2':float(np.linalg.norm(density[ch]-op[ch])/np.linalg.norm(op[ch])),
                    'transmission_max_absolute_error':float(max(abs(t-Td))),
                    'current_same_H_occupation_split_a':float(I),'full_current_same_H_occupation_split_a':float(Id),
                    'current_relative_error_to_full':float(abs(I-Id)/max(abs(Id),1e-30)),
                    'lowest_eigenvalue_error_ev':float(e[0]-er[0])}
                rows.append(row)
                atomic_npz(out/f'dx{dx:g}_L{length:g}_{name}.npz',xyz_nm=grid.xyz,closed_density=density,
                    full_closed_density=truth,open_density=op,energy_ev=Es,transmission=t,full_transmission=Td)
            write_json(out/'summary.json',{'source_sha256':source_digest(),'cases':rows,'complete':False,
                'scope':'reduced abrupt wire frozen fields; occupation-split current diagnostic; no GAA grid plateau'})
            print('constriction',dx,length,flush=True)
    fullrows=[r for r in rows if r['basis']=='full']
    write_json(out/'summary.json',{'source_sha256':source_digest(),'cases':rows,'complete':True,
        'full_basis_identity_pass':all(r['projection_error_ev']<1e-10 and r['closed_density_relative_to_full']<1e-10 and r['transmission_max_absolute_error']<1e-9 for r in fullrows),
        'scope':'reduced abrupt wire frozen fields; occupation-split current diagnostic; no GAA grid plateau'})

if __name__=='__main__':run(sys.argv[1] if len(sys.argv)>1 else 'results/constriction')
