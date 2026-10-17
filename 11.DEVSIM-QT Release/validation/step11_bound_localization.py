"""Energy-dependent pole residue and localization in a frozen projected model.

Classification as a real device state requires grid and contact convergence;
this script does not replace that requirement or invent biased occupations.
"""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.sparse import bmat,csc_matrix,eye
from scipy.sparse.linalg import eigsh
from config import Config
from modespace.hamiltonian import Grid
from modespace.junction_basis import JunctionCompatibleBuilder
from device.definition import DeviceDefinition
from transport.contact_model import solve_reference
from transport.negf_reference import device_blocks,lead_matrices
from transport.bound_states import audit

def run(state,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);data=np.load(state);c=Config(**json.loads(str(data['config_json'])));g=Grid.build(c);valleys=DeviceDefinition().frame.valleys();ref=solve_reference(c,g,valleys);rows=[]
    for v in valleys:
        b=JunctionCompatibleBuilder(c,g,v);b.block_only=True;b.contacts=ref;s=b.build(data['phi_v']);mu=ref.mu-c.potential_reference_v
        aud=audit(s,c,mu,0);ons,hop,t,lead=device_blocks(s,c,(c.potential_reference_v,)*2);n=len(ons)
        for root in aud['bound_poles_below_both_leads']:
            energy=root['energy_ev'];sl,sr=lead_matrices(np.array([energy+0j]),t,lead)
            blocks=[[None]*n for _ in range(n)]
            for i,h in enumerate(ons):blocks[i][i]=csc_matrix(h+(sl[0].real if i==0 else 0)+(sr[0].real if i==n-1 else 0))
            for i,h in enumerate(hop):blocks[i][i+1]=csc_matrix(h);blocks[i+1][i]=csc_matrix(h.T)
            H=bmat(blocks,format='csc');e,U=eigsh(H,k=1,sigma=energy+1e-9,which='LM',tol=1e-11);u=U[:,0]
            delta=1e-6;sp=lead_matrices(np.array([energy+delta+0j]),t,lead);sm=lead_matrices(np.array([energy-delta+0j]),t,lead)
            derivative=[(p[0].real-m[0].real)/(2*delta) for p,m in zip(sp,sm)]
            norm=1-u[s.sl(0)]@derivative[0]@u[s.sl(0)]-u[s.sl(n-1)]@derivative[1]@u[s.sl(n-1)]
            weight=np.array([np.sum(abs(u[s.sl(i)])**2)/norm for i in range(n)])
            gamma=[1j*(x[0]-x[0].conj().T) for x in (sl,sr)]
            row={'valley':v.name,'energy_ev':energy,'eigen_residual_ev':float(np.linalg.norm(H@u-energy*u)),
                'projected_device_spectral_weight':float(1/norm),'lead_tail_weight':float(1-1/norm),
                'gamma_coupling_ev':float(abs(u[s.sl(0)].conj()@gamma[0]@u[s.sl(0)]+u[s.sl(n-1)].conj()@gamma[1]@u[s.sl(n-1)])),
                'source_weight':float(weight[g.x<c.ls].sum()),'central_weight':float(weight[(g.x>=c.ls)&(g.x<=c.drain_start)].sum()),'drain_weight':float(weight[g.x>c.drain_start].sum()),
                'classification':'projected isolated pole; physical vs contact/grid artifact INCONCLUSIVE'}
            rows.append(row);np.savez_compressed(out/f'{v.name}_{len(rows)}.npz',x_nm=g.x,spectral_weight_per_slice=weight,energy_ev=energy)
    (out/'summary.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
if __name__=='__main__':run(*sys.argv[1:])
