"""Sparse diagnostic of narrow/dark candidates in the actual biased Hamiltonian.

Finite sampling of nonlinear real-part pole branches is not an exhaustive dark
state audit. Width/rank diagnostics never create an occupation or certify one.
"""
import sys,json,shutil,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.sparse import bmat,csc_matrix
from scipy.sparse.linalg import eigsh
from config import Config
from run import prepare
from transport.negf_reference import device_blocks,lead_matrices
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'results/embedded_candidates';out.mkdir(exist_ok=True)
if not (out/'mesh').exists():shutil.copytree(ROOT/'results/raised_global_bias/mesh',out/'mesh')
c=Config(**json.loads((ROOT/'configs/release_raised.json').read_text()));p,g,mapping,builders,mu,backend=prepare(c,out)
seed=np.load(ROOT/'results/current_seed8.npz');p.set_charge(seed['poisson_density_cm3']);phi_node=p.solve(.4,.01);phi=mapping.potential(phi_node)
np.savez_compressed(out/'frozen_initial_trial.npz',node_xyz_nm=p.xyz['silicon'],node_phi_v=phi_node,xyz_nm=g.xyz,phi_v=phi,mu_l_ev=mu,vd_v=.01,numerically_converged=False)
rows=[]
try:
 for builder in builders:
    tic=time.time();s=builder.build(phi);ons,hop,t,lead=device_blocks(s,c,(0.,.01));n=len(ons);blocks=[[None]*n for _ in range(n)]
    for i,h in enumerate(ons):blocks[i][i]=csc_matrix(h)
    for i,h in enumerate(hop):blocks[i][i+1]=csc_matrix(h);blocks[i+1][i]=csc_matrix(h.T)
    H=bmat(blocks,format='csc');N=H.shape[0];off=s.offsets
    def embedded(energy):
        sl,sr=lead_matrices(np.array([energy+0j]),t,lead)
        rr,cc=np.indices(sl[0].shape);r2,c2=np.indices(sr[0].shape)
        B=csc_matrix((np.r_[sl[0].real.ravel(),sr[0].real.ravel()],(np.r_[rr.ravel(),(N-len(sr[0])+r2).ravel()],np.r_[cc.ravel(),(N-len(sr[0])+c2).ravel()])),shape=(N,N))
        return H+B,sl[0],sr[0]
    seeds=[]
    for energy in (mu-.02,mu+.03,mu+.07):
        e,U=eigsh(H,k=6,sigma=energy,which='LM',tol=1e-10)
        seeds.extend(zip(e,U.T))
    found=[]
    for energy,vector in seeds:
        converged=False
        try:
            for iteration in range(15):
                A,sl,sr=embedded(energy);e,U=eigsh(A,k=3,sigma=energy+1e-10,which='LM',tol=1e-10,v0=vector)
                j=int(np.argmax(abs(U.T@vector)));new=float(e[j]);vector=U[:,j]
                if abs(new-energy)<1e-9:energy=new;converged=True;break
                energy=new
            if any(abs(energy+x['energy_ev']*-1)<1e-7 for x in found):continue
            A,sl,sr=embedded(energy);vl=vector[s.sl(0)];vr=vector[s.sl(n-1)];gl=1j*(sl-sl.conj().T);gr=1j*(sr-sr.conj().T)
            width=float(np.real(vl.conj()@gl@vl+vr.conj()@gr@vr));residual=float(np.linalg.norm(A@vector-energy*vector))
            probability=np.array([np.sum(abs(vector[s.sl(i)])**2) for i in range(n)])
            row={'energy_ev':energy,'gamma_expectation_ev':width,'gamma_vector_norm_ev':float(np.hypot(np.linalg.norm(gl@vl),np.linalg.norm(gr@vr))),
                'real_part_pole_residual_ev':residual,'root_iteration_converged':converged,
                'below_global_lead_edge':bool(energy<min(a.min()-2*t for a,w in lead)),
                'source_probability':float(probability[g.x<c.ls].sum()),'channel_extension_probability':float(probability[(g.x>=c.ls)&(g.x<=c.drain_start)].sum()),'drain_probability':float(probability[g.x>c.drain_start].sum()),
                'classification':'narrow_or_dark_candidate' if converged and width<1e-10 else 'finite_width_or_unresolved_candidate',
                'occupation':'NOT_ASSIGNED','device_physical_origin':'INCONCLUSIVE'}
            found.append(row)
        except Exception as exc:found.append({'energy_ev':float(energy),'numerical_failure':repr(exc),'occupation':'NOT_ASSIGNED'})
    rows.append({'valley':s.valley.name,'seconds':time.time()-tic,'candidates':found});(out/'summary.json').write_text(json.dumps({'scope':'sampled real-part pole candidates of actual initial biased trial; NOT an SCF result or exhaustive exclusion','rows':rows},indent=2));print(s.valley.name,len(found),'candidates',time.time()-tic,flush=True)
finally:p.close()
