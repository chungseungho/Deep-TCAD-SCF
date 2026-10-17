"""Nonequilibrium open-density checks with independent dense real-axis integration."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import replace
import numpy as np
from scipy.integrate import quad_vec
from config import Q,HBAR_EV_S
from validation.analytical import chain
from multispace.occupation import fermi,solve_states,occupations
from transport.negf_reference import solve_negf,device_blocks,lead_matrices
from transport.ms_postprocess import transmission
from transport.finite_spectral import solve_finite_spectral
from transport.reservoir import bond_currents
from solver import write_json,atomic_npz
from provenance import source_digest


def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[]
    for vd in (0.,.001,.01,.05,.1):
        c,s,edge=chain(n=13,barrier=.03);mu=edge+.12
        phi=np.linspace(0,vd,13);s.H-=np.diag(phi)
        ons,hop,t,leads=device_blocks(s,c,(0.,vd));Ho=s.H.copy();Ho[0,0]-=t;Ho[-1,-1]-=t
        edges=np.array([edge-vd,edge,edge+4*t-vd,edge+4*t]);eye=np.eye(len(Ho))
        def integrand(e):
            sl,sr=lead_matrices(np.array([complex(e)]),t,leads)
            A=e*eye-Ho.astype(complex);A[0,0]-=sl[0,0,0];A[-1,-1]-=sr[0,0,0]
            G=np.linalg.solve(A,eye);inj=np.zeros_like(Ho)
            inj[0,0]=-2*sl[0,0,0].imag*fermi(e,mu,300);inj[-1,-1]=-2*sr[0,0,0].imag*fermi(e,mu-vd,300)
            return G@inj@G.conj().T/(2*np.pi)
        direct,error=quad_vec(integrand,edges.min(),edges.max(),points=edges,epsabs=2e-10,epsrel=2e-9)
        I=-2*Q/HBAR_EV_S*np.imag(-t*np.diag(direct,-1))
        for panel in (.02,.01,.005):
            cc=replace(c,energy_step=panel,energy_rtol=1e-6,energy_atol=1e-9)
            blocks,bonds,diag,_=solve_negf(s,cc,mu,mu-vd,vd)
            density=np.array(blocks)[:,0,0].real
            row={'vd_v':vd,'energy_panel_ev':panel,'density_max_abs_error':float(max(abs(density-np.diag(direct).real))),
                'current_a_per_degeneracy':float(bonds.mean()),'direct_current_a_per_degeneracy':float(I.mean()),
                'current_abs_error_a':float(max(abs(bonds-I))),
                'current_conservation_relative':float(max(abs(bonds-bonds.mean()))/max(abs(bonds.mean()),1e-15)),
                'direct_density_matrix_min_eigenvalue':float(min(np.linalg.eigvalsh(direct))),
                'gf_nodes':diag['real_energy_nodes'],'dense_integrator_error_estimate':float(error)}
            rows.append(row)
        atomic_npz(out/f'vd{vd:g}.npz',full_density_matrix=direct,density=density,phi_v=phi,current_per_bond_a=bonds)
    # Standing-state counterexample at fixed uniform H: spatial probability is
    # symmetric while the open wire transmits. No tracking can infer direction.
    c,s,edge=chain(n=31);states=occupations(solve_states(s.H,s.regions,c),edge+.15,edge+.10,300,strict=False)
    counter={'active_states':int(states['active'].sum()),'ambiguous_active_states':int(states['ambiguous'].sum()),
        'max_absolute_LR_weight_difference':float(max(abs(states['weights'][:,0]-states['weights'][:,2]))),
        'open_transmission_at_mu':float(transmission(s,c,[edge+.15],0)[0]),
        'closed_H_time_reversal_real':bool(np.max(abs(s.H.imag))==0)}
    # The generic chain partition has an asymmetric number of cells for n=31;
    # use exactly equal outer regions for this analytical counterexample.
    regions=np.full(31,1);regions[:10]=0;regions[-10:]=2
    states=occupations(solve_states(s.H,regions,c),edge+.15,edge+.10,300,strict=False)
    counter.update(active_states=int(states['active'].sum()),ambiguous_active_states=int(states['ambiguous'].sum()),
                   max_absolute_LR_weight_difference=float(max(abs(states['weights'][:,0]-states['weights'][:,2]))))
    bath=[]
    for L in (8,16,32,64):
        c,s,edge=chain(n=2*L+5);mu=edge+.12;vd=.05;n=len(s.H)
        phi=np.r_[np.zeros(L),np.linspace(0,vd,7)[1:-1],np.full(L,vd)]
        s.H-=np.diag(phi);s.H[L:L+5,L:L+5]+=.03*np.eye(5)
        s.regions=np.r_[np.zeros(L,int),np.ones(5,int),np.full(L,2,int)]
        blocks,flux,_,_=solve_negf(s,c,mu,mu-vd,vd)
        refn=np.array(blocks)[:,0,0].real;refI=flux.mean()*2
        for gamma in (.0025,.005,.01,.02):
            rho,diag,_=solve_finite_spectral(s.H,s.regions,gamma,mu,mu-vd,300)
            bonds=bond_currents(s,rho,2);I=float(bonds[L:L+4].mean())
            bath.append({'electrode_cells':L,'gamma_ev':gamma,'current_a':I,'open_current_a':float(refI),
                'current_relative_to_open':float(abs(I-refI)/abs(refI)),
                'channel_density_relative_to_open':float(np.linalg.norm(np.diag(rho).real[L:L+5]-refn[L:L+5])/np.linalg.norm(refn[L:L+5])),
                'current_conservation_relative':float(max(abs(bonds[L:L+4]-I))/abs(I))})
    result={'source_sha256':source_digest(),'open_chain':rows,'standing_state_counterexample':counter,
            'finite_bath_contact_sweep':bath,'scope':'independent small frozen chains, not fully SCF GAA I-V',
            'pass':all(r['density_max_abs_error']<1e-7 and r['current_abs_error_a']<1e-12 for r in rows)}
    write_json(out/'summary.json',result);print('open chain pass',result['pass'],flush=True)

if __name__=='__main__':run(sys.argv[1] if len(sys.argv)>1 else 'results/open_chain')
