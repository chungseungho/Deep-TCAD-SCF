#!/usr/bin/env python3
"""Executed reduced references; not a proxy for transistor SCF convergence."""
import csv,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh
from config import KIN,VALLEYS,Q,HBAR_EV_S,KB_EV
from modespace.hamiltonian import transverse_matrix
from transport.finite_spectral import solve_finite_spectral
from transport.negf_reference import direct_finite_gf,surface_sigma,rgf,fermi_poles
from multispace.occupation import fermi


def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);tic=time.perf_counter()
    H=np.diag([.2,.3,.25])+np.diag([-.07,-.1],1)+np.diag([-.07,-.1],-1)
    regions=np.array([0,1,2]);rho,diag,_=solve_finite_spectral(H,regions,.03,.3,.2,300)
    reference,error=direct_finite_gf(H,regions,.03,.3,.2,300)
    bond=-2*Q/HBAR_EV_S*np.imag(H[0,1]*rho[1,0])
    def landauer(e):
        G=np.linalg.inv(e*np.eye(3)-H+.5j*np.diag([.03,0,.03]))
        return Q/HBAR_EV_S/(2*np.pi)*.03**2*abs(G[0,2])**2*(fermi(e,.3,300)-fermi(e,.2,300))
    current=quad(landauer,-2,2,epsabs=1e-18,epsrel=1e-10)[0]
    subbands=[]
    for factor in [1,2,4]:
        ny,nz=12*factor,5*factor;dy=dz=1/factor
        for v in VALLEYS:
            from scipy.sparse.linalg import eigsh
            matrix=transverse_matrix(np.ones((ny,nz),bool),KIN/(v.my*dy**2),KIN/(v.mz*dz**2))
            numerical=float(eigsh(matrix,k=1,which='SA',return_eigenvectors=False,tol=1e-12)[0])
            exact=np.pi**2*KIN*(1/(v.my*12**2)+1/(v.mz*5**2))
            subbands.append({'valley':v.name,'dy_dz_nm':dy,'subband_ev':numerical,'continuum_ev':exact,'relative_error':abs(numerical/exact-1)})
    with (out/'subband_convergence.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(subbands[0]));w.writeheader();w.writerows(subbands)
    t=.3;ml=.3;mr=.2;ideal=2*Q/HBAR_EV_S/(2*np.pi)*quad(lambda E:fermi(E,ml,300)-fermi(E,mr,300),0,4*t,epsabs=1e-12)[0]
    contact=[]
    for length in [8,24,64]:
        N=2*length+9;chain=np.diag(np.full(N,2*t))+np.diag(np.full(N-1,-t),1)+np.diag(np.full(N-1,-t),-1)
        reg=np.r_[np.zeros(length,int),np.ones(9,int),np.full(length,2)]
        for gamma in [.005,.01,.03,.06]:
            C,d,_=solve_finite_spectral(chain,reg,gamma,ml,mr,300)
            flux=-4*Q/HBAR_EV_S*np.imag(-t*np.diag(C,-1));value=float(flux[length-1:length+9].mean())
            contact.append({'lead_sites':length,'gamma_ev':gamma,'current_a':value,'semiinfinite_ideal_a':ideal,'relative_error':abs(value/ideal-1)})
    with (out/'contact_convergence.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(contact[0]));w.writeheader();w.writerows(contact)
    summary={'finite_spectral_vs_direct_gf_max_abs':float(max(abs(rho-reference).ravel())),
             'direct_gf_quadrature_error_estimate':float(error),
             'bond_current_a':float(bond),'direct_transmission_current_a':float(current),
             'bond_vs_landauer_relative':float(abs(bond/current-1)),
             'finite_density_matrix_pauli_min':diag['pauli_min'],'finite_density_matrix_pauli_max':diag['pauli_max'],
             'best_chain_finite_vs_semiinfinite_error':min(r['relative_error'] for r in contact),
             'runtime_s':time.perf_counter()-tic,'scope':'reduced systems only; no actual GAA contact plateau certified'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2));return summary


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/benchmarks');a=p.parse_args()
    print(json.dumps(run(a.out),indent=2))
