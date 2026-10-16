#!/usr/bin/env python3
"""Actual independent finite/open benchmarks; saves numeric evidence."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh
from config import Config,KIN,Q,HBAR_EV_S,Valley
from modespace.hamiltonian import System,Grid
from transport.ms_postprocess import transmission,direct_transmission,postprocess
from transport.negf_reference import solve_negf
from solver import write_json,atomic_npz
from provenance import source_digest


def chain(n=15,barrier=0.,weak=None):
    c=Config(dx=1,dy=1,dz=1,energy_max_nodes=180000)
    mass=KIN/.2;valley=Valley('analytic_chain',mass,mass,mass,1);tx=.2
    x=np.arange(n)+.5;mask=np.ones((n,1,1),bool);xyz=np.column_stack([x,x*0,x*0])
    grid=Grid(x,np.array([0.]),np.array([0.]),mask,xyz,np.ones(n,int),np.arange(n+1),1.)
    transverse=8*tx;onsite=np.full(n,transverse+2*tx);onsite[[0,-1]]+=tx
    onsite[n//3:2*n//3]+=barrier
    hopping=np.full(n-1,-tx)
    if weak is not None:hopping[:]= -weak
    H=np.diag(onsite)+np.diag(hopping,1)+np.diag(hopping,-1)
    sys=System(H,[np.ones((1,1,1)) for _ in x],[np.array([transverse]) for _ in x],
        np.arange(n+1),np.where(x<n/3,0,np.where(x>2*n/3,2,1)),grid,
        [np.array([[h]]) for h in hopping],valley)
    return c,sys,transverse


def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);report={'source_sha256':source_digest(),'cases':{}}
    for barrier in (0.,.15,.35):
        c,s,edge=chain(barrier=barrier);e=np.linspace(edge+.001,edge+.79,161)
        r=transmission(s,c,e,0);direct=direct_transmission(s,c,e,0)
        maxerr=float(max(abs(r-direct)))
        result={'dense_rgf_absolute_error':maxerr,'transmission_min':float(min(r)),
            'transmission_max':float(max(r)),'config':c.to_dict()}
        if barrier==0:
            result['uniform_ballistic_max_error']=float(max(abs(r-1)))
            mu=edge+.2;kbt=c.temperature*8.617333262145e-5
            blocks,flux,diag,_=solve_negf(s,c,mu,mu,0)
            # Infinite-chain analytic k integral; every site must have same density.
            expected=quad(lambda k:1/(1+np.exp((edge+2*.2*(1-np.cos(k))-mu)/kbt))/np.pi,0,np.pi,epsabs=1e-12)[0]
            result['equilibrium_population_error']=float(max(abs(np.array(blocks)[:,0,0].real-expected)))
            result['equilibrium_current_max_a']=float(max(abs(flux)))
            # A deliberately fixed frozen H, not a biased Poisson solution.
            post=postprocess(s,c,mu,.05)
            current,err=quad(lambda E:transmission(s,c,[E],.05)[0]*(1/(1+np.exp((E-mu)/kbt))-1/(1+np.exp((E-mu+.05)/kbt))),
                             edge,edge+.8,epsabs=1e-12,epsrel=1e-10,points=[edge+.05])
            current*=Q/HBAR_EV_S/(2*np.pi)*2
            result['current_relative_error']=float(abs(current-post['current_a'])/abs(current))
            result['current_conservation_relative']=post['conservation_relative']
        report['cases'][f'barrier_{barrier:g}']=result
        atomic_npz(out/f'barrier_{barrier:g}.npz',energy_ev=e,transmission_rgf=r,transmission_dense=direct,
                   source_sha256=source_digest(),config_json=json.dumps(c.to_dict()))
    # Quantum-well eigenspectrum against the exact discrete Dirichlet formula.
    n=31;t=.2;H=np.diag(np.full(n,2*t))+np.diag(np.full(n-1,-t),1)+np.diag(np.full(n-1,-t),-1)
    computed=eigh(H,eigvals_only=True);exact=2*t*(1-np.cos(np.arange(1,n+1)*np.pi/(n+1)))
    report['quantum_well_max_error_ev']=float(max(abs(computed-exact)))
    checks=[v['dense_rgf_absolute_error']<1e-10 for v in report['cases'].values()]
    uniform=report['cases']['barrier_0']
    checks.extend([uniform['uniform_ballistic_max_error']<1e-10,uniform['equilibrium_population_error']<1e-8,
        uniform['current_relative_error']<1e-5,report['quantum_well_max_error_ev']<1e-12])
    report['pass']=bool(all(checks));write_json(out/'summary.json',report)
    print(json.dumps(report,indent=2));return report


if __name__=='__main__':run(Path(sys.argv[1]) if len(sys.argv)>1 else Path('results/analytical'))
