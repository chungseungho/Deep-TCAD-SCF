#!/usr/bin/env python3
"""Independent finite-bath RGF transmission AFTER analytic-spectral SCF.

Eliminate the actual finite broadened S/D blocks by Schur complements.
This uses the same bath model, not a semi-infinite contact substitution and
not the complex eigensystem/digamma formula used for SCF density.
"""
import argparse,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config,VALLEYS,KB_EV,Q,HBAR_EV_S
from modespace.hamiltonian import Grid,ModeBuilder
from transport.negf_reference import rgf
from transport.quadrature import integrate_vector
from multispace.occupation import fermi


def finite_rgf_transmission(system,gamma,energies):
    e=np.asarray(energies,dtype=complex);n=len(system.grid.x)
    h=[system.H[system.sl(i),system.sl(i)] for i in range(n)];t=system.hop
    regions=[system.regions[system.offsets[i]] for i in range(n)]
    central=np.where(np.array(regions)==1)[0];a,b=central[0],central[-1]
    def inv(i,extra=0): return np.linalg.inv(e[:,None,None]*np.eye(len(h[i]))-h[i]+.5j*gamma*np.eye(len(h[i]))-extra)
    gl=inv(0)
    for i in range(1,a): gl=inv(i,t[i-1].conj().T@gl@t[i-1])
    sl=t[a-1].conj().T@gl@t[a-1]
    gr=inv(n-1)
    for i in range(n-2,b,-1): gr=inv(i,t[i]@gr@t[i].conj().T)
    sr=t[b]@gr@t[b].conj().T
    _,amp=rgf(e,h[a:b+1],t[a:b],sl,sr,False,0)
    gam=1j*(sr-sr.conj().transpose(0,2,1));A=amp[-1]
    transmission=np.einsum('eij,ejk,eik->e',gam,A,A.conj(),optimize=True).real
    if transmission.min()<-1e-10: raise RuntimeError('Nonpositive Schur transmission')
    return np.maximum(transmission,0)


def run(directory,key,out):
    directory=Path(directory);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    c=Config(**json.loads((directory/'settings.json').read_text())['config'])
    if c.backend!='spectral': raise ValueError('This independent check is for the finite spectral bath model')
    s=np.load(directory/'states'/(key+'.npz'));grid=Grid.build(c);systems=[ModeBuilder(c,grid,v).build(s['phi_v']) for v in VALLEYS]
    ml=float(s['mu_l_ev']);mr=float(s['mu_r_ev']);kt=KB_EV*c.temperature
    low=min(ml,mr)-20*kt;high=max(ml,mr)+20*kt;start=time.perf_counter()
    def integrand(e):
        df=fermi(e,ml,c.temperature)-fermi(e,mr,c.temperature)
        return np.column_stack([2*sys.valley.multiplicity*finite_rgf_transmission(sys,c.gamma,e)*df for sys in systems])
    values=[]
    for panel in [.02,.01,.005]:
        breaks=np.r_[np.arange(low,high,panel),high]
        integral,diag=integrate_vector(integrand,breaks,8,2e-7,1e-16,60000)
        current=float(integral.real.sum()*Q/HBAR_EV_S/(2*np.pi));bond=float(np.mean(s['bond_current_a'][(grid.x[:-1]>=c.ls-c.dx/2)&(grid.x[1:]<=c.drain_start+c.dx/2)]))
        values.append({'initial_panel_ev':panel,'current_a':current,'bond_current_a':bond,'relative_mismatch':abs(current-bond)/max(abs(bond),1e-18),**diag})
    E=np.linspace(low,high,600);T=sum(2*sys.valley.multiplicity*finite_rgf_transmission(sys,c.gamma,E) for sys in systems)
    np.savez_compressed(out/'transmission.npz',energy_ev=E,spin_valley_summed_transmission=T)
    result={'same_hamiltonian_same_finite_bath':True,'gf_only_after_scf':True,'checks':values,'runtime_s':time.perf_counter()-start,
            'pass':bool(all(r['relative_mismatch']<2e-5 for r in values)),
            'spectrum_note':'600-point plot samples; independent current uses adaptive quadrature, not this sampled curve'}
    (out/'summary.json').write_text(json.dumps(result,indent=2));return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory');p.add_argument('key');p.add_argument('--out',default='results/post_scf_current');a=p.parse_args()
    print(json.dumps(run(a.directory,a.key,a.out),indent=2))
