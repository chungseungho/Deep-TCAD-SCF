"""Independent separable uniform-wire hard-wall/open convergence and pole audit."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
from config import VALLEYS,KIN,KB_EV
from multispace.occupation import fermi
from solver import write_json,atomic_npz
from provenance import source_digest


def transverse(v,dy=1.,dz=1.):
    ny,nz=round(20/dy),round(15/dz)
    ey=4*KIN/(v.my*dy**2)*np.sin(np.arange(1,ny+1)*np.pi/(2*ny))**2
    ez=4*KIN/(v.mz*dz**2)*np.sin(np.arange(1,nz+1)*np.pi/(2*nz))**2
    return np.sort((ey[:,None]+ez).ravel())


def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);dx=.5;T=300.;kbt=KB_EV*T
    x,w=leggauss(512);k=(x+1)*np.pi/2;w=w/2
    spectra=[transverse(v) for v in VALLEYS];tx=[KIN/v.mx/dx**2 for v in VALLEYS]
    def bulk(mu,m=300):
        return sum(4*np.sum(fermi(e[:m,None]+2*t*(1-np.cos(k)),mu,T)@w)/dx for e,t in zip(spectra,tx))
    mu=brentq(lambda mu:bulk(mu)-3.,-1,1,xtol=1e-14)
    rows=[];levels=[4,8,12,16,20,32,48,64,100,300]
    for L in (12.,24.,48.,96.):
        N=round(L/dx);pos=(np.arange(N)+.5)*dx
        eig=[]
        for t in tx:
            onsite=np.full(N,2*t);onsite[[0,-1]]+=t
            e,u=eigh_tridiagonal(onsite,np.full(N-1,-t));eig.append((e,u*u))
        for m in levels:
            def finite_count(mu):return sum(4*fermi(ex[:,None]+et[None,:m],mu,T).sum() for (ex,u),et in zip(eig,spectra))
            mN=brentq(lambda mm:finite_count(mm)-3*L,-1,3,xtol=1e-14)
            for prescription,mm in [('fixed_mu',mu),('fixed_N',mN)]:
                density=np.zeros(N);thermal=np.zeros(N)
                for (ex,u),et in zip(eig,spectra):
                    f=fermi(ex[:,None]+et[None,:m],mm,T)
                    density+=4*(u@f.sum(axis=1))/dx
                    thermal+=4*(u@(f*(1-f)/kbt).sum(axis=1))/dx
                interior=abs(pos-L/2)<3.;outer=pos<3.
                row={'length_nm':L,'modes':m,'prescription':prescription,'mu_ev':mm,
                     'electrons':float(density.sum()*dx),'donor_electrons':3*L,
                     'open_line_density_same_mu':bulk(mm,m),'open_line_density_full_basis_same_mu':bulk(mm),
                     'center_line_density':float(density[interior].mean()),
                     'outer_3nm_line_density':float(density[outer].mean()),
                     'center_relative_error_to_open':float(abs(density[interior].mean()/bulk(mm,m)-1)),
                     'end_region_neutrality':float(density[outer].mean()/3-1),
                     'closed_total_neutrality':float(density.sum()*dx/(3*L)-1),
                     'center_thermal_ldos_per_ev_nm':float(thermal[interior].mean()),
                     'potential_gradient_v_nm':0.,'potential_is_imposed_uniform':True,'scf':False}
                rows.append(row)
                if m in (4,20,300):atomic_npz(out/f'L{L:g}_M{m}_{prescription}.npz',x_nm=pos,line_density_per_nm=density,
                    thermal_ldos_per_ev_nm=thermal,mu_ev=mm,source_sha256=source_digest())
    bulkrows=[{'modes':m,'line_density_fixed_mu':bulk(mu,m),'relative_missing_density':1-bulk(mu,m)/3} for m in levels]
    gridrows=[]
    for h in (1.,.5,.25):
        es=[transverse(v,h,h) for v in VALLEYS]
        def num(mm):return sum(4*np.sum(fermi(e[:,None]+2*t*(1-np.cos(k)),mm,T)@w)/dx for e,t in zip(es,tx))
        gridrows.append({'dy_dz_nm':h,'full_modes':len(es[0]),'fixed_mu_line_density':num(mu),
            'neutral_mu_ev':brentq(lambda mm:num(mm)-3,-1,1)})
    # Finite-chain self-energy is a positive-residue pole expansion. Evaluate
    # only off the real axis: this is not a steady biased density solver.
    from transport.negf_reference import surface_sigma
    aux=[];z=np.linspace(-.1,1.,201)+.025j;t=.2;eps=.4
    truth=surface_sigma(z,np.array([eps]),t)[:,0]
    for n in (8,16,32,64,128):
        e,u=eigh_tridiagonal(np.full(n,eps),np.full(n-1,-t))
        value=np.sum((t*t*u[0]**2)[None,:]/(z[:,None]-e),axis=1)
        aux.append({'poles':n,'max_self_energy_error_ev':float(max(abs(value-truth))),'imaginary_argument_ev':.025})
    write_json(out/'summary.json',{'source_sha256':source_digest(),'temperature_k':T,'dx_nm':dx,
        'full_basis_neutral_mu_ev':mu,'bulk_mode_sweep':bulkrows,'cases':rows,'transverse_grid':gridrows,
        'pole_embedding_complex_axis_only':aux,'biased_auxiliary_density_implemented':False,
        'ldos_definition':'thermal convolution of delta DOS with -df/dE at the specified mu, not artificial GF broadening',
        'scope':'separable uniform 20x15 nm wire, frozen zero potential; no GAA screening claim'})
    print('uniform completed',mu,flush=True)

if __name__=='__main__':run(sys.argv[1] if len(sys.argv)>1 else 'results/uniform')
