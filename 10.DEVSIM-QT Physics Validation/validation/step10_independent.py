"""Independent scattering, dense-GF and full-stencil transport benchmarks."""
import sys,json,csv
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.integrate import quad_vec
from config import Config,KIN,Q,HBAR_EV_S,VALLEYS
from modespace.hamiltonian import Grid,ModeBuilder
from validation.analytical import chain
from validation.reduced_realspace import full_operator
from multispace.occupation import fermi
from transport.negf_reference import surface_sigma,device_blocks,lead_matrices,rgf,solve_negf
from transport.ms_postprocess import transmission
from transport.bound_states import audit


def scattering_transmission(energy,potentials,t,left=0.,right=0.):
    """Independent wave-matching solve of the finite-difference equation.

    Unknowns are reflection amplitude, all interior wavefunction values and
    transmission amplitude. No Green function or contact self-energy is used.
    """
    if not(left<energy<left+4*t and right<energy<right+4*t):return 0.
    kl=np.arccos(1-(energy-left)/(2*t));kr=np.arccos(1-(energy-right)/(2*t))
    v=np.asarray(potentials);n=len(v);A=np.zeros((n+2,n+2),complex);b=np.zeros(n+2,complex)
    a=left+2*t-energy
    A[0,0]=a-t*np.exp(1j*kl);A[0,1]=-t;b[0]=-a+t*np.exp(-1j*kl)
    for j in range(n):
        row=j+1;A[row,j+1]=v[j]+2*t-energy
        if j==0:A[row,0]=-t;b[row]=t
        else:A[row,j]=-t
        if j==n-1:A[row,-1]=-t*np.exp(1j*kr*(n+1))
        else:A[row,j+2]=-t
    A[-1,-2]=-t
    A[-1,-1]=(right+2*t-energy)*np.exp(1j*kr*(n+1))-t*np.exp(1j*kr*(n+2))
    amplitude=np.linalg.solve(A,b)[-1]
    return float(abs(amplitude)**2*np.sin(kr)/np.sin(kl))


def chain_checks(out):
    rows=[];curves={}
    for barrier in (0.,.08,.25):
        c,s,edge=chain(n=13,barrier=barrier)
        energies=edge+np.linspace(.001,.799,120)
        actual=transmission(s,c,energies,0.)
        potentials=np.diag(s.H).copy()-edge-.4;potentials[[0,-1]]-=.2
        expected=np.array([scattering_transmission(e-edge,potentials,.2) for e in energies])
        err=float(max(abs(actual-expected)))
        rows.append({'case':f'barrier_{barrier:g}','maximum_transmission_error':err,'pass':err<2e-10})
        curves[f'barrier_{barrier:g}']=(energies,actual,expected)
    # A real electrostatic step includes the same shift in the right lead.
    c,s,edge=chain(n=13);vd=.08;s.H-=np.diag(np.r_[np.zeros(6),np.full(7,vd)])
    energies=edge+np.linspace(.002,.70,90)
    actual=transmission(s,c,energies,vd)
    expected=np.array([scattering_transmission(e-edge,np.r_[np.zeros(6),np.full(7,-vd)],.2,0.,-vd) for e in energies])
    err=float(max(abs(actual-expected)))
    rows.append({'case':'step','maximum_transmission_error':err,'pass':err<2e-10})
    # Independent dense lesser-GF energy integration, not the contour formula.
    for vd in (0.,.01,.05,.1):
        c,s,edge=chain(n=13,barrier=.03);mu=edge+.12;s.H-=np.diag(np.linspace(0,vd,13))
        ons,hop,t,leads=device_blocks(s,c,(0.,vd));H=s.H.copy();H[0,0]-=t;H[-1,-1]-=t
        edges=np.array([edge-vd,edge,edge+4*t-vd,edge+4*t]);eye=np.eye(len(H))
        def integrand(e):
            sl,sr=lead_matrices(np.array([complex(e)]),t,leads)
            A=e*eye-H.astype(complex);A[0,0]-=sl[0,0,0];A[-1,-1]-=sr[0,0,0]
            G=np.linalg.solve(A,eye);injection=np.zeros_like(H)
            injection[0,0]=-2*sl[0,0,0].imag*fermi(e,mu,300)
            injection[-1,-1]=-2*sr[0,0,0].imag*fermi(e,mu-vd,300)
            return G@injection@G.conj().T/(2*np.pi)
        density,error=quad_vec(integrand,edges.min(),edges.max(),points=edges,epsabs=2e-10,epsrel=2e-9)
        blocks,bonds,diagnostic,_=solve_negf(s,c,mu,mu-vd,vd)
        direct_i=-2*Q/HBAR_EV_S*np.imag(-t*np.diag(density,-1))
        ne=float(max(abs(np.asarray(blocks)[:,0,0].real-np.diag(density).real)))
        ie=float(max(abs(bonds-direct_i)))
        rows.append({'case':f'dense_lesser_vd{vd:g}','density_max_abs_error':ne,
            'current_max_abs_error_a':ie,'dense_integrator_error':float(error),
            'pass':ne<1e-7 and ie<1e-12})
        np.savez_compressed(out/f'chain_vd{vd:g}.npz',density_dense=density,
            density_rgf=np.asarray(blocks),bond_current_a=bonds,dense_bond_current_a=direct_i)
    # Surface quadratic residual, causality and band-gap transmission.
    z=np.r_[np.linspace(1.01,2.19,100)+1e-4j,np.array([1.,2.4])+1e-4j]
    sig=surface_sigma(z,np.array([2.]),.2)[:,0]
    residual=float(max(abs(sig*(z-2.-sig)-.2**2)))
    rows.append({'case':'surface_retarded','quadratic_residual_ev2':residual,
        'max_imaginary_sigma_ev':float(max(sig.imag)),
        'pass':residual<1e-12 and max(sig.imag)<=0})
    c,s,edge=chain(n=13)
    outside=transmission(s,c,[edge-.1,edge+.9],0.)
    rows.append({'case':'outside_lead_bands','transmission_max':float(max(abs(outside))),
        'pass':max(abs(outside))<1e-12})
    return rows,curves


def full3d_checks():
    c=Config(ls=2.,lg=2.,ld=2.,spacer=1.,dx=1.,dy=1.,dz=1.,modes=6,sd_modes=6)
    nx,ny,nz=8,3,2;x=np.arange(nx)+.5;y=np.arange(ny)-1.;z=np.arange(nz)-.5
    X,Y,Z=np.meshgrid(x,y,z,indexing='ij');xyz=np.column_stack([X.ravel(),Y.ravel(),Z.ravel()])
    mask=np.ones((nx,ny,nz),bool);g=Grid(x,y,z,mask,xyz,np.full(nx,ny*nz),np.arange(nx+1)*ny*nz,1.)
    phi=np.zeros(len(xyz));phi[(xyz[:,0]>3)&(xyz[:,0]<5)]=-.04
    v=VALLEYS[0];system=ModeBuilder(c,g,v).build(phi);U=system.realspace_transform()
    real=full_operator(g,c,v,phi);projection=float(max(abs((system.H-U.T@real@U).ravel())))
    ons,hop,t,leads=device_blocks(system,c,(0.,0.));m=ny*nz
    real[:m,:m]-=t*np.eye(m);real[-m:,-m:]-=t*np.eye(m)
    full_basis=system.bases[0][mask[0]];edge=float(min(leads[0][0])-2*t)
    errors=[];thresholds=[]
    # Uniform multichannel thresholds must count all propagating subbands.
    uniform=ModeBuilder(c,g,v).build(np.zeros(len(xyz)))
    for e in np.linspace(edge+.01,edge+1.0,45):
        sig_l,sig_r=lead_matrices(np.array([complex(e)]),t,leads)
        Sl=full_basis@sig_l[0]@full_basis.T;Sr=full_basis@sig_r[0]@full_basis.T
        A=e*np.eye(len(real))-real.astype(complex);A[:m,:m]-=Sl;A[-m:,-m:]-=Sr
        G=np.linalg.inv(A);gl=np.zeros_like(A);gr=gl.copy()
        gl[:m,:m]=1j*(Sl-Sl.conj().T);gr[-m:,-m:]=1j*(Sr-Sr.conj().T)
        dense=(G@(gl*.8+gr*.2)@G.conj().T)/(2*np.pi)
        _,left=rgf(np.array([e]),ons,hop,sig_l,sig_r,False,0)
        _,right=rgf(np.array([e]),ons,hop,sig_l,sig_r,False,1)
        reconstructed=[]
        for i,(a,b) in enumerate(zip(left,right)):
            block=(.8*a[0]@a[0].conj().T+.2*b[0]@b[0].conj().T)/(2*np.pi)
            u=system.bases[i][mask[i]];reconstructed.extend(np.diag(u@block@u.T).real)
        errors.append(float(max(abs(np.array(reconstructed)-np.diag(dense).real))))
        count=int(sum((e>leads[0][0]-2*t)&(e<leads[0][0]+2*t)))
        thresholds.append(abs(float(transmission(uniform,c,[e],0.)[0])-count))
    return [{'case':'full_3d_projection','max_operator_error_ev':projection,'pass':projection<1e-11},
        {'case':'full_3d_open_dense_lesser','max_density_spectral_error_per_ev':max(errors),'pass':max(errors)<1e-8},
        {'case':'uniform_multimode','max_transmission_channel_count_error':max(thresholds),'pass':max(thresholds)<1e-9}]


def main(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rows,curves=chain_checks(out);rows+=full3d_checks()
    result={'scope':'independent frozen-H open quantum transport; not a GAA SCF certification',
        'cases':rows,'pass':bool(all(r['pass'] for r in rows))}
    (out/'summary.json').write_text(json.dumps(result,indent=2,default=lambda v:v.item())+'\n')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4))
    for name,(e,t,ref) in curves.items():
        ax.plot(e,t,label=name);ax.plot(e[::6],ref[::6],'k.',ms=2)
    ax.set(xlabel='Energy (eV)',ylabel='Transmission per spin/valley');ax.legend();fig.tight_layout()
    fig.savefig(out/'independent_transmission.png',dpi=180);plt.close(fig)
    print(json.dumps(result,indent=2,default=lambda v:v.item()),flush=True)
    return 0 if result['pass'] else 1

if __name__=='__main__':raise SystemExit(main(sys.argv[1] if len(sys.argv)>1 else 'results/independent'))
