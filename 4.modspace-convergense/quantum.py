"""Fixed-basis coupled mode-space ballistic NEGF for rectangular mono-Si nanosheet.

Fixed transverse EMA basis (same masses/constants as STEP01), projected
3D potential, full off-diagonal intermode electrostatics. Source/drain
semi-infinite tight-binding lead self-energies. Independent of DEVSIM.

No DD and no fitted effective mass. Valley multiplicity and spin explicit.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import eigsh

Q=1.602176634e-19
H_PLANCK=6.62607015e-34
HB2_OVER_2M0=0.0380998212  # same STEP01 constant, eV nm^2
KB_EV=8.617333262145e-5


def fermi(e, mu, kt):
    return 1/(1+np.exp(np.clip((e-mu)/kt,-150,150)))


def transverse_sparse(cfg, valley):
    ny,nz=cfg.ny,cfg.nz
    dy=cfg.width_nm/(ny+1); dz=cfg.thickness_nm/(nz+1)
    ty=HB2_OVER_2M0/(valley.my*dy*dy)
    tz=HB2_OVER_2M0/(valley.mz*dz*dz)
    ly=diags([-np.ones(ny-1),2*np.ones(ny),-np.ones(ny-1)],[-1,0,1],format='csr')*ty
    lz=diags([-np.ones(nz-1),2*np.ones(nz),-np.ones(nz-1)],[-1,0,1],format='csr')*tz
    return kron(ly,eye(nz,format='csr'),format='csr')+kron(eye(ny,format='csr'),lz,format='csr')


@dataclass
class Basis:
    valley: object
    levels: np.ndarray
    U: np.ndarray
    hopping_ev: float


def make_bases(cfg):
    result=[]
    for v in cfg.valleys:
        mat=transverse_sparse(cfg,v)
        if cfg.modes >= mat.shape[0]-1:
            raise ValueError('Need fewer retained modes than transverse lattice sites-1')
        levels,U=eigsh(mat,k=cfg.modes,which='SA',tol=1e-10,maxiter=20000)
        inds=np.argsort(levels)
        levels,U=levels[inds],U[:,inds]
        # Set deterministic signs for reproducibility across jobs.
        for j in range(U.shape[1]):
            q=np.argmax(np.abs(U[:,j]))
            if U[q,j]<0: U[:,j]*=-1
        dx=cfg.length_nm/(cfg.nx-1)
        hopping=HB2_OVER_2M0/(v.mx*dx*dx)
        result.append(Basis(v,levels,U,hopping))
    return result


def get_onsite(phi,b,ec0):
    nx,n_sites=phi.shape[0],phi.shape[1]*phi.shape[2]
    flattened=phi.reshape(nx,n_sites)
    U=b.U
    shift=np.einsum('pm,xp,pn->xmn',U,flattened,U,optimize=True)
    d=b.levels+ec0+2*b.hopping_ev
    hs=-shift.astype(complex)
    hs[:,np.arange(len(d)),np.arange(len(d))] += d
    # Leads must inherit the FULL endpoint transverse matrix: diagonal-only
    # approximation breaks when transverse contact potential couples modes.
    return hs,hs[0].copy(),hs[-1].copy()


def lead_self_energy(E, e_sub, t, eta):
    """Retarded semi-infinite homogeneous 1D lead, vectorized [energy, mode]."""
    z=E[:,None]+1j*eta
    a=z-(e_sub[None,:]+2*t)
    root=np.sqrt(a*a-4*t*t)
    s1,s2=(a-root)*.5,(a+root)*.5
    s=np.where(np.abs(s1)<=np.abs(s2),s1,s2)
    if np.any(s.imag>1e-8):
        raise ArithmeticError('Nonretarded contact self-energy')
    return s


def lead_matrix_self_energy(E, lead_onsite, t, eta):
    """Consistent semi-infinite contact for H0 with inter-mode mixing.

    Hlead is the full onsite matrix INCLUDING +2t. Because the longitudinal
    coupling is -t*I, it commutes with the unitary transverse lead basis.
    Hence diagonalize Hlead and apply the exact 1D retarded surface solution
    separately to its eigenchannels, then rotate the self energy back.
    """
    mat=np.asarray(lead_onsite)
    if mat.ndim==1:  # backward compatibility with STEP03 diagonal subband energies
        mat=np.diag(mat+2*t)
    if mat.shape[0]!=mat.shape[1]:raise ValueError('Lead onsite must be square')
    if np.max(np.abs(mat-mat.conj().T))>1e-8:raise ValueError('Non-Hermitian lead onsite')
    eig,U=np.linalg.eigh(mat)
    se=lead_self_energy(np.asarray(E,float),eig-2*t,t,eta)
    full=np.einsum('am,em,bm->eab',U,se,U.conj(),optimize=True)
    gamma=1j*(full-full.conj().swapaxes(-1,-2))
    # Gamma must be positive semidefinite for a passive contact.
    if np.min(np.linalg.eigvalsh(gamma)) < -1e-8:
        raise ArithmeticError('Lead broadening is not positive semidefinite')
    return full,gamma


def spectral_batched(onsite,t,es,ed,E,eta,with_density=True):
    """Batched-energy RGF transmission and injected diagonal density matrices.

    Returns (T[E], AL[E,x,m,n], AR[E,x,m,n]). Either density may be None.
    Matrix operations are batched across energies to avoid per-energy Python loops.
    """
    ne=len(E); nx,nm=onsite.shape[:2]
    ident=np.eye(nm,dtype=complex)
    z=E[:,None,None]+1j*eta
    sleft,GL=lead_matrix_self_energy(E,es,t,eta)
    sright,GR=lead_matrix_self_energy(E,ed,t,eta)
    # H_{i,i+1}=-t I; coupling self-energy = t² G.
    c=t*t
    right_conn=np.empty((ne,nx,nm,nm),complex)
    right_conn[:,-1]=np.linalg.inv(z*ident-onsite[-1][None,:,:]-sright)
    for i in range(nx-2,-1,-1):
        right_conn[:,i]=np.linalg.inv(z*ident-onsite[i][None,:,:]-c*right_conn[:,i+1])
    G00=np.linalg.inv(z*ident-onsite[0][None,:,:]-sleft-c*right_conn[:,1])
    G_i0=np.empty_like(right_conn)
    G_i0[:,0]=G00
    for i in range(1,nx):
        G_i0[:,i]=(-t)*np.matmul(right_conn[:,i],G_i0[:,i-1])
    # Right-to-left injection requires a left-connected recursion.
    left_conn=np.empty_like(right_conn)
    left_conn[:,0]=np.linalg.inv(z*ident-onsite[0][None,:,:]-sleft)
    for i in range(1,nx):
        left_conn[:,i]=np.linalg.inv(z*ident-onsite[i][None,:,:]-c*left_conn[:,i-1])
    G_iN=np.empty_like(right_conn)
    G_iN[:,-1]=np.linalg.inv(z*ident-onsite[-1][None,:,:]-sright-c*left_conn[:,-2])
    for i in range(nx-2,-1,-1):
        G_iN[:,i]=(-t)*np.matmul(left_conn[:,i],G_iN[:,i+1])
    G0N=G_iN[:,0]
    prod=np.matmul(np.matmul(GL,G0N),GR)
    transmission=np.maximum(0.0,np.einsum('eij,eji->e',prod,G0N.conj().swapaxes(-1,-2)).real)
    if not with_density: return transmission,None,None
    AL=np.matmul(np.matmul(G_i0,GL[:,None,:,:]),G_i0.conj().swapaxes(-1,-2))
    AR=np.matmul(np.matmul(G_iN,GR[:,None,:,:]),G_iN.conj().swapaxes(-1,-2))
    return transmission,AL,AR


def energy_grid(cfg, vd, basis, es, ed, de):
    """Energy range for *charge*: UNION of occupied lead continua.

    This is essential for high Vds: drain-only injected states below source
    band are charged even though they do not transmit any current. Using the
    intersection of the bands gives a severely wrong Poisson charge at 6 V.
    """
    t=basis.hopping_ev
    def extrema(lead):
        arr=np.asarray(lead)
        if arr.ndim==2:
            ev=np.linalg.eigvalsh(arr).real-2*t
        else:ev=arr
        return float(np.min(ev)),float(np.max(ev))
    lo_s,hi_s=extrema(es);lo_d,hi_d=extrema(ed)
    lower=min(lo_s,lo_d)
    upper=max(hi_s+4*t,hi_d+4*t)
    upper=min(upper,max(cfg.efermi_source_ev,cfg.efermi_source_ev-vd)
              +12*KB_EV*cfg.temperature_k)
    if upper<=lower: return np.array([]),np.array([])
    # Anchored mesh: identical physical energy gets the SAME quadrature node
    # as phi evolves; continuously shifting linspace caused SCF 2-cycle noise.
    kmin=int(np.floor(lower/de));kmax=int(np.ceil(upper/de))
    energy=np.arange(kmin,kmax+1,dtype=float)*de
    w=np.full(len(energy),de,dtype=float)
    w[[0,-1]]*=.5
    return energy,w


def solve_quantum(cfg,phi,bases,vd,energy_step_ev=0.03,density=True):
    if phi.shape != (cfg.nx,cfg.ny,cfg.nz): raise ValueError('Potential grid shape mismatch')
    kt=KB_EV*cfg.temperature_k
    per_valley={}; total_current=0.0
    total_density=np.zeros_like(phi,float)
    transmissions={}; energies_dict={}
    dy=cfg.width_nm/(cfg.ny+1); dz=cfg.thickness_nm/(cfg.nz+1)
    dx=cfg.length_nm/(cfg.nx-1)
    volume_nm3=dx*dy*dz
    for b in bases:
        H,es,ed=get_onsite(phi,b,cfg.ec0_ev)
        E,weights=energy_grid(cfg,vd,b,es,ed,energy_step_ev)
        if len(E)==0:
            per_valley[b.valley.name]=0.0
            transmissions[b.valley.name]=np.zeros(0)
            energies_dict[b.valley.name]=E
            continue
        T,AL,AR=spectral_batched(H,b.hopping_ev,es,ed,E,cfg.eta_ev,density)
        fl=fermi(E,cfg.efermi_source_ev,kt)
        fr=fermi(E,cfg.efermi_source_ev-vd,kt)
        dI=2*b.valley.multiplicity*Q*Q/H_PLANCK*np.sum(weights*T*(fl-fr))
        total_current+=dI
        per_valley[b.valley.name]=float(dI)
        transmissions[b.valley.name]=T
        energies_dict[b.valley.name]=E
        if density:
            coeff=weights/(2*np.pi)
            RHO=np.einsum('e,exmn->xmn',coeff*fl,AL,optimize=True)
            RHO+=np.einsum('e,exmn->xmn',coeff*fr,AR,optimize=True)
            n=np.einsum('am,xmn,an->xa',b.U,RHO,b.U,optimize=True).real
            total_density+=2*b.valley.multiplicity*n.reshape(phi.shape)/volume_nm3*1e21
    if density and (not np.all(np.isfinite(total_density)) or np.min(total_density)<-1e7):
        raise ArithmeticError('Invalid NEGF electron density')
    return float(total_current),np.maximum(total_density,0),per_valley,energies_dict,transmissions


def min_bandwidth_ev(cfg):
    dx=cfg.length_nm/(cfg.nx-1)
    return min(4*HB2_OVER_2M0/(v.mx*dx*dx) for v in cfg.valleys)
