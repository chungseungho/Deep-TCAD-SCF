"""Semi-infinite coherent leads, independent of finite-bath digamma filling.

Real-energy RGF injectivity + continued-fraction Fermi contour equilibrium.
No artificial eta, scattering or fitted broadening. The current path and all
interior lesser diagonal blocks are obtained from the same lead injectivity.
References: Ozaki, PRB 75,035123 (2007); Ozaki et al. PRB 81,035116 (2010).
This backend DOES use GF evaluations inside SCF; it is not the GF-free route.
"""
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh_tridiagonal,eigh
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
from config import KIN,KB_EV,Q,HBAR_EV_S
from multispace.occupation import fermi


@lru_cache(maxsize=16)
def fermi_poles(n):
    j=np.arange(1,2*n)
    b=1/(2*np.sqrt((2*j-1)*(2*j+1)))
    eigen,u=eigh_tridiagonal(np.zeros(2*n),b)
    keep=eigen>0;lam=eigen[keep]
    poles=1/lam;weight=u[0,keep]**2/(4*lam**2)
    order=np.argsort(poles)
    return poles[order],weight[order]


def surface_sigma(z,onsite,hopping):
    """Analytic retarded t^2*g_s for a uniform 1D FD channel (eV)."""
    a=np.asarray(z,dtype=complex)[...,None]-np.asarray(onsite)
    root=np.sqrt(a-2*hopping+0j)*np.sqrt(a+2*hopping+0j)
    return 2*hopping**2/(a+root)


def device_blocks(system,c,lead_phi=None):
    if lead_phi is None:lead_phi=(c.potential_reference_v,c.potential_reference_v)
    tx=KIN/(system.valley.mx*c.dx**2)
    onsite=[system.H[system.sl(i),system.sl(i)].copy() for i in range(len(system.grid.x))]
    # Replace the two closed half-cell x hard walls by continuation into leads.
    onsite[0]-=tx*np.eye(len(onsite[0]));onsite[-1]-=tx*np.eye(len(onsite[-1]))
    from modespace.hamiltonian import transverse_matrix
    lead=[]
    for i,phi in [(0,lead_phi[0]),(-1,lead_phi[1])]:
        mask=system.grid.mask[i];u=system.bases[i][mask]
        ht=transverse_matrix(mask,KIN/(system.valley.my*c.dy**2),KIN/(system.valley.mz*c.dz**2))
        e,w=eigh(u.T@ht@u-phi*np.eye(u.shape[1])+2*tx*np.eye(u.shape[1]))
        lead.append((e,w))
    return onsite,system.hop,tx,lead


def lead_matrices(z,tx,lead):
    out=[]
    for e,u in lead:
        s=surface_sigma(z,e,tx)
        out.append((u[None,:,:]*s[:,None,:])@u.T)
    return out


def rgf(z,onsite,hop,sigma_l,sigma_r,diagonal=True,injection=None):
    """Batched block recursion. No reuse of the finite-spectral eigensystem.

    injection=None: diagonal retarded G. 0/1: G_iL sqrt(Gamma_L) or G_iR
    sqrt(Gamma_R). Only a lead's diagonal gamma square root is diagonalized.
    """
    z=np.atleast_1d(z).astype(complex);n=len(onsite)
    if injection==1:
        reversed_hop=[h.conj().T for h in hop[::-1]]
        ans=rgf(z,onsite[::-1],reversed_hop,sigma_r,sigma_l,diagonal,0)
        if diagonal: return ans[0][::-1],ans[1][::-1]
        return None,ans[1][::-1]
    right=[None]*n
    def invblock(i,extra):
        return np.linalg.inv(z[:,None,None]*np.eye(len(onsite[i]))-onsite[i]-extra)
    right[-1]=invblock(n-1,sigma_r)
    for i in range(n-2,-1,-1): right[i]=invblock(i,hop[i]@right[i+1]@hop[i].conj().T)
    diag=None
    if diagonal:
        diag=[];left=sigma_l
        for i in range(n):
            sr=sigma_r if i==n-1 else hop[i]@right[i+1]@hop[i].conj().T
            diag.append(invblock(i,left+sr))
            if i<n-1:
                gl=invblock(i,left);left=hop[i].conj().T@gl@hop[i]
    amplitudes=None
    if injection is not None:
        gam=1j*(sigma_l-sigma_l.conj().transpose(0,2,1))
        e,u=np.linalg.eigh(gam)
        if e.min()<-1e-10: raise RuntimeError('Noncausal lead spectral function')
        sqrtgam=(u*np.sqrt(np.maximum(e,0))[:,None,:])@u.conj().transpose(0,2,1)
        g0=diag[0] if diag is not None else invblock(0,sigma_l+hop[0]@right[1]@hop[0].conj().T)
        amplitudes=[g0@sqrtgam]
        for i in range(1,n): amplitudes.append(right[i]@hop[i-1].conj().T@amplitudes[-1])
    return diag,amplitudes


def energy_quadrature(low,high,panel,order,edges=()):
    breaks=np.r_[np.arange(low,high,panel),high]
    breaks=np.unique(np.r_[breaks,[v for v in edges if low<v<high]])
    x,w=leggauss(order);half=np.diff(breaks)/2;middle=(breaks[1:]+breaks[:-1])/2
    return (middle[:,None]+half[:,None]*x).ravel(),(half[:,None]*w).ravel()


def solve_negf(system,c,mu_l,mu_r,vd,detail=False):
    ref=c.potential_reference_v
    onsite,hop,tx,leads=device_blocks(system,c,(ref,ref+vd))
    kbt=KB_EV*c.temperature
    from transport.bound_states import audit
    bound_audit=audit(system,c,mu_l,vd) if abs(vd)>1e-14 else None
    if bound_audit is not None and bound_audit['bias_dependent_bound_occupation_unresolved']:
        from transport.bound_states import BoundOccupationError
        raise BoundOccupationError(bound_audit)
    # Verify the rational Fermi approximation over a conservative H range.
    H=system.H;radius=np.sum(abs(H),axis=1)-abs(np.diag(H))
    extremum=max(abs(np.diag(H)-radius-min(mu_l,mu_r)).max(),abs(np.diag(H)+radius-min(mu_l,mu_r)).max())/kbt
    sample=np.linspace(-extremum,extremum,2001);used=c.fermi_poles
    from scipy.special import expit
    while True:
        pole,weight=fermi_poles(used)
        approximation=.5-2*np.sum(weight[None,:]*sample[:,None]/(sample[:,None]**2+pole[None,:]**2),axis=1)
        fermi_error=float(max(abs(approximation-expit(-sample))))
        if fermi_error<min(1e-10,c.energy_atol*.01):break
        used*=2
        if used>768:raise RuntimeError('Fermi contour could not resolve Hamiltonian energy range')
    base=min(mu_l,mu_r);z=base+1j*pole*kbt
    # Each pole contributes linearly. Batching changes memory, not the contour,
    # contact physics, number of poles, or integration tolerance.
    equilibrium=[.5*np.eye(len(h)) for h in onsite]
    for start in range(0,len(z),c.contour_batch_size):
        end=start+c.contour_batch_size
        sl,sr=lead_matrices(z[start:end],tx,leads)
        diagonal,_=rgf(z[start:end],onsite,hop,sl,sr)
        for i,g in enumerate(diagonal):
            equilibrium[i]+=2*kbt*np.einsum('e,eij->ij',weight[start:end],g.real)
        del diagonal

    excess=[np.zeros_like(h,dtype=complex) for h in onsite]
    bonds=np.zeros(len(hop));current_integral=0.;points=0
    # Error from tails is separately controlled by fermi_tail, not by eta.
    if abs(mu_l-mu_r)>1e-14:
        origin=0 if mu_l>mu_r else 1
        low=base-c.fermi_tail*kbt;high=max(mu_l,mu_r)+c.fermi_tail*kbt
        edges=np.concatenate([np.r_[e-2*tx,e+2*tx] for e,u in leads])
        from transport.resonances import resonance_breaks
        breaks,pole_diagnostic=resonance_breaks(onsite,hop,tx,leads,low,high,c.energy_step)
        sizes=[len(h)**2 for h in onsite];offsets=np.r_[0,np.cumsum(sizes)]
        def integrand(e):
            sl,sr=lead_matrices(e.astype(complex),tx,leads)
            _,amp=rgf(e,onsite,hop,sl,sr,False,origin)
            df=abs(fermi(e,mu_l,c.temperature)-fermi(e,mu_r,c.temperature))/(2*np.pi)
            pieces=[(A@A.conj().transpose(0,2,1)).reshape(len(e),-1) for A in amp]
            flux=[]
            for i,h in enumerate(hop):
                flux.append(-2*np.einsum('ab,ebk,eak->e',h,amp[i+1],amp[i].conj(),optimize=True).imag)
            pieces.append(np.stack(flux,axis=1))
            target=sr if origin==0 else sl;A=amp[-1] if origin==0 else amp[0]
            gamma=1j*(target-target.conj().transpose(0,2,1))
            T=np.einsum('eij,ejk,eik->e',gamma,A,A.conj(),optimize=True).real
            pieces.append(((1 if origin==0 else -1)*T)[:,None])
            return np.concatenate(pieces,axis=1)*df[:,None]
        from transport.quadrature import integrate_vector
        absolute=np.r_[np.full(offsets[-1],c.energy_atol),np.full(len(hop)+1,2e-15)]
        integral,integration=integrate_vector(integrand,breaks,c.energy_order,c.energy_rtol,absolute,c.energy_max_nodes)
        points=integration['gf_nodes']
        for i,h in enumerate(onsite): excess[i]=integral[offsets[i]:offsets[i+1]].reshape(h.shape)
        bonds=integral[offsets[-1]:offsets[-1]+len(hop)].real*Q/HBAR_EV_S
        current_integral=float(integral[-1].real*Q/HBAR_EV_S)
    blocks=[(eq+ex+eq.T+ex.conj().T)/2 for eq,ex in zip(equilibrium,excess)]
    eigen=[np.linalg.eigvalsh(r) for r in blocks]
    minimum=min(e.min() for e in eigen);maximum=max(e.max() for e in eigen)
    # These are necessary local Pauli tests, not a full-matrix PSD proof.
    if minimum< -2e-7 or maximum>1+2e-7:
        raise RuntimeError(f'NEGF local Pauli bound failed: {minimum:g},{maximum:g}; refine Fermi/energy quadrature')
    scale=max(abs(current_integral),1e-15)
    diagnostics={'density_algorithm':'semi-infinite GF contour + real-axis injectivity',
                 'gf_inside_scf':True,'contour_batch_size':c.contour_batch_size,'fermi_poles':used,'requested_fermi_poles':c.fermi_poles,
                 'fermi_range_max_error':fermi_error,'bound_state_audit':bound_audit,'real_energy_nodes':points,
                 'energy_panel_ev':c.energy_step,'gauss_order':c.energy_order,
                 'integration_error_estimator':integration if points else None,
                 'pole_guidance':pole_diagnostic if points else None,
                 'local_pauli_min':float(minimum),'local_pauli_max':float(maximum),
                 'full_density_matrix_psd_checked':False,
                 'independent_transmission_current_per_degeneracy_a':float(current_integral),
                 'bond_vs_transmission_relative':float(max(abs(bonds-current_integral))/scale),
                 'bound_state_occupation':'equilibrium base; below-lead poles with material reservoir occupation difference rejected; embedded dark states not globally excluded',
                 'contact_model':'analytic semi-infinite uniform hard-wall nanowire leads; no phenomenological gamma'}
    return blocks,bonds,diagnostics,{'equilibrium':equilibrium,'excess':excess}


def blocks_density(system,blocks,degeneracy):
    out=[]
    for i,mask in enumerate(system.grid.mask):
        u=system.bases[i][mask]
        n=np.einsum('im,mn,in->i',u,blocks[i],u.conj()).real
        if n.min()<-2e-7: raise RuntimeError('Negative real-space diagonal')
        out.append(np.maximum(n,0)*degeneracy/system.grid.cell_volume*1e21)
    return np.concatenate(out)


def neutral_lead_mu(builders,c,grid):
    spectra=[];x,w=leggauss(256);k=(x+1)*np.pi/2;wk=w/2
    for builder in builders:
        system=builder.build(np.zeros(len(grid.xyz)))
        _,_,tx,leads=device_blocks(system,c,(0.,0.))
        sub=leads[0][0]-2*tx
        spectra.append((sub[:,None]+2*tx*(1-np.cos(k)),2*builder.v.multiplicity))
    target=c.donor_sd*c.sd_width*c.sd_thickness*1e-21
    def number(mu):
        return sum(g*np.dot(fermi(e,mu,c.temperature),wk).sum()/c.dx for e,g in spectra)
    mu=brentq(lambda mu:number(mu)-target,-2.,2.,xtol=2e-14)
    physical_mu=float(mu-c.potential_reference_v)
    return physical_mu,{'mu_source_ev':physical_mu,'zero_reference_mu_ev':float(mu),'potential_reference_v':c.potential_reference_v,'donor_electrons_per_nm':float(target),
                     'lead_electrons_per_nm':float(number(mu)),
                     'neutrality_relative_error':float(abs(number(mu)/target-1)),
                     'calibration':'infinite nanowire FD dispersion integral at prescribed donor line density; no current fit',
                     'k_quadrature_nodes':256,'connected_contact_neutrality':'separate finite-device diagnostic'}


def direct_finite_gf(H,regions,gamma,mu_l,mu_r,temperature):
    """Small-system independent real-axis integration of the finite-bath model."""
    from scipy.integrate import quad_vec
    gl=np.diag((regions==0)*gamma);gr=np.diag((regions==2)*gamma);eye=np.eye(len(H))
    def integrand(e):
        G=np.linalg.solve(e*eye-H+.5j*(gl+gr),eye)
        return G@(gl*fermi(e,mu_l,temperature)+gr*fermi(e,mu_r,temperature))@G.conj().T/(2*np.pi)
    return quad_vec(integrand,-np.inf,np.inf,epsabs=2e-11,epsrel=2e-10)
