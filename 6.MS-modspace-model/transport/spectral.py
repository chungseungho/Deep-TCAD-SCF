"""POST-SCF finite-reservoir spectral diagnostic, no semi-infinite surface GF.

Lindblad current uses f(e_k), not f(E); T(E)*(fL(E)-fR(E)) is NOT substituted.
The separately integrated matrix injection kernel must agree with the Lyapunov
bond current upon energy-window/quadrature convergence.
"""
import numpy as np
from scipy.linalg import solve,eig
from config import Q,HBAR_EV_S


def finite_transmission(H,regions,gamma,energies):
    """Post-SCF G_LR evaluated by complex-pole expansion, in small energy batches.

    Exact for the specified finite H and constant gamma projectors. This is a
    GF postprocessor, not an extra density approximation or semi-infinite lead.
    """
    G=np.diag((regions!=1).astype(float)*gamma)
    z,V=eig(H-.5j*G,check_finite=False)
    inverse=solve(V,np.eye(len(H)),assume_a='gen',check_finite=False)
    L=V[regions==0];R=inverse[:,regions==2]
    energies=np.asarray(energies);out=np.empty(len(energies))
    for i in range(0,len(energies),64):
        e=energies[i:i+64]
        glr=(L[None,:,:]/(e[:,None,None]-z[None,None,:]))@R
        out[i:i+len(e)]=gamma**2*np.sum(abs(glr)**2,axis=(1,2))
    return out


def spectrum(H,baths,energies,degeneracy,backend='lindblad',mu_l=None,mu_r=None,temperature=300):
    if backend in ('spectral','ms'):
        from multispace.occupation import fermi
        gamma=float(np.max(np.diag(baths.gamma)))
        T=finite_transmission(H,baths.regions,gamma,energies)
        J=degeneracy*Q/HBAR_EV_S/(2*np.pi)*T*(fermi(energies,mu_l,temperature)-fermi(energies,mu_r,temperature))
        return T,J
    n=len(H);ind_l=np.where(np.diag(baths.gamma)>0)[0]
    # Injection supports identify reservoir projectors even for vanishing Fermi factors.
    # Caller supplies L/R partition through the zero cross-block structure of Pi.
    # Explicit gamma components are derived from indices in the passed attributes.
    if not hasattr(baths,'regions'): raise ValueError('Spectral analysis needs L/R region labels')
    l=np.where(baths.regions==0)[0];r=np.where(baths.regions==2)[0]
    GL=np.zeros_like(H);GR=np.zeros_like(H)
    GL[l,l]=np.diag(baths.gamma)[l];GR[r,r]=np.diag(baths.gamma)[r]
    T=[];J=[]
    for e in energies:
        G=solve(e*np.eye(n)-H+.5j*baths.gamma,np.eye(n),assume_a='gen',check_finite=False)
        T.append(float(np.trace(GR@G@GL@G.conj().T).real))
        if backend=='spectral':
            from multispace.occupation import fermi
            kernel=T[-1]*(fermi(e,mu_l,temperature)-fermi(e,mu_r,temperature))
        else:
            kernel=np.trace(GR@G@baths.injection_l@G.conj().T-GL@G@baths.injection_r@G.conj().T).real
        J.append(degeneracy*Q/HBAR_EV_S/(2*np.pi)*kernel)
    return np.array(T),np.array(J)
