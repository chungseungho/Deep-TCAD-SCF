"""Translationally invariant EMA wire with insulating transverse electrostatics.

Cell-centred carriers use hard-wall confinement; electrostatics has zero normal
field on the exposed S/D side walls and zero mean gauge. The same full-basis
reference supplies neutrality, surface Green functions and Poisson end data.
This is an ideal infinitely long isolated doped wire, not a metal spill-in model.
"""
from dataclasses import dataclass
import time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.sparse import lil_matrix, bmat, csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.interpolate import RegularGridInterpolator
from numpy.polynomial.legendre import leggauss
from config import KIN,KB_EV,Q
from modespace.hamiltonian import transverse_matrix
from multispace.occupation import fermi

@dataclass
class LeadReference:
    y: np.ndarray
    z: np.ndarray
    phi: np.ndarray
    density: np.ndarray
    mu: float
    spectra: dict
    diagnostics: dict

    def profile(self,xyz):
        yz=np.asarray(xyz)[:,1:3].copy()
        yz[:,0]=np.clip(yz[:,0],self.y[0],self.y[-1]);yz[:,1]=np.clip(yz[:,1],self.z[0],self.z[-1])
        return RegularGridInterpolator((self.y,self.z),self.phi.reshape(len(self.y),len(self.z)),bounds_error=False,fill_value=None)(yz)

    def embedding(self,system,c,lead_phi):
        e,u,tx=self.spectra[system.valley.name];out=[]
        for i,shift in [(0,lead_phi[0]),(-1,lead_phi[1])]:
            device_basis=system.bases[i][system.grid.mask[i]]
            # V=-tx U_device^T U_lead; sigma=V g_surface V^dagger.
            overlap=device_basis.conj().T@u
            out.append((e+2*tx-shift,overlap))
        return out


def transverse_laplacian(ny,nz,dy,dz):
    K=lil_matrix((ny*nz,ny*nz))
    for j in range(ny):
        for k in range(nz):
            a=j*nz+k
            for jj,kk,h in [(j+1,k,dy),(j-1,k,dy),(j,k+1,dz),(j,k-1,dz)]:
                if 0<=jj<ny and 0<=kk<nz:
                    K[a,a]+=1/h**2;K[a,jj*nz+kk]-=1/h**2
    return K.tocsc()


def solve_reference(c,grid,valleys):
    start=time.perf_counter();ny,nz=len(grid.y),len(grid.z);size=ny*nz
    mask=np.ones((ny,nz),bool);K=transverse_laplacian(ny,nz,c.dy,c.dz)
    alpha=Q*1e-14/(8.8541878128e-14*c.eps_si)
    x,w=leggauss(c.lead_k_nodes);k=(x+1)*np.pi/2;weights=w/2
    transverse={v.name:transverse_matrix(mask,KIN/(v.my*c.dy**2),KIN/(v.mz*c.dz**2)).toarray() for v in valleys}
    def density(phi):
        spectra={};terms=[]
        for v in valleys:
            e,u=eigh(transverse[v.name]-np.diag(phi),check_finite=False)
            tx=KIN/(v.mx*c.dx**2);spectra[v.name]=(e,u,tx)
            terms.append((v,e[:,None]+2*tx*(1-np.cos(k)),u))
        target=c.donor_sd*size*c.dx*c.dy*c.dz*1e-21
        def populations(mu):return [(2*v.multiplicity)*(fermi(E,mu,c.temperature)@weights) for v,E,u in terms]
        mu=brentq(lambda m:sum(np.sum(o) for o in populations(m))-target,-5,5,xtol=2e-14)
        n=sum((abs(u)**2@o for (v,E,u),o in zip(terms,populations(mu))),np.zeros(size))/(c.dx*c.dy*c.dz)*1e21
        return n,mu,spectra
    phi=np.zeros(size);flat_n,flat_mu,_=density(phi);history=[]
    for iteration in range(c.lead_max_iterations):
        n,mu,spectra=density(phi);res=alpha*(c.donor_sd-n)-K@phi
        # Approximate compressibility preconditioner; neutrality is solved exactly
        # at every image. Gauge constraint removes the Neumann nullspace.
        chi=alpha*n/(KB_EV*c.temperature)
        M=K+csc_matrix(np.diag(chi));one=csc_matrix(np.ones((size,1)))
        step=spsolve(bmat([[M,one],[one.T,None]],format='csc'),np.r_[res,0.])[:-1]
        history.append({'iteration':iteration,'potential_correction_v':float(max(abs(step))),'poisson_residual_v_nm2':float(max(abs(res))),'mu_ev':float(mu)})
        if max(abs(res))<c.lead_residual_tol:break
        phi+=.7*step;phi-=phi.mean()
    else:raise RuntimeError('NUMERICAL_FAILURE: transverse lead SP did not converge')
    # Flat-profile inconsistency measured by a zero-mean Poisson inverse, not
    # by altering the finite device length or fitting its current.
    one=csc_matrix(np.ones((size,1)))
    flat_defect=spsolve(bmat([[K,one],[one.T,None]],format='csc'),np.r_[alpha*(c.donor_sd-flat_n),0.])[:-1]
    diag={'model':'transverse SP; full transverse EMA basis; Neumann exposed side walls; zero mean phi',
        'converged':True,'iterations':iteration+1,'runtime_s':time.perf_counter()-start,
        'mu_source_ev':float(mu-c.potential_reference_v),'zero_reference_mu_ev':float(mu),
        'neutrality_relative_error':float(abs(n.mean()/c.donor_sd-1)),
        'phi_range_v':float(np.ptp(phi)),'flat_poisson_correction_range_v':float(np.ptp(flat_defect)),
        'flat_vs_sp_density_relative_l2':float(np.linalg.norm(n-flat_n)/np.linalg.norm(n)),
        'flat_mu_ev':float(flat_mu),'flat_vs_sp_mu_change_ev':float(mu-flat_mu),
        'full_transverse_dimension':size,'k_nodes':c.lead_k_nodes,'history':history,
        'gauge':'phi -> phi+V, H -> H-V eV, mu -> mu-V eV; electron charge -q',
        'boundary_assumption':'zero normal electric field on exposed transverse S/D side faces'}
    return LeadReference(grid.y,grid.z,phi,n,float(mu),spectra,diag)
