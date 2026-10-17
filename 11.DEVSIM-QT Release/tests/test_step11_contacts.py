import numpy as np
from numpy.testing import assert_allclose
from dataclasses import replace
from config import Config,Valley
from modespace.hamiltonian import Grid
from modespace.junction_basis import JunctionCompatibleBuilder
from transport.contact_model import solve_reference,transverse_laplacian
from transport.negf_reference import lead_matrices,rgf,device_blocks,solve_negf
from validation.analytical import chain


def test_reference_neutrality_gauge_and_actual_rotated_coupling():
    c=Config(geometry_profile='uniform',width=2,thickness=2,sd_width=2,sd_thickness=2,dy=1,dz=1,modes=2,sd_modes=2)
    g=Grid.build(c);v=Valley('test',.19,.19,.916,2);r=solve_reference(c,g,[v])
    assert r.diagnostics['neutrality_relative_error']<1e-10
    assert max(abs(transverse_laplacian(2,2,1,1)@r.phi-1.602176634e-19*1e-14/(8.8541878128e-14*c.eps_si)*(c.donor_sd-r.density)))<1e-9
    b=JunctionCompatibleBuilder(c,g,v);b.contacts=r;b.block_only=True;s=b.build(np.tile(r.phi,len(g.x)))
    ons,hop,t,lead=device_blocks(s,c,(0,0));energy=.8+.02j
    actual=lead_matrices(np.array([energy]),t,lead)[0][0]
    e,u,tx=r.spectra[v.name];U=s.bases[0][g.mask[0]]
    from transport.negf_reference import surface_sigma
    V=-tx*U.T@u;gs=np.diag(surface_sigma([energy],e+2*tx,tx)[0]/tx**2)
    assert_allclose(actual,V@gs@V.T,atol=1e-14)
    assert s.H is None


def test_cut_invariance_rotated_multichannel_density_and_transmission():
    # Same open two-channel Hamiltonian; add ideal layers on either end, then
    # compare the identical central region. Complex z tests spectral density.
    t=.3;eps=np.array([.2,.45]);ang=.41;U=np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang),np.cos(ang)]])
    lead=[(eps,U),(eps,U)];E=np.array([.29,.37,.61]);z=E+.001j
    H=U@np.diag(eps)@U.T;barrier=np.array([[.07,.01],[.01,.02]])
    values=[]
    for padding in (0,1,4):
        ons=[H.copy() for _ in range(5+2*padding)];ons[padding+2]+=barrier
        hop=[-t*np.eye(2) for _ in range(len(ons)-1)];sl,sr=lead_matrices(z,t,lead)
        G,_=rgf(z,ons,hop,sl,sr)
        sl,sr=lead_matrices(E,t,lead);_,A=rgf(E,ons,hop,sl,sr,False,0)
        gam=1j*(sr-sr.conj().transpose(0,2,1));T=np.einsum('eij,ejk,eik->e',gam,A[-1],A[-1].conj()).real
        values.append((np.array(G[padding:padding+5]),np.array(A[padding:padding+5])@np.array(A[padding:padding+5]).conj().transpose(0,1,3,2),T))
    for v in values[1:]:
        for a,b in zip(values[0],v):assert_allclose(a,b,rtol=1e-10,atol=1e-11)


def test_left_right_contour_same_density_without_isolated_states():
    c,s,e=chain(n=9,barrier=.05);mu=e+.14;vd=.04
    left=solve_negf(s,replace(c,density_reference='left'),mu,mu-vd,vd)
    right=solve_negf(s,replace(c,density_reference='right'),mu,mu-vd,vd)
    assert_allclose(left[0],right[0],atol=2e-8,rtol=2e-7)
    assert_allclose(left[1],right[1],atol=1e-13,rtol=1e-7)


def test_rectangular_device_contact_cut_invariance():
    # Two device states coupled to four lead modes. Moving the cut adds a
    # complete ideal lead layer; it must not turn V into a same-size identity.
    t=.24;eps=np.array([.1,.3,.6,.9]);rng=np.random.default_rng(9)
    U=np.linalg.qr(rng.normal(size=(4,4)))[0];W=U[:2];lead=[(eps,W),(eps,W)]
    device=[np.array([[.4,.03],[.03,.5]])]*3;hop=[-.15*np.eye(2)]*2;E=np.array([.23,.44,.67])+.003j
    sl,sr=lead_matrices(E,t,lead);G,_=rgf(E,device,hop,sl,sr)
    complete=[(eps,np.eye(4)),(eps,np.eye(4))];sl2,sr2=lead_matrices(E,t,complete)
    ons=[np.diag(eps),*device,np.diag(eps)];hops=[-t*W.T,*hop,-t*W]
    G2,_=rgf(E,ons,hops,sl2,sr2)
    assert_allclose(np.array(G),np.array(G2[1:-1]),atol=2e-12,rtol=2e-12)


def test_small_full_realspace_sp_contact_density_matches_projected_rgf():
    from validation.reduced_realspace import full_operator
    from modespace.hamiltonian import ModeBuilder,transverse_matrix
    from config import KIN
    from scipy.linalg import eigh
    c=Config(geometry_profile='uniform',ls=2,ld=2,lg=2,spacer=1,width=4,thickness=2,sd_width=4,sd_thickness=2,dx=1,dy=1,dz=1,modes=8,sd_modes=8)
    g=Grid.build(c);v=Valley('small',.19,.19,.916,2);ref=solve_reference(c,g,[v]);b=ModeBuilder(c,g,v);b.contacts=ref;b.block_only=True
    phi=np.tile(ref.phi,len(g.x));phi[(g.xyz[:,0]>3)&(g.xyz[:,0]<5)]-=.03;s=b.build(phi);ons,hop,t,lead=device_blocks(s,c,(0,0))
    H=full_operator(g,c,v,phi);m=g.counts[0];H[:m,:m]-=t*np.eye(m);H[-m:,-m:]-=t*np.eye(m)
    ht=transverse_matrix(g.mask[0],KIN/(v.my*c.dy**2),KIN/(v.mz*c.dz**2)).toarray()-np.diag(ref.phi)+2*t*np.eye(m)
    eps,U=eigh(ht);transform=s.realspace_transform()
    for energy in (.25,.35,.5):
        # Independent dense open real-space lesser density; analytic scalar
        # surface roots reconstructed separately from the block implementation.
        a=energy-eps;root=np.sqrt(a-2*t+0j)*np.sqrt(a+2*t+0j);sig=U@np.diag((a-root)/2)@U.T
        S=np.zeros_like(H,dtype=complex);S[:m,:m]+=sig;S[-m:,-m:]+=sig;G=np.linalg.inv(energy*np.eye(len(H))-H-S)
        gamma=np.zeros_like(S);gamma[:m,:m]=1j*(sig-sig.conj().T)*.7;gamma[-m:,-m:]+=1j*(sig-sig.conj().T)*.2;rho=G@gamma@G.conj().T/(2*np.pi)
        sl,sr=lead_matrices(np.array([energy]),t,lead);_,AL=rgf([energy],ons,hop,sl,sr,False,0);_,AR=rgf([energy],ons,hop,sl,sr,False,1)
        diagonal=np.concatenate([np.diag(s.bases[i][g.mask[i]]@(.7*L[0]@L[0].conj().T+.2*R[0]@R[0].conj().T)@s.bases[i][g.mask[i]].T).real/(2*np.pi) for i,(L,R) in enumerate(zip(AL,AR))])
        assert_allclose(diagonal,np.diag(rho).real,atol=3e-12,rtol=3e-12)
