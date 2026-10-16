"""Independent numerical identities; these do not certify the full device."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import block_diag,eigh
from scipy.sparse import diags
import pytest
from config import Config,Valley,VALLEYS,KIN,KB_EV,Q,HBAR_EV_S
from modespace.hamiltonian import Grid,ModeBuilder,transverse_matrix,density_from_matrix
from transport.reservoir import reservoirs,solve_steady,bond_currents,contact_currents,neutral_mu
from multispace.occupation import fermi,closed_ms_density,partition_eigenstates
from coupling.mapping import affine_matrix,balanced_matrix
from coupling.mixing import Anderson


def chain(n=9,t=.2):
    H=np.diag(np.ones(n)*2*t)+np.diag(np.ones(n-1)*-t,1)+np.diag(np.ones(n-1)*-t,-1)
    regions=np.r_[np.zeros(3,int),np.ones(n-6,int),np.full(3,2)]
    return H,regions


def test_unchanged_masses():
    assert [(v.mx,v.my,v.mz,v.multiplicity) for v in VALLEYS]==[(.19,.19,.916,2),(.916,.19,.19,2),(.19,.916,.19,2)]


def test_sparse_start_captures_odd_and_even_transverse_modes():
    c=Config(dx=3,modes=2,sd_modes=4);grid=Grid.build(c)
    for v in VALLEYS:
        b=ModeBuilder(c,grid,v);s=b.build(np.zeros(len(grid.xyz)))
        for i in [0,len(grid.x)//2]:
            direct=eigh(b.cache[grid.mask[i].tobytes()].toarray(),eigvals_only=True)
            assert_allclose(s.subbands[i],direct[:len(s.subbands[i])],atol=1e-10)


def test_transverse_cellcentred_box_analytic():
    ny,nz=16,10;w,t=12.,5.;v=VALLEYS[0];dy,dz=w/ny,t/nz
    e=eigh(transverse_matrix(np.ones((ny,nz),bool),KIN/(v.my*dy**2),KIN/(v.mz*dz**2)).toarray(),eigvals_only=True)
    truth=np.pi**2*KIN*(1/(v.my*w*w)+1/(v.mz*t*t))
    assert abs(e[0]/truth-1)<.01
    exact=4*KIN/(v.my*dy**2)*np.sin(np.pi/(2*ny))**2+4*KIN/(v.mz*dz**2)*np.sin(np.pi/(2*nz))**2
    assert_allclose(e[0],exact,atol=1e-12)


def test_full_basis_galerkin_variable_cross_section():
    c=Config(ls=2,ld=2,lg=2,spacer=1,width=2,thickness=1,sd_width=4,sd_thickness=3,tox=.5,
             dx=1,dy=1,dz=1,modes=2,sd_modes=12)
    grid=Grid.build(c);b=ModeBuilder(c,grid,VALLEYS[0]);s=b.build(np.zeros(len(grid.xyz)))
    # Independent 3D real-space finite-volume assembly, all directions at once.
    tx,ty,tz=[KIN/m for m in (VALLEYS[0].mx,VALLEYS[0].my,VALLEYS[0].mz)]
    mask=grid.mask;ids=np.full(mask.shape,-1,int);ids[mask]=np.arange(mask.sum());R=np.zeros((mask.sum(),mask.sum()))
    for ix,iy,iz in zip(*np.where(mask)):
        a=ids[ix,iy,iz]
        for axis,hop in enumerate((tx,ty,tz)):
            for direction in (-1,1):
                pos=[ix,iy,iz];pos[axis]+=direction
                inside=all(0<=p<n for p,n in zip(pos,mask.shape)) and mask[tuple(pos)]
                R[a,a]+=hop if inside else 2*hop
                if inside: R[a,ids[tuple(pos)]]=-hop
    U=s.realspace_transform();assert_allclose(s.H,U.T@R@U,atol=1e-12)
    assert_allclose(eigh(s.H,eigvals_only=True),eigh(R,eigvals_only=True),atol=1e-11)
    assert_allclose(s.H,s.H.T,atol=1e-14)


def test_lyapunov_pauli_continuity_and_polarity():
    H,regions=chain();bath=reservoirs(H,regions,.02,.3,.2,300)
    rho,diag,_=solve_steady(H,bath)
    assert diag['lyapunov_residual']<1e-12
    assert np.linalg.eigvalsh(rho).min()>-1e-12
    assert np.linalg.eigvalsh(rho).max()<1+1e-12
    currents=-2*Q/HBAR_EV_S*np.imag(H[np.arange(8),np.arange(1,9)]*rho[np.arange(1,9),np.arange(8)])
    central=currents[2:6]
    assert central.mean()>0
    assert_allclose(central,central.mean(),rtol=1e-11,atol=1e-18)
    l=np.where(regions==0)[0];r=np.where(regions==2)[0]
    il=Q/HBAR_EV_S*np.trace(bath.injection_l[np.ix_(l,l)]-.02*rho[np.ix_(l,l)]).real
    ir=Q/HBAR_EV_S*np.trace(bath.injection_r[np.ix_(r,r)]-.02*rho[np.ix_(r,r)]).real
    assert_allclose(il,central.mean(),rtol=1e-11,atol=1e-18)
    assert_allclose(il,-ir,rtol=1e-11,atol=1e-18)


def test_equilibrium_zero_current_and_bias_reversal():
    H,regions=chain()
    r0,_,_=solve_steady(H,reservoirs(H,regions,.02,.25,.25,300))
    rp,_,_=solve_steady(H,reservoirs(H,regions,.02,.3,.2,300))
    rm,_,_=solve_steady(H,reservoirs(H,regions,.02,.2,.3,300))
    # Nonlocal mode baths can have reservoir coherences even at equal mu;
    # the UNDAMPED CHANNEL current must vanish, and reverse under mu reversal.
    flux=lambda r:np.imag(H[np.arange(2,6),np.arange(3,7)]*r[np.arange(3,7),np.arange(2,6)])
    assert max(abs(flux(r0)))<1e-13
    assert_allclose(flux(rp),-flux(rm),atol=1e-12)


def test_gauge_invariance():
    H,regions=chain();a=reservoirs(H,regions,.02,.3,.2,300)
    r,_,_=solve_steady(H,a)
    H2=H+.4*np.eye(len(H));r2,_,_=solve_steady(H2,reservoirs(H2,regions,.02,.7,.6,300))
    assert_allclose(r,r2,atol=2e-13)


def test_injection_decomposition_positive_and_additive():
    H,regions=chain();rho,_,parts=solve_steady(H,reservoirs(H,regions,.02,.3,.2,300),split=True)
    assert_allclose(parts[0]+parts[1],rho,atol=1e-13)
    for p in parts: assert np.linalg.eigvalsh(p).min()>-1e-12


def test_ms_standing_wave_is_rejected_and_equilibrium_is_fermi():
    H,regions=chain()
    with pytest.raises(RuntimeError,match='origin unresolved'): closed_ms_density(H,regions,.3,.2,300)
    r,s=closed_ms_density(H,regions,.25,.25,300)
    assert_allclose(np.diag(s['u'].T@r@s['u']),fermi(s['e'],.25,300),atol=1e-13)
    # Real standing-state occupations cannot produce DC bond current.
    assert np.max(abs(r.imag))==0


def test_original_tetra_affine_map():
    vertices=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1.]])
    sites=np.array([[.1,.2,.3],[.05,.05,.05],[.3,.3,.3]])
    P=affine_matrix(vertices,np.array([[0,1,2,3]]),vertices,sites)
    f=lambda p:2+.1*p[:,0]-.2*p[:,1]+.3*p[:,2]
    assert_allclose(P@f(vertices),f(sites),atol=1e-14)
    assert P.data.min()>=0


def test_balanced_charge_constant_positive_and_conservative():
    rng=np.random.default_rng(2);nodes=rng.random((80,3));sites=rng.random((70,3))
    nv=rng.random(80)+.2;qv=rng.random(70)+.2;qv*=nv.sum()/qv.sum()
    C,report=balanced_matrix(nodes,sites,nv,qv)
    assert_allclose(np.asarray(C.sum(axis=1)).ravel(),nv,rtol=3e-9)
    assert_allclose(np.asarray(C.sum(axis=0)).ravel(),qv,rtol=3e-9)
    n=rng.random(70);mapped=C@n/nv
    assert mapped.min()>=0
    assert_allclose(np.dot(mapped,nv),np.dot(n,qv),rtol=1e-12)


def test_anderson_fixedpoint_and_positivity():
    mixer=Anderson(np.ones(4),mix=.3,history=3);x=np.ones(4)*2
    for i in range(25): x,_,_=mixer.update(x,.2*x+.8)
    assert_allclose(x,1,atol=1e-8)
    assert x.min()>=0


def test_quantum_cell_volume_matches_actual_raised_geometry():
    c=Config();g=Grid.build(c)
    truth=(c.ls+c.ld)*c.sd_width*c.sd_thickness+(c.lg+2*c.spacer)*c.width*c.thickness
    assert_allclose(len(g.xyz)*g.cell_volume,truth,rtol=1e-14)


def test_insufficient_alignment_fails():
    with pytest.raises(ValueError,match='align'): Config(dx=2).validate()
