import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
from scipy.special import expit
from transport.negf_reference import fermi_poles,surface_sigma,rgf,energy_quadrature


def test_continued_fraction_fermi_over_populated_energy_range():
    x=np.linspace(-500,500,1001);pole,weight=fermi_poles(48)
    estimate=.5-2*np.sum(weight[None,:]*x[:,None]/(x[:,None]**2+pole[None,:]**2),axis=1)
    assert_allclose(estimate,expit(-x),atol=3e-13)


def test_retarded_lead_branch_and_asymptote():
    z=np.array([-10.,-.1,.5,1.,10.])+1e-8j
    sigma=surface_sigma(z,np.array([.4]),.2)[:,0]
    assert np.max(sigma.imag)<0
    assert_allclose(sigma,.2**2/(z-.4-sigma),atol=1e-13)
    assert_allclose(surface_sigma(np.array([1e5j]),[.4],.2)[0,0],.2**2/(1e5j-.4),rtol=1e-10)


def test_rgf_diagonals_and_injectivities_match_dense_unequal_blocks():
    rng=np.random.default_rng(13);sizes=[2,3,2,2]
    h=[]
    for m in sizes:
        a=rng.normal(size=(m,m));h.append((a+a.T)*.03+np.eye(m)*.4)
    t=[rng.normal(size=(a,b))*.02 for a,b in zip(sizes,sizes[1:])]
    offsets=np.r_[0,np.cumsum(sizes)];H=np.zeros((sum(sizes),sum(sizes)))
    slices=[slice(offsets[i],offsets[i+1]) for i in range(len(sizes))]
    for s,hh in zip(slices,h): H[s,s]=hh
    for i,tt in enumerate(t): H[slices[i],slices[i+1]]=tt;H[slices[i+1],slices[i]]=tt.T
    for z in [np.array([.2,.4,.6])+1e-7j,np.array([.4+1j,.5+3j])]:
        sl=np.broadcast_to(-.03j*np.eye(2),(len(z),2,2));sr=sl.copy()
        diag,amp=rgf(z,h,t,sl,sr,True,0)
        _,amp_r=rgf(z,h,t,sl,sr,False,1)
        for j,e in enumerate(z):
            M=e*np.eye(len(H))-H.astype(complex);M[slices[0],slices[0]]-=sl[j];M[slices[-1],slices[-1]]-=sr[j]
            G=np.linalg.inv(M)
            for i,s in enumerate(slices):
                assert_allclose(diag[i][j],G[s,s],rtol=1e-12,atol=1e-13)
                assert_allclose(amp[i][j],G[s,slices[0]]*np.sqrt(.06),rtol=1e-11,atol=1e-12)
                assert_allclose(amp_r[i][j],G[s,slices[-1]]*np.sqrt(.06),rtol=1e-11,atol=1e-12)


def test_uniform_infinite_chain_contour_density_matches_k_integration():
    from numpy.polynomial.legendre import leggauss
    t=.2;h=[np.array([[2*t]])]*9;hop=[np.array([[-t]])]*8
    kt=.0258519998;mu=.12;p,w=fermi_poles(48);z=mu+1j*p*kt
    sigma=surface_sigma(z,[2*t],t)[:,:,None]
    diag,_=rgf(z,h,hop,sigma,sigma)
    density=np.array([.5+2*kt*np.dot(w,g[:,0,0].real) for g in diag])
    x,wk=leggauss(512);k=(x+1)*np.pi/2
    truth=np.dot(wk,expit((mu-2*t*(1-np.cos(k)))/kt))/2
    assert_allclose(density,truth,atol=3e-14)


def test_uniform_chain_transmission_one_without_eta():
    e=np.linspace(.001,.799,401);t=.2
    h=[np.array([[2*t]])]*11;hop=[np.array([[-t]])]*10
    sigma=surface_sigma(e,[2*t],t)[:,:,None]
    _,amp=rgf(e,h,hop,sigma,sigma,False,0)
    T=(-2*sigma[:,0,0].imag)*abs(amp[-1][:,0,0])**2
    assert_allclose(T,1,atol=2e-14)


def test_quadrature_splits_band_edges_and_integrates_polynomial():
    x,w=energy_quadrature(-.2,.7,.031,6,[-.12,.01,.26])
    assert_allclose(np.dot(w,x**5),(.7**6-(-.2)**6)/6,rtol=2e-14)


def test_uniform_barrier_matches_analytical_single_site_transmission():
    t=.2;barrier=.1;e=np.linspace(.01,.79,200)
    sigma=surface_sigma(e,[2*t],t)[:,:,None]
    h=[np.array([[2*t]]),np.array([[2*t+barrier]]),np.array([[2*t]])]
    _,amp=rgf(e,h,[np.array([[-t]])]*2,sigma,sigma,False,0)
    T=(-2*sigma[:,0,0].imag)*abs(amp[-1][:,0,0])**2
    k=np.arccos(1-e/(2*t));truth=1/(1+(barrier/(2*t*np.sin(k)))**2)
    assert_allclose(T,truth,rtol=3e-14)


def test_complete_negf_solver_density_and_current_match_dense_open_gf():
    from types import SimpleNamespace
    from scipy.integrate import quad_vec
    from config import Config,VALLEYS,KIN,KB_EV,Q,HBAR_EV_S
    from modespace.hamiltonian import System,transverse_matrix
    from transport.negf_reference import solve_negf,device_blocks,lead_matrices
    c=Config(dx=1.,dy=1.,dz=1.,energy_step=.02,energy_rtol=1e-6)
    valley=VALLEYS[0];nx=7;vd=.05;t=KIN/valley.mx
    transverse=float(transverse_matrix(np.ones((1,1),bool),KIN/valley.my,KIN/valley.mz)[0,0])
    phi=np.linspace(0,vd,nx);H=np.diag(transverse+2*t-phi)+np.diag(np.full(nx-1,-t),1)+np.diag(np.full(nx-1,-t),-1)
    H[0,0]+=t;H[-1,-1]+=t
    grid=SimpleNamespace(x=np.arange(nx),mask=np.ones((nx,1,1),bool),cell_volume=1.)
    system=System(H,[np.ones((1,1,1))]*nx,[np.array([transverse-v]) for v in phi],np.arange(nx+1),np.array([0,0,1,1,1,2,2]),grid,[np.array([[-t]])]*(nx-1),valley)
    ml=transverse+.1;mr=ml-vd
    blocks,bonds,diag,_=solve_negf(system,c,ml,mr,vd)
    h,hop,tx,leads=device_blocks(system,c,(0.,vd));openH=H.copy();openH[0,0]-=t;openH[-1,-1]-=t
    def integrand(e):
        sl,sr=lead_matrices(np.array([e+0j]),tx,leads)
        M=e*np.eye(nx)-openH.astype(complex);M[0,0]-=sl[0,0,0];M[-1,-1]-=sr[0,0,0]
        G=np.linalg.inv(M);injection=np.zeros((nx,nx));injection[0,0]=-2*sl[0,0,0].imag*fermi(e,ml,300);injection[-1,-1]=-2*sr[0,0,0].imag*fermi(e,mr,300)
        return G@injection@G.conj().T/(2*np.pi)
    from multispace.occupation import fermi
    edges=np.r_[transverse-vd,transverse,transverse+4*t-vd,transverse+4*t]
    dense,error=quad_vec(integrand,edges.min(),edges.max(),points=edges,epsabs=1e-9,epsrel=1e-8)
    assert_allclose([b[0,0].real for b in blocks],np.diag(dense).real,atol=2e-8)
    direct=-2*Q/HBAR_EV_S*np.imag(-t*np.diag(dense,-1))
    assert_allclose(bonds,direct,rtol=2e-7,atol=1e-17)
    assert np.linalg.eigvalsh(dense).min()>-1e-9
