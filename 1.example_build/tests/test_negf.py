"""No Gmsh/DEVSIM needed for these elementary transport sanity checks."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config, Valley
from negf.mode_space import spectral_rgf, evaluate_negf, surface_self_energy


def test_uniform_1d_perfect_chain():
    N=19
    t=.2
    onsite=np.full((N,1,1), 2*t, complex)
    hop=np.full((N-1,1,1),-t,complex)
    val, *_=spectral_rgf(onsite,hop,t,np.array([0.0]),np.array([0.0]),
                         .4,1e-8,with_density=False)
    assert abs(val-1)<1e-5, f'perfect channel transmission {val}'


def test_surface_self_energy_branch():
    for energy in [-1,0,.1,.4,.8,1.2]:
        s=surface_self_energy(energy,np.array([0.0]),.2,1e-6)[0]
        assert s.imag<=1e-12, (energy,s)


def test_mode_space_3d_density_and_zero_current():
    cfg=Config(nx=9,ny=4,nz=3,modes=2,energy_points=70,
               valleys=(Valley('test',.19,.19,.916,2),))
    cfg.validate()
    potential=np.zeros((cfg.nx,cfg.ny,cfg.nz))
    res=evaluate_negf(cfg,potential,vg=0,vd=0)
    assert abs(res.current_a)<1e-15
    assert np.isfinite(res.electron_density_cm3).all()
    assert np.min(res.electron_density_cm3)>-1e-2
    assert np.max(res.transmission)>0


if __name__=='__main__':
    test_uniform_1d_perfect_chain()
    test_surface_self_energy_branch()
    test_mode_space_3d_density_and_zero_current()
    print('PASS: 1D unit transmission, retarded self-energy, transverse modes, I(V=0)=0')
