"""Runnable on development machine, no gmsh or devsim required.

These are controlled physics unit tests, not a claim of a calibrated GAA IV.
"""
import sys
import unittest
from dataclasses import replace
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parent/'1.example_build'))
from config import Config
from quantum import (make_bases,transverse_sparse,get_onsite,
                     spectral_batched,solve_quantum,energy_grid,min_bandwidth_ev)
from analysis import voltages,extract_vth,choose_vgs,extract_ss


class TestPhysics(unittest.TestCase):
    def cfg(self):
        c=replace(Config(),nx=21,ny=9,nz=7,modes=3)
        c.valleys=c.valleys[:1]
        return c

    def test_transverse_orthogonality(self):
        c=self.cfg();b=make_bases(c)[0]
        self.assertLess(np.max(np.abs(b.U.T@b.U-np.eye(c.modes))),1e-10)
        self.assertLess(np.linalg.norm(transverse_sparse(c,b.valley)@b.U-b.U*b.levels),1e-7)

    def test_uniform_channel_transmission(self):
        n=11;m=3;t=.3
        H=np.tile(np.diag(np.array([.03,.12,.2])+2*t),(n,1,1)).astype(complex)
        E=np.array([.37,.48,.59])
        T,al,ar=spectral_batched(H,t,np.array([.03,.12,.2]),np.array([.03,.12,.2]),E,1e-9)
        self.assertTrue(np.allclose(T,3,atol=1e-5))
        self.assertTrue(np.all(np.linalg.eigvalsh(al[:,n//2]).real>-1e-8))
        self.assertTrue(np.all(np.linalg.eigvalsh(ar[:,n//2]).real>-1e-8))

    def test_zero_bias_current(self):
        c=self.cfg();phi=np.zeros((c.nx,c.ny,c.nz))
        I,n,_,_,_=solve_quantum(c,phi,make_bases(c),0,0.03)
        self.assertLess(abs(I),1e-17)
        self.assertTrue(np.isfinite(n).all())
        self.assertTrue((n>=0).all())

    def test_gate_barrier_suppresses_current(self):
        c=replace(self.cfg(),nx=61)
        x=np.linspace(0,c.length_nm,c.nx)
        w=np.where((x>c.ls_nm)&(x<c.ls_nm+c.lg_nm),1.,0.)
        b=make_bases(c)
        phioff=(-.6*w)[:,None,None]*np.ones((1,c.ny,c.nz))
        phion=(.3*w)[:,None,None]*np.ones((1,c.ny,c.nz))
        ioff,*_=solve_quantum(c,phioff,b,.1,.01,False)
        ion,*_=solve_quantum(c,phion,b,.1,.01,False)
        self.assertGreater(ion,ioff*10)

    def test_hermitian_projected_hamiltonian(self):
        c=self.cfg();b=make_bases(c)[0]
        phi=np.random.default_rng(1).normal(size=(c.nx,c.ny,c.nz))*.1
        H,_,_=get_onsite(phi,b,c.ec0_ev)
        self.assertLess(np.max(np.abs(H-H.conj().swapaxes(-1,-2))),1e-12)

    def test_biased_density_energy_union(self):
        c=replace(self.cfg(),nx=345)
        b=make_bases(c)[0]
        phi=np.zeros((c.nx,c.ny,c.nz));phi[-1]=6.
        _,s,d=get_onsite(phi,b,c.ec0_ev)
        E,_=energy_grid(c,6,b,s,d,.02)
        # Drain continuum extends 6eV below the source: must integrate it.
        self.assertLess(E[0],-5.5)
        self.assertGreater(E[-1],0.3)

    def test_biased_dx_bandwidth_guard(self):
        c=replace(self.cfg(),nx=345)
        self.assertGreater(min_bandwidth_ev(c),7.0)

    def test_voltage_grid_vth_triplet(self):
        v=voltages();self.assertEqual(len(v),31)
        self.assertEqual(v[-1],6)
        x=voltages(2,.2);I=1e-12*10**(4*x)
        t=extract_vth(x,I,12)
        self.assertEqual(t['status'],'OK')
        a=choose_vgs(t['vth_v'],.2,6.)
        self.assertLess(a[0],t['vth_v']);self.assertGreater(a[1],t['vth_v'])
        self.assertGreater(a[2],a[1])
        self.assertIsNone(extract_vth(x,np.ones_like(x)*1e-3,12)['vth_v'])


if __name__=='__main__': unittest.main()
