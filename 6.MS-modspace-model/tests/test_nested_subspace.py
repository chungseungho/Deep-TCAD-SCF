import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
from config import Config,VALLEYS
from modespace.hamiltonian import Grid,ModeBuilder
from validation.junction_compatible import JunctionCompatibleBuilder


def test_enriched_raised_basis_contains_narrow_trace_without_hopping_fit():
    c=Config(dx=3.,modes=3,sd_modes=6);g=Grid.build(c)
    for v in VALLEYS:
        s=JunctionCompatibleBuilder(c,g,v).build(np.zeros(len(g.xyz)));i=np.where(g.x<c.ls)[0][-1]
        a=s.bases[i].reshape(-1,s.bases[i].shape[-1]);b=s.bases[i+1].reshape(-1,s.bases[i+1].shape[-1])
        assert_allclose((a.T@b).T@(a.T@b),np.eye(c.modes),atol=2e-14)
        assert_allclose(s.H,s.H.T,atol=1e-14)
        assert s.bases[0].shape[-1]>c.sd_modes


def test_full_basis_enrichment_does_not_change_real_space_operator():
    c=Config(ls=2,ld=2,lg=2,spacer=1,width=2,thickness=1,sd_width=4,sd_thickness=3,tox=.5,
             dx=1,dy=1,dz=1,modes=2,sd_modes=12)
    g=Grid.build(c);phi=np.sin(g.xyz[:,0])*.01
    a=ModeBuilder(c,g,VALLEYS[0]).build(phi);b=JunctionCompatibleBuilder(c,g,VALLEYS[0]).build(phi)
    Ua=a.realspace_transform();Ub=b.realspace_transform()
    assert_allclose(Ua@a.H@Ua.T,Ub@b.H@Ub.T,atol=5e-14)
