"""Reduced-device numerical golden regression against preserved STEP08."""
from pathlib import Path
import numpy as np


def test_all_builtin_backends_match_step08_quantum_outputs():
    from config import Config
    from device.definition import Frame
    from modespace.hamiltonian import Grid,ModeBuilder
    from transport.registry import create_transport_backend

    root=Path(__file__).parents[1]
    reference=np.load(root/'validation/modularization_reference.npz')
    for name in ('ms','negf','spectral'):
        c=Config(backend=name,ls=4.,lg=4.,ld=4.,spacer=2.,width=4.,thickness=4.,
            sd_width=8.,sd_thickness=8.,tox=2.,dx=1.,dy=1.,dz=1.,modes=2,sd_modes=2,
            ms_boundary_shift_slices=0)
        grid=Grid.build(c)
        builders=[ModeBuilder(c,grid,v) for v in Frame(c.transport_axis).valleys()]
        backend=create_transport_backend(name,c)
        phi=np.zeros(len(grid.xyz))
        result=backend.solve(phi,0.,{'config':c,'builders':builders,
            'grid':grid,'fixed_mu':-.1,'detail':True,'bias':0.})
        backend.finalize(result,0.,grid,True)
        for key,value in (('n',result['n']),('bonds',result['bonds']),
                          ('mu',np.asarray(result['mu'])),('tracking',result['tracking']),
                          ('current',np.asarray(result['id_a'])),('post_bonds',result['bonds']),
                          ('phi',phi)):
            np.testing.assert_allclose(value,reference[f'{name}_{key}'],rtol=1e-12,atol=1e-14,
                err_msg=f'{name} {key} differs from STEP08 baseline')
