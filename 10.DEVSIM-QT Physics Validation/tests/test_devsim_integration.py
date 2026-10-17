"""Actual 3D DEVSIM heterogeneous dielectric benchmark; requires native libraries."""
import sys,uuid
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
import pytest
try:
    import gmsh
    import devsim as d
except (ImportError,OSError,RuntimeError) as exc:
    pytest.skip('Gmsh/DEVSIM native dependencies unavailable: '+str(exc),allow_module_level=True)
from poisson.discretization import configure
from geometry.build import scale_mesh


def test_3d_dielectric_interface_potential_and_displacement(tmp_path):
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal',0);gmsh.model.add('slab')
        o=gmsh.model.occ;a=o.addBox(0,0,0,1,2,2);b=o.addBox(1,0,0,1,2,2)
        _,mapping=o.fragment([(3,a)],[(3,b)]);o.synchronize()
        groups={'A':[v for dim,v in mapping[0] if dim==3],'B':[v for dim,v in mapping[1] if dim==3]}
        faces={'Left':[],'Right':[],'Interface':[]}
        for _,s in gmsh.model.getEntities(2):
            up,_=gmsh.model.getAdjacencies(2,s);x,_,_=o.getCenterOfMass(2,s)
            if len(up)==2: faces['Interface'].append(s)
            elif abs(x)<1e-8: faces['Left'].append(s)
            elif abs(x-2)<1e-8: faces['Right'].append(s)
        for dim,gs in [(3,groups),(2,faces)]:
            for name,tags in gs.items():
                g=gmsh.model.addPhysicalGroup(dim,tags);gmsh.model.setPhysicalName(dim,g,name)
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0),.45)
        gmsh.option.setNumber('Mesh.MshFileVersion',2.2);gmsh.model.mesh.generate(3)
        nm=tmp_path/'slab_nm.msh';cm=tmp_path/'slab_cm.msh';gmsh.write(str(nm));scale_mesh(nm,cm)
    finally: gmsh.finalize()
    mesh='slab_'+uuid.uuid4().hex;device=mesh+'_device'
    d.create_gmsh_mesh(mesh=mesh,file=str(cm))
    for r in ['A','B']: d.add_gmsh_region(mesh=mesh,gmsh_name=r,region=r,material=r)
    d.add_gmsh_contact(mesh=mesh,gmsh_name='Left',region='A',name='left',material='Metal')
    d.add_gmsh_contact(mesh=mesh,gmsh_name='Right',region='B',name='right',material='Metal')
    d.add_gmsh_interface(mesh=mesh,gmsh_name='Interface',region0='A',region1='B',name='AB')
    d.finalize_mesh(mesh=mesh);d.create_device(mesh=mesh,device=device)
    try:
        report=configure(d,device,cm,{'A':'A','B':'B'})
        for r,eps in [('A',1.),('B',2.)]:
            args=dict(device=device,region=r)
            d.set_parameter(**args,name='epsilon',value=eps)
            d.node_solution(**args,name='Potential');d.edge_from_node_model(**args,node_model='Potential')
            for name,eq in [('Flux','epsilon*(Potential@n0-Potential@n1)*EdgeInverseLength'),
                            ('Flux:Potential@n0','epsilon*EdgeInverseLength'),('Flux:Potential@n1','-epsilon*EdgeInverseLength')]:
                d.edge_model(**args,name=name,equation=eq)
            d.equation(**args,name='Poisson',variable_name='Potential',edge_model='Flux')
        args=dict(device=device,interface='AB')
        for name,eq in [('continuity','Potential@r0-Potential@r1'),('continuity:Potential@r0','1'),('continuity:Potential@r1','-1')]:
            d.interface_model(**args,name=name,equation=eq)
        d.interface_equation(**args,name='Poisson',interface_model='continuity',type='continuous')
        for contact,bias in [('left',0),('right',1)]:
            args=dict(device=device,contact=contact)
            d.contact_node_model(**args,name='bc',equation=f'Potential-{bias}')
            d.contact_node_model(**args,name='bc:Potential',equation='1')
            d.contact_equation(**args,name='Poisson',node_model='bc',edge_charge_model='Flux')
        d.solve(type='dc',absolute_error=1e-12,relative_error=1e-8,maximum_iterations=10)
        slopes=[]
        for r in ['A','B']:
            x=np.asarray(d.get_node_model_values(device=device,region=r,name='x'))*1e7
            phi=np.asarray(d.get_node_model_values(device=device,region=r,name='Potential'))
            truth=np.where(x<=1,2*x/3,2/3+(x-1)/3)
            assert_allclose(phi,truth,atol=3e-13)
            slopes.append(np.polyfit(x,phi,1)[0])
            assert_allclose(report[r]['physical_volume_nm3'],4.,rtol=1e-12)
        assert_allclose(slopes[0],2*slopes[1],rtol=1e-12)
        left=d.get_contact_charge(device=device,contact='left',equation='Poisson')
        right=d.get_contact_charge(device=device,contact='right',equation='Poisson')
        assert_allclose(left,-right,rtol=1e-12)
    finally:
        d.delete_device(device=device);d.delete_mesh(mesh=mesh)
