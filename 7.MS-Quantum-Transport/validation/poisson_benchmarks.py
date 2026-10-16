"""Independent shared-DOF DEVSIM capacitors and charged-slab mesh refinement."""
import json,sys,uuid
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from config import Q
from geometry.build import scale_mesh
from poisson.device import merged_mesh,EPS0
from poisson.discretization import mesh_regions,assemble


def slab(directory,h,eps_a=11.7,eps_b=3.9,density=0.):
    import gmsh,devsim as d
    out=Path(directory);out.mkdir(parents=True,exist_ok=True)
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal',0);gmsh.model.add('capacitor')
        o=gmsh.model.occ;a=o.addBox(0,0,0,1,2,2);b=o.addBox(1,0,0,1,2,2)
        _,mapping=o.fragment([(3,a)],[(3,b)]);o.synchronize()
        for label,items in zip(['A','B'],mapping):
            tag=gmsh.model.addPhysicalGroup(3,[t for dim,t in items if dim==3]);gmsh.model.setPhysicalName(3,tag,label)
        for label,side in [('Left',0),('Right',2)]:
            faces=[s for dim,s in gmsh.model.getEntities(2) if abs(o.getCenterOfMass(2,s)[0]-side)<1e-8]
            tag=gmsh.model.addPhysicalGroup(2,faces);gmsh.model.setPhysicalName(2,tag,label)
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0),h);gmsh.option.setNumber('Mesh.MshFileVersion',2.2)
        gmsh.model.mesh.generate(3);gmsh.write(str(out/'slab_nm.msh'));scale_mesh(out/'slab_nm.msh',out/'slab_cm.msh')
    finally: gmsh.finalize()
    mesh='slab_'+uuid.uuid4().hex;device=mesh+'_device';path=out/'slab_cm.msh'
    d.create_gmsh_mesh(mesh=mesh,file=str(merged_mesh(path)))
    d.add_gmsh_region(mesh=mesh,gmsh_name='Electrostatic',region='all',material='ExplicitTetra')
    for label in ['Left','Right']: d.add_gmsh_contact(mesh=mesh,gmsh_name=label,region='all',name=label.lower(),material='Metal')
    d.finalize_mesh(mesh=mesh);d.create_device(mesh=mesh,device=device);args=dict(device=device,region='all')
    try:
        xyz=np.column_stack([d.get_node_model_values(**args,name=a) for a in ['x','y','z']]);points,tets=mesh_regions(path)
        ids=cKDTree(xyz).query(points)[1];K=csr_matrix((len(xyz),len(xyz)));vol=np.zeros(len(xyz))
        for name,eps in [('A',eps_a),('B',eps_b)]:
            k,v=assemble(xyz,ids[tets[name]]);K+=k*eps;vol+=v
        for name in ['Potential','Volume']:d.node_solution(**args,name=name)
        d.set_node_values(**args,name='Volume',values=vol.tolist());d.set_parameter(**args,name='node_volume_model',value='Volume')
        d.edge_from_node_model(**args,node_model='node_index');i=np.array(d.get_edge_model_values(**args,name='node_index@n0'),int);j=np.array(d.get_edge_model_values(**args,name='node_index@n1'),int)
        coupling=-np.asarray(K[i,j]).ravel()*np.linalg.norm(xyz[i]-xyz[j],axis=1)
        d.edge_solution(**args,name='Coupling');d.set_edge_values(**args,name='Coupling',values=coupling.tolist())
        d.set_parameter(**args,name='edge_couple_model',value='Coupling');d.set_parameter(**args,name='eps0',value=EPS0)
        d.edge_from_node_model(**args,node_model='Potential')
        for name,eq in [('Flux','eps0*(Potential@n0-Potential@n1)*EdgeInverseLength'),('Flux:Potential@n0','eps0*EdgeInverseLength'),('Flux:Potential@n1','-eps0*EdgeInverseLength')]:d.edge_model(**args,name=name,equation=eq)
        d.node_model(**args,name='Rho',equation=str(-Q*density));d.equation(**args,name='Poisson',variable_name='Potential',edge_model='Flux',node_model='Rho')
        for label,bias in [('left',0),('right',0 if density else 1)]:
            d.contact_node_model(device=device,contact=label,name=label+'_bc',equation=f'Potential-{bias}')
            d.contact_node_model(device=device,contact=label,name=label+'_bc:Potential',equation='1')
            d.contact_equation(device=device,contact=label,name='Poisson',node_model=label+'_bc',edge_charge_model='Flux')
        d.solve(type='dc',absolute_error=1e-13,relative_error=1e-8,maximum_iterations=10)
        phi=np.array(d.get_node_model_values(**args,name='Potential'));x=xyz[:,0];length=2e-7
        if density:
            if eps_a!=eps_b:raise ValueError('Charged analytical slab uses uniform epsilon')
            truth=Q*density*x*(length-x)/(2*eps_a*EPS0)
        else:
            slope=1/(1e-7/eps_a+1e-7/eps_b)
            truth=np.where(x<=1e-7,slope*x/eps_a,slope*(1e-7/eps_a+(x-1e-7)/eps_b))
        boundary=(abs(x)<1e-14)|(abs(x-length)<1e-14)
        flux=(K*EPS0)@phi;rho=Q*density*vol
        scale=max(abs(flux[boundary]).sum(),abs(rho).sum(),1e-25)
        charge={label:float(d.get_contact_charge(device=device,contact=label,equation='Poisson')) for label in ['left','right']}
        return {'mesh_h_nm':h,'eps_a':eps_a,'eps_b':eps_b,'uniform_charge_cm3':density,'nodes':len(x),
                'maximum_potential_error_v':float(max(abs(phi-truth))),
                'weighted_l2_error_v':float(np.sqrt(np.dot(vol,(phi-truth)**2)/vol.sum())),
                'gauss_relative_error':float(abs(flux[boundary].sum()+rho[~boundary].sum())/scale),
                'contact_charge_c':charge,'silicon_oxide_displacement_continuity':'analytical piecewise-linear exact field; shared interface DOF'}
    finally:d.delete_device(device=device);d.delete_mesh(mesh=mesh)


def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[]
    for h in [.7,.45,.3]:
        rows.append(slab(out/f'layered_{h}',h))
        rows.append(slab(out/f'uniform_{h}',h,11.7,11.7))
        rows.append(slab(out/f'charged_{h}',h,11.7,11.7,1e18))
    report={'results':rows,'scope':'actual shared-DOF DEVSIM analytical capacitors and uniform charged-slab refinement'}
    (out/'summary.json').write_text(json.dumps(report,indent=2));return report


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/poisson_benchmarks');a=p.parse_args();print(json.dumps(run(a.out),indent=2))
