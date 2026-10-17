"""DEVSIM on one conformal P1 potential space with explicit material stiffness.

Every physical mesh vertex has ONE electrostatic DOF, including triple lines.
Material permittivities enter the original tetra stiffness before assembly.
Silicon loads use V_Si/V_all times the positive global lumped volume. This is
P1 FEM expressed as signed edge balance, not native circumcentric FVM.
"""
import uuid,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from config import Q
from poisson.discretization import mesh_regions,assemble

EPS0=8.8541878128e-14


def merged_mesh(path):
    path=Path(path);lines=path.read_text().splitlines()
    i=lines.index('$PhysicalNames');n=int(lines[i+1])
    tag=max(int(row.split()[1]) for row in lines[i+2:i+2+n])+1
    lines[i+1]=str(n+1);lines.insert(i+2+n,f'3 {tag} "Electrostatic"')
    j=lines.index('$Elements');count=int(lines[j+1])
    for k in range(j+2,j+2+count):
        f=lines[k].split()
        if f[1]=='4': f[3]=str(tag);lines[k]=' '.join(f)
    out=path.with_name('gaa_poisson_global_cm.msh')
    out.write_text('\n'.join(lines)+'\n');return out


class Poisson:
    def __init__(self,c,mesh):
        import devsim as d
        self.d,self.c=d,c;self.mesh_path=str(mesh)
        self.device='gaa_'+uuid.uuid4().hex[:8];self.mesh=self.device+'_mesh'
        points,tets=mesh_regions(mesh)
        d.create_gmsh_mesh(mesh=self.mesh,file=str(merged_mesh(mesh)))
        d.add_gmsh_region(mesh=self.mesh,gmsh_name='Electrostatic',region='electrostatic',material='ExplicitTetraMaterials')
        for contact in ['Source','Drain','Gate']:
            d.add_gmsh_contact(mesh=self.mesh,gmsh_name=contact,region='electrostatic',name=contact.lower(),material='Metal')
        d.finalize_mesh(mesh=self.mesh);d.create_device(mesh=self.mesh,device=self.device)
        self.args=dict(device=self.device,region='electrostatic')
        coords=np.column_stack([d.get_node_model_values(**self.args,name=a) for a in ['x','y','z']])
        self.allxyz=coords*1e7;dist,index=cKDTree(coords).query(points)
        if max(dist)>1e-12: raise RuntimeError('Global mesh nodes mismatch')
        K=csr_matrix((len(coords),len(coords)));total=np.zeros(len(coords));si=np.zeros(len(coords))
        self.indices={};self.xyz={};self.material_volume={};self.discretization={}
        for region,group,eps in [('silicon','Silicon',c.eps_si),('oxide','GateOxide',c.eps_ox),('spacer','Spacer',c.eps_spacer)]:
            tet=index[tets[group]];kr,volume=assemble(coords,tet)
            K=K+eps*kr;total+=volume
            if region=='silicon': si=volume
            self.indices[region]=np.unique(tet);self.xyz[region]=self.allxyz[self.indices[region]]
            self.material_volume[region]=volume
            self.discretization[region]={'physical_volume_nm3':float(volume.sum()*1e21),'relative_permittivity':eps}
        if np.any(total<=0): raise RuntimeError('Nonpositive shared Poisson volume')
        self.K=K*EPS0;self.global_volume=total;self.si_volume=si
        self.volume=si[self.indices['silicon']]*1e21
        self.donors=c.doping(self.xyz['silicon'])
        for name in ['PhysicalControlVolume','SiliconVolumeFraction','Donors','Electrons','Susceptibility','PhiReference','ContactReference','Potential']:
            d.node_solution(**self.args,name=name)
        self._set('PhysicalControlVolume',total);self._set('SiliconVolumeFraction',si/total)
        d.edge_from_node_model(**self.args,node_model='node_index')
        e0=np.array(d.get_edge_model_values(**self.args,name='node_index@n0'),int)
        e1=np.array(d.get_edge_model_values(**self.args,name='node_index@n1'),int)
        coupling=-np.asarray(K[e0,e1]).ravel()*np.linalg.norm(coords[e0]-coords[e1],axis=1)
        d.edge_solution(**self.args,name='SignedMaterialCouple')
        d.set_edge_values(**self.args,name='SignedMaterialCouple',values=coupling.tolist())
        d.set_parameter(**self.args,name='node_volume_model',value='PhysicalControlVolume')
        d.set_parameter(**self.args,name='edge_couple_model',value='SignedMaterialCouple')
        d.set_parameter(**self.args,name='epsilon0',value=EPS0);d.set_parameter(**self.args,name='q',value=Q)
        d.edge_from_node_model(**self.args,node_model='Potential')
        for name,equation in [('DField','epsilon0*(Potential@n0-Potential@n1)*EdgeInverseLength'),
                ('DField:Potential@n0','epsilon0*EdgeInverseLength'),('DField:Potential@n1','-epsilon0*EdgeInverseLength')]:
            d.edge_model(**self.args,name=name,equation=equation)
        d.node_model(**self.args,name='SpaceCharge',equation='-q*SiliconVolumeFraction*(Donors-Electrons-Susceptibility*(Potential-PhiReference))')
        d.node_model(**self.args,name='SpaceCharge:Potential',equation='q*SiliconVolumeFraction*Susceptibility')
        d.equation(**self.args,name='PotentialEquation',variable_name='Potential',edge_model='DField',node_model='SpaceCharge')
        self.boundary=np.zeros(len(coords),bool)
        for name in ['source','drain','gate']:
            d.set_parameter(device=self.device,name=name+'_bias',value=0.)
            d.contact_node_model(device=self.device,contact=name,name=name+'_bc',equation='Potential-ContactReference-'+name+'_bias')
            d.contact_node_model(device=self.device,contact=name,name=name+'_bc:Potential',equation='1')
            d.contact_equation(device=self.device,contact=name,name='PotentialEquation',node_model=name+'_bc',edge_charge_model='DField')
            elements=d.get_element_node_list(device=self.device,region='electrostatic',contact=name)
            self.boundary[np.unique(elements)]=True
        self.set_values('Donors',self.donors);self.set_charge(self.donors)

    def _get(self,name):
        return np.array(self.d.get_node_model_values(**self.args,name=name),float)
    def _set(self,name,value):
        self.d.set_node_values(**self.args,name=name,values=np.asarray(value,float).tolist())
    def values(self,name,region='silicon'):
        if name=='PhysicalControlVolume': return self.material_volume[region][self.indices[region]]
        return self._get(name)[self.indices[region]]
    def coords(self,region): return self.xyz[region]
    def set_values(self,name,value,region='silicon'):
        allvalues=np.zeros(len(self.allxyz));allvalues[self.indices[region]]=np.asarray(value)
        self._set(name,allvalues)
    def set_charge(self,n):
        n=np.asarray(n)
        if n.shape!=self.donors.shape or np.any(n<0) or not np.isfinite(n).all(): raise ValueError('Invalid Poisson density')
        self.set_values('Electrons',n)
    def set_boundary_profile(self,provider):
        # Only end-face boundary data change. Stiffness, loads and mapping stay
        # identical; gate reference remains zero.
        profile=np.zeros(len(self.allxyz))
        end=np.isclose(self.allxyz[:,0],0,atol=1e-8)|np.isclose(self.allxyz[:,0],self.c.length,atol=1e-8)
        profile[end]=provider(self.allxyz[end]);self._set('ContactReference',profile)

    def solve(self,vg,vd,screening=None):
        if screening is None: self._set('Susceptibility',np.zeros(len(self.allxyz)))
        else:
            reference,chi=screening;self.set_values('PhiReference',reference);self.set_values('Susceptibility',chi)
        ref=self.c.potential_reference_v
        for name,bias in [('source',ref),('drain',ref+vd),('gate',ref+vg-self.c.flatband)]:
            self.d.set_parameter(device=self.device,name=name+'_bias',value=float(bias))
        self.d.solve(type='dc',absolute_error=1e-14,relative_error=1e-8,maximum_iterations=30,solver_type='direct')
        return self.values('Potential')
    def gauss_balance(self):
        phi=self._get('Potential');source=Q*(self._get('Donors')-self._get('Electrons'))*self.si_volume
        flux=self.K@phi;interior=~self.boundary;residual=flux-source
        scale=max(float(abs(source).sum()),float(abs(flux[self.boundary]).sum()),1e-25)
        charges={name:float(self.d.get_contact_charge(device=self.device,contact=name,equation='PotentialEquation')) for name in ['source','drain','gate']}
        return {'contact_displacement_charge_c':charges,'interior_free_charge_c':float(source[interior].sum()),
                'integrated_boundary_reaction_c':float(flux[self.boundary].sum()),
                'gauss_relative_error':float(abs(flux[self.boundary].sum()+source[interior].sum())/scale),
                'maximum_interior_equation_residual_c':float(max(abs(residual[interior]))),
                'contact_sum_relative_error':float(abs(sum(charges.values())+source[interior].sum())/scale),
                'contact_node_charge_treatment':'Dirichlet rows replaced; only interior volume sources enter balance'}
    def diagnostics(self):
        return {'regions':list(self.d.get_region_list(device=self.device)),
                'physical_materials':['Silicon','GateOxide','Spacer'],
                'material_interface_treatment':'one shared conformal P1 potential DOF; element permittivity assembled before shared edge flux',
                'interfaces':list(self.d.get_interface_list(device=self.device)),
                'contacts':list(self.d.get_contact_list(device=self.device)),
                'silicon_node_count':len(self.volume),'nodevolume_min_nm3':float(self.volume.min()),
                'nodevolume_sum_nm3':float(self.volume.sum()),'discretization':self.discretization,
                'silicon_oxide_potential_jump_v':0.,'triple_junction_potential_jump_v':0.}
    def export(self,path): self.d.write_devices(file=str(path),type='vtk')
    def close(self):
        self.d.delete_device(device=self.device);self.d.delete_mesh(mesh=self.mesh)
