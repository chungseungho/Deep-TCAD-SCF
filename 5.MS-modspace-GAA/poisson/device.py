"""Independent STEP05 DEVSIM Poisson with actual NodeVolume integration."""
import uuid
import numpy as np
from config import Q


class Poisson:
    def __init__(self,c,mesh):
        import devsim as d
        self.d,self.c=d,c;self.device='gaa_'+uuid.uuid4().hex[:8];m=self.device+'_mesh'
        d.create_gmsh_mesh(mesh=m,file=str(mesh))
        regions={'silicon':('Silicon',c.eps_si),'oxide':('GateOxide',c.eps_ox),
                 'spacer':('Spacer',c.eps_spacer)}
        for r,(group,eps) in regions.items():
            d.add_gmsh_region(mesh=m,gmsh_name=group,region=r,material=group)
        for contact,region in [('Source','silicon'),('Drain','silicon'),('Gate','oxide')]:
            d.add_gmsh_contact(mesh=m,gmsh_name=contact,region=region,name=contact.lower(),material='Metal')
        ints=[('SiOx','silicon','oxide'),('SiSpacer','silicon','spacer'),('OxSpacer','oxide','spacer')]
        for name,r0,r1 in ints: d.add_gmsh_interface(mesh=m,gmsh_name=name,region0=r0,region1=r1,name=name)
        d.finalize_mesh(mesh=m);d.create_device(mesh=m,device=self.device)
        from poisson.discretization import configure
        self.discretization=configure(d,self.device,mesh,{r:g for r,(g,_) in regions.items()})
        for r,(_,eps) in regions.items():
            args={'device':self.device,'region':r}
            d.set_parameter(**args,name='Permittivity',value=eps*8.8541878128e-14)
            d.node_solution(**args,name='Potential');d.edge_from_node_model(**args,node_model='Potential')
            for name,eq in [('DField','Permittivity*(Potential@n0-Potential@n1)*EdgeInverseLength'),
                            ('DField:Potential@n0','Permittivity*EdgeInverseLength'),
                            ('DField:Potential@n1','-Permittivity*EdgeInverseLength')]:
                d.edge_model(**args,name=name,equation=eq)
            kwargs={}
            if r=='silicon':
                for name in ['Donors','Electrons','Susceptibility','PhiReference']: d.node_solution(**args,name=name)
                d.set_parameter(**args,name='q',value=Q)
                d.node_model(**args,name='SpaceCharge',equation='-q*(Donors-Electrons-Susceptibility*(Potential-PhiReference))')
                d.node_model(**args,name='SpaceCharge:Potential',equation='q*Susceptibility')
                kwargs['node_model']='SpaceCharge'
            d.equation(**args,name='PotentialEquation',variable_name='Potential',edge_model='DField',
                       variable_update='default',**kwargs)
        for name,r0,r1 in ints:
            args={'device':self.device,'interface':name}
            for n,eq in [('continuous','Potential@r0-Potential@r1'),
                         ('continuous:Potential@r0','1'),('continuous:Potential@r1','-1')]:
                d.interface_model(**args,name=n,equation=eq)
            d.interface_equation(**args,name='PotentialEquation',interface_model='continuous',type='continuous')
        for name in ['source','drain','gate']:
            d.set_parameter(device=self.device,name=name+'_bias',value=0.)
            d.contact_node_model(device=self.device,contact=name,name=name+'_bc',equation='Potential-'+name+'_bias')
            d.contact_node_model(device=self.device,contact=name,name=name+'_bc:Potential',equation='1')
            d.contact_equation(device=self.device,contact=name,name='PotentialEquation',node_model=name+'_bc',edge_charge_model='DField')
        self.xyz={r:self.coords(r) for r in regions}
        self.volume=self.values('PhysicalControlVolume')*1e21
        if np.any(self.volume<=0) or not np.isfinite(self.volume).all():
            raise RuntimeError('DEVSIM NodeVolume is not positive; refine/rebuild mesh')
        self.donors=c.doping(self.xyz['silicon'])
        self.set_values('Donors',self.donors);self.set_charge(self.donors)
        self.zero_laplace_check={}

    def values(self,name,region='silicon'):
        return np.asarray(self.d.get_node_model_values(device=self.device,region=region,name=name),float)
    def coords(self,r): return np.column_stack([self.values(a,r)*1e7 for a in ['x','y','z']])
    def set_values(self,name,v,region='silicon'):
        self.d.set_node_values(device=self.device,region=region,name=name,values=np.asarray(v,float).tolist())
    def set_charge(self,n):
        n=np.asarray(n)
        if n.shape!=self.donors.shape or np.any(n<0) or not np.isfinite(n).all(): raise ValueError('Invalid Poisson density')
        self.set_values('Electrons',n)
    def solve(self,vg,vd,screening=None):
        # Optional quasi-Newton preconditioner only. Final residual is ALWAYS
        # computed with Susceptibility=0, the original physical Poisson equation.
        if screening is None:
            self.set_values('Susceptibility',np.zeros(len(self.volume)))
        else:
            reference,chi=screening
            self.set_values('PhiReference',reference);self.set_values('Susceptibility',chi)
        for name,bias in [('source',0.),('drain',vd),('gate',vg-self.c.flatband)]:
            self.d.set_parameter(device=self.device,name=name+'_bias',value=float(bias))
        # Relative update denominators near zero can amplify roundoff; the
        # absolute update remains <1e-14 V, far below the outer 1e-4 V criterion.
        self.d.solve(type='dc',absolute_error=1e-14,relative_error=1e-8,maximum_iterations=30,solver_type='direct')
        return self.values('Potential')
    def export(self,path): self.d.write_devices(file=str(path),type='vtk')
    def close(self):
        self.d.delete_device(device=self.device)
        self.d.delete_mesh(mesh=self.device+'_mesh')

    def diagnostics(self):
        return {'regions':list(self.d.get_region_list(device=self.device)),
                'interfaces':list(self.d.get_interface_list(device=self.device)),
                'contacts':list(self.d.get_contact_list(device=self.device)),
                'silicon_node_count':len(self.volume),'nodevolume_min_nm3':float(self.volume.min()),
                'nodevolume_sum_nm3':float(self.volume.sum()),'discretization':self.discretization}
