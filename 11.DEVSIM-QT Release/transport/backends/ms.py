from transport.backend_common import initialize, finish
from transport.results import TransportResult
from transport import backend_hooks

class MSBackend:
    name='ms'
    def __init__(self,config): self.config=config
    def prepare(self,context): return None
    def calibrate_contacts(self,context): return backend_hooks.ms_calibrate(self.config,context)
    def solve(self,potential,bias,context):
        from device.regions import mode_regions
        from multispace.occupation import solve_states
        from multispace.constraints import solve_number
        from modespace.hamiltonian import density_from_matrix
        from transport.reservoir import bond_currents
        d=initialize(dict(context,config=self.config,potential=potential)); c=d['c'];grid=d['grid'];detail=d['detail']
        for s in d['systems']:
            labels=mode_regions(s,c); alternatives=[]; step=c.ms_boundary_shift_slices*c.dx
            if step:
                for sign in (-1,1): alternatives.append(mode_regions(s,c,sign*step))
            d['states'].append(solve_states(s.H,labels,c,alternatives))
        target=float(c.doping(grid.xyz).sum()*grid.cell_volume*1e-21)
        d['mu'],d['states'],d['constraint']=solve_number(d['states'],d['deg'],target,bias,c,context['fixed_mu'])
        for s,g,state in zip(d['systems'],d['deg'],d['states']):
            u,f=state['u'],state['occupation']; rho=(u*f)@u.conj().T
            d['density']+=density_from_matrix(s,rho,g)
            for origin,out in ((0,d['left']),(2,d['right']),(-1,d['unresolved'])):
                use=state['origin']==origin; part=(u[:,use]*f[use])@u[:,use].conj().T
                out+=density_from_matrix(s,part,g)
            closed=bond_currents(s,rho,g)
            d['diagnostics'].append({'valley':s.valley.name,'dimension':s.dimension,
                'orthogonality_error':float(__import__('numpy').max(abs(u.conj().T@u-__import__('numpy').eye(len(u))))),
                'projector_sum_error':float(max(abs(state['weights'].sum(axis=1)-1))),
                'closed_stationary_max_current_a':float(max(abs(closed))),
                'pauli_min':float(min(f)),'pauli_max':float(max(f)),
                'origin_counts':{str(o):int(sum(state['origin']==o)) for o in (-1,0,2)},
                'active_states':int(sum(state['active'])),'degenerate_clusters':state['clusters'],'gf_evaluations_inside_scf':0})
            d['tracking'].extend(f)
            if detail:d['matrices'].append(rho)
        d['current_is_undefined']=True
        return finish(d)
    def convergence_metrics(self,previous,current): return backend_hooks.metrics(self,previous,current,self.config)
    def initial_convergence_metrics(self): return {'occupation_max_change':float('inf'),'classification_stable':False,'current_relative_change':0.}
    def current_relative_change(self,previous,current): return 0.
    def state_stability(self,old,new): return backend_hooks.ms_state_stability(old,new,self.config)
    def converged(self,iteration,metrics,phi_residual,density_residual,mu_change):
        return bool(iteration>=3 and phi_residual<self.config.phi_tol and density_residual<self.config.charge_tol and mu_change<self.config.phi_tol and metrics['classification_stable'] and metrics['current_relative_change']<self.config.current_tol)
    def finalize(self,result,bias,grid,accepted): return backend_hooks.ms_finalize(self.config,result,bias,grid,accepted)
    def physical_validity(self,result,fields):return backend_hooks.physical_validity(self.config,result,fields)
    def export_fields(self,result): return backend_hooks.ms_export(result)
    def result_metadata(self): return {'gf_inside_scf':False,
        'density_current_distinction':'closed stationary MS density; independent open scattering current after SCF'}
    def validate(self,result): return backend_hooks.validate_density(result)
