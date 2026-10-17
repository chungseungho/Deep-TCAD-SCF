from transport.backend_common import initialize, finish
from transport.results import TransportResult
from transport import backend_hooks

class NEGFBackend:
    name='negf'
    def __init__(self,config):
        from transport.contracts import RGFDensitySolver
        self.config=config;self.density_solver=RGFDensitySolver()
    def prepare(self,context):
        if 'poisson' in context and getattr(self,'contact',None) is not None:
            context['poisson'].set_boundary_profile(self.contact.profile)
    def calibrate_contacts(self,context):
        if self.config.contact_model=='legacy_flat':return backend_hooks.negf_calibrate(self.config,context)
        from transport.contact_model import solve_reference
        self.contact=solve_reference(self.config,context['grid'],[b.v for b in context['builders']])
        for b in context['builders']:b.contacts=self.contact;b.block_only=True
        return self.contact.mu-self.config.potential_reference_v,self.contact.diagnostics
    def solve(self,potential,bias,context):
        import numpy as np
        from transport.negf import solve_negf,blocks_density
        for builder in context['builders']:builder.block_only=True
        d=initialize(dict(context,config=self.config,potential=potential)); c=d['c']; detail=d['detail']
        if abs(bias)>1e-14:
            from transport.bound_states import audit,BoundOccupationError
            for system in d['systems']:
                check=audit(system,c,d['mu'],bias)
                if check['bias_dependent_bound_occupation_unresolved']:raise BoundOccupationError(check)
        for s,g in zip(d['systems'],d['deg']):
            blocks,flux,diag,components=self.density_solver.solve(s,c,d['mu'],d['mu']-bias,bias,detail)
            d['density']+=blocks_density(s,blocks,g); d['bonds']+=g*flux
            d['tracking'].extend(np.concatenate([np.linalg.eigvalsh(x) for x in blocks]))
            diag['valley']=s.valley.name;d['diagnostics'].append(diag)
            if detail:d['matrices'].append(blocks)
        d['left']=d['right']=d['unresolved']=None
        return finish(d)
    def convergence_metrics(self,previous,current): return backend_hooks.metrics(self,previous,current,self.config)
    def initial_convergence_metrics(self): return {'occupation_max_change':float('inf'),'classification_stable':True,'current_relative_change':float('inf')}
    def current_relative_change(self,previous,current): return abs(current['id_a']-previous['id_a'])/max(abs(current['id_a']),1e-15)
    def state_stability(self,old,new): return True,[]
    def converged(self,iteration,metrics,phi_residual,density_residual,mu_change):
        return bool(iteration>=3 and phi_residual<self.config.phi_tol and density_residual<self.config.charge_tol and mu_change<self.config.phi_tol and metrics['classification_stable'] and metrics['current_relative_change']<self.config.current_tol)
    def finalize(self,result,bias,grid,accepted): return backend_hooks.open_finalize(self.config,result,bias,grid,accepted)
    def physical_validity(self,result,fields):return backend_hooks.physical_validity(self.config,result,fields)
    def export_fields(self,result): return backend_hooks.open_export(result)
    def result_metadata(self): return {'gf_inside_scf':True,'source_drain_density':'UNSUPPORTED','contact_model':self.config.contact_model,'bound_occupation_policy':'equilibrium Fermi; biased isolated states with unresolved population rejected','training_ready':False}
    def validate(self,result): return backend_hooks.validate_density(result)
