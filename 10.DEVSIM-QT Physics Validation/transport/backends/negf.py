from transport.backend_common import initialize, finish
from transport.results import TransportResult
from transport import backend_hooks

class NEGFBackend:
    name='negf'
    def __init__(self,config): self.config=config
    def prepare(self,context): return None
    def calibrate_contacts(self,context): return backend_hooks.negf_calibrate(self.config,context)
    def solve(self,potential,bias,context):
        import numpy as np
        from transport.negf import solve_negf,blocks_density
        d=initialize(dict(context,config=self.config,potential=potential)); c=d['c']; detail=d['detail']
        for s,g in zip(d['systems'],d['deg']):
            blocks,flux,diag,components=solve_negf(s,c,d['mu'],d['mu']-bias,bias,detail)
            d['density']+=blocks_density(s,blocks,g); d['bonds']+=g*flux
            d['tracking'].extend(np.concatenate([np.linalg.eigvalsh(x) for x in blocks]))
            diag['valley']=s.valley.name;d['diagnostics'].append(diag)
            if detail:d['matrices'].append(blocks)
        return finish(d)
    def convergence_metrics(self,previous,current): return backend_hooks.metrics(self,previous,current,self.config)
    def initial_convergence_metrics(self): return {'occupation_max_change':float('inf'),'classification_stable':True,'current_relative_change':float('inf')}
    def current_relative_change(self,previous,current): return abs(current['id_a']-previous['id_a'])/max(abs(current['id_a']),1e-15)
    def state_stability(self,old,new): return True,[]
    def converged(self,iteration,metrics,phi_residual,density_residual,mu_change):
        return bool(iteration>=3 and phi_residual<self.config.phi_tol and density_residual<self.config.charge_tol and mu_change<self.config.phi_tol and metrics['classification_stable'] and metrics['current_relative_change']<self.config.current_tol)
    def finalize(self,result,bias,grid,accepted): return backend_hooks.open_finalize(self.config,result,bias,grid,accepted)
    def export_fields(self,result): return backend_hooks.open_export(result)
    def result_metadata(self): return {'gf_inside_scf':True}
    def validate(self,result): return backend_hooks.validate_density(result)
