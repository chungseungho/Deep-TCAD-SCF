from transport.backend_common import initialize, finish
from transport.results import TransportResult
from transport import backend_hooks

class SpectralBackend:
    name='spectral'
    def __init__(self,config): self.config=config
    def prepare(self,context): return None
    def calibrate_contacts(self,context): return backend_hooks.spectral_calibrate(self.config,context)
    def solve(self,potential,bias,context):
        import numpy as np
        from modespace.hamiltonian import density_from_matrix
        from transport.finite_spectral import solve_finite_spectral
        from transport.reservoir import bond_currents
        d=initialize(dict(context,config=self.config,potential=potential));c=d['c'];detail=d['detail']
        for s,g in zip(d['systems'],d['deg']):
            rho,diag,parts=solve_finite_spectral(s.H,s.regions,c.gamma,d['mu'],d['mu']-bias,c.temperature,split=detail,diagnostic_detail=detail)
            d['tracking'].extend(np.linalg.eigvalsh(rho) if detail else diag.pop('_correlation_spectrum'))
            d['density']+=density_from_matrix(s,rho,g);d['bonds']+=bond_currents(s,rho,g)
            if detail:
                d['left']+=density_from_matrix(s,parts[0],g);d['right']+=density_from_matrix(s,parts[1],g);d['matrices'].append(rho)
            diag['valley']=s.valley.name;d['diagnostics'].append(diag)
        return finish(d)
    def convergence_metrics(self,previous,current): return backend_hooks.metrics(self,previous,current,self.config)
    def initial_convergence_metrics(self): return {'occupation_max_change':float('inf'),'classification_stable':True,'current_relative_change':float('inf')}
    def current_relative_change(self,previous,current): return abs(current['id_a']-previous['id_a'])/max(abs(current['id_a']),1e-15)
    def state_stability(self,old,new): return True,[]
    def converged(self,iteration,metrics,phi_residual,density_residual,mu_change):
        return bool(iteration>=3 and phi_residual<self.config.phi_tol and density_residual<self.config.charge_tol and mu_change<self.config.phi_tol and metrics['classification_stable'] and metrics['current_relative_change']<self.config.current_tol)
    def finalize(self,result,bias,grid,accepted): return backend_hooks.open_finalize(self.config,result,bias,grid,accepted)
    def export_fields(self,result): return backend_hooks.open_export(result)
    def result_metadata(self): return {'gf_inside_scf':False}
    def validate(self,result): return backend_hooks.validate_density(result)
