"""Deterministic backend used only by architecture and SCF wiring tests."""
import numpy as np
from transport.results import TransportResult


class MockQuantumBackend:
    name='mock_quantum'

    def __init__(self,config): self.config=config; self.solve_calls=0
    def prepare(self,context): self.prepared=True
    def calibrate_contacts(self,context): return 0.2,{'status':'mock_calibrated'}
    def solve(self,potential,bias,context):
        self.solve_calls+=1
        phi=np.asarray(potential,dtype=float)
        density=np.ones(len(phi),dtype=float)
        current=float(bias)*1e-9
        result=TransportResult(density,current,np.full(max(len(density)-1,1),current),
            0.2,0.2-float(bias),diagnostics={'backend':'mock_quantum'},
            optional={'tracking':np.array([float(phi.mean())]),'systems':[],'states':[],
                      'valleys':[],'constraint':{}})
        return result
    def convergence_metrics(self,previous,current):
        if previous is None:return self.initial_convergence_metrics()
        return {'occupation_max_change':float(abs(current['tracking'][0]-previous['tracking'][0])),
            'classification_stable':True,'current_relative_change':0.0}
    def initial_convergence_metrics(self):return {'occupation_max_change':float('inf'),
        'classification_stable':True,'current_relative_change':float('inf')}
    def current_relative_change(self,previous,current):return 0.
    def converged(self,iteration,metrics,phi_residual,density_residual,mu_change):
        return bool(iteration>=1 and phi_residual<1e-8 and density_residual<1e-8)
    def state_stability(self,old,new):return True,[]
    def finalize(self,result,bias,grid,accepted):return {'postprocess':[],'postprocess_error':None,
        'current_conservation_abs_a':0.0,'current_conservation_relative':0.0,
        'current_ok':True,'transport_status':'CALCULATED' if accepted else 'UNVALIDATED'}
    def export_fields(self,result):return {}
    def result_metadata(self):return {'gf_inside_scf':False}
    def validate(self,result):
        if not np.all(np.isfinite(result.electron_density_cm3)):
            raise RuntimeError('mock density invalid')


def register_mock_backend():
    from transport.registry import register_backend
    register_backend('mock_quantum',MockQuantumBackend,replace=True)
