"""Adapter preserving STEP08 calculations while exposing the new contract."""
from transport.results import TransportResult


class LegacyBackendAdapter:
    name = ''

    def __init__(self, config):
        self.config = config

    def prepare(self, context):
        return None

    def calibrate_contacts(self, context):
        # STEP08 contact calibration remains in its existing transport utilities.
        return {'status': 'delegated_to_step08_transport_utilities'}

    def solve(self, potential, bias, context):
        values=context['legacy_evaluator'](self.config,context['builders'],context['grid'],
            potential,context['fixed_mu'],bias,context['detail'])
        current=values['id_a']
        result=TransportResult(values['n'],current,values['bonds'],values['mu'],
            values['mu']-bias, diagnostics={'valleys':values['valleys']}, optional=values)
        self.validate(result)
        return values

    def convergence_metrics(self, previous, current):
        if previous is None:
            return {'density_relative_change': float('inf'), 'current_relative_change': float('inf')}
        dn=current.electron_density_cm3-previous.electron_density_cm3
        scale=max(float(abs(current.electron_density_cm3).max()),1e-30)
        di=abs(current.drain_current_a-previous.drain_current_a)/max(abs(current.drain_current_a),1e-15)
        return {'density_relative_change':float(abs(dn).max()/scale),
                'current_relative_change':float(di)}

    def validate(self, result):
        import numpy as np
        if result.electron_density_cm3.ndim != 1 or not np.all(np.isfinite(result.electron_density_cm3)):
            raise RuntimeError(f'{self.name} backend returned invalid electron density')

