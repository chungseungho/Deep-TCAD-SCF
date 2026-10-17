"""Transport-agnostic coupling of Poisson, mapping, and a backend contract."""
from dataclasses import dataclass
import numpy as np


@dataclass
class SCFStep:
    poisson_potential: np.ndarray
    quantum_potential: np.ndarray
    transport_result: object
    mapped_density: np.ndarray
    poisson_image: np.ndarray
    potential_residual: float


def evaluate_transport_step(poisson,mapping,backend,builders,grid,charge,
                            fixed_mu,gate_bias,drain_bias,detail=False):
    """One physical-map evaluation; no backend algorithm is inspected here."""
    poisson.set_charge(charge)
    phi_p=poisson.solve(gate_bias,drain_bias)
    phi_q=mapping.potential(phi_p)
    context={'config':backend.config,'builders':builders,'grid':grid,
             'fixed_mu':fixed_mu,'detail':detail,'bias':drain_bias}
    result=backend.solve(phi_q,drain_bias,context)
    backend.validate(result)
    mapped=mapping.charge(result['n'])
    poisson.set_charge(mapped)
    image=poisson.solve(gate_bias,drain_bias)
    return SCFStep(phi_p,phi_q,result,mapped,image,
                   float(max(abs(image-phi_p))))
