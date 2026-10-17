"""Extension contracts: density and contact kernels are independent of SCF."""
from typing import Protocol,Any

class ContactModel(Protocol):
    def profile(self,xyz:Any)->Any: ...
    def embedding(self,system:Any,config:Any,lead_phi:Any)->Any: ...

class HamiltonianProvider(Protocol):
    def build(self,potential:Any)->Any: ...

class DensitySolver(Protocol):
    def solve(self,system:Any,config:Any,mu_left:float,mu_right:float,bias:float,detail:bool=False)->Any: ...

class RGFDensitySolver:
    def solve(self,system,config,mu_left,mu_right,bias,detail=False):
        from transport.negf_reference import solve_negf
        return solve_negf(system,config,mu_left,mu_right,bias,detail)

# Reserved data-only contracts; no DFT correction or learned predictor is run.
FUTURE_INPUT_CONTRACT={'coordinates':'nm','potential':'V','hamiltonian':'eV',
    'electron_density':'cm^-3','bias_gauge':'electron energy shifts by -potential',
    'predictor_role':'initial guess only; raw physical residual required',
    'kinetic_correction_role':'Hermitian block operators with declared basis/units/provenance'}
