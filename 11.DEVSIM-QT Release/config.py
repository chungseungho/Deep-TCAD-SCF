"""Units: geometry nm, Hamiltonian eV, Poisson cm/V, carrier density cm^-3.

Masses and +/- valley multiplicities are copied unchanged from STEP01.
No numerical parameter is a fit to an I-V curve.
"""
from dataclasses import asdict, dataclass
import numpy as np

Q = 1.602176634e-19
HBAR_EV_S = 6.582119569e-16
KB_EV = 8.617333262145e-5
KIN = 0.0380998212  # hbar^2/(2m0), eV nm^2 (STEP01)


@dataclass(frozen=True)
class Valley:
    name: str
    mx: float
    my: float
    mz: float
    multiplicity: int = 2


VALLEYS = (Valley('Delta_z', .19, .19, .916),
           Valley('Delta_x', .916, .19, .19),
           Valley('Delta_y', .19, .916, .19))


@dataclass
class Config:
    geometry_profile: str = 'raised'
    ls: float = 12.
    lg: float = 24.
    ld: float = 12.
    spacer: float = 3.
    width: float = 12.
    thickness: float = 5.
    sd_width: float = 20.
    sd_thickness: float = 15.
    tox: float = 2.
    eps_si: float = 11.7
    eps_ox: float = 3.9
    eps_spacer: float = 7.5
    donor_sd: float = 1e19
    donor_channel: float = 0.
    transition: float = 1.
    temperature: float = 300.
    flatband: float = .5  # explicit inherited assumption; NOT calibrated
    dx: float = .5
    dy: float = 1.
    dz: float = 1.
    modes: int = 2
    sd_modes: int = 4
    mesh_h: float = 3.5
    channel_h: float = 1.8
    junction_h: float = 1.3
    gamma: float = .01  # hbar*relaxation rate, eV, finite RESERVOIRS only
    backend: str = 'ms'
    mapping: str = 'integrated'
    fermi_poles: int = 48
    contour_batch_size: int = 8
    transmission_batch_size: int = 16
    energy_step: float = .01  # composite GL panel width, NOT node spacing (eV)
    energy_order: int = 8
    fermi_tail: float = 16.
    energy_rtol: float = 2e-5
    energy_atol: float = 1e-7
    energy_max_nodes: int = 120000
    basis: str = 'fixed'
    preconditioner: str = 'screened'
    mix: float = .15
    history: int = 5
    max_iterations: int = 160
    charge_tol: float = 2e-3
    phi_tol: float = 1e-4
    current_tol: float = 2e-3
    ms_tie_tol: float = 1e-7
    junction_basis: str = 'nested'
    transport_axis: int = 0
    reservoir_left_nm: float | None = None
    reservoir_right_nm: float | None = None
    ms_contrast_min: float = .10
    ms_reach_min: float = 1e-5
    ms_occupation_uncertainty: float = 1e-7
    ms_degenerate_ev: float = 1e-7
    ms_boundary_shift_slices: int = 1
    ms_number_constraint: str = 'finite_total'
    mu_source_ev: float | None = None
    resonance_order: int = 32
    resonance_max_iterations: int = 12
    high_field_v_nm: float = .5
    dispersion_tolerance: float = .02
    newton_start: int = 25
    newton_trust_radius_v: float = .1
    newton_max_backtracks: int = 16
    potential_reference_v: float = 0.

    density_reference: str = 'lower'
    contact_model: str = 'legacy_flat'
    lead_k_nodes: int = 256
    lead_max_iterations: int = 240
    lead_residual_tol: float = 1e-9

    def __post_init__(self):
        if self.reservoir_left_nm is None:self.reservoir_left_nm=self.ls/2
        if self.reservoir_right_nm is None:self.reservoir_right_nm=self.ld/2

    @property
    def gate_start(self): return self.ls + self.spacer
    @property
    def gate_end(self): return self.gate_start + self.lg
    @property
    def drain_start(self): return self.gate_end + self.spacer
    @property
    def length(self): return self.drain_start + self.ld

    def validate(self):
        if not all(np.isfinite(v) for v in self.to_dict().values() if isinstance(v,(int,float))):
            raise ValueError('All numerical inputs must be finite')
        if min(self.donor_sd,self.eps_si,self.eps_ox,self.eps_spacer)<=0 or self.donor_channel<0 or self.transition<0 or self.history<0:
            raise ValueError('Invalid doping, permittivity, transition length or mixing history')
        from transport.registry import registered_backends
        if self.backend not in registered_backends(): raise ValueError(f'Unknown backend: {self.backend}')
        if self.density_reference not in ('lower','left','right'):raise ValueError('Invalid contour reference')
        if self.contact_model not in ('legacy_flat','transverse_sp'):raise ValueError('Unknown contact model')
        if self.lead_k_nodes<32 or self.lead_residual_tol<=0:raise ValueError('Invalid lead controls')
        if self.mapping not in ('integrated','adjoint','balanced_adjoint','legacy'): raise ValueError('Unknown mapping')
        if min(self.contour_batch_size,self.transmission_batch_size)<1:raise ValueError('Invalid GF batch size')
        if self.fermi_poles<4 or self.energy_order<2 or self.energy_step<=0 or self.fermi_tail<8 or min(self.energy_rtol,self.energy_atol,self.energy_max_nodes)<=0:
            raise ValueError('Invalid quantum integration settings')
        if self.basis not in ('fixed','local'): raise ValueError('basis must be fixed or local')
        if self.junction_basis not in ('nested','plain'): raise ValueError('Unknown junction basis')
        if self.transport_axis not in (0,1,2): raise ValueError('Invalid transport axis')
        if self.newton_trust_radius_v<=0 or self.newton_max_backtracks<1:raise ValueError('Invalid domain Newton controls')
        if self.newton_start<1:raise ValueError('Newton start must be positive')
        if self.ms_number_constraint not in ('finite_total','fixed_contact_mu'): raise ValueError('Invalid number constraint')
        if not (0<self.reservoir_left_nm<=self.ls and 0<self.reservoir_right_nm<=self.ld): raise ValueError('Reservoir must lie inside its electrode')
        if not (0<=self.ms_contrast_min<1 and 0<self.ms_reach_min<1 and self.ms_occupation_uncertainty>0 and self.ms_degenerate_ev>=0):raise ValueError('Invalid MS rejection controls')
        for r in (self.reservoir_left_nm,self.reservoir_right_nm):
            if abs(r/self.dx-round(r/self.dx))>1e-8:raise ValueError('Reservoir boundary must align with quantum cell face')
        if self.resonance_order<8 or self.resonance_max_iterations<2: raise ValueError('Invalid resonance controls')
        if self.preconditioner not in ('screened','none'): raise ValueError('Unknown preconditioner')
        if min(self.ls, self.ld, self.lg, self.spacer, self.width, self.thickness,
               self.sd_width, self.sd_thickness, self.tox, self.dx, self.dy,
               self.dz, self.mesh_h, self.channel_h, self.junction_h,
               self.gamma, self.temperature, self.charge_tol, self.phi_tol,
               self.current_tol) <= 0: raise ValueError('Parameters must be positive')
        if not 0 < self.mix <= 1: raise ValueError('mix must be in (0,1]')
        if self.geometry_profile not in ('raised','uniform'):raise ValueError('Unknown geometry profile')
        if self.geometry_profile=='uniform' and (self.sd_width!=self.width or self.sd_thickness!=self.thickness):
            raise ValueError('Uniform wire must have identical source/channel/drain cross sections')
        if self.geometry_profile=='raised' and (self.sd_width < self.width+2*self.tox or self.sd_thickness < self.thickness+2*self.tox):
            raise ValueError('Raised S/D must enclose nanosheet and oxide envelope')
        if min(self.modes, self.sd_modes, self.max_iterations) < 1: raise ValueError('Invalid count')
        # One common cell grid: geometric faces MUST coincide with cell faces.
        for length, step in ((self.length,self.dx),(self.ls,self.dx),(self.lg,self.dx),
                             (self.spacer,self.dx),(self.ld,self.dx),
                             (self.sd_width,self.dy),(self.width,self.dy),
                             ((self.sd_width-self.width)/2,self.dy),
                             (self.sd_thickness,self.dz),(self.thickness,self.dz),
                             ((self.sd_thickness-self.thickness)/2,self.dz)):
            if abs(length/step-round(length/step)) > 1e-8:
                raise ValueError('Quantum cell faces must align with geometry; change grid spacing')
        return self

    def to_dict(self): return asdict(self)

    def doping(self, xyz):
        x = np.asarray(xyz)[...,0]
        if self.transition == 0:
            f = ((x < self.ls) | (x > self.drain_start)).astype(float)
        else:
            # Smooth inside the nanosheet extensions, while raised blocks stay n+.
            s = np.clip((x-self.ls)/self.transition,0,1)
            d = np.clip((self.drain_start-x)/self.transition,0,1)
            f = 1-(s*s*(3-2*s))*(d*d*(3-2*d))
        return self.donor_channel+(self.donor_sd-self.donor_channel)*f
