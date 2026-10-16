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
    dx: float = .25
    dy: float = 1.
    dz: float = 1.
    modes: int = 6
    sd_modes: int = 10
    mesh_h: float = 2.5
    channel_h: float = 1.2
    junction_h: float = .9
    gamma: float = .01  # hbar*relaxation rate, eV, finite RESERVOIRS only
    backend: str = 'spectral'
    basis: str = 'fixed'
    preconditioner: str = 'screened'
    mix: float = .15
    history: int = 5
    max_iterations: int = 160
    charge_tol: float = 2e-3
    phi_tol: float = 1e-4
    current_tol: float = 2e-3
    ms_tie_tol: float = 1e-7

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
        if self.backend not in ('spectral', 'lindblad', 'ms'): raise ValueError('Unknown backend')
        if self.basis not in ('fixed','local'): raise ValueError('basis must be fixed or local')
        if self.preconditioner not in ('screened','none'): raise ValueError('Unknown preconditioner')
        if min(self.ls, self.ld, self.lg, self.spacer, self.width, self.thickness,
               self.sd_width, self.sd_thickness, self.tox, self.dx, self.dy,
               self.dz, self.mesh_h, self.channel_h, self.junction_h,
               self.gamma, self.temperature, self.charge_tol, self.phi_tol,
               self.current_tol) <= 0: raise ValueError('Parameters must be positive')
        if not 0 < self.mix <= 1: raise ValueError('mix must be in (0,1]')
        if self.sd_width < self.width+2*self.tox or self.sd_thickness < self.thickness+2*self.tox:
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
