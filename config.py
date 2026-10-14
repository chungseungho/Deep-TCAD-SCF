"""Mono-nanosheet GAA proof-of-concept. Geometric dimensions: nm.

EMA VALUES ARE EXPLICIT INPUTS, NOT FITTED. Replace the example silicon valley
mass configuration with the user's independently validated EMA when available.
"""
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Valley:
    name: str
    mx: float
    my: float
    mz: float
    multiplicity: int = 2  # +/- valleys, excludes spin (counted separately)


@dataclass
class Config:
    ls_nm: float = 12.0
    lg_nm: float = 24.0
    ld_nm: float = 12.0
    width_nm: float = 12.0
    thickness_nm: float = 5.0
    tox_nm: float = 2.0
    mesh_h_nm: float = 2.0
    nx: int = 49                    # dx=1 nm (48 nm overall)
    ny: int = 7                     # transverse interior finite difference nodes
    nz: int = 5
    modes: int = 2
    temperature_k: float = 300.0
    efermi_source_ev: float = 0.15 # illustrative μ relative to flat-band Si E_c
    ec0_ev: float = 0.0
    donor_sd_cm3: float = 1.0e19
    donor_channel_cm3: float = 0.0
    eps_si: float = 11.7
    eps_ox: float = 3.9
    energy_min_ev: float = -0.45
    energy_max_ev: float = 0.85
    energy_points: int = 220
    eta_ev: float = 1.0e-5
    mix: float = 0.15
    tol_v: float = 2.0e-3
    valleys: tuple = field(default_factory=lambda: (
        Valley('Delta_z', mx=0.19, my=0.19, mz=0.916),
        Valley('Delta_x', mx=0.916, my=0.19, mz=0.19),
        Valley('Delta_y', mx=0.19, my=0.916, mz=0.19),
    ))

    @property
    def length_nm(self):
        return self.ls_nm + self.lg_nm + self.ld_nm

    def validate(self):
        assert self.nx >= 3 and self.ny >= 2 and self.nz >= 2
        assert 1 <= self.modes < self.ny * self.nz
        assert 0.0 < self.mix <= 1.0
        assert self.energy_points >= 16
        assert all(v.mx > 0 and v.my > 0 and v.mz > 0 for v in self.valleys)
        assert self.mesh_h_nm > 0 and self.tox_nm > 0
