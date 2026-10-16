"""Electrochemical energies in eV; positive drain voltage lowers electron μ."""
from dataclasses import dataclass

@dataclass(frozen=True)
class Contacts:
    mu_left: float
    drain_voltage: float
    source_potential: float = 0.

    @property
    def mu_right(self): return self.mu_left-self.drain_voltage

    @property
    def lead_potentials(self): return (self.source_potential,self.source_potential+self.drain_voltage)
