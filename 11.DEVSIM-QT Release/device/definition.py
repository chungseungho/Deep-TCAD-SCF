"""Device coordinates and model interfaces; transport uses ordered blocks."""
from dataclasses import dataclass
from typing import Protocol
import numpy as np
from config import VALLEYS, Valley


@dataclass(frozen=True)
class Frame:
    axis: int = 0

    @property
    def order(self): return ((0,1,2),(1,2,0),(2,0,1))[self.axis]

    def to_world(self,canonical):
        a=np.asarray(canonical);out=np.empty_like(a);out[...,self.order]=a
        return out

    def to_canonical(self,world): return np.asarray(world)[...,self.order]

    def valleys(self):
        return tuple(Valley(v.name,*np.array([v.mx,v.my,v.mz])[list(self.order)],v.multiplicity) for v in VALLEYS)


class QuantumModel(Protocol):
    def build(self,potential): ...


class ElectrostaticModel(Protocol):
    def set_charge(self,density): ...
    def solve(self,vg,vd): ...


class TransportModel(Protocol):
    def density(self,system,contacts): ...
    def postprocess(self,system,contacts): ...


@dataclass(frozen=True)
class DeviceDefinition:
    name: str = 'raised_source_drain_nanosheet_GAA'
    transport_axis: int = 0

    @property
    def frame(self): return Frame(self.transport_axis)

    def describe(self):
        return {'name':self.name,'world_transport_axis':self.transport_axis,
                'coordinates':'quantum/Poisson calculation in canonical longitudinal/transverse frame',
                'other_geometry_validation':'interface only; no vertical transistor physics certification'}
