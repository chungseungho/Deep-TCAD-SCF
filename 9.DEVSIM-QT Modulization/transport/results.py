"""Backend-neutral quantum transport result contract."""
from dataclasses import dataclass, field
from typing import Any
import numpy as np


@dataclass
class TransportResult:
    electron_density_cm3: np.ndarray
    drain_current_a: float
    bond_current_a: np.ndarray
    mu_left_ev: float
    mu_right_ev: float
    diagnostics: dict[str, Any] = field(default_factory=dict)
    optional: dict[str, Any] = field(default_factory=dict)

    # Temporary STEP08 output-key compatibility. New code should use named fields.
    _ALIASES = {'n':'electron_density_cm3','id_a':'drain_current_a',
                'bonds':'bond_current_a','mu':'mu_left_ev'}

    def __getitem__(self, key):
        if key in self._ALIASES:
            return getattr(self,self._ALIASES[key])
        if hasattr(self,key):
            return getattr(self,key)
        return self.optional[key]

    def __setitem__(self,key,value):
        if key in self._ALIASES:
            setattr(self,self._ALIASES[key],value)
        elif hasattr(self,key):
            setattr(self,key,value)
        else:
            self.optional[key]=value

    def get(self,key,default=None):
        try:return self[key]
        except KeyError:return default

    def keys(self):
        return tuple(self.optional)

    def items(self):
        return self.optional.items()

