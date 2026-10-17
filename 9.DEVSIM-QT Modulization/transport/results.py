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

