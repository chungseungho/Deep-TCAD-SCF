# STEP07 status

This is a tested research prototype, **not a validated nonequilibrium MS TCAD**.
No physical biased MS I–V has been reproduced on this GAA under the stated
origin, reservoir and neutrality conditions.

| Level | Evidence-based result |
|---|---|
| 1 — Multi-Space implemented | Achieved for the admissible-state contract; 56 tests pass |
| 2 — Actual GAA MS–Poisson | Partial: equilibrium converges; no accepted biased MS point |
| 3 — MS transport reproduced | Not achieved for biased GAA; independent post-SCF/analytical checks pass separately |
| 4 — Physics/numerical reference | Not achieved: density/contact/basis errors, no joint GAA plateau |
| 5 — General quantum TCAD | Not achieved: frame/interfaces implemented, other device physics unvalidated |

## Actual outcomes

- Finite-N MS at Vg=.8 V, Vd=0: 28 iterations,
  phi residual 6.218e-08 V, density residual 8.826e-08.
- Transfer: 25 attempted, **zero accepted biased points**.
- Output: 75 attempted, **three accepted equilibrium points, 72 rejected biased points**.
- Frozen MS/NEGF density error at finite-N MS phi and mu:
  6.065% relative volume L2.
- Same neutral-lead mu, separately converged coupled fixed points:
  24.066% density error and
  9.522 mV maximum nodal potential difference.
- Finite-spectral regression: 4.4922834576e-06 A, independent finite-bath
  Landauer mismatch below 8.175e-09 relative.
- Repaired frozen semi-infinite density integration agrees across two panel
  settings within 1.011e-13, with no added eta.
  This is not a coupled biased reference.
- Fully coupled NEGF equilibrium points converge. Biased .8/.1 V and .4/.01 V
  runs stop when a Newton line-search trial has an undetermined bound-state
  filling. This does **not** prove no feasible fixed point exists; a certified
  biased solution and robust treatment of that domain remain unfinished.
- Gate/drain stress points through 6 V execute with rejected densities,
  available-iterate electric fields and explicit physical-scope flags.

## Blocking issues

Finite closed ends deplete reservoir density. A real stationary eigenstate has
no inherent injection direction. Active C-state origin changes under a one-slice
boundary shift. Rotating states or assigning a mean mu would conceal missing
evidence. Larger subspaces change low-energy density substantially. Joint
mesh/mode/contact and high-bias dispersion limits are not established.

The SI downloads were blocked and are not claimed as inspected. The supplied
SIESTA archive has MS extensions, but is not authenticated as the exact 2020
implementation. Source paths and access limits are documented separately.
