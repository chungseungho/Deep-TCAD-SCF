"""Build evidence summaries from recorded results, with no inferred successes."""
import sys,json,csv,platform,importlib.metadata as im
from pathlib import Path
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from provenance import source_digest
R=ROOT/'results'
def read(p):return json.loads((R/p).read_text())


def run():
    counts={}
    for name,file in [('transfer','stage_a_transfer/idvgs.csv'),('output','stage_a_output/idvds.csv'),
                      ('stress_gate','stress_gate/stress.csv'),('stress_drain','stress_drain/stress.csv')]:
        rows=list(csv.DictReader((R/file).open()))
        counts[name]={'points':len(rows),'numerical_status':dict(Counter(r['numerical_status'] for r in rows))}
    a=read('gaa_ms_lowbias/states/vg0.800000_vd0.000000.json')
    c=read('gaa_spectral_regression/states/vg0.800000_vd0.100000.json')
    f=read('ms_equilibrium_frozen_negf/comparison.json')
    m=read('matched_mu_coupled_comparison/summary.json');mf=read('matched_mu_frozen_comparison/comparison.json')
    repair=read('repaired_frozen_negf/summary.json');finite=read('finite_current_check/summary.json')
    analytic=read('analytical/summary.json');tests=read('test_report/summary.json')
    s=np.load(R/'gaa_ms_lowbias/states/vg0.800000_vd0.000000.npz')
    donor_p=float(np.dot(s['donors_cm3'],s['node_volume_nm3'])*1e-21)
    donor_q=float(s['doping_quantum_cm3'].sum()*float(s['cell_volume_nm3'])*1e-21)
    env={'python':platform.python_version(),'packages':{n:im.version(n) for n in ['numpy','scipy','matplotlib','gmsh','devsim','pytest']},
         'source_sha256':source_digest(),'threads':1,'external_native_libraries':'private OpenBLAS/UMFPACK and Gmsh dependencies; not included'}
    (R/'environment.json').write_text(json.dumps(env,indent=2))
    status=f'''# STEP07 status

This is a tested research prototype, **not a validated nonequilibrium MS TCAD**.
No physical biased MS I–V has been reproduced on this GAA under the stated
origin, reservoir and neutrality conditions.

| Level | Evidence-based result |
|---|---|
| 1 — Multi-Space implemented | Achieved for the admissible-state contract; {tests['total_tests']} tests pass |
| 2 — Actual GAA MS–Poisson | Partial: equilibrium converges; no accepted biased MS point |
| 3 — MS transport reproduced | Not achieved for biased GAA; independent post-SCF/analytical checks pass separately |
| 4 — Physics/numerical reference | Not achieved: density/contact/basis errors, no joint GAA plateau |
| 5 — General quantum TCAD | Not achieved: frame/interfaces implemented, other device physics unvalidated |

## Actual outcomes

- Finite-N MS at Vg=.8 V, Vd=0: {a['iterations']} iterations,
  phi residual {a['phi_residual_v']:.3e} V, density residual {a['density_relative_residual']:.3e}.
- Transfer: 25 attempted, **zero accepted biased points**.
- Output: 75 attempted, **three accepted equilibrium points, 72 rejected biased points**.
- Frozen MS/NEGF density error at finite-N MS phi and mu:
  {100*f['density_relative_volume_l2']:.3f}% relative volume L2.
- Same neutral-lead mu, separately converged coupled fixed points:
  {100*m['density_relative_volume_l2']:.3f}% density error and
  {1e3*m['maximum_node_phi_difference_v']:.3f} mV maximum nodal potential difference.
- Finite-spectral regression: {c['id_a']:.12g} A, independent finite-bath
  Landauer mismatch below {max(x['relative_mismatch'] for x in finite['checks']):.3e} relative.
- Repaired frozen semi-infinite density integration agrees across two panel
  settings within {repair['density_relative_volume_l2']:.3e}, with no added eta.
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
'''
    (ROOT/'STATUS.md').write_text(status)
    uniform=analytic['cases']['barrier_0']
    validation=f'''# Executed validation

## Test completion

{tests['total_tests']} tests pass: 52 quantum/numerical and four actual native
DEVSIM/Gmsh tests in isolated processes, each with a completed JUnit XML.
`results/test_report/summary.json` is the completion authority. An earlier mixed
run exited before its summary; its zero exit code was not counted as a pass.
Native logs retain a NumPy long-double warning; model arrays use float64.

Coverage includes positive projector sums, electron count, C-dominant origin,
weak/bound-state rejection, degenerate rotations, equilibrium/gauge limits,
Pauli bounds, masses, full 3D projection, nested trace, independent density/GF,
RGF/dense agreement, Poisson/Gauss/interface/capacitor checks and the screened
Newton inverse. The native/quantum test isolation is also the recommended
portable test command, not a suppressed failure.

## Independent small models

Chain, barrier and well checks pass (`results/analytical/`). Uniform T error is
{uniform['uniform_ballistic_max_error']:.3e}; equilibrium population error
{uniform['equilibrium_population_error']:.3e}; current relative error
{uniform['current_relative_error']:.3e}. The discrete well spectrum error is
{analytic['quantum_well_max_error_ev']:.3e} eV.

`results/reduced_realspace/` executes three grids (1,.5,.25 nm), with full,
fixed-low, nested-low and local-nested-low bases, against an independent 3D
stencil. Full-basis density and projection identities pass. This is not a
continuum or self-consistent GAA convergence plateau.

## Actual GAA and conservation

| Case | Recorded result | Scope |
|---|---|---|
| MS finite-N .8/0 V | phiR {a['phi_residual_v']:.3e} V; nR {a['density_relative_residual']:.3e}; Gauss {a['gauss_relative_error']:.3e} | Equilibrium only |
| MS transfer | 25/25 origins rejected | No current curve |
| MS output | 3 equilibrium converged, 72 biased rejected | No biased current curve |
| Finite spectral .8/.1 V | I {c['id_a']:.12g} A; conservation {c['current_conservation_relative']:.3e} | Separate finite-bath model |
| Independent finite-bath T integral | max bond mismatch {max(x['relative_mismatch'] for x in finite['checks']):.3e} | Same frozen H, three panel settings |
| Repaired semi-infinite integral | density difference {repair['density_relative_volume_l2']:.3e}; I difference {repair['current_relative_difference']:.3e} | Two panels, no eta; not coupled SCF |

All accepted-state Gauss checks and raw unmixed residuals are recorded per
state. One potential DOF spans material interfaces; the operator is the
configured tetrahedral P1 discretization, not native circumcentric FVM.

## Density and ensemble comparisons

At the finite-N MS field and mu: N_MS={f['ms_electrons']:.10g},
N_NEGF={f['negf_electrons']:.10g}; relative volume density L2 error
{100*f['density_relative_volume_l2']:.5f}%. That mu differs from the independently
neutral infinite-lead mu. These are frozen-field results, not same-contact
coupled equivalence.

At the **same neutral-lead mu**, separately coupled MS and NEGF fixed points
differ by {100*m['density_relative_volume_l2']:.5f}% in density and
{1e3*m['maximum_node_phi_difference_v']:.6g} mV in maximum nodal potential.
The frozen open-density error on this MS field is
{100*mf['density_relative_volume_l2']:.5f}%. Regional populations and all-valley
lowest-subband channel barriers are in `matched_mu_coupled_comparison/`.
Equilibrium I=0 does not establish density equivalence.

Finite-N MS reservoir deficits at Vg=.8 are
{100*a['regions']['left_reservoir']['neutrality_relative']:.3f}% (L) and
{100*a['regions']['right_reservoir']['neutrality_relative']:.3f}% (R).
Screening targets fail. Boundary placement is a sensitivity test, not an IV fit.

## Numerical limits

Actual frozen GAA subspaces 2/4, 4/8 and 6/12, plus plain/local comparisons,
record dimensions, orthogonality, junction missing norm, density and population
above mu+8 kBT. Large added-mode density is not called a negligible high-energy
tail. No mode plateau is certified.

The donor field has different discrete measures: quantum cells contain
{donor_q:.12g} donor electrons, while Poisson nodal lumping gives
{donor_p:.12g}, a {100*(donor_p/donor_q-1):.5f}% quadrature discrepancy.
The finite-N constraint uses the declared quantum-cell count; Poisson Gauss
uses its actual nodal load. This error is reported, not hidden by donor fitting
or a claim of exact global neutrality.

`results/numerical_plan/` contains 22 **unexecuted** configurations for three
levels of longitudinal/transverse/mode/Poisson/reservoir/contact studies and
fixed/local/plain/nested comparisons. A plan is not a result. A joint coupled
GAA current/density plateau and a bias-domain contact limit remain unfinished.
STEP06's two-grid 0.47% difference is not treated as a plateau.

## Stress, gauge and restart

Six gate points through 6 V are all origin-rejected. Six drain points through
6 V give one equilibrium solve and five rejections. Numerical status and
available-iterate field/scope flags are separate. Fields during a rejected
Newton trial can be the last available Poisson field, not the unevaluated
trial Hamiltonian; they are not SCF predictions.

A .4 V common potential shift is tested in actual GAA SCF and in frozen open
density/current/transmission. The finite-tolerance SCF differences are in
`gauge_scf_comparison.json`; exact frozen gauge tests are separate. Resume
preserves accepted NPZ bytes; regenerated CSV physical values match that state.

## Not executed or not established

- No accepted biased MS GAA fixed point or reproduced MS transport IV.
- No certified biased GAA NEGF fixed point under the bound-preparation contract.
- No constrained-search functional or variational energy proof.
- No full GAA joint mesh/mode/contact convergence certificate.
- No same-geometry DD solve: supplied U-channel topology, dimensions, channel
  acceptor doping and oxide differ (see DD_COMPARISON).
- No SIESTA binary compilation/run or atomistic published-device reproduction.
- No validated vertical/planar device family or missing high-field physics.

Every quantitative_validation_pass flag remains false. Failures, source/config
IDs, source snapshots and state residuals remain available for review.
'''
    (ROOT/'VALIDATION.md').write_text(validation)
    (R/'campaign_summary.json').write_text(json.dumps({'source_sha256':source_digest(),'counts':counts,
        'levels':{'1':'achieved','2':'equilibrium_only','3':'not_achieved_biased_MS','4':False,'5':False},
        'tests':tests['total_tests'],'no_AI_training':True,'no_git_mutations_or_hpc_jobs':True},indent=2))
    print(json.dumps({'tests':tests['total_tests'],'counts':counts},indent=2))


if __name__=='__main__':run()
