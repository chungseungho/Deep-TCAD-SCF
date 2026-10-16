# STEP07 — Multi-Space-Inspired EMA Quantum Transport

**구현과 실제 계산은 수행했지만, 검증된 비평형 MS Quantum TCAD에는 아직 도달하지 못했다.**
평형 3D GAA MS–Poisson SCF는 수렴했다. 기본 MS 스윕의 비평형 97점은
모두 점유 기원 모호성으로 거부됐다. 같은 화학퍼텐셜로 비교해도 MS/NEGF
밀도 차이가 커서 Level4·5 또는 물리적 기준해를 주장하지 않는다.

This is a standalone CPU research implementation. **MS means Multi-Space**;
mode space is a separate dimension reduction. A finite EMA/Hartree occupation
solver is not the original constrained-search MS-DFT functional.

## Implemented paths

| Backend | Density during electronic SCF | Current |
|---|---|---|
| A: `ms` | Finite Hermitian eigenstates, positive spatial projectors, published channel-origin rule with explicit rejection; finite-N or fixed-contact-μ ensemble | Semi-infinite GF/RGF **after** convergence |
| B: `negf` | Open semi-infinite contour/injectivity density, resonance-guided real integration and a below-lead bound-pole audit | Same open scattering model inside SCF |
| C: `spectral` | Analytic finite-wide-band spectral density | Finite-bath bond current, independently checked by Schur/RGF Landauer integration |

The GAA geometry, shared-potential 3D DEVSIM P1 electrostatics, conservative
adjoint cell/tetra exchange, silicon valley masses and unrescaled Hamiltonian
projection descend from STEP06. The main new work is the actual A density path,
five-region diagnostics, number/occupation constraints, degeneracy and boundary
rejection, post-SCF separation, open-density comparison and narrow-integral
regression. The nested basis fixes a missing junction trace; it does not prove
mode completeness. No AI, surrogate or DeepSCF training is included.

## Read the evidence

- [STATUS.md](STATUS.md): achievement levels and blocking issues.
- [VALIDATION.md](VALIDATION.md): executed cases, exact metrics and unexecuted work.
- [MS_IMPLEMENTATION_STUDY.md](MS_IMPLEMENTATION_STUDY.md): practical published workflow, actual source occupation loop and limits.
- [SIESTA_TRANSPORT_CODE_STUDY.md](SIESTA_TRANSPORT_CODE_STUDY.md): actual supplied-source paths; standard and modified routes distinguished.
- [THEORY.md](THEORY.md): units, boundary conditions, physical and numerical contracts.
- [HOW_TO_RUN.md](HOW_TO_RUN.md): private environment, CLI, testing, continuation, stress and user-submitted job template.
- [DD_COMPARISON.md](DD_COMPARISON.md): why the supplied U-channel DD device is not the same geometry/doping benchmark.

`results/` contains actual CSV/NPZ/VTK/mesh/log/figure outputs, including rejected
points. Missing currents stay NaN. The transfer/output status figures show
failure locations, not invented I–V curves. Source snapshots reproduce the
executed core versions; the package manifest hashes every delivered payload.

## Minimal execution

```bash
.venv/bin/python validation/run_tests.py
.venv/bin/python run.py point --config configs/gaa_ms.json --vg .8 --vd 0 --out results/my_ms
.venv/bin/python run.py point --config configs/gaa_negf.json --vg .8 --vd .1 --out results/my_negf
```

See HOW_TO_RUN for native dependencies. Both gate and drain range parsers
support 0–6 V. Numerical completion and physical validity are separate fields;
the stress results are not validated high-voltage characteristics.

## Architecture and boundaries

`device/` defines coordinates, contacts and five regions; `modespace/` builds
the projected Hamiltonian; `multispace/` classifies and fills finite states;
`poisson/` and `coupling/` solve electrostatics and exchange density;
`transport/` contains the three contact/density models; `validation/` holds
independent checks and campaign drivers. Ordered transport blocks and a world
coordinate frame support extension work. Only the nanosheet family was run;
a physical vertical transistor is not validated.

STEP01–06 were preserved byte-for-byte against the recorded baseline. No git
commit/push/branch/PR, school server access or batch-job submission was made.
