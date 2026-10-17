# Minimal STEP09 interface changes

## Checkpoint metadata

STEP09 checkpoints stored density under `n` but omitted node coordinates. The existing `run.py --initial-state` loader expects `node_xyz_nm` and `poisson_density_cm3` to validate mesh identity and read density. A timed-out run therefore could not directly use its own checkpoint through that loader.

STEP10 adds those two fields to the two checkpoint-writing calls in `solver.py`. They contain the existing Poisson coordinates and exactly the same saved density; old `n`, potential, iteration, bias, config and source fields remain. Poisson equations, SCF evaluation order, mixing, Newton logic, residual thresholds, backend calls and acceptance equations were not changed. The full mock-SCF test checks both aliases against the actual stored arrays. Native tests were rerun after the addition.

For the earlier Vd=0.05 V timeout generated before this addition, `results/final_uniform_vd005_resume/convert_legacy_checkpoint.py` checks that the reference and run meshes have identical SHA256 and transfers only coordinate metadata from that reference; the saved density is not modified. The conversion inputs, iteration and mesh hash are retained. New checkpoints work directly with `--initial-state`.

## Uniform diagnostic geometry

`Config.geometry_profile` defaults to `raised`. The explicit `uniform` option accepts identical source/channel/drain silicon cross sections. The same geometry builder, quantum stencil, mapping and Poisson code then run on that declared geometry. Existing Raised S/D validation remains the default.

## Bound-pole audit

The audit's dense eigenvalue root search was replaced by block-LDL inertia with real bisection. This is backend numerical linear algebra, not a new occupation model. The same pole equation and rejection policy remain. Tests compare inertia with dense eigenvalue counts and the roots with analytic and frozen-dense references. No eta is introduced.
