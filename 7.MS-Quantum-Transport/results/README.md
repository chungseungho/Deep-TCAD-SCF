# Results index

- `campaign_summary.json`: executed counts and achievement status.
- `test_report/`: 56 completed tests with isolated native/quantum JUnit XML.
- `stage_a_transfer/`, `stage_a_output/`: all 100 requested MS Stage-A attempts, including NaN currents for rejected points.
- `gaa_ms_lowbias/`: finite-N equilibrium and explicit biased occupation failures.
- `gaa_ms_fixed_mu_equilibrium/`, `gaa_negf_equilibrium/`: matched-mu independently coupled equilibrium models.
- `matched_mu_coupled_comparison/`, `matched_mu_frozen_comparison/`: density, charge and barrier comparisons; not equivalent models.
- `gaa_negf_preconditioned/`, `gaa_negf_small_bias/`: actual unsuccessful biased Newton trials; bound-filling ambiguity retained.
- `gaa_negf_lowbias/`: earlier interrupted unconverged run, with explicit status and seed.
- `gaa_spectral_regression/`, `finite_current_check/`: separate finite-bath density/current and independent Landauer regression.
- `repaired_frozen_negf/`: no-eta semi-infinite density integration at the difficult STEP06 frozen field, two panel settings.
- `frozen_subspace/`: executed GAA frozen mode/basis study; not a coupled plateau.
- `reduced_realspace/`, `analytical/`: independent small-model numerical identities and transport tests.
- `stress_gate/`, `stress_drain/`: actual 0–6 V trial statuses and available-iterate fields, not high-voltage device predictions.
- `final_gauge_regression/`, `gauge_scf_comparison.json`, `resume_regression.json`: final-core gauge and immutable-state restart checks.
- `numerical_plan/`: 22 configurations, explicitly unexecuted.
- `source_snapshots/`: five core snapshots with verified source hashes; historical results retain their original version.
- `inputs/`: copied prior-state inputs and their identity; no import from STEP06 is needed.
- `figures/`: scientific figures generated only from these saved data.
- `dataset_catalog.json`, `failed_artifact_catalog.json`, `package_audit.json`: raw artifact context and integrity audit.

Source SHA is the hash of sorted production Python path/content pairs, excluding
results/tests/validation/environments. The package manifest additionally hashes
all delivered source, validation, documentation and result payloads. Historical
core code reconstructed from exact edits is accepted only when its recomputed
SHA equals the originally recorded execution SHA. No private SIESTA source is
redistributed. Development-error logs are retained but not counted as passes.
