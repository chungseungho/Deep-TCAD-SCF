# Calculated results

These are actual finite-bath EMA-Poisson calculations. No analytic I-V fit,
invented data or failed-bias current is used as an accepted I-V point.

* `idvgs.csv`: separate screened-Anderson transfer run.
* `idvds.csv`: unscreened-Anderson output run with Newton-Krylov fallback.
* Full physical/numerical settings and actual command lines are saved separately.
* `summary.json` contains exact attempted/converged counts and failed biases.
* `states/` contains representative 3D snapshots, not every bias snapshot.
* `field_vtk/` is the transfer endpoint's actual DEVSIM 3D export.
* `mesh_iv/` is the mesh used for the I-V calculations; `mesh_default/` is the
  independently checked, finer default Poisson mesh.
* `convergence.csv` is a real single-bias SCF sensitivity campaign.
  Agreement of all parameters has NOT been achieved.
* `learning_example/` is explicitly PROTOTYPE training data and an RBF baseline,
  not the lab DeepSCF network or quantitatively certified teacher data.

All I-V here uses dx=1 nm, dy=dz=1 nm, 3 channel / 6 S/D modes per valley.
In particular the high-Vds output fall is affected by the narrow discrete FD
lead bandwidth and cannot be interpreted as physical Si NDR/saturation.
The code default dx=.25 nm has a wider band, but its FULL quantum I-V has not
been executed here. Its real Poisson/mapping checks were executed successfully.

Source/Drain neutrality, longitudinal/mode/basis/contact convergence remain
unresolved. `quantitative_validation_pass` is false throughout.
