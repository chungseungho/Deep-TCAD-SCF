> Historical baseline document retained for research provenance. Its results and paths are not STEP11 validation evidence. See RESULTS.md and RELEASE_CERTIFICATE.json for this candidate.

# Supplied drift-diffusion source audit

The 19 supplied Python/GEO/README files were inspected read-only. They describe
a **U-shaped fin channel**, not STEP07's rectangular nanosheet. For example,
`13-u_channel_sim.geo` uses W=100 nm, H=80 nm, T=20 nm, Lg=80 nm,
Lsource=Ldrain=20 nm, spacers=8 nm and oxide=5 nm, with transport along z.
The three boxes forming the silicon U are explicitly joined by BooleanUnion.
The no-metal variant applies the gate contact directly to the oxide exterior.
The documented mesh is generated in nm and scaled to cm for DEVSIM.

In `02-u_channel_fin_idvgs_vds0p1_step0p05.py`, ND_SD=1e19 cm^-3 and
NA_CH=1e16 cm^-3; the channel net doping is **negative**. It uses DEVSIM
semiconductor continuity equations for electrons/holes and small-bias
transfer/output sweeps. Other files vary doping and sweep intervals.

STEP07 instead uses a 12×5 nm rectangular sheet, 24 nm gate, 2 nm oxide,
raised 20×15 nm source/drain, 3 nm spacers and an undoped channel by default.
It has electron-only hard-wall EMA and ballistic contacts. Geometry, channel
doping, oxide thickness and contact physics therefore do not match. No numeric
DD/quantum current comparison is presented, no unexecuted DD result is
inferred, and no EMA parameter is tuned to DD current. A same-geometry DD solve
is **not executed** in this deliverable. Extending the geometry interface to a
U-section and supporting acceptor/hole charge would be necessary for a direct
comparison; an axis permutation alone does not supply those physical changes.

The relevant qualitative questions for that future comparison are gate control,
sign and extent of accumulated/depleted charge, electrostatic boundary
conditions and profiles at identical bias. A DD versus ballistic current
difference would not itself diagnose an implementation error.
