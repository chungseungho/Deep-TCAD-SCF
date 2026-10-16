"""Conservative memory planning, never a hardware-specific runtime guarantee."""
def estimate(c):
    nx=round(c.length/c.dx);ns=round((c.ls+c.ld)/c.dx)
    m=min(c.modes,round(c.width/c.dy)*round(c.thickness/c.dz))
    sd=min(c.sd_modes,round(c.sd_width/c.dy)*round(c.sd_thickness/c.dz))
    dimension=ns*sd+(nx-ns)*m
    dense=dimension**2*16
    batch=256*nx*sd**2*16
    return {'backend':c.backend,'slices':nx,'maximum_valley_dimension':dimension,
            'one_dense_complex_matrix_mib':dense/2**20,
            'estimated_quantum_peak_mib':(32*dense if c.backend!='negf' else 5*dense+4*batch)/2**20,
            'dense_runtime_scaling':'O(N^3) per valley per evaluation' if c.backend!='negf' else 'RGF O(Nx*Nm^3*NE); no dense diagonalization in SCF',
            'runtime_estimate':'measure one point on target CPU; allow 100-160 electronic evaluations; walltime depends on energy/mode refinement',
            'hpc_submission':'not performed; specify partition/resources according to local rules'}
