"""Numerical convergence and model scope are independent judgments."""
import numpy as np
from config import KIN,KB_EV
from poisson.discretization import mesh_regions


def element_fields(mesh,p):
    xyz,tets=mesh_regions(mesh);phi=p._get('Potential')
    from scipy.spatial import cKDTree
    distance,indices=cKDTree(p.allxyz*1e-7).query(xyz)
    if max(distance)>1e-12:raise ValueError('Field mesh mismatch')
    out={}
    for name,tt in tets.items():
        v=xyz[tt]*1e7;values=phi[indices[tt]]
        M=np.stack([v[:,1]-v[:,0],v[:,2]-v[:,0],v[:,3]-v[:,0]],axis=1)
        grad=np.linalg.solve(M,(values[:,1:]-values[:,0,None])[...,None])[...,0]
        field=np.linalg.norm(grad,axis=1)
        out[name]={'maximum_v_nm':float(max(field)),'p99_v_nm':float(np.quantile(field,.99))}
    return out


def evaluate(c,systems,mu,vd,fields):
    dispersion=[]
    for s in systems:
        bottom=min(float(x[0]) for x in s.subbands)
        kinetic=max(0.,max(mu,mu-vd)+c.fermi_tail*KB_EV*c.temperature-bottom)
        bandwidth=4*KIN/(s.valley.mx*c.dx**2);fraction=kinetic/bandwidth
        dispersion.append({'valley':s.valley.name,'kinetic_window_ev':kinetic,
            'fd_bandwidth_ev':bandwidth,'fraction':fraction,
            'velocity_error_estimate':1-np.sqrt(max(0.,1-fraction))})
    high=any(x['maximum_v_nm']>c.high_field_v_nm for x in fields.values())
    band=all(x['fraction']<1 for x in dispersion)
    accurate=all(x['velocity_error_estimate']<=c.dispersion_tolerance for x in dispersion)
    flags=[]
    if high:flags.append('HIGH_FIELD_WARNING')
    if not accurate:flags.append('FD_DISPERSION_WARNING')
    if high or not band:flags.append('MODEL_OUT_OF_SCOPE')
    return {'physical_status':'MODEL_OUT_OF_SCOPE' if high or not band else 'PROVISIONAL_EMA',
        'flags':flags,'dispersion':dispersion,'element_fields':fields,
        'dispersion_tolerance_pass':accurate,'quantitative_physics_validated':False,
        'missing_physics':['oxide_tunneling','dielectric_breakdown','impact_ionization',
            'high_field_scattering','self_heating','nonparabolic_bands'],
        'field_threshold_v_nm':c.high_field_v_nm,'field_threshold_is':'screening diagnostic, not a breakdown prediction'}
