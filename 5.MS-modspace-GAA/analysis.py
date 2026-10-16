"""Explicit definitions for prototype I-V metrics; no fit to experimental data."""
import numpy as np


def transfer_metrics(rows,c,threshold_a_per_um=1e-7):
    good=[r for r in rows if r['status']=='SCF_CONVERGED_PROTOTYPE' and float(r['id_a'])>0]
    good.sort(key=lambda r:float(r['vg_v']))
    if len(good)<4: return {'status':'insufficient_converged_points'}
    v=np.array([float(r['vg_v']) for r in good]);i=np.array([float(r['id_a']) for r in good])
    gm=np.gradient(i,v);effective_width_um=2*(c.width+c.thickness)*1e-3
    target=threshold_a_per_um*effective_width_um
    crossings=np.where((i[:-1]<=target)&(i[1:]>target))[0]
    vth=None
    if len(crossings):
        j=crossings[0];vth=float(np.interp(np.log10(target),np.log10(i[j:j+2]),v[j:j+2]))
    slope=np.diff(np.log10(i))/np.diff(v)
    sub=(i[:-1]>1e-16)&(i[1:]<target)&(slope>0)
    ss=float(1000/np.max(slope[sub])) if sub.any() else None
    return {'status':'prototype_metrics_only','vth_constant_current_v':vth,
            'threshold_definition_a_per_um':threshold_a_per_um,'threshold_current_a':target,
            'normalizing_perimeter_um':effective_width_um,'ss_min_mv_dec':ss,
            'ss_definition':'smallest adjacent-sample swing below threshold and above 1e-16 A',
            'gm_max_a_per_v':float(gm.max()),'ion_a':float(i[-1]),'ioff_a':float(i[0]),
            'ion_ioff':float(i[-1]/i[0]),'vg_on_v':float(v[-1]),'vg_off_v':float(v[0]),
            'limitation':'uncalibrated finite-bath EMA; mesh/modes/contact convergence required'}
