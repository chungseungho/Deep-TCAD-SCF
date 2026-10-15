"""Electrical characterization / defensible extraction rules for STEP03."""
from __future__ import annotations
import numpy as np


def voltages(vmax=6.0,step=0.2):
    if step<=0 or vmax<=0: raise ValueError('vmax and step must be positive')
    arr=np.arange(int(round(vmax/step))+1,dtype=float)*step
    if abs(arr[-1]-vmax)>1e-8: raise ValueError('vmax must be multiple of step')
    return np.round(arr,10)


def extract_vth(vg,ids,width_nm,iref_a_per_um=1e-7):
    """Constant-current extraction, with actual geometry-based current reference.

    Does not invent Vth when gate sweep never crosses Iref.
    """
    x=np.asarray(vg,float); y=np.abs(np.asarray(ids,float))
    if len(x)<3 or np.any(~np.isfinite(y)) or np.any(np.diff(x)<=0):
        raise ValueError('Invalid Id-Vgs points')
    iref=iref_a_per_um * (width_nm / 1000.0)
    ratio=np.log10(np.maximum(y,1e-30)/iref)
    candidates=np.flatnonzero((ratio[:-1]<=0)&(ratio[1:]>=0))
    if len(candidates)==0:
        return {'vth_v':None,'iref_a':float(iref),'method':'constant_current_log_interp',
                'status':'NO_REFERENCE_CURRENT_CROSSING',
                'reason':'I_D(V_G) never crosses reference from below.'}
    ix=int(candidates[0]); y0,y1=ratio[ix:ix+2]
    vt=float(x[ix]+(x[ix+1]-x[ix])*(-y0)/(y1-y0)) if y1!=y0 else float(x[ix])
    return {'vth_v':vt,'iref_a':float(iref),'method':'constant_current_log_interp',
            'status':'OK','bracket_v':[float(x[ix]),float(x[ix+1])]}


def choose_vgs(vth,step=0.2,vmax=6.0,below=0.4,above=0.4,strong=1.4):
    """Return [below Vth, above Vth, well above Vth] at user-selected step.

    Fail rather than silently claim that no available voltage is 'below'.
    """
    if vth is None or not np.isfinite(vth):
        raise ValueError('No defensible Vth: cannot choose below/above/strong-on Vgs.')
    picks=np.array([vth-below,vth+above,vth+strong])
    picks=np.clip(np.rint(picks/step)*step,0,vmax)
    if len(set(np.round(picks,7)))!=3 or not (picks[0]<vth<picks[1]<picks[2]):
        raise ValueError(f'Vth={vth:.3f} V too close to sweep bounds for 3 distinct Vgs; '
                         'use --manual-vth only if physically independently justified.')
    return np.round(picks,10).tolist()


def extract_ss(vg,ids,min_decades=1.5):
    """Find strongest linear semilog interval spanning >= min_decades.

    Avoid reporting a numeric SS when the data have no credible subthreshold
    exponential regime (common with uncalibrated contact/workfunction models).
    """
    x=np.asarray(vg,float); y=np.log10(np.maximum(np.abs(ids),1e-30))
    best=None
    for i in range(len(x)-3):
        for j in range(i+3,len(x)):
            if y[j]-y[i]<min_decades: continue
            vals=y[i:j+1]
            slopes=np.diff(vals)
            if np.any(slopes<=0): continue
            p=np.polyfit(x[i:j+1],vals,1)
            residual=vals-np.polyval(p,x[i:j+1])
            score=float(np.sqrt(np.mean(residual**2)))
            if p[0]>0 and score<0.08 and (best is None or score<best['fit_rms_dec']):
                best={'ss_mv_dec':float(1000/p[0]),'vg_window_v':[float(x[i]),float(x[j])],
                      'fit_rms_dec':score}
    return best


def summarize_iv(vg,ids,width_nm):
    x=np.asarray(vg,float); y=np.asarray(ids,float)
    vth=extract_vth(x,y,width_nm)
    ss=extract_ss(x,y)
    gm=np.gradient(y,x)
    return {'vth':vth,'ss':ss,'ion_a':float(np.max(np.abs(y))),
            'ioff_a':float(abs(y[0])),'max_gm_a_per_v':float(np.max(gm))}
