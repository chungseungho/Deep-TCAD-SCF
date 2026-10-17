"""Backend-owned SCF diagnostics, final current checks, and exports."""
import numpy as np


def validate_density(result):
    if result.electron_density_cm3.ndim != 1 or not np.all(np.isfinite(result.electron_density_cm3)):
        raise RuntimeError('Transport backend returned invalid electron density')


def metrics(backend,previous,current,config):
    if previous is None:
        return backend.initial_convergence_metrics()
    a=np.asarray(previous['tracking']);b=np.asarray(current['tracking'])
    occ=float(max(abs(b-a))) if a.shape==b.shape and a.size else float('inf')
    stable=backend.state_stability(previous,current)[0]
    ir=backend.current_relative_change(previous,current)
    return {'occupation_max_change':occ,'classification_stable':stable,
            'current_relative_change':float(ir)}


def ms_state_stability(old,new,config):
    from multispace.tracking import match_states
    from scipy.linalg import block_diag
    stable=True;reports=[]
    for a,b,sa,sb in zip(old['states'],new['states'],old['systems'],new['systems']):
        overlaps=[ua.reshape(-1,ua.shape[-1]).T@ub.reshape(-1,ub.shape[-1])
                  for ua,ub in zip(sa.bases,sb.bases)]
        match=match_states(a,b,config.ms_degenerate_ev,basis_overlap=block_diag(*overlaps))
        i,j=match['old_indices'],match['new_indices'];active=a['active'][i]|b['active'][j]
        stable=stable and bool(np.all(a['origin'][i][active]==b['origin'][j][active]))
        reports.append(match)
    return stable,reports


def ms_calibrate(config,context):
    from transport.negf import neutral_lead_mu
    return neutral_lead_mu(context['builders'],config,context['grid'])

def negf_calibrate(config,context):
    from transport.negf import neutral_lead_mu
    return neutral_lead_mu(context['builders'],config,context['grid'])

def spectral_calibrate(config,context):
    from transport.reservoir import neutral_finite_mu
    return neutral_finite_mu(context['builders'],config,context['grid'],spectral=True)


def ms_finalize(config,result,bias,grid,accepted):
    posts=[];error=None
    if accepted:
        from transport.ms_postprocess import postprocess
        try:
            posts=[postprocess(s,config,result['mu'],bias) for s in result['systems']]
            result['id_a']=sum(x['current_a'] for x in posts)
            result['bonds']=sum((x['bonds_a'] for x in posts),np.zeros(len(grid.x)-1))
        except Exception as exc:error=repr(exc)
    defect=max(abs(result['bonds']-result['id_a'])) if posts else float('nan')
    return {'postprocess':posts,'postprocess_error':error,'current_conservation_abs_a':float(defect),
            'current_conservation_relative':float(defect/max(abs(result['id_a']),1e-15)),
            'current_ok':bool(np.isfinite(result['id_a']) and np.isfinite(defect) and defect<max(1e-15,abs(result['id_a'])*1e-6)),
            'transport_status':'POST_SCF_FAILED' if error else 'CALCULATED' if accepted and np.isfinite(defect) and defect<max(1e-15,abs(result['id_a'])*1e-6) else 'UNVALIDATED'}


def open_finalize(config,result,bias,grid,accepted):
    central=(grid.x[:-1]>=config.ls)&(grid.x[1:]<=config.drain_start)
    defect=max(abs(result['bonds'][central]-result['id_a']))
    ok=bool(np.isfinite(result['id_a']) and defect<max(1e-15,abs(result['id_a'])*1e-6))
    return {'postprocess':[],'postprocess_error':None,'current_conservation_abs_a':float(defect),
        'current_conservation_relative':float(defect/max(abs(result['id_a']),1e-15)),
        'current_ok':ok,'transport_status':'CALCULATED' if accepted and ok else 'UNVALIDATED'}


def ms_export(result):
    out={'density_from_l_cm3':result['density_l'],'density_from_r_cm3':result['density_r'],
         'density_unresolved_origin_cm3':result['density_unresolved']}
    for i,s in enumerate(result['systems']):
        out[f'valley{i}_subbands_ev']=np.stack([a[:min(result['config'].modes,len(a))] for a in s.subbands])
        for key in ('e','weights','origin','largest_region','occupation','ambiguous','reach','contrast','boundary_unstable','cluster_unstable'):
            out[f'valley{i}_{key}']=result['states'][i][key]
    return out


def open_export(result):
    out={'density_from_l_cm3':result['density_l'],'density_from_r_cm3':result['density_r'],
         'density_unresolved_origin_cm3':result['density_unresolved']}
    for i,s in enumerate(result['systems']):
        out[f'valley{i}_subbands_ev']=np.stack([a[:min(result['config'].modes,len(a))] for a in s.subbands])
    return out
