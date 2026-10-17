"""Overlap tracking of states and degenerate subspaces, never injection labels."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from multispace.projection import degenerate_groups


def match_states(old,new,tolerance=1e-7,old_transform=None,new_transform=None,basis_overlap=None):
    # U columns must share the same physical inner product. For changed mode
    # bases supply the isometries to common real space; raw coefficients do not.
    u=old['u'] if old_transform is None else old_transform@old['u']
    v=new['u'] if new_transform is None else new_transform@new['u']
    if basis_overlap is None and u.shape[0]!=v.shape[0]:raise ValueError('Tracking requires a common physical grid')
    metric_v=v if basis_overlap is None else basis_overlap@v
    overlap=abs(u.conj().T@metric_v)**2
    i,j=linear_sum_assignment(-overlap)
    groups_old=degenerate_groups(old['e'],tolerance);groups_new=degenerate_groups(new['e'],tolerance)
    clusters=[]
    for a in groups_old:
        if a.stop-a.start<2:continue
        candidates=[]
        for b in groups_new:
            sv=np.linalg.svd(u[:,a].conj().T@metric_v[:,b],compute_uv=False)
            score=float(np.sum(sv**2)/max(a.stop-a.start,b.stop-b.start))
            candidates.append((score,b,sv))
        score,b,sv=max(candidates,key=lambda x:x[0])
        clusters.append({'old_range':[a.start,a.stop],'new_range':[b.start,b.stop],
            'normalized_projector_overlap':score,'singular_values':sv.tolist(),
            'individual_labels_defined':False})
    return {'old_indices':i,'new_indices':j,'overlap_squared':overlap[i,j],
        'energy_change_ev':new['e'][j]-old['e'][i],
        'degenerate_subspaces':clusters,'defines_injection_origin':False}
