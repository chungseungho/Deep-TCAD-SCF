"""Pole-guided real integration: no lifetime or eta is added to the model.

Closed-device levels seed panels. Endpoint broadening estimates seed several
scales around their real embedding shifts. This reduces missed narrow peaks;
the independent adaptive budget must still pass. It is not a pole proof.
"""
import numpy as np
from scipy.linalg import eigh
from transport.negf_reference import surface_sigma


def resonance_breaks(onsite,hop,tx,leads,low,high,panel):
    sizes=[len(h) for h in onsite];off=np.r_[0,np.cumsum(sizes)]
    H=np.zeros((off[-1],off[-1]))
    for i,h in enumerate(onsite):H[off[i]:off[i+1],off[i]:off[i+1]]=h
    for i,t in enumerate(hop):
        H[off[i]:off[i+1],off[i+1]:off[i+2]]=t
        H[off[i+1]:off[i+2],off[i]:off[i+1]]=t.conj().T
    # Only levels relevant to the real bias integral need eigenvectors.
    e,u=eigh(H,subset_by_value=(low-2*tx,high+2*tx),check_finite=False)
    centres=e.copy();width=np.zeros(len(e));shift=np.zeros(len(e))
    for idx,(spectrum,w) in enumerate(leads):
        endpoint=u[:sizes[0]] if idx==0 else u[-sizes[-1]:]
        amplitudes=w.conj().T@endpoint
        sig=surface_sigma(e,spectrum,tx)
        shift+=(abs(amplitudes.T)**2*sig.real).sum(axis=1)
        width+=(abs(amplitudes.T)**2*(-2*sig.imag)).sum(axis=1)
    centres+=shift
    peaks=(centres>low)&(centres<high)&(width<panel)
    additions=[e[(e>low)&(e<high)],centres[peaks]]
    # Include distance scales much larger than a first-order width: shifts are
    # only estimates. Adaptive refinement, not this estimate, accepts the result.
    for scale in (.125,.5,2.,8.,32.,128.):
        additions.extend([centres[peaks]-scale*width[peaks],centres[peaks]+scale*width[peaks]])
    edges=np.concatenate([np.r_[a-2*tx,a+2*tx] for a,w in leads])
    breaks=np.unique(np.r_[low,np.arange(low,high,panel),high,*additions,edges])
    breaks=breaks[(breaks>=low)&(breaks<=high)]
    # Remove sub-ulp intervals rather than evaluating duplicate floating nodes.
    distance=np.diff(breaks);keep=np.r_[True,distance>32*np.finfo(float).eps*np.maximum(1,abs(breaks[1:]))]
    breaks=breaks[keep]
    if breaks[-1]<high:breaks=np.r_[breaks,high]
    return breaks,{'seed_levels':len(e),'narrow_peak_seeds':int(peaks.sum()),
        'minimum_seed_width_ev':float(width[peaks].min()) if peaks.any() else None,
        'algorithm':'closed-H pole/embedding-shift and endpoint-width guided panels; adaptive independent acceptance',
        'eta_added_ev':0.}
