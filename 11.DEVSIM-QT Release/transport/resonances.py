"""Pole-guided real integration: no lifetime or eta is added to the model.

Closed-device levels seed panels. Endpoint broadening estimates seed several
scales around their real embedding shifts. This reduces missed narrow peaks;
the independent adaptive budget must still pass. It is not a pole proof.
"""
import numpy as np
from scipy.linalg import eig_banded
from transport.negf_reference import surface_sigma


def resonance_breaks(onsite,hop,tx,leads,low,high,panel,order=8):
    sizes=[len(h) for h in onsite];off=np.r_[0,np.cumsum(sizes)]
    # Hermitian lower band storage: O(N * block width), no global dense H.
    bandwidth=max([len(h)-1 for h in onsite]+[len(onsite[i])+len(onsite[i+1])-1 for i in range(len(hop))])
    band=np.zeros((bandwidth+1,off[-1]))
    for i,h in enumerate(onsite):
        for j in range(len(h)):
            for k in range(j,len(h)):band[k-j,off[i]+j]=h[k,j]
    for i,t in enumerate(hop):
        for j in range(t.shape[0]):
            for k in range(t.shape[1]):band[off[i+1]+k-off[i]-j,off[i]+j]=t[j,k]
    # On the real axis ||Sigma|| <= tx for orthonormal rectangular overlaps.
    # Hence a closed level farther than tx cannot shift into this window.
    e,u=eig_banded(band,lower=True,select='v',select_range=(low-tx,high+tx),check_finite=False)
    centres=e.copy();width=np.zeros(len(e));shift=np.zeros(len(e))
    for idx,(spectrum,w) in enumerate(leads):
        endpoint=u[:sizes[0]] if idx==0 else u[-sizes[-1]:]
        amplitudes=w.conj().T@endpoint
        sig=surface_sigma(e,spectrum,tx)
        shift+=(abs(amplitudes.T)**2*sig.real).sum(axis=1)
        width+=(abs(amplitudes.T)**2*(-2*sig.imag)).sum(axis=1)
    centres+=shift
    # Broad peaks are already sampled by the closed/shifted centres and GL
    # panels. Reserve multiscale seeding for widths below the node scale.
    peaks=(centres>low)&(centres<high)&(width<panel/max(2,order))
    additions=[e[(e>low)&(e<high)],centres[(centres>low)&(centres<high)]]
    # Include distance scales much larger than a first-order width: shifts are
    # only estimates. Adaptive refinement, not this estimate, accepts the result.
    for scale in (.25,1.,4.,16.):
        additions.extend([centres[peaks]-scale*width[peaks],centres[peaks]+scale*width[peaks]])
    edges=np.concatenate([np.r_[a-2*tx,a+2*tx] for a,w in leads])
    breaks=np.unique(np.r_[low,np.arange(low,high,panel),high,*additions,edges])
    breaks=breaks[(breaks>=low)&(breaks<=high)]
    # Remove sub-ulp intervals rather than evaluating duplicate floating nodes.
    distance=np.diff(breaks);keep=np.r_[True,distance>32*np.finfo(float).eps*np.maximum(1,abs(breaks[1:]))]
    breaks=breaks[keep]
    if breaks[-1]<high:breaks=np.r_[breaks,high]
    return breaks,{'storage':'Hermitian block-band; no dense Hamiltonian','band_storage_bytes':band.nbytes,'embedding_shift_bound_ev':tx,'multiscale_width_threshold_ev':panel/max(2,order),'seed_levels':len(e),'narrow_peak_seeds':int(peaks.sum()),
        'minimum_seed_width_ev':float(width[peaks].min()) if peaks.any() else None,
        'algorithm':'closed-H pole/embedding-shift and endpoint-width guided panels; adaptive independent acceptance',
        'eta_added_ev':0.}
