"""Positive orthogonal spatial projectors, including subspace diagnostics."""
import numpy as np


def weights(vectors,regions):
    u=np.asarray(vectors);regions=np.asarray(regions)
    if u.shape[0]!=len(regions) or not np.isin(regions,[0,1,2]).all():
        raise ValueError('Invalid spatial partition')
    w=np.stack([(abs(u[regions==a])**2).sum(axis=0) for a in range(3)],axis=1)
    if max(abs(w.sum(axis=1)-1))>1e-9: raise ValueError('Projectors do not resolve normalized states')
    return w


def degenerate_groups(energies,tolerance):
    start=0
    for i in range(1,len(energies)+1):
        if i==len(energies) or energies[i]-energies[start]>tolerance:
            yield slice(start,i);start=i


def origin_spectrum(vectors,regions):
    d=(np.asarray(regions)==0).astype(float)-(np.asarray(regions)==2)
    return np.linalg.eigvalsh(vectors.conj().T@(d[:,None]*vectors))
