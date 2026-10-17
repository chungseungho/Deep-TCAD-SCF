"""Batched adaptive Gauss integration with a componentwise global error budget.

The q/2q difference is an estimator, not a rigorous missed-resonance bound.
Independent changes in initial panel width/order remain mandatory.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss


def integrate_vector(fun,breaks,order,rtol,atol,max_nodes=60000):
    lower=np.asarray(breaks[:-1]);upper=np.asarray(breaks[1:]);levels=np.zeros(len(lower),int)
    nodes=0
    def evaluate(a,b,q):
        nonlocal nodes
        x,w=leggauss(q)
        # Each lead band edge is a panel boundary. A sin^2 coordinate
        # removes its integrable 1/sqrt(E-Eedge) DOS singularity.
        angle=np.pi*(x+1)/4
        e=(a[:,None]+(b-a)[:,None]*np.sin(angle)**2).ravel();nodes+=len(e)
        jac=(b-a)[:,None]*(np.pi/4)*np.sin(2*angle)
        if nodes>max_nodes: raise RuntimeError(f'Energy integration exhausted {max_nodes} GF nodes; increase integration budget, not SCF tolerance')
        pieces=[]
        for start in range(0,len(e),192): pieces.append(fun(e[start:start+192]))
        values=np.concatenate(pieces).reshape(len(a),q,-1)
        return np.einsum('pqk,q,pq->pk',values,w,jac,optimize=True)
    low=evaluate(lower,upper,order);high=evaluate(lower,upper,2*order)
    errors=abs(high-low)
    for iteration in range(45):
        total=high.sum(axis=0);target=np.asarray(atol)+rtol*abs(total)
        global_error=errors.sum(axis=0)
        if np.all(global_error<=target):
            return total,{'gf_nodes':nodes,'adaptive_panels':len(lower),'maximum_normalized_error_estimate':float(np.max(global_error/target))}
        unresolved=global_error>target
        bad=np.any(errors[:,unresolved]>target[unresolved][None,:]/(2*len(lower)),axis=1)
        if not bad.any(): bad[np.argmax(np.max(errors/target,axis=1))]=True
        if np.max(levels[bad])>=44:
            worst=int(np.argmax(global_error/target))
            raise RuntimeError(f'Unresolved energy integral after 44 subdivisions: component={worst}, error estimate={global_error[worst]:g}, target={target[worst]:g}')
        middle=(lower[bad]+upper[bad])/2
        if np.any((middle==lower[bad])|(middle==upper[bad])):
            raise RuntimeError('UNRESOLVED_RESONANCE: energy floating-point resolution exhausted; no artificial broadening')
        a=np.r_[lower[bad],middle];b=np.r_[middle,upper[bad]];depth=np.tile(levels[bad]+1,2)
        ll=evaluate(a,b,order);hh=evaluate(a,b,2*order)
        lower=np.r_[lower[~bad],a];upper=np.r_[upper[~bad],b];levels=np.r_[levels[~bad],depth]
        high=np.concatenate([high[~bad],hh]);errors=np.concatenate([errors[~bad],abs(hh-ll)])
    raise RuntimeError('Energy quadrature did not meet its independent integration error budget')
