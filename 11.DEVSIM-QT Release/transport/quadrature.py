"""Constant-memory global-error adaptive Gauss integration.

The accepted componentwise error sum has the same global interpretation as the
legacy adaptive method. Parent panels are recomputed on subdivision to avoid
retaining panel x density-component arrays. No artificial broadening is added.
The q/2q difference is an estimator; independent resonance checks remain required.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss


def integrate_vector(fun,breaks,order,rtol,atol,max_nodes=60000):
    breaks=np.asarray(breaks);nodes=0;rules={}
    if len(breaks)<2 or np.any(np.diff(breaks)<=0):raise ValueError('Strictly increasing energy breaks required')
    minimum_nodes=4*order*(len(breaks)-1)
    if minimum_nodes>max_nodes:raise RuntimeError(f'Energy integration guidance requires at least {minimum_nodes} GF nodes, budget is {max_nodes}; increase budget or validate a more economical guidance rule')
    for q in (order,2*order):
        x,w=leggauss(q);a=np.pi*(x+1)/4;rules[q]=(np.sin(a)**2,w*np.pi/4*np.sin(2*a))
    def panel(a,b,q):
        nonlocal nodes
        if nodes+q>max_nodes:raise RuntimeError(f'Energy integration exhausted {max_nodes} GF nodes; increase integration budget, not SCF tolerance')
        x,w=rules[q];e=a+(b-a)*x;nodes+=q;value=None
        for first in range(0,q,16):
            last=min(first+16,q);part=np.einsum('e,ek->k',w[first:last]*(b-a),fun(e[first:last]),optimize=True)
            value=part if value is None else value+part
        return value
    estimate=None
    for a,b in zip(breaks[:-1],breaks[1:]):
        value=panel(a,b,order);estimate=value if estimate is None else estimate+value
    target=np.asarray(atol)+rtol*abs(estimate)
    total=np.zeros_like(estimate);errors=np.zeros_like(abs(estimate));leaves=[]
    def pair(a,b):
        low=panel(a,b,order);high=panel(a,b,2*order);return high,abs(high-low)
    for a,b in zip(breaks[:-1],breaks[1:]):
        high,err=pair(a,b);total+=high;errors+=err;leaves.append((a,b,0,float(max(err/target))))
    for iteration in range(90):
        final_target=np.asarray(atol)+rtol*abs(total)
        if np.all(errors<=final_target):
            return total,{'gf_nodes':nodes,'adaptive_panels':len(leaves),'initial_panels':len(breaks)-1,
                'maximum_normalized_error_estimate':float(max(errors/final_target)),
                'algorithm':'constant-memory globally budgeted q/2q; recompute parent before replacement',
                'energy_batch_limit':16,'stored_panel_density_vectors':0}
        # A scalar conservative envelope identifies panels needing attention;
        # acceptance always uses the actual componentwise GLOBAL error sum.
        scale=float(max(target/final_target));target=final_target
        leaves=[(a,b,d,score*scale) for a,b,d,score in leaves]
        bad=[i for i,(a,b,d,score) in enumerate(leaves) if score>1/(2*len(leaves)) and d<44 and a<(a+b)/2<b]
        if not bad:
            candidates=[(x[3],i) for i,x in enumerate(leaves) if x[2]<44 and x[0]<(x[0]+x[1])/2<x[1]]
            if not candidates:
                raise RuntimeError(f'UNRESOLVED_RESONANCE: global estimated error ratio {max(errors/target):g}; energy resolution exhausted')
            bad=[max(candidates)[1]]
        marked=set(bad);updated=[x for i,x in enumerate(leaves) if i not in marked]
        for i in bad:
            a,b,depth,score=leaves[i];old,old_error=pair(a,b);total-=old;errors-=old_error;mid=(a+b)/2
            for l,h in ((a,mid),(mid,b)):
                value,error=pair(l,h);total+=value;errors+=error;updated.append((l,h,depth+1,float(max(error/target))))
        leaves=updated
    raise RuntimeError('Energy quadrature failed its global error budget')
