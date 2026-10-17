"""Inexact Newton with explicit admissible-domain line search.

Only declared physical-domain exceptions are backtracked. Programming and
quadrature errors propagate. Invalid derivatives are retried on smaller,
one-sided probes. A stalled search is not proof that no fixed point exists.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator,gmres


class DomainSearchFailure(RuntimeError):
    def __init__(self,x,residual,events):
        self.x=np.array(x);self.residual=np.array(residual);self.events=events
        super().__init__('DOMAIN_SEARCH_STALLED: no admissible decreasing trial; existence of a fixed point is undecided')


def solve(fun,x0,*,invalid=(),preconditioner=None,tol=1e-5,maxiter=30,
          trust_radius=.1,max_backtracks=16,callback=None,log=None,rdiff=1e-5):
    x=np.array(x0,dtype=float,copy=True);events=[]
    def record(event):
        events.append(event)
        if log:log(event)
    def evaluate(z,kind,iteration,scale):
        try:
            value=np.asarray(fun(z),dtype=float)
            if not np.isfinite(value).all():raise FloatingPointError('Nonfinite nonlinear residual')
            return value
        except invalid as e:
            event={'kind':kind,'iteration':iteration,'step_scale':float(scale),'valid':False,'reason':str(e)}
            if hasattr(e,'diagnostics'):event['transport_diagnostics']=e.diagnostics
            record(event);return None
    # An invalid initial iterate is reported as such, never as a failed solution.
    f=fun(x)
    for iteration in range(maxiter):
        if np.max(abs(f))<tol:
            if callback:callback(x,f)
            return x,{'converged':True,'iterations':iteration,'events':events}
        def jv(v):
            norm=np.linalg.norm(v)
            if norm==0:return np.zeros_like(v)
            h=rdiff*max(1.,np.linalg.norm(x))/norm
            for k in range(10):
                for sign in (1.,-1.):
                    ff=evaluate(x+sign*h*v,'derivative',iteration,sign*h)
                    if ff is not None:return (ff-f)/(sign*h)
                h*=.25
            raise DomainSearchFailure(x,f,events)
        J=LinearOperator((len(x),len(x)),matvec=jv,dtype=float)
        step,info=gmres(J,-f,M=preconditioner,rtol=.08,atol=0.,restart=12,maxiter=1)
        directions=[('newton',step),('preconditioned_residual',-(preconditioner@f) if preconditioner is not None else f)]
        accepted=False;old_merit=np.linalg.norm(f)
        for method,d in directions:
            scale=max(1.,np.max(abs(d))/trust_radius);d=d/scale
            for k in range(max_backtracks+1):
                alpha=2.**(-k);trial=x+alpha*d
                ft=evaluate(trial,'line_search',iteration,alpha)
                if ft is None:continue
                merit=np.linalg.norm(ft)
                good=merit<old_merit*(1-1e-4*alpha)
                record({'kind':'line_search','iteration':iteration,'step_scale':alpha,'valid':True,
                    'accepted':bool(good),'method':method,'gmres_info':int(info),
                    'residual_inf':float(max(abs(ft))),'residual_l2':float(merit)})
                if good:x,f=trial,ft;accepted=True;break
            if accepted:break
        if not accepted:raise DomainSearchFailure(x,f,events)
        if callback:callback(x,f)
    return x,{'converged':False,'iterations':maxiter,'events':events}
