"""Anderson type-II with volume weighted least squares and positivity guard."""
import numpy as np


def screened_inverse(K,capacitance,indices,boundary):
    """Approximate inverse of -I-K^-1 C on the silicon potential space.

    This is only a Krylov preconditioner. The nonlinear residual always uses
    the original quantum density and unscreened physical Poisson operator.
    """
    from scipy.sparse import diags
    from scipy.sparse.linalg import splu,LinearOperator
    free=np.flatnonzero(~np.asarray(boundary));indices=np.asarray(indices)
    C=np.asarray(capacitance);factor=splu((K+diags(C))[free][:,free].tocsc())
    def apply(v):
        rhs=np.zeros(K.shape[0]);rhs[indices]=C[indices]*v
        response=np.zeros(K.shape[0]);response[free]=factor.solve(rhs[free])
        return -np.asarray(v)+response[indices]
    return LinearOperator((len(indices),len(indices)),matvec=apply,dtype=float)


class Anderson:
    def __init__(self,volume,mix=.15,history=5):
        self.weight=np.sqrt(volume/np.mean(volume));self.mix=mix;self.base_mix=mix
        self.history=history;self.x=[];self.r=[];self.previous=None

    def update(self,x,F):
        r=F-x;norm=np.linalg.norm(self.weight*r)
        if self.previous is not None and norm>self.previous*1.4:
            self.mix=max(self.base_mix/8,self.mix*.5);self.x=[];self.r=[]
        elif self.previous is not None and norm<self.previous*.8:
            self.mix=min(self.base_mix,self.mix*1.1)
        self.previous=norm;self.x.append(x.copy());self.r.append(r.copy())
        self.x=self.x[-self.history-1:];self.r=self.r[-self.history-1:]
        proposal=x+self.mix*r;method='linear'
        if len(self.x)>1 and self.history>0:
            dx=np.column_stack(np.diff(self.x,axis=0));dr=np.column_stack(np.diff(self.r,axis=0))
            A=dr*self.weight[:,None];b=r*self.weight
            # Small Tikhonov term limits ill-conditioned multi-secant directions.
            scale=max(np.linalg.norm(A),1e-20)
            coef=np.linalg.lstsq(np.vstack([A,np.eye(A.shape[1])*scale*1e-5]),
                                 np.r_[b,np.zeros(A.shape[1])],rcond=1e-10)[0]
            candidate=proposal-(dx+self.mix*dr)@coef
            if np.isfinite(candidate).all() and candidate.min()>=0 and np.linalg.norm((candidate-x)*self.weight)<5*norm:
                proposal=candidate;method='anderson'
        return np.maximum(proposal,0),method,self.mix
