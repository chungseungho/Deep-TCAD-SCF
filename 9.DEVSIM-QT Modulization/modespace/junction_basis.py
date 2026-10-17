"""Nested junction trace enrichment from STEP06; independently tested here."""
import numpy as np
from scipy.linalg import orth,eigh
from config import KIN
from modespace.hamiltonian import ModeBuilder,System,transverse_matrix

class JunctionCompatibleBuilder(ModeBuilder):
    def build(self,phi):
        original=super().build(phi);c,g,v=self.c,self.g,self.v
        n=len(g.x);channel=np.where((g.x>=c.ls)&(g.x<=c.drain_start))[0]
        targets={0:original.bases[channel[0]],2:original.bases[channel[-1]]}
        bases=[];onsite=[];subbands=[];regions=[];tx=KIN/(v.mx*c.dx**2)
        for i,x in enumerate(g.x):
            mask=g.mask[i];region=0 if x<c.ls else 2 if x>c.drain_start else 1
            U=original.bases[i].reshape(-1,original.bases[i].shape[-1])
            if region!=1:
                narrow=targets[region].reshape(-1,targets[region].shape[-1])
                U=orth(np.column_stack([U,narrow]),rcond=1e-11)
            u=U[mask.ravel()];m=u.shape[1]
            transverse=transverse_matrix(mask,KIN/(v.my*c.dy**2),KIN/(v.mz*c.dz**2))
            ht=u.T@transverse@u-u.T@(g.slice_values(phi,i)[:,None]*u)
            left=g.mask[i-1] if i else np.zeros_like(mask);right=g.mask[i+1] if i+1<n else np.zeros_like(mask)
            wall=(~left).astype(float)+(~right).astype(float)
            h=ht+2*tx*np.eye(m)+tx*(u.T@(wall[mask,None]*u))
            bases.append(U.reshape(*mask.shape,m));onsite.append(h)
            subbands.append(eigh(ht,eigvals_only=True,check_finite=False));regions.extend([region]*m)
        off=np.r_[0,np.cumsum([len(h) for h in onsite])];H=np.zeros((off[-1],off[-1]));hop=[]
        for i,h in enumerate(onsite):H[off[i]:off[i+1],off[i]:off[i+1]]=h
        for i in range(n-1):
            a=bases[i].reshape(-1,bases[i].shape[-1]);b=bases[i+1].reshape(-1,bases[i+1].shape[-1]);t=-tx*a.T@b;hop.append(t)
            H[off[i]:off[i+1],off[i+1]:off[i+2]]=t;H[off[i+1]:off[i+2],off[i]:off[i+1]]=t.T
        return System(H,bases,subbands,off,np.array(regions),g,hop,v)

