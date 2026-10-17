"""Cell-centred hard-wall EMA and its Galerkin projection on slice modes.

All slices share one physical transverse cell grid. Missing silicon is zero
extension, NOT a same-size identity coupling across a raised-S/D junction.
At hard-wall faces, half-cell Dirichlet distances add one hopping to diagonal.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh, block_diag
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh
from config import KIN


@dataclass
class Grid:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    mask: np.ndarray
    xyz: np.ndarray
    counts: np.ndarray
    offsets: np.ndarray
    cell_volume: float

    @classmethod
    def build(cls,c):
        c.validate()
        x=(np.arange(round(c.length/c.dx))+.5)*c.dx
        y=-c.sd_width/2+(np.arange(round(c.sd_width/c.dy))+.5)*c.dy
        z=-c.sd_thickness/2+(np.arange(round(c.sd_thickness/c.dz))+.5)*c.dz
        X,Y,Z=np.meshgrid(x,y,z,indexing='ij')
        mask=(X<c.ls)|(X>c.drain_start)|((abs(Y)<c.width/2)&(abs(Z)<c.thickness/2))
        xyz=np.column_stack([X[mask],Y[mask],Z[mask]])
        counts=mask.sum(axis=(1,2));offsets=np.r_[0,np.cumsum(counts)]
        return cls(x,y,z,mask,xyz,counts,offsets,c.dx*c.dy*c.dz)

    def slice_values(self,v,i): return np.asarray(v)[self.offsets[i]:self.offsets[i+1]]


def transverse_matrix(mask,ty,tz):
    idx=np.full(mask.shape,-1,int);idx[mask]=np.arange(mask.sum())
    rows=[];cols=[];values=[];diag=np.zeros(mask.sum())
    for j,k in zip(*np.where(mask)):
        a=idx[j,k]
        for dj,dk,t in [(1,0,ty),(-1,0,ty),(0,1,tz),(0,-1,tz)]:
            jj,kk=j+dj,k+dk
            neighbor=0<=jj<mask.shape[0] and 0<=kk<mask.shape[1] and mask[jj,kk]
            diag[a]+=t if neighbor else 2*t
            if neighbor: rows.append(a);cols.append(idx[jj,kk]);values.append(-t)
    rows.extend(range(len(diag)));cols.extend(range(len(diag)));values.extend(diag)
    return csr_matrix((values,(rows,cols)),shape=(len(diag),len(diag)))


@dataclass
class System:
    H: np.ndarray | None
    bases: list
    subbands: list
    offsets: np.ndarray
    regions: np.ndarray
    grid: Grid
    hop: list
    valley: object

    onsite: list | None = None
    contacts: object = None

    @property
    def dimension(self): return int(self.offsets[-1])
    def diagonal_blocks(self):
        return [self.H[self.sl(i),self.sl(i)] for i in range(len(self.grid.x))] if self.H is not None else self.onsite
    def dense(self):
        """Explicit opt-in for finite MS/Spectral or small reference calculations."""
        if self.H is None:
            H=np.zeros((self.dimension,self.dimension))
            for i,h in enumerate(self.onsite):H[self.sl(i),self.sl(i)]=h
            for i,t in enumerate(self.hop):
                H[self.sl(i),self.sl(i+1)]=t;H[self.sl(i+1),self.sl(i)]=t.conj().T
            self.H=H
        return self.H
    def sl(self,i): return slice(self.offsets[i],self.offsets[i+1])
    def realspace_transform(self):
        """Only for independent small-mesh verification, not production."""
        return block_diag(*[u[m] for u,m in zip(self.bases,self.grid.mask)])


class ModeBuilder:
    def __init__(self,c,grid,valley):
        self.c,self.g,self.v=c,grid,valley;self.cache={};self.fixed={};self.config_signature=c.to_dict()

    def build(self,phi):
        c,g,v=self.c,self.g,self.v
        if c.to_dict()!=self.config_signature:raise ValueError('Configuration changed; create a new Hamiltonian provider to invalidate caches')
        if np.asarray(phi).shape!=(len(g.xyz),): raise ValueError('Invalid potential shape')
        tx=KIN/(v.mx*c.dx**2);ty=KIN/(v.my*c.dy**2);tz=KIN/(v.mz*c.dz**2)
        bases=[];subbands=[];onsite=[];regs=[]
        for i,x in enumerate(g.x):
            mask=g.mask[i];key=mask.tobytes()
            if key not in self.cache: self.cache[key]=transverse_matrix(mask,ty,tz)
            ht=self.cache[key]-diags(g.slice_values(phi,i))
            m=min(c.sd_modes if x<c.ls or x>c.drain_start else c.modes,mask.sum())
            reference=self.cache[key] if c.basis=='fixed' else ht
            if c.basis=='fixed' and (key,m) in self.fixed:
                e,u=self.fixed[(key,m)]
            elif m>=mask.sum()-1:
                e,u=eigh(reference.toarray(),subset_by_index=(0,m-1),check_finite=False)
            else:
                # Deterministic start; sparse eigenproblem includes ALL retained modes.
                # A constant start can lie entirely in the even-parity sector
                # and miss low odd modes in a symmetric rectangular slice.
                e,u=eigsh(reference,k=m,which='SA',tol=1e-11,
                          v0=np.random.default_rng(195).normal(size=mask.sum()))
                order=np.argsort(e);e,u=e[order],u[:,order]
            if c.basis=='fixed':
                self.fixed[(key,m)]=(e,u)
                transverse=np.diag(e)-u.T@(g.slice_values(phi,i)[:,None]*u)
                energies=eigh(transverse,eigvals_only=True,check_finite=False)
            else: transverse=np.diag(e);energies=e
            full=np.zeros((*mask.shape,m));full[mask]=u
            # Arbitrary signs/degenerate rotations cancel between hopping and density.
            bases.append(full);subbands.append(energies)
            left=g.mask[i-1] if i else np.zeros_like(mask)
            right=g.mask[i+1] if i+1<len(g.x) else np.zeros_like(mask)
            wall=(~left).astype(float)+(~right).astype(float)
            h=transverse+2*tx*np.eye(m)+tx*(u.T@(wall[mask,None]*u))
            onsite.append(h);regs.extend([0 if x<c.ls else 2 if x>c.drain_start else 1]*m)
        off=np.r_[0,np.cumsum([len(e) for e in subbands])]
        hop=[]
        for i in range(len(g.x)-1):
            u=bases[i].reshape(-1,bases[i].shape[-1]);w=bases[i+1].reshape(-1,bases[i+1].shape[-1])
            hop.append(-tx*(u.T@w))
        system=System(None,bases,subbands,off,np.array(regs),g,hop,v,onsite,getattr(self,'contacts',None))
        if not getattr(self,'block_only',False):system.dense()
        return system



def density_from_matrix(sys,rho,degeneracy):
    n=[]
    for i,mask in enumerate(sys.grid.mask):
        u=sys.bases[i][mask];r=rho[sys.sl(i),sys.sl(i)]
        values=np.einsum('am,mn,an->a',u,r,u.conj()).real
        if values.min() < -1e-9: raise RuntimeError('Non-positive carrier density')
        n.append(np.maximum(values,0)*degeneracy/sys.grid.cell_volume*1e21)
    return np.concatenate(n)
