"""Original-tetra affine potential interpolation + positive balanced charge map.

The charge matrix C has marginals (actual DEVSIM NodeVolume, quantum cell
volume). Therefore n_node = C n_cell / V_node preserves BOTH constants and
total charge. It is not the adjoint of the affine interpolation P; that work
defect is measured. No claim of a variational discrete energy is made.
"""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, coo_matrix, diags
from scipy.spatial import cKDTree


def load_tetra(path):
    lines=Path(path).read_text().splitlines();i=lines.index('$PhysicalNames');n=int(lines[i+1]);physical={}
    for row in lines[i+2:i+2+n]:
        dim,tag,name=row.split(maxsplit=2);physical[(int(dim),int(tag))]=name.strip('"')
    tag=next(t for (d,t),name in physical.items() if d==3 and name=='Silicon')
    i=lines.index('$Nodes');n=int(lines[i+1]);ns=[v.split() for v in lines[i+2:i+2+n]]
    ids={int(r[0]):i for i,r in enumerate(ns)};xyz=np.array([[float(x) for x in r[1:]] for r in ns])
    i=lines.index('$Elements');n=int(lines[i+1]);tet=[]
    for row in lines[i+2:i+2+n]:
        a=list(map(int,row.split()));nt=a[2]
        if a[1]==4 and a[3]==tag: tet.append([ids[v] for v in a[3+nt:]])
    return xyz,np.array(tet,int)


def affine_matrix(vertices,tetra,nodes,sites):
    used=np.unique(tetra);dist,indices=cKDTree(nodes).query(vertices[used])
    if max(dist)>1e-5 or len(np.unique(indices))!=len(used): raise ValueError('Gmsh/DEVSIM node mismatch')
    convert=np.full(len(vertices),-1,int);convert[used]=indices;tet=convert[tetra]
    v=nodes[tet];M=np.stack([v[:,1]-v[:,0],v[:,2]-v[:,0],v[:,3]-v[:,0]],axis=-1)
    inv=np.linalg.inv(M);tree=cKDTree(v.mean(axis=1));rows=[];cols=[];data=[]
    for start in range(0,len(sites),1000):
        p=sites[start:start+1000];unresolved=np.arange(len(p));chosen=np.full(len(p),-1);weights=np.zeros((len(p),4))
        for k in (16,64,256,1024):
            if not len(unresolved): break
            kk=min(k,len(tet));cand=tree.query(p[unresolved],k=kk)[1].reshape(len(unresolved),kk)
            for j in range(kk):
                if not len(unresolved): break
                t=cand[:,j];uvw=np.einsum('nij,nj->ni',inv[t],p[unresolved]-v[t,0])
                lam=np.column_stack([1-uvw.sum(axis=1),uvw]);ok=(lam.min(axis=1)>-1e-8)&(lam.max(axis=1)<1+1e-8)
                found=unresolved[ok];chosen[found]=t[ok];b=np.maximum(lam[ok],0)
                weights[found]=b/b.sum(axis=1,keepdims=True)
                unresolved=unresolved[~ok];cand=cand[~ok]
        if len(unresolved): raise ValueError(f'Quantum point not in original Si tetra: {p[unresolved[0]]}')
        rows.extend(np.repeat(np.arange(start,start+len(p)),4));cols.extend(tet[chosen].ravel());data.extend(weights.ravel())
    return csr_matrix((data,(rows,cols)),shape=(len(sites),len(nodes)))


def balanced_matrix(nodes,sites,nodevol,sitevol,max_iterations=40000,tol=2e-9):
    if not np.isclose(nodevol.sum(),sitevol.sum(),rtol=1e-7): raise ValueError('Geometry and quantum integration volumes differ')
    # Sparse Gaussian support, guaranteed at least one neighbor for every row/column.
    # The width is deliberately larger than one mesh cell to allow feasible marginals.
    k=min(48,len(sites));dist,cols=cKDTree(sites).query(nodes,k=k);dist=dist.reshape(len(nodes),k);cols=cols.reshape(len(nodes),k)
    rr=np.repeat(np.arange(len(nodes)),k);cc=cols.ravel();dd=dist.ravel()
    k2=min(12,len(nodes));d2,r2=cKDTree(nodes).query(sites,k=k2)
    rr=np.r_[rr,r2.ravel()];cc=np.r_[cc,np.repeat(np.arange(len(sites)),k2)];dd=np.r_[dd,d2.ravel()]
    sigma=max(float(np.median(dist[:,-1]))*.6,.5)
    K=coo_matrix((np.exp(-.5*(dd/sigma)**2),(rr,cc)),shape=(len(nodes),len(sites))).tocsr()
    u=np.ones(len(nodes));v=np.ones(len(sites));err=np.inf
    # Last scale removes only floating point CAD-volume mismatch, never physical charge.
    nv=nodevol*(sitevol.sum()/nodevol.sum())
    for it in range(max_iterations):
        u=nv/(K@v);v=sitevol/(K.T@u)
        if it%20==0:
            err=max(np.max(abs(u*(K@v)/nv-1)),np.max(abs(v*(K.T@u)/sitevol-1)))
            if err<tol: break
    if err>=tol: raise RuntimeError(f'Charge map balancing did not converge: {err:g}')
    return diags(u)@K@diags(v),{'balancing_error':float(err),'iterations':it+1,'kernel_sigma_nm':sigma}


@dataclass
class Mapping:
    P: object
    C: object
    nodevol: np.ndarray
    sitevol: np.ndarray
    nodes: np.ndarray
    sites: np.ndarray
    balance: dict

    @classmethod
    def build(cls,path,nodes,nodevol,grid,method='integrated'):
        v,t=load_tetra(path);P=affine_matrix(v,t,nodes,grid.xyz)
        sitevol=np.full(len(grid.xyz),grid.cell_volume)
        if method=='integrated':
            from coupling.cell_integrals import integrated_matrix
            C,report=integrated_matrix(v,t,nodes,nodevol,grid)
            P=diags(1/sitevol)@C.T
        elif method=='adjoint':
            C=P.T@diags(sitevol)
            report={'method':method,'balancing_error':None,
                    'limitation':'positive affine FE load quadrature; nodal constant-density reproduction is not exact'}
        else:
            C,report=balanced_matrix(nodes,grid.xyz,nodevol,sitevol)
            report['method']=method
            if method=='balanced_adjoint':
                P=diags(1/sitevol)@C.T
                report['limitation']='exact discrete exchange work and constants; finite spatial averaging, affine reproduction must be refined'
            else:
                report['limitation']='STEP05 non-adjoint legacy map; work defect remains'
        return cls(P,C,nodevol,sitevol,nodes,grid.xyz,report)
    def potential(self,phi): return np.asarray(self.P@phi)
    def charge(self,n):
        n=np.asarray(n)
        if np.any(n<0) or not np.isfinite(n).all(): raise ValueError('Invalid quantum density')
        return np.asarray(self.C@n)/self.nodevol
    def diagnostics(self):
        f=lambda x: 1.3+.03*x[:,0]-.04*x[:,1]+.01*x[:,2]
        rng=np.random.default_rng(73);n=rng.random(len(self.sites))*1e19
        mapped=self.charge(n);qa=np.dot(n,self.sitevol);qb=np.dot(mapped,self.nodevol)
        p=f(self.nodes);work_q=np.dot(self.potential(p)*n,self.sitevol);work_p=np.dot(p*mapped,self.nodevol)
        return {**self.balance,'linear_potential_error_v':float(max(abs(self.potential(p)-f(self.sites)))),
                'constant_density_relative_error':float(max(abs(self.charge(np.ones(len(n)))-1))),
                'charge_relative_error':float(abs(qa-qb)/qa),'minimum_mapped_density_cm3':float(mapped.min()),
                'nodevolume_nm3':float(self.nodevol.sum()),'quantum_volume_nm3':float(self.sitevol.sum()),
                'electrostatic_work_relative_defect':float(abs(work_q-work_p)/max(abs(work_q),1e-30)),
                'charge_measure':'DEVSIM configured node_volume_model=PhysicalControlVolume','potential_measure':'original silicon tetrahedra',
                'limitation':self.balance.get('limitation'),
                'work_adjoint_pass':bool(abs(work_q-work_p)<1e-12*max(abs(work_q),1e-30)),
                'constant_density_exact':bool(max(abs(self.charge(np.ones(len(n)))-1))<1e-7)}
