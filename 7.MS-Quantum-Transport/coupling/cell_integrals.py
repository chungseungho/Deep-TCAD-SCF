"""Exact geometric intersection of P1 tetrahedra and rectangular FD cells.

C_ij = integral_{tet mesh intersect quantum cell j} N_i(r) dV.
Each clipped convex polyhedron is integrated using its exact volume/centroid.
This yields positive loads, both volume marginals, linear precision and the
potential/charge work adjoint without a fitted kernel or Sinkhorn smoothing.
"""
from itertools import combinations
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.spatial import ConvexHull,cKDTree


TRIPLES=np.array(list(combinations(range(10),3)))
CORNERS=np.array(list(__import__('itertools').product([0.,1.],repeat=3)))


def polyhedron_volume_centroid(points):
    p=np.unique(np.round(points,12),axis=0)
    if len(p)<4: return 0.,None
    try: hull=ConvexHull(p)
    except Exception as exc:
        # A plane/edge contact has zero 3D measure. A genuine failed 3D hull
        # must not silently remove a physical intersection.
        if np.linalg.matrix_rank(p-p.mean(axis=0),tol=1e-10)<3: return 0.,None
        raise RuntimeError('Failed full-dimensional cell/tetra intersection') from exc
    centre=p.mean(axis=0);tri=p[hull.simplices]
    matrices=np.stack([tri[:,0]-centre,tri[:,1]-centre,tri[:,2]-centre],axis=-1)
    volume=abs(np.linalg.det(matrices))/6
    total=volume.sum()
    if total<1e-17: return 0.,None
    centroid=np.sum(volume[:,None]*(tri.sum(axis=1)+centre)/4,axis=0)/total
    return float(total),centroid


def integrated_matrix(vertices,tetra,nodes,nodevol,grid):
    tic=time.perf_counter();used=np.unique(tetra)
    dist,ids=cKDTree(nodes).query(vertices[used])
    if max(dist)>1e-7: raise ValueError('Tetra/Poisson nodes mismatch')
    convert=np.full(len(vertices),-1,int);convert[used]=ids;tet=convert[tetra]
    step=np.array([np.diff(grid.x)[0],np.diff(grid.y)[0],np.diff(grid.z)[0]])
    origin=np.array([grid.x[0],grid.y[0],grid.z[0]])-step/2
    shape=np.array(grid.mask.shape);index=np.full(tuple(shape),-1,int)
    index[grid.mask]=np.arange(grid.mask.sum())
    rr=[];cc=[];vv=[];intersections=0
    for ordinal,t in enumerate(tet):
        p=nodes[t];v0=p[0]
        M=np.column_stack([p[1]-v0,p[2]-v0,p[3]-v0]);inverse=np.linalg.inv(M)
        gradients=np.vstack([-inverse.sum(axis=0),inverse])
        tetvol=abs(np.linalg.det(M))/6
        lo=np.maximum(0,np.floor((p.min(axis=0)-origin)/step+1e-10).astype(int))
        hi=np.minimum(shape-1,np.ceil((p.max(axis=0)-origin)/step-1e-10).astype(int)-1)
        # Include every rectangular cell having a positive-measure intersection.
        if np.any(hi<lo): continue
        cells=np.array(np.meshgrid(*[np.arange(a,b+1) for a,b in zip(lo,hi)],indexing='ij')).reshape(3,-1).T
        js=index[tuple(cells.T)];keep=js>=0;cells=cells[keep];js=js[keep]
        if len(js)==1 and np.all(p>=origin+cells[0]*step-1e-10) and np.all(p<=origin+(cells[0]+1)*step+1e-10):
            rr.extend(t);cc.extend([js[0]]*4);vv.extend([tetvol/4]*4);intersections+=1;continue
        if not len(js): continue
        A=np.vstack([-gradients,np.eye(3),-np.eye(3)])
        norm=np.linalg.norm(A,axis=1);A=A/norm[:,None]
        bc=np.zeros((len(cells),10));bc[:,0]=1
        bc[:,4:7]=origin+(cells+1)*step-v0
        bc[:,7:10]=v0-origin-cells*step
        bc/=norm
        matrices=A[TRIPLES];valid=abs(np.linalg.det(matrices))>1e-12
        combos=TRIPLES[valid];inverses=np.linalg.inv(matrices[valid])
        candidates=np.einsum('tij,ctj->cti',inverses,bc[:,combos])
        inside=np.max(np.einsum('cti,ai->cta',candidates,A)-bc[:,None,:],axis=2)<3e-10
        for j,points in enumerate(candidates):
            volume,centre=polyhedron_volume_centroid(points[inside[j]])
            if volume==0: continue
            lam=np.r_[1.,0.,0.,0.]+gradients@centre
            if lam.min()<-1e-8: raise RuntimeError('Intersection centroid outside tetrahedron')
            lam=np.maximum(lam,0);lam/=lam.sum()
            rr.extend(t);cc.extend([js[j]]*4);vv.extend(volume*lam);intersections+=1
        if ordinal and ordinal%5000==0:
            print(f'Cell/tetra integration {ordinal}/{len(tet)} tetra, {time.perf_counter()-tic:.1f}s',flush=True)
    C=coo_matrix((vv,(rr,cc)),shape=(len(nodes),len(grid.xyz))).tocsr()
    sitevol=np.full(len(grid.xyz),grid.cell_volume)
    row_error=float(max(abs(np.asarray(C.sum(axis=1)).ravel()/nodevol-1)))
    col_error=float(max(abs(np.asarray(C.sum(axis=0)).ravel()/sitevol-1)))
    if max(row_error,col_error)>2e-7:
        raise RuntimeError(f'Geometric integration did not cover the actual Si volume: {row_error:g},{col_error:g}')
    return C,{'method':'integrated','balancing_error':max(row_error,col_error),
              'node_marginal_error':row_error,'cell_marginal_error':col_error,
              'intersections':intersections,'integration_runtime_s':time.perf_counter()-tic,
              'limitation':'exact P1 potential cell average / piecewise-constant quantum charge exchange; not a nonequilibrium variational theorem'}
