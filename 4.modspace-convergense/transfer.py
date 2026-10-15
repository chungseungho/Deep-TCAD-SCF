"""Tetrahedral (Gmsh MSH2) -> Cartesian quantum-grid conservative transfer.

Potential: P @ phi_tet, using *original Gmsh tetrahedron* barycentric weights.
Charge:    M_tet^{-1} P.T M_grid @ n_grid (mass-lumped weak projection).
Conserves the same geometric discrete total charge on the two sides. This does
NOT prove equality to DEVSIM internal control-volume weights (verify separately).
The finite difference grid does not refine the source Poisson solution itself.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree


def load_si_tetra_msh2(path):
    """Return (node_xyz_nm, tetra_indices) for original physical Silicon tets."""
    path=Path(path)
    with path.open() as f:
        s=f.read().splitlines()
    try:
        i=s.index('$PhysicalNames'); n=int(s[i+1]); names={}
        for line in s[i+2:i+2+n]:
            spl=line.split(maxsplit=2)
            names[(int(spl[0]),int(spl[1]))]=spl[2].strip('"')
        si_id=next(tag for (dim,tag),name in names.items() if dim==3 and name=='Silicon')
        i=s.index('$Nodes'); num=int(s[i+1]); lines=s[i+2:i+2+num]
        node_ids=np.fromiter((int(v.split()[0]) for v in lines),dtype=np.int64,count=num)
        coord=np.array([[float(a) for a in v.split()[1:4]] for v in lines],float)
        i=s.index('$Elements');num=int(s[i+1]);lines=s[i+2:i+2+num]
        tets=[]
        for line in lines:
            vals=line.split(); et=int(vals[1]);ntags=int(vals[2]);
            if et==4 and ntags>=1 and int(vals[3])==si_id:
                tets.append([int(v) for v in vals[3+ntags:3+ntags+4]])
        if not tets:raise ValueError('No physical Silicon tetrahedra in mesh')
        byid={int(k):i for i,k in enumerate(node_ids)}
        tets=np.array([[byid[k] for k in t] for t in tets],np.int32)
        if coord.max()<1e-3:coord=coord*1e7 # cm -> nm; ordinary device coordinates are 48 nm
        return coord,tets
    except (ValueError,IndexError,KeyError,StopIteration) as e:
        raise ValueError(f'Unable to parse Silicon MSH2 original connectivity in {path}: {e}') from e


@dataclass
class Transfer:
    P: csr_matrix
    node_volume_nm3: np.ndarray
    site_volume_nm3: np.ndarray
    axes: tuple
    nodes_nm: np.ndarray
    tetra_count: int
    @classmethod
    def from_tetra(cls,all_points_nm,tets,node_xyz_nm,axes,chunk_size=2000):
        vertices=np.asarray(all_points_nm,float)
        tets=np.asarray(tets,int)
        node_xyz_nm=np.asarray(node_xyz_nm,float)
        if len(vertices)==0 or tets.shape[1]!=4:raise ValueError('Expected tetrahedral mesh')
        # Restrict volumes to original Si tetra vertices; map into DEVSIM Si node order.
        relevant=np.unique(tets.ravel())
        kd_nodes=cKDTree(node_xyz_nm)
        distances,indices=kd_nodes.query(vertices[relevant])
        if np.max(distances)>1e-5:
            raise ValueError('Gmsh tetrahedral vertices do not match DEVSIM silicon nodes '
                             f'(max mismatch {np.max(distances):.4g} nm)')
        idmap=np.full(len(vertices),-1,dtype=int);idmap[relevant]=indices
        if len(np.unique(indices))!=len(indices):
            raise ValueError('More than one Gmsh silicon vertex mapped to a single DEVSIM node')
        si_tets=idmap[tets]
        vv=node_xyz_nm[si_tets]
        # Inverse of tetra affine map from barycentric coordinate (u,v,w) to xyz.
        M=np.stack([vv[:,1]-vv[:,0],vv[:,2]-vv[:,0],vv[:,3]-vv[:,0]],axis=-1)
        det=np.linalg.det(M)
        if np.any(np.abs(det)<1e-12): raise ValueError('Degenerate source tetrahedra found')
        inv=np.linalg.inv(M)
        volumes=np.abs(det)/6.0
        nodevol=np.bincount(si_tets.ravel(),weights=np.repeat(volumes/4,4),minlength=len(node_xyz_nm))
        if np.any(nodevol<=0): raise ValueError('Silicon mesh has nodes without positive tetra volume')
        mid=vv.mean(axis=1); centroids=cKDTree(mid)
        x,y,z=[np.asarray(v,float) for v in axes]
        shp=(len(x),len(y),len(z))
        dx=np.mean(np.diff(x));dy=np.mean(np.diff(y));dz=np.mean(np.diff(z))
        wx=np.ones(len(x));wx[0]=wx[-1]=.5
        sitevol=np.broadcast_to(wx[:,None,None]*dx*dy*dz,shp).copy().ravel()
        # Cartesian grid is ordered [x,y,z], with explicit interval subdivisions.
        rows=[];cols=[];weights=[]
        npts=int(np.prod(shp));n_tet=len(tets)
        for a in range(0,npts,chunk_size):
            ix=np.arange(a,min(a+chunk_size,npts))
            xyz=np.column_stack([x[ix//(len(y)*len(z))],
                     y[(ix//len(z))%len(y)],z[ix%len(z)]])
            chosen=np.full(len(ix),-1,int)
            bary=np.zeros((len(ix),4))
            unresolved=np.arange(len(ix))
            for k in (16,64,256):
                if not len(unresolved):break
                kk=min(k,n_tet)
                _,cand=centroids.query(xyz[unresolved],k=kk)
                cand=np.asarray(cand).reshape(len(unresolved),kk)
                # Boundary points can belong to >1 simplex. First valid original tetra wins.
                for j in range(kk):
                    if not len(unresolved):break
                    tet_idx=cand[:,j]
                    uvw=np.einsum('nij,nj->ni',inv[tet_idx],xyz[unresolved]-vv[tet_idx,0],optimize=True)
                    lam=np.column_stack([1-uvw.sum(axis=1),uvw])
                    hit=np.all(lam>=-3e-8,axis=1)&np.all(lam<=1+3e-8,axis=1)
                    if np.any(hit):
                        original=unresolved[hit]
                        chosen[original]=tet_idx[hit]
                        keep=np.maximum(lam[hit],0)
                        bary[original]=keep/keep.sum(axis=1,keepdims=True)
                    unresolved=unresolved[~hit]
                    cand=cand[~hit]
            if len(unresolved):
                raise ValueError(f'{len(unresolved)} FD sites outside ORIGINAL tetrahedra or unlocated. '
                    f'First: {xyz[unresolved[0]]} nm. Verify Si mesh connectivity/bounds.')
            rows.extend(np.repeat(ix,4).tolist())
            cols.extend(si_tets[chosen].ravel().tolist())
            weights.extend(bary.ravel().tolist())
        P=csr_matrix((weights,(rows,cols)),shape=(npts,len(node_xyz_nm)))
        obj=cls(P,nodevol,sitevol,(x,y,z),node_xyz_nm,n_tet)
        obj.check_structure()
        return obj
    @classmethod
    def from_msh(cls,mesh_path,si_node_positions_nm,axes):
        points,tetra=load_si_tetra_msh2(mesh_path)
        return cls.from_tetra(points,tetra,si_node_positions_nm,axes)
    @property
    def shape(self):return tuple(len(a) for a in self.axes)
    @property
    def quantum_charge_volume_nm3(self):return float(np.sum(self.site_volume_nm3))
    @property
    def mesh_volume_nm3(self):return float(np.sum(self.node_volume_nm3))
    def check_structure(self):
        row_sums=np.asarray(self.P.sum(axis=1)).ravel()
        err=np.max(np.abs(row_sums-1))
        if err>1e-9 or self.P.data.min() < -1e-12:
            raise RuntimeError(f'Invalid tetra interpolation: row-sum error={err:g}')
    def to_regular(self,phi_nodes):
        phi=np.asarray(phi_nodes,float)
        if phi.shape!=(len(self.node_volume_nm3),):raise ValueError('Invalid Si potential shape')
        return (self.P@phi).reshape(self.shape)
    def to_nodes(self,density_grid):
        arr=np.asarray(density_grid,float)
        if arr.shape!=self.shape:raise ValueError('Invalid quantum density shape')
        if not np.all(np.isfinite(arr)) or arr.min()<-1e-7:raise ValueError('Invalid density values')
        return np.asarray(self.P.T@(arr.ravel()*self.site_volume_nm3)).ravel()/self.node_volume_nm3
    def totals(self,density_grid,density_nodes):
        a=np.sum(np.asarray(density_grid).ravel()*self.site_volume_nm3)
        b=np.sum(np.asarray(density_nodes).ravel()*self.node_volume_nm3)
        return {'grid_particles':float(a*1e-21),'tet_particles':float(b*1e-21),
                'relative_charge_mismatch':float(abs(a-b)/max(abs(a),1e-30))}
    def diagnostics(self):
        xyz=self.nodes_nm
        phi=1.3 + .03*xyz[:,0]-.04*xyz[:,1]+.01*xyz[:,2]
        x,y,z=self.axes
        X,Y,Z=np.meshgrid(x,y,z,indexing='ij')
        truth=1.3+.03*X-.04*Y+.01*Z
        err=np.max(np.abs(self.to_regular(phi)-truth))
        v=np.ones(self.shape)*2e17
        totals=self.totals(v,self.to_nodes(v))
        return {'linear_phi_error_V':float(err),**totals,
                'tetra_count':self.tetra_count,'fd_grid_points':int(self.P.shape[0]),
                'silicon_node_count':int(self.P.shape[1]),
                'tet_mesh_volume_nm3':self.mesh_volume_nm3,
                'interior_fd_coverage_nm3':self.quantum_charge_volume_nm3,
                'caveat':'Geometric P.T mass lumping; not proven equal to DEVSIM internal dual-volume integration.'}
