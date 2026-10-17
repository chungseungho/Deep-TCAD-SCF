"""Signed element-integrated Poisson flux and positive lumped cell volumes.

DEVSIM supports user-selected edge_couple_model and node_volume_model. We
assemble these from original conformal tetrahedra. The resulting operator is
the linear tetrahedron Galerkin stiffness (equivalently its conservative edge
balance), with positive barycentric control volumes V_T/4. This avoids the
over-counted absolute circumcentric NodeVolume of obtuse coarse tetrahedra.
Fluxes are antisymmetric on edges; the assembled stiffness is symmetric PSD.
Off-diagonal conductances may be negative on non-Delaunay meshes, so a discrete
maximum principle is NOT asserted. Dirichlet conditions remove the nullspace.
"""
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree


def mesh_regions(path):
    s=Path(path).read_text().splitlines();i=s.index('$PhysicalNames');n=int(s[i+1]);names={}
    for row in s[i+2:i+2+n]:
        dim,tag,name=row.split(maxsplit=2)
        if int(dim)==3: names[int(tag)]=name.strip('"')
    i=s.index('$Nodes');n=int(s[i+1]);raw=[r.split() for r in s[i+2:i+2+n]]
    ids={int(r[0]):j for j,r in enumerate(raw)};xyz=np.array([[float(v) for v in r[1:]] for r in raw])
    tets={name:[] for name in names.values()};i=s.index('$Elements');n=int(s[i+1])
    for row in s[i+2:i+2+n]:
        a=list(map(int,row.split()));nt=a[2]
        if a[1]==4: tets[names[a[3]]].append([ids[v] for v in a[3+nt:]])
    return xyz,{k:np.array(v,int) for k,v in tets.items()}


def assemble(nodes,tets):
    p=nodes[tets];M=np.stack([p[:,1]-p[:,0],p[:,2]-p[:,0],p[:,3]-p[:,0]],axis=-1)
    inv=np.linalg.inv(M);volume=abs(np.linalg.det(M))/6
    grads=np.concatenate([-inv.sum(axis=1)[:,None,:],inv],axis=1)
    local=np.einsum('tik,tjk,t->tij',grads,grads,volume)
    rows=np.repeat(tets,4,axis=1).ravel();cols=np.tile(tets,(1,4)).ravel()
    K=coo_matrix((local.ravel(),(rows,cols)),shape=(len(nodes),len(nodes))).tocsr()
    vol=np.bincount(tets.ravel(),weights=np.repeat(volume/4,4),minlength=len(nodes))
    return K,vol


def configure(d,device,mesh,regions):
    points,tets=mesh_regions(mesh);report={}
    for region,group in regions.items():
        args={'device':device,'region':region}
        nodes=np.column_stack([d.get_node_model_values(**args,name=a) for a in ['x','y','z']])
        used=np.unique(tets[group]);dist,idx=cKDTree(nodes).query(points[used])
        if max(dist)>1e-12: raise RuntimeError('Original tetra/DEVSIM region nodes mismatch')
        conversion=np.full(len(points),-1,int);conversion[used]=idx;tet=conversion[tets[group]]
        K,vol=assemble(nodes,tet)
        native=np.asarray(d.get_node_model_values(**args,name='NodeVolume'))
        d.node_solution(**args,name='PhysicalControlVolume');d.set_node_values(**args,name='PhysicalControlVolume',values=vol.tolist())
        d.edge_from_node_model(**args,node_model='node_index')
        i=np.asarray(d.get_edge_model_values(**args,name='node_index@n0'),int)
        j=np.asarray(d.get_edge_model_values(**args,name='node_index@n1'),int)
        conductance=-np.asarray(K[i,j]).ravel()
        coupling=conductance*np.linalg.norm(nodes[i]-nodes[j],axis=1)
        d.edge_solution(**args,name='SignedElementCouple')
        d.set_edge_values(**args,name='SignedElementCouple',values=coupling.tolist())
        d.set_parameter(**args,name='node_volume_model',value='PhysicalControlVolume')
        d.set_parameter(**args,name='edge_couple_model',value='SignedElementCouple')
        report[region]={'physical_volume_nm3':float(vol.sum()*1e21),
                        'native_NodeVolume_nm3':float(native.sum()*1e21),
                        'native_relative_volume_excess':float(native.sum()/vol.sum()-1),
                        'stiffness_row_sum_max_cm':float(max(abs(np.asarray(K.sum(axis=1)).ravel()))),
                        'signed_edge_couple_negative_fraction':float(np.mean(coupling<0))}
    return report
