"""Independent checks on exported tetrahedra and shared physical triangles."""
from collections import defaultdict
from pathlib import Path
import json
import numpy as np
from poisson.discretization import mesh_regions,assemble


def validate_mesh(path,c):
    points,tets=mesh_regions(path)
    from plotting import parse_mesh
    _,faces,_=parse_mesh(path)
    ring=(c.width+2*c.tox)*(c.thickness+2*c.tox)-c.width*c.thickness
    expected={'Silicon':(c.ls+c.ld)*c.sd_width*c.sd_thickness+(c.lg+2*c.spacer)*c.width*c.thickness,
              'GateOxide':c.lg*ring,'Spacer':2*c.spacer*ring}
    volume={k:float(assemble(points,t)[1].sum()) for k,t in tets.items()}
    incident=defaultdict(list)
    for material,tt in tets.items():
        for t in tt:
            for q in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]:
                incident[tuple(sorted(t[list(q)]))].append((material,points[t].mean(axis=0)))
    expected_materials={'SiOx':{'Silicon','GateOxide'},'SiSpacer':{'Silicon','Spacer'},'OxSpacer':{'GateOxide','Spacer'}}
    checks={};areas={}
    for group in ['Source','Drain','Gate','SiOx','SiSpacer','OxSpacer']:
        ts=np.array(faces[group],int);p=points[ts]
        areas[group]=float(np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1).sum()/2)
        incidents=[incident[tuple(sorted(t))] for t in ts]
        if group in expected_materials:
            checks[group]=all(len(v)==2 and {p[0] for p in v}==expected_materials[group] for v in incidents)
        else: checks[group]=all(len(v)==1 for v in incidents)
    gate_normals=[]
    for t in faces['Gate']:
        p=points[t];centre=p.mean(axis=0);normal=np.cross(p[1]-p[0],p[2]-p[0]);normal/=np.linalg.norm(normal)
        inside=incident[tuple(sorted(t))][0][1]
        if np.dot(normal,centre-inside)<0: normal=-normal
        axis=np.argmax(abs(normal));gate_normals.append((int(axis),int(np.sign(normal[axis]))))
    sides=sorted(set(gate_normals));gate_pass=sides==[(1,-1),(1,1),(2,-1),(2,1)]
    subregions={};tt=tets['Silicon'];x=points[tt].mean(axis=1)[:,0]
    for label,mask in [('source',x<c.ls),('channel_and_extensions',(x>c.ls)&(x<c.drain_start)),('drain',x>c.drain_start)]:
        subregions[label]=float(assemble(points,tt[mask])[1].sum())
    error=max(abs(volume[k]/expected[k]-1) for k in expected)
    report={'tetra_volumes_nm3':volume,'cad_expected_nm3':expected,'maximum_volume_relative_error':error,
            'interface_conformity':checks,'surface_areas_nm2':areas,'four_sided_outward_gate_normals':sides,
            'silicon_spatial_partitions_nm3':subregions,'nonmanifold_tetra_faces':sum(len(v)>2 for v in incident.values()),
            'gate_area_relative_error':abs(areas['Gate']/(2*(c.width+c.thickness+4*c.tox)*c.lg)-1),
            'gap_check':'conformal shared faces and CAD Boolean volume partition; no exhaustive floating-point pairwise overlap test',
            'pass':bool(error<1e-9 and all(checks.values()) and gate_pass and not any(len(v)>2 for v in incident.values()))}
    if not report['pass']: raise RuntimeError('Geometry validation failed: '+json.dumps(report))
    return report
