"""Conformal OCC raised S/D + suspended silicon + oxide ring + SiN spacers.

Metal is the continuous four-sided outer oxide Dirichlet surface, not a volume.
Physical groups and analytic CAD volumes are checked before meshing.
"""
from pathlib import Path
import json
import numpy as np


def scale_mesh(source, target):
    s=Path(source).read_text().splitlines(); i=s.index('$Nodes'); n=int(s[i+1])
    for j in range(i+2,i+2+n):
        f=s[j].split(); s[j]=f[0]+' '+' '.join(f'{float(v)*1e-7:.17g}' for v in f[1:])
    Path(target).write_text('\n'.join(s)+'\n')


def build_mesh(c, directory, quiet=True):
    import gmsh
    c.validate(); out=Path(directory); out.mkdir(parents=True, exist_ok=True)
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal', int(not quiet))
        gmsh.option.setNumber('General.NumThreads', 1)
        gmsh.option.setNumber('Mesh.MshFileVersion',2.2)
        gmsh.option.setNumber('Mesh.Binary',0)
        gmsh.model.add('RaisedSD_GAA'); o=gmsh.model.occ
        def box(x,l,w,t): return o.addBox(x,-w/2,-t/2,l,w,t)
        # Boolean cuts produce disjoint materials; fragment makes shared faces conformal.
        left=box(0,c.ls,c.sd_width,c.sd_thickness)
        bridge=box(c.ls,2*c.spacer+c.lg,c.width,c.thickness)
        right=box(c.drain_start,c.ld,c.sd_width,c.sd_thickness)
        oxide_outer=box(c.gate_start,c.lg,c.width+2*c.tox,c.thickness+2*c.tox)
        oxide_inner=box(c.gate_start,c.lg,c.width,c.thickness)
        oxide=o.cut([(3,oxide_outer)],[(3,oxide_inner)])[0]
        spacers=[]
        for x in (c.ls,c.gate_end):
            outer=box(x,c.spacer,c.width+2*c.tox,c.thickness+2*c.tox)
            inner=box(x,c.spacer,c.width,c.thickness)
            spacers.extend(o.cut([(3,outer)],[(3,inner)])[0])
        vols=[(3,left),(3,bridge),(3,right)]+oxide+spacers
        _,fragments=o.fragment(vols[:1],vols[1:]); o.synchronize()
        labels={}; groups={'Silicon':[],'GateOxide':[],'Spacer':[]}
        # An annulus centroid lies inside its HOLE: never classify by centroid.
        input_labels=['Silicon']*3+['GateOxide']*len(oxide)+['Spacer']*len(spacers)
        for name,pieces in zip(input_labels,fragments):
            for dim,v in pieces:
                if dim!=3: continue
                if v in labels and labels[v]!=name: raise RuntimeError('Materials overlap')
                labels[v]=name
        for v,label in labels.items(): groups[label].append(v)
        ring=(c.width+2*c.tox)*(c.thickness+2*c.tox)-c.width*c.thickness
        expected={'Silicon':(c.ls+c.ld)*c.sd_width*c.sd_thickness+(c.lg+2*c.spacer)*c.width*c.thickness,
                  'GateOxide':c.lg*ring,'Spacer':2*c.spacer*ring}
        actual={k:sum(o.getMass(3,v) for v in vs) for k,vs in groups.items()}
        for k in expected:
            if not np.isclose(expected[k],actual[k],rtol=1e-8): raise RuntimeError(f'CAD volume mismatch {k}')
        faces={k:[] for k in ['Source','Drain','Gate','SiOx','SiSpacer','OxSpacer','SiExterior','SpacerExterior']}
        exterior=[];gate_sides=set()
        for _,s in gmsh.model.getEntities(2):
            up,_=gmsh.model.getAdjacencies(2,s); mats={labels[int(v)] for v in up}
            x,y,z=o.getCenterOfMass(2,s)
            bb=gmsh.model.getBoundingBox(2,s)
            if mats=={'Silicon','GateOxide'}: faces['SiOx'].append(s)
            elif mats=={'Silicon','Spacer'}: faces['SiSpacer'].append(s)
            elif mats=={'GateOxide','Spacer'}: faces['OxSpacer'].append(s)
            elif len(up)==1:
                exterior.append(s)
                if mats=={'Silicon'} and abs(x)<1e-5 and bb[3]-bb[0]<1e-5: faces['Source'].append(s)
                elif mats=={'Silicon'} and abs(x-c.length)<1e-5 and bb[3]-bb[0]<1e-5: faces['Drain'].append(s)
                elif mats=={'Silicon'}: faces['SiExterior'].append(s)
                elif mats=={'Spacer'}: faces['SpacerExterior'].append(s)
                elif mats=={'GateOxide'}:
                    for axis,coord,dim in [(1,y,c.width),(2,z,c.thickness)]:
                        if abs(abs(coord)-(dim/2+c.tox))<1e-5 and bb[axis+3]-bb[axis]<1e-5:
                            faces['Gate'].append(s);gate_sides.add((axis,int(np.sign(coord))))
        if len(gate_sides)!=4 or any(not v for v in faces.values()):
            raise RuntimeError(f'Incomplete GAA physical surfaces: {faces}, sides={gate_sides}')
        tags={}
        for dim, gs in ((3,groups),(2,faces)):
            for name,vs in gs.items():
                tag=gmsh.model.addPhysicalGroup(dim,sorted(set(vs)))
                gmsh.model.setPhysicalName(dim,tag,name);tags[name]={'dimension':dim,'tag':tag,'entities':vs}
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0),c.mesh_h)
        fields=[]
        def refine(x0,x1,y0,y1,z0,z1,h):
            f=gmsh.model.mesh.field.add('Box')
            for k,v in dict(VIn=h,VOut=c.mesh_h,XMin=x0,XMax=x1,YMin=y0,YMax=y1,ZMin=z0,ZMax=z1,Thickness=1.).items():
                gmsh.model.mesh.field.setNumber(f,k,v)
            fields.append(f)
        refine(c.ls-1,c.drain_start+1,-c.width/2-c.tox,c.width/2+c.tox,
               -c.thickness/2-c.tox,c.thickness/2+c.tox,c.channel_h)
        for x in (c.ls,c.gate_start,c.gate_end,c.drain_start):
            refine(x-1,x+1,-c.width/2-c.tox,c.width/2+c.tox,
                   -c.thickness/2-c.tox,c.thickness/2+c.tox,c.junction_h)
        f=gmsh.model.mesh.field.add('Min');gmsh.model.mesh.field.setNumbers(f,'FieldsList',fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(f)
        gmsh.option.setNumber('Mesh.MeshSizeMin',min(c.channel_h,c.junction_h))
        gmsh.option.setNumber('Mesh.MeshSizeMax',c.mesh_h)
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0)
        gmsh.model.mesh.generate(3)
        elem=gmsh.model.mesh.getElements(3)[1];nt=sum(len(e) for e in elem)
        quality=gmsh.model.mesh.getElementQualities(np.concatenate(elem).tolist(),'minSICN')
        if np.min(quality)<=0: raise RuntimeError('Inverted tetrahedron')
        # CAD + msh retained for geometry inspection, cm copy ONLY for DEVSIM.
        gmsh.write(str(out/'gaa.brep'));gmsh.write(str(out/'gaa_nm.msh'))
        scale_mesh(out/'gaa_nm.msh',out/'gaa_cm.msh')
        report={'physical_groups':tags,'cad_volume_nm3':actual,'expected_volume_nm3':expected,
                'gate_sides':sorted(gate_sides),'tetrahedra':nt,'minimum_minSICN':float(min(quality)),
                'external_surface_count':len(exterior),'gate_model':'continuous outer oxide surface; no metal volume'}
        (out/'geometry.json').write_text(json.dumps(report,indent=2))
        return out/'gaa_cm.msh'
    finally: gmsh.finalize()
