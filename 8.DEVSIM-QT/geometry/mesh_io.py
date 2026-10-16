from pathlib import Path
import numpy as np

def parse_mesh(path):
    lines=Path(path).read_text().splitlines();i=lines.index('$PhysicalNames');n=int(lines[i+1]);names={}
    for row in lines[i+2:i+2+n]:
        dim,tag,name=row.split(maxsplit=2);names[(int(dim),int(tag))]=name.strip('"')
    i=lines.index('$Nodes');n=int(lines[i+1]);raw=[r.split() for r in lines[i+2:i+2+n]]
    byid={int(r[0]):j for j,r in enumerate(raw)};xyz=np.array([[float(v) for v in r[1:]] for r in raw])
    faces={};tets={};i=lines.index('$Elements');n=int(lines[i+1])
    for row in lines[i+2:i+2+n]:
        a=list(map(int,row.split()));typ,nt=a[1:3];dim={2:2,4:3}.get(typ)
        if dim is None: continue
        name=names.get((dim,a[3]),'unlabelled');ind=[byid[v] for v in a[3+nt:]]
        (faces if dim==2 else tets).setdefault(name,[]).append(ind)
    return xyz,faces,tets

