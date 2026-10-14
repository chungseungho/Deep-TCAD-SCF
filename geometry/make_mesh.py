"""Build conformal Gmsh 3D mesh of one rectangular silicon nanosheet.

Transport x, sheet width y, sheet thickness z. Gate dielectric is a rectangular
annulus covering only the channel section; gate is a Dirichlet outer boundary,
NOT a volumetric metal mesh. Physical groups: Silicon, GateOxide,
Source, Drain, Gate, SiOx_Interface. Gmsh uses nm, ASCII msh2 scaled to cm
only AFTER all physical tags have been generated.
"""
from pathlib import Path


def scale_msh2_nm_to_cm(source: Path, target: Path):
    lines = source.read_text().splitlines()
    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        out.append(line)
        if line.strip() == '$Nodes':
            i += 1
            n = int(lines[i].strip())
            out.append(lines[i])
            for _ in range(n):
                i += 1
                fields = lines[i].split()
                if len(fields) != 4:
                    raise ValueError('Only ASCII Gmsh MSH 2.2 node lines supported')
                out.append(fields[0] + ' ' + ' '.join(f'{float(v)*1e-7:.17g}' for v in fields[1:]))
        i += 1
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text('\n'.join(out) + '\n')


def build_mesh(cfg, build_dir='build', quiet=False):
    import gmsh
    directory = Path(build_dir)
    directory.mkdir(parents=True, exist_ok=True)
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal', 0 if quiet else 1)
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.option.setNumber('Mesh.Binary', 0)
        gmsh.model.add('MonoNanosheetGAA')
        o = gmsh.model.occ
        ls, lg, ld = cfg.ls_nm, cfg.lg_nm, cfg.ld_nm
        w, t, ox = cfg.width_nm, cfg.thickness_nm, cfg.tox_nm
        xg = ls
        # Silica and Si are separate volumes, but all silicon belongs to one region.
        boxes = [
            o.addBox(0, -w/2, -t/2, ls, w, t),
            o.addBox(ls, -w/2, -t/2, lg, w, t),
            o.addBox(ls+lg, -w/2, -t/2, ld, w, t),
            o.addBox(xg, -w/2-ox, -t/2-ox, lg, ox, t+2*ox),
            o.addBox(xg, w/2, -t/2-ox, lg, ox, t+2*ox),
            o.addBox(xg, -w/2, -t/2-ox, lg, w, ox),
            o.addBox(xg, -w/2, t/2, lg, w, ox),
        ]
        # Fragment all touching solids so every Si/oxide interface is conformal.
        o.fragment([(3, boxes[0])], [(3, b) for b in boxes[1:]])
        o.synchronize()
        silicon, oxide = [], []
        for _, tag in gmsh.model.getEntities(3):
            _, y, z = o.getCenterOfMass(3, tag)
            if abs(y) < w/2 - 1e-5 and abs(z) < t/2 - 1e-5:
                silicon.append(tag)
            else:
                oxide.append(tag)
        assert silicon and oxide, 'Volume classification failed'
        labels = {}
        for v in silicon:
            labels[v] = 'Silicon'
        for v in oxide:
            labels[v] = 'GateOxide'
        contacts = {'Source': [], 'Drain': [], 'Gate': [], 'SiOx_Interface': []}
        for _, s in gmsh.model.getEntities(2):
            upward, _ = gmsh.model.getAdjacencies(2, s)
            types = {labels.get(int(v)) for v in upward}
            x, y, z = o.getCenterOfMass(2, s)
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, s)
            tol = 1.0e-4
            if 'Silicon' in types and 'GateOxide' in types:
                contacts['SiOx_Interface'].append(s)
            elif types == {'Silicon'}:
                if abs(x) < tol and abs(xmax-xmin) < tol:
                    contacts['Source'].append(s)
                elif abs(x-cfg.length_nm) < tol and abs(xmax-xmin) < tol:
                    contacts['Drain'].append(s)
            elif types == {'GateOxide'}:
                is_outer_y = abs(abs(y)-(w/2+ox)) < tol and abs(ymax-ymin) < tol
                is_outer_z = abs(abs(z)-(t/2+ox)) < tol and abs(zmax-zmin) < tol
                if is_outer_y or is_outer_z:
                    contacts['Gate'].append(s)
        if any(not s for s in contacts.values()):
            raise RuntimeError('Missing physical surfaces: ' + str({k:len(v) for k,v in contacts.items()}))
        def physical(dim, tags, name):
            group = gmsh.model.addPhysicalGroup(dim, sorted(set(tags)))
            gmsh.model.setPhysicalName(dim, group, name)
        physical(3, silicon, 'Silicon')
        physical(3, oxide, 'GateOxide')
        for name, tags in contacts.items():
            physical(2, tags, name)
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), cfg.mesh_h_nm)
        gmsh.option.setNumber('Mesh.MeshSizeMin', cfg.mesh_h_nm * 0.5)
        gmsh.option.setNumber('Mesh.MeshSizeMax', cfg.mesh_h_nm)
        gmsh.model.mesh.generate(3)
        nm_path = directory/'mono_nanosheet_nm.msh'
        cm_path = directory/'mono_nanosheet_cm.msh'
        gmsh.write(str(nm_path))
        scale_msh2_nm_to_cm(nm_path, cm_path)
        print('Physical groups:', {k:len(v) for k,v in contacts.items()})
        print('Mesh saved:', cm_path)
        return cm_path
    finally:
        gmsh.finalize()
