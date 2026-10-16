"""ParaView VTU voxels with cell-centred quantum fields and explicit units."""
from pathlib import Path
from itertools import product
import json
import numpy as np


def export_state(state,destination):
    s=np.load(state);xyz=s['xyz_nm'];x,y,z=s['x_nm'],s['y_nm'],s['z_nm']
    step=np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    # VTK hexahedron ordering (not the different voxel ordering).
    corners=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    nodes=(xyz[:,None,:]+(corners-.5)*step).reshape(-1,3)
    points,conn=np.unique(np.round(nodes,10),axis=0,return_inverse=True);conn=conn.reshape(-1,8)
    fields={'Potential_V':s['phi_v'],'Ec_eV':-s['phi_v'],'Electrons_cm3':s['electron_density_cm3'],
            'Donors_cm3':s['doping_quantum_cm3']}
    for name in ['density_from_l_cm3','density_from_r_cm3','equilibrium_density_cm3','bias_injected_density_cm3']:
        if name in s: fields[name]=s[name]
    path=Path(destination);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w') as f:
        f.write('<?xml version="1.0"?><VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian"><UnstructuredGrid>')
        f.write(f'<Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(xyz)}"><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">')
        np.savetxt(f,points,fmt='%.10g');f.write('</DataArray></Points><Cells>')
        for name,data,typ in [('connectivity',conn,'Int64'),('offsets',np.arange(1,len(xyz)+1)*8,'Int64'),('types',np.full(len(xyz),12),'UInt8')]:
            f.write(f'<DataArray type="{typ}" Name="{name}" format="ascii">');np.savetxt(f,np.atleast_2d(data),fmt='%d');f.write('</DataArray>')
        f.write('</Cells><CellData>')
        for name,data in fields.items():
            f.write(f'<DataArray type="Float64" Name="{name}" format="ascii">');np.savetxt(f,data,fmt='%.10g');f.write('</DataArray>')
        f.write('</CellData></Piece></UnstructuredGrid></VTKFile>')
    return path


def export_last(out):
    out=Path(out);files=[p for p in sorted((out/'states').glob('vg*.npz')) if not p.name.endswith('_spectra.npz')]
    valid=[p for p in files if bool(np.load(p)['converged'])]
    if valid: return export_state(valid[-1],out/'quantum_fields_nm.vtu')


def restore_last_poisson(p,out):
    """Restore the same accepted state used by quantum export, including resume.

    This solves only the unchanged linear Poisson problem for its saved input
    charge. The silicon potential must reproduce the saved field before any
    native artifact is exported. DEVSIM's mesh coordinates remain in cm.
    """
    out=Path(out);valid=[]
    for file in sorted((out/'states').glob('vg*.npz')):
        if file.name.endswith('_spectra.npz'):continue
        with np.load(file) as state:
            if bool(state['converged']):valid.append(file)
    if not valid:return None
    file=valid[-1]
    with np.load(file) as state:
        nodes=state['node_xyz_nm']
        if nodes.shape!=p.xyz['silicon'].shape or not np.allclose(nodes,p.xyz['silicon'],atol=1e-9,rtol=0):
            raise ValueError('Saved export state uses different Poisson nodes')
        p.set_charge(state['poisson_density_cm3'])
        restored=p.solve(float(state['vg_v']),float(state['vd_v']))
        error=float(np.max(abs(restored-state['node_phi_v'])))
        if error>1e-10:raise RuntimeError(f'Saved electrostatic export mismatch: {error:g} V')
        report={'state':file.name,'vg_v':float(state['vg_v']),'vd_v':float(state['vd_v']),
            'saved_state_source_sha256':str(state['source_sha256']),
            'configuration_fingerprint':str(state['configuration_fingerprint']),
            'maximum_saved_node_phi_error_v':error,'native_devsim_coordinates':'cm',
            'quantum_vtu_coordinates':'nm','potential_unit':'V','carrier_density_unit':'cm^-3',
            'method':'Linear Poisson restoration from accepted saved input charge; quantum data unchanged'}
    (out/'field_export.json').write_text(json.dumps(report,indent=2));return report
