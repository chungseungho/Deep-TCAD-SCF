#!/usr/bin/env python3
"""Re-export accepted saved states; never recalculate or replace quantum data."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from config import Config
from export_fields import restore_last_poisson,export_last
from poisson.device import Poisson


EXCLUDE={'development_point','negf_pilot','spectral_pilot'}


def reexport(root):
    root=Path(root);reports=[]
    for file in sorted(root.rglob('settings.json')):
        directory=file.parent
        if set(directory.relative_to(root).parts)&EXCLUDE:continue
        if not list((directory/'states').glob('vg*.npz')):continue
        mesh=directory/'mesh/gaa_cm.msh'
        if not mesh.exists():continue
        settings=json.loads(file.read_text());p=Poisson(Config(**settings['config']),mesh)
        try:
            restored=restore_last_poisson(p,directory)
            if restored:
                p.export(directory/'electrostatics.vtk');export_last(directory)
                reports.append({'run':str(directory.relative_to(root)),**restored})
        finally:p.close()
    result={'runs':reports,'count':len(reports),
            'maximum_saved_phi_reproduction_error_v':max((r['maximum_saved_node_phi_error_v'] for r in reports),default=None),
            'scope':'Only the linear Poisson field is restored from the exact saved input charge; NPZ/CSV quantum results are unchanged. Native cm, quantum VTU nm.'}
    (root/'field_export_audit.json').write_text(json.dumps(result,indent=2));return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--root',default='results');args=parser.parse_args()
    print(json.dumps(reexport(args.root),indent=2))
