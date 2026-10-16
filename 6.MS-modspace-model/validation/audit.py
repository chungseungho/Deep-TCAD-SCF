#!/usr/bin/env python3
"""Audit the preserved STEP05 raw files, not numbers transcribed from prose."""
from pathlib import Path
import csv,json,hashlib


def run(old,destination):
    old=Path(old);out=Path(destination);out.parent.mkdir(parents=True,exist_ok=True)
    table=list(csv.DictReader((old/'results/convergence.csv').open()))
    baseline=next(float(r['id_a']) for r in table if r['case']=='baseline')
    variations={r['case']:float(r['id_a'])/baseline-1 for r in table}
    output=list(csv.DictReader((old/'results/idvds.csv').open()));transfer=list(csv.DictReader((old/'results/idvgs.csv').open()))
    failed=[[float(r['vg_v']),float(r['vd_v'])] for r in output if r['status']!='SCF_CONVERGED_PROTOTYPE']
    report={'source_step05':str(old.name),'baseline_id_a':baseline,'relative_current_sensitivity':variations,
            'idvgs_scf':[sum(r['status']=='SCF_CONVERGED_PROTOTYPE' for r in transfer),len(transfer)],
            'idvds_scf':[len(output)-len(failed),len(output)],'failed_idvds_biases':failed,
            'classified_code_issues':['non-adjoint potential/density map; work defect','coarse FD band exceeds applied-bias validity','saved final fixed-point residual needed re-evaluation'],
            'classified_model_limits':['wide-band finite contact width and length are model parameters','truncated transverse basis convergence absent','isolated neutral calibration does not imply connected neutral contacts','closed stationary spatial occupation does not define a DC steady state'],
            'step06_remedies':['exact tetra/cell intersection exchange integration','energy-window dispersion diagnostics per valley','independent semi-infinite reference with explicit GF cost','adaptive bias integration with failure on unresolved error','unchanged physical tolerances and separately labelled finite-bath outputs'],
            'remaining_validation':'No fitting, gamma plateau or AI-teacher readiness inferred from numerical identities'}
    out.write_text(json.dumps(report,indent=2));return report


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--old',default='../5.MS-modspace-GAA');p.add_argument('--out',default='results/step05_audit.json');a=p.parse_args()
    print(json.dumps(run(a.old,a.out),indent=2))
