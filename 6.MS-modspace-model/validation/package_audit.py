#!/usr/bin/env python3
"""Read-only cross-check of saved physical states and delivery provenance."""
import csv
import hashlib
import json
import sys
import zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from provenance import source_digest


EXCLUDE={'development_point','negf_pilot','spectral_pilot'}


def audit(root):
    root=Path(root).resolve();result=root/'results';issues=[];runs=[];count=0;accepted=0
    digest=source_digest();origin_errors=[];known_sources={digest};archives=[]
    index=result/'source_snapshots/index.json'
    if index.exists():
        for expected,record in json.loads(index.read_text()).items():
            archive=index.parent/record['archive'];check=hashlib.sha256()
            if hashlib.sha256(archive.read_bytes()).hexdigest()!=record['archive_sha256']:
                issues.append('Source archive file checksum differs: '+record['archive']);continue
            with zipfile.ZipFile(archive) as bundle:
                for name in sorted(bundle.namelist()):
                    check.update(name.encode());check.update(bundle.read(name))
            valid=check.hexdigest()==expected
            if valid:known_sources.add(expected)
            else:issues.append('Source archive core identity differs: '+record['archive'])
            archives.append({'source_sha256':expected,'archive':record['archive'],'verified':valid})
    for path in sorted(result.rglob('settings.json')):
        directory=path.parent
        if any(part in EXCLUDE for part in directory.relative_to(result).parts):continue
        settings=json.loads(path.read_text());conf=settings['config'];cases=[]
        if settings['source_sha256'] not in known_sources:issues.append(str(path)+': source identity is neither current nor exactly archived')
        for phase in ['point','idvgs','idvds']:
            table=directory/(phase+'.csv')
            if table.exists():cases+=list(csv.DictReader(table.open()))
        for row in cases:
            good=row['status']=='SCF_CONVERGED_PROTOTYPE'
            state=directory/'states'/(row['key']+'.npz')
            if not state.exists():
                if good:issues.append(str(state)+': accepted state missing')
                continue
            with np.load(state) as data:
                count+=1;accepted+=int(good)
                def check(condition,message):
                    if not condition:issues.append(str(state)+': '+message)
                check(str(data['source_sha256'])==settings['source_sha256'],'state/run source identity differs')
                identity=settings.get('extension_configuration_fingerprint',settings['fingerprint'])
                check(str(data['configuration_fingerprint'])==identity,'configuration identity differs')
                check(json.loads(str(data['config_json']))==conf,'configuration JSON differs')
                check(bool(data['converged'])==good,'CSV and NPZ status differ')
                check(not bool(data['quantitative_validation_pass']),'unjustified reference label')
                for field in ['phi_v','electron_density_cm3','node_phi_v','poisson_density_cm3','bond_current_a']:
                    check(np.isfinite(data[field]).all(),'nonfinite '+field)
                check(data['electron_density_cm3'].min()>=-1e-10,'negative density')
                check(data['poisson_density_cm3'].min()>=-1e-10,'negative Poisson input density')
                for name in ['vg_v','vd_v']:
                    check(float(data[name])==float(row[name]),'bias differs: '+name)
                if good:
                    check(float(data['phi_residual_v'])<conf['phi_tol'],'unaccepted potential residual')
                    check(float(data['charge_residual'])<conf['charge_tol'],'unaccepted charge residual')
                    meta=json.loads(state.with_suffix('.json').read_text())
                    check(meta['gauss_relative_error']<1e-7,'Gauss balance fails')
                    check(meta['mapping_work_relative_defect']<1e-8,'work adjoint fails')
                    check(np.isclose(meta['id_a'],float(row['id_a']),rtol=1e-13,atol=1e-25),'current CSV differs')
                    if 'density_from_l_cm3' in data:
                        total=data['density_from_l_cm3']+data['density_from_r_cm3']
                        density=data['electron_density_cm3'];difference=total-density
                        # Matrix contractions and the individual nonnegative
                        # roundoff floors need a global density scale. A local
                        # relative test at a depleted cell overstates this error.
                        peak=float(abs(difference).max()/max(density.max(),1.))
                        norm=float(np.linalg.norm(difference)/max(np.linalg.norm(density),1.))
                        origin_errors.append({'state':str(state.relative_to(result)),
                                              'peak_normalized_error':peak,'relative_l2_error':norm})
                        check(peak<2e-11 and norm<2e-11,'origin density additivity exceeds contraction roundoff budget')
        runs.append({'directory':str(directory.relative_to(result)),'completed_rows':len(cases)})
    old=result/'step01_05_sha256.json';old_count=0;old_bad=[];old_status='not available'
    if old.exists():
        original=json.loads(old.read_text())
        old_roots={name.split('/')[0] for name in original}
        if all((root.parent/name).is_dir() for name in old_roots):
            old_count=len(original);old_status='checked against the preserved pre-work inventory'
            for name,expected in original.items():
                file=root.parent/name
                if not file.exists() or hashlib.sha256(file.read_bytes()).hexdigest()!=expected:old_bad.append(name)
            issues += ['STEP01-05 changed: '+name for name in old_bad]
        else:old_status='skipped: standalone STEP06 extraction has no full STEP01-05 tree'
    return {'pass':not issues,'issues':issues,'source_sha256':digest,
            'verified_source_archives':archives,
            'saved_states_checked':count,'accepted_states_checked':accepted,
            'maximum_origin_additivity_peak_error':max((r['peak_normalized_error'] for r in origin_errors),default=None),
            'maximum_origin_additivity_l2_error':max((r['relative_l2_error'] for r in origin_errors),default=None),
            'old_files_checked':old_count,'old_files_changed':old_bad,'old_files_check_status':old_status,'runs':runs,
            'scope':'File/state/provenance integrity and recorded numerical acceptance; not contact or continuum validation'}


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];report=audit(root)
    (root/'results/package_audit.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2));raise SystemExit(0 if report['pass'] else 1)
