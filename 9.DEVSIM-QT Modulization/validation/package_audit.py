"""Verify saved-state identity, invariant checks, and prior-stage protection."""
import sys,json,hashlib,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from provenance import source_digest


def run():
    results=ROOT/'results';source_index=[];available=set()
    for archive in sorted((results/'source_snapshots').glob('*.zip')):
        h=hashlib.sha256()
        with zipfile.ZipFile(archive) as z:
            for n in sorted(z.namelist()):h.update(n.encode());h.update(z.read(n))
        exact=h.hexdigest()==archive.stem
        if not exact:raise AssertionError('Source snapshot mismatch '+archive.name)
        available.add(archive.stem)
        source_index.append({'source_sha256':archive.stem,'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
                             'path':str(archive.relative_to(ROOT)),'payload_identity_verified':exact})
    (results/'source_snapshots/index.json').write_text(json.dumps(source_index,indent=2))
    states=[];maximum_parts_error=0.
    for f in sorted(results.glob('*/states/*.npz')):
        with np.load(f) as s:
            config=json.loads(str(s['config_json']))
            fp=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
            assert fp==str(s['configuration_fingerprint']),f
            assert str(s['source_sha256']) in available,f
            assert not bool(s['quantitative_validation_pass']),f
            density=s['electron_density_cm3'];assert np.isfinite(density).all() and min(density)>=0,f
            accepted=bool(s['converged'])
            if accepted:
                assert float(s['phi_residual_v'])<config['phi_tol'],f
                assert float(s['density_relative_residual'])<config['charge_tol'],f
                assert float(s['gauss_relative_error'])<1e-7,f
            if str(s['backend'])=='ms':
                parts=s['density_from_l_cm3']+s['density_from_r_cm3']+s['density_unresolved_origin_cm3']
                defect=float(np.linalg.norm(parts-density)/max(np.linalg.norm(density),1e-300))
                maximum_parts_error=max(maximum_parts_error,defect);assert defect<1e-10,f
                for j in range(3):
                    assert max(abs(s[f'valley{j}_weights'].sum(axis=1)-1))<1e-9,f
                    assert not s[f'valley{j}_ambiguous'].any(),f
                    assert min(s[f'valley{j}_occupation'])>=0 and max(s[f'valley{j}_occupation'])<=1,f
            states.append({'path':str(f.relative_to(ROOT)),'source_sha256':str(s['source_sha256']),
                'configuration_fingerprint':fp,'backend':str(s['backend']),'numerically_converged':accepted,
                'quantitative_validation_pass':False})
    (results/'dataset_catalog.json').write_text(json.dumps(states,indent=2))
    prior=json.loads((ROOT/'research/prior_stages_readonly_manifest.json').read_text());changed=[];unavailable=[]
    for name,expected in prior.items():
        file=ROOT.parent/name
        if not file.exists():unavailable.append(name)
        elif hashlib.sha256(file.read_bytes()).hexdigest()!=expected:changed.append(name)
    assert not changed,changed
    tests=json.loads((results/'test_report/summary.json').read_text());assert tests['pass'] and tests['total_tests']==56
    report={'pass':True,'source_sha256':source_digest(),'prior_stage_files_unchanged':len(prior)-len(unavailable),
            'prior_stage_files_unavailable':len(unavailable),
            'state_files':len(states),'numerically_converged_states':sum(x['numerically_converged'] for x in states),
            'verified_core_snapshots':len(source_index),'maximum_ms_origin_density_additivity_error':maximum_parts_error,
            'note':'packaging/numerical integrity audit; not a Level4 physical validation'}
    (results/'package_audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));return report


if __name__=='__main__':run()
