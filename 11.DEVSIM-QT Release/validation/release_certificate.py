"""Fail-closed release certificate. Missing evidence never becomes a pass."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ALLOWED={'PASS','FAIL','NOT_RUN','INCONCLUSIVE','OUT_OF_SCOPE'}

def all_required_pass(items):
    return bool(items) and all(x.get('status')=='PASS' for x in items.values() if x.get('required',True))

def certify(root):
    root=Path(root);results=root/'results';items={}
    from provenance import source_digest
    current_source=source_digest()
    def put(name,status,evidence,reason=''):
        if status not in ALLOWED:raise ValueError(status)
        items[name]={'status':status,'required':True,'evidence':evidence,'reason':reason}
    def execution(name):
        p=results/name/'execution.json'
        return json.loads(p.read_text()) if p.exists() else {}
    xml=results/'pytest_final/junit.xml'
    if xml.exists():
        tree=ET.parse(xml);cases=tree.findall('.//testcase');bad=tree.findall('.//failure')+tree.findall('.//error');skips=tree.findall('.//skipped')
        native=[x for x in cases if any(k in x.attrib.get('classname','') for k in ('devsim','triple_junction','shared_poisson','restart_boundary'))]
        put('native_integration','PASS' if native and not bad and not skips else 'FAIL',str(xml.relative_to(root)),f'{len(cases)} tests; {len(bad)} failures/errors; {len(skips)} skipped')
        mock=[x for x in cases if 'mock_quantum' in x.attrib.get('name','')]
        put('backend_replacement','PASS' if len(mock)>=2 and not bad else 'FAIL',str(xml.relative_to(root)))
    else:
        put('native_integration','NOT_RUN','results/pytest_final/junit.xml');put('backend_replacement','NOT_RUN','results/pytest_final/junit.xml')
    independent=results/'independent_final/summary.json'
    if independent.exists():
        d=json.loads(independent.read_text());ok=d.get('pass',d.get('all_pass',False));put('independent_quantum_reference','PASS' if ok else 'FAIL',str(independent.relative_to(root)))
    else:put('independent_quantum_reference','NOT_RUN','results/independent_final/summary.json')
    statefiles=list(results.glob('*/states/*.json'));states=[]
    for f in statefiles:
        row=json.loads(f.read_text());npz=f.with_suffix('.npz')
        if not npz.exists():continue
        import numpy as np
        with np.load(npz,allow_pickle=False) as z:config=json.loads(str(z['config_json']))
        states.append((row,config,str(f.relative_to(root))))
    target={'geometry_profile':'raised','backend':'negf','contact_model':'transverse_sp','ls':12.,'ld':12.,'lg':24.,'spacer':3.,'width':12.,'thickness':5.,'sd_width':20.,'sd_thickness':15.,'tox':2.,'donor_sd':1e19,'donor_channel':0.,'temperature':300.,'flatband':.5}
    raised=[(r,c,f) for r,c,f in states if all(c.get(k)==v for k,v in target.items()) and r.get('numerical_status')=='SCF_CONVERGED' and r.get('source_sha256')==current_source]
    nonzero=[f for r,c,f in raised if r['vd_v']>0]
    attempted=list(results.glob('raised*/failures/*.json'))
    put('raised_nonequilibrium_scf','PASS' if nonzero else 'FAIL' if attempted else 'NOT_RUN',nonzero or [str(f.relative_to(root)) for f in attempted], 'Numerical SCF alone is not a physical or training certificate.')
    have={(round(r['vg_v'],6),round(r['vd_v'],6)) for r,c,f in raised}
    required={(round(i*.05,6),vd) for vd in (.05,.1) for i in range(25)}|{(vg,round(i*.05,6)) for vg in (.4,.8,1.2) for i in range(25)}
    missing=sorted(required-have)
    put('gate_drain_iv','PASS' if not missing else 'NOT_RUN','results/assessment/requested_sweep_status.json',f'{len(required-have)} of {len(required)} unique requested bias points missing; no interpolation')
    assessment=results/'assessment/gates.json';gates=json.loads(assessment.read_text()) if assessment.exists() else {}
    for name in ('lead_device_matching_and_cut','bound_state_completeness','mode_convergence','grid_convergence','electrostatic_buffer_convergence','energy_convergence','clean_install_and_resume','extracted_archive_execution'):
        item=gates.get(name,{'status':'NOT_RUN','evidence':None,'reason':'No completed evidence'})
        put(name,item['status'],item.get('evidence'),item.get('reason',''))
    ready=all_required_pass(items)
    cert={'product':'DEVSIM-QT','version':'1.0' if ready else '1.0-rc1','classification':'RELEASE' if ready else 'RELEASE_CANDIDATE',
        'RELEASE_READY':ready,'TEACHER_READY':ready,'source_sha256':current_source,'mandatory':items,
        'scope':'six-valley spin-degenerate silicon EMA ballistic GAA, 300 K; no scattering or learned correction',
        'status_values':sorted(ALLOWED),'numerically_converged_states':[{'path':f,'vg_v':r['vg_v'],'vd_v':r['vd_v'],'geometry':c.get('geometry_profile','raised'),'teacher_ready':False,'source_matches_current':r.get('source_sha256')==current_source} for r,c,f in states if r['numerical_status']=='SCF_CONVERGED'],
        'missing_requested_raised_biases':missing,
        'out_of_scope':{'AI_training':'OUT_OF_SCOPE','DFT_correction':'OUT_OF_SCOPE','phonon_scattering':'OUT_OF_SCOPE','drift_diffusion':'OUT_OF_SCOPE'},
        'policy':'Every mandatory item must be PASS. A toy-kernel pass never substitutes for device convergence. All exported states are non-teacher data unless separately certified within validated scope.'}
    (root/'RELEASE_CERTIFICATE.json').write_text(json.dumps(cert,indent=2,allow_nan=False));return cert

if __name__=='__main__':
    import sys
    c=certify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]);print(c['classification'],c['RELEASE_READY'])
