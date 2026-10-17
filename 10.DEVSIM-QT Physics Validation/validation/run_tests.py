#!/usr/bin/env python3
"""Run scientific and native tests in isolated processes; require completion XML.

A native runtime can terminate a mixed-process test run with exit status zero.
Exit status alone is therefore never sufficient evidence of completed tests.
"""
import sys,json,subprocess,xml.etree.ElementTree as ET
from pathlib import Path


def run(out='results/test_report'):
    root=Path(__file__).resolve().parents[1];out=(root/out).resolve();out.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(root))
    native=['tests/test_devsim_integration.py','tests/test_shared_poisson_capacitor.py','tests/test_triple_junction.py']
    groups=[('quantum',['tests']+['--ignore='+p for p in native]),('native',native)]
    results=[]
    for name,args in groups:
        xml=out/f'{name}.xml';log=out/f'{name}.log'
        if xml.exists():xml.unlink()
        command=[sys.executable,'-m','pytest','-q',*args,'--junitxml='+str(xml)]
        with log.open('w') as f:process=subprocess.run(command,cwd=root,stdout=f,stderr=subprocess.STDOUT)
        if process.returncode!=0 or not xml.exists():raise RuntimeError(f'{name} test group did not complete; inspect {log}')
        suites=ET.parse(xml).getroot().findall('.//testsuite')
        tests=sum(int(s.attrib.get('tests',0)) for s in suites)
        failures=sum(int(s.attrib.get('failures',0))+int(s.attrib.get('errors',0)) for s in suites)
        skipped=sum(int(s.attrib.get('skipped',0)) for s in suites)
        if not tests or failures:raise RuntimeError(f'{name} failed or empty test report')
        results.append({'group':name,'tests':tests,'skipped':skipped,'failures':failures,'exit_code':process.returncode,'command':command})
    from provenance import source_digest
    report={'source_sha256':source_digest(),'groups':results,'total_tests':sum(x['tests'] for x in results),
        'pass':True,'process_isolation':'native and quantum groups isolated; each requires JUnit completion marker',
        'monolithic_run':'not used as pass evidence; earlier mixed run terminated before its completion report'}
    (out/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));return report


if __name__=='__main__':run(sys.argv[1] if len(sys.argv)>1 else 'results/test_report')
