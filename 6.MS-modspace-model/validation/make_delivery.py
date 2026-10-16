#!/usr/bin/env python3
"""Build a reproducible STEP06 source/result ZIP after all solvers stop."""
import argparse
import hashlib
import json
import os
import zipfile
from pathlib import Path


EXCLUDED_PARTS={'__pycache__','.pytest_cache','.venv','.git',
                'development_point','negf_pilot','spectral_pilot'}
EXCLUDED_NAMES={'PACKAGE_MANIFEST.json','tests_verbose.txt','tests_native.log'}


def make(root,destination):
    root=Path(root).resolve();destination=Path(destination).resolve()
    if root in destination.parents:raise ValueError('Put the delivery ZIP outside the source tree')
    files=[p for p in sorted(root.rglob('*')) if p.is_file()
           and not set(p.relative_to(root).parts)&EXCLUDED_PARTS
           and p.name not in EXCLUDED_NAMES and not p.name.endswith('.tmp.npz')]
    inventory=[{'path':str(p.relative_to(root)),'bytes':p.stat().st_size,
                'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
    summary=json.loads((root/'results/step06_summary.json').read_text())
    manifest={'files':inventory,'file_count':len(inventory),
              'payload_bytes':sum(item['bytes'] for item in inventory),
              'validation_levels':summary['levels'],
              'quantitative_device_validation':False,
              'excluded':'caches, private environments, superseded development pilots and intermediate test-capture logs; final tests, actual failed reference attempts and incomplete checkpoints are included',
              'note':'This manifest hashes every delivered payload file except itself. Saved state source identities refer to the simulator core, with validation extensions separately fingerprinted.'}
    file=root/'PACKAGE_MANIFEST.json';file.write_text(json.dumps(manifest,indent=2))
    destination.parent.mkdir(parents=True,exist_ok=True);temporary=destination.with_suffix('.tmp.zip')
    with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as archive:
        for path in files+[file]:archive.write(path,str(Path(root.name)/path.relative_to(root)))
    with zipfile.ZipFile(temporary) as archive:
        bad=archive.testzip()
        if bad:raise RuntimeError('ZIP CRC verification failed: '+bad)
        entries=len(archive.infolist())
    os.replace(temporary,destination)
    digest=hashlib.sha256(destination.read_bytes()).hexdigest()
    destination.with_suffix('.sha256').write_text(digest+'  '+destination.name+'\n')
    return {'zip':str(destination),'bytes':destination.stat().st_size,'sha256':digest,
            'entries':entries,'crc_check':'pass'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',required=True);args=parser.parse_args()
    print(json.dumps(make(Path(__file__).resolve().parents[1],args.out),indent=2))
