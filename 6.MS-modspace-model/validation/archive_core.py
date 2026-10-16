#!/usr/bin/env python3
"""Keep the exact simulator core behind already computed source identities."""
import hashlib
import json
import sys
import zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from provenance import source_digest


def archive(root):
    root=Path(root).resolve();digest=source_digest()
    excluded={'output','results','tests','validation','__pycache__','.venv','venv','.git','.pytest_cache'}
    files=[p for p in sorted(root.rglob('*.py')) if not set(p.relative_to(root).parts)&excluded]
    check=hashlib.sha256()
    for file in files:check.update(str(file.relative_to(root)).encode());check.update(file.read_bytes())
    if check.hexdigest()!=digest:raise RuntimeError('Archive selection differs from active core digest')
    folder=root/'results/source_snapshots';folder.mkdir(exist_ok=True)
    destination=folder/(digest+'.zip')
    with zipfile.ZipFile(destination,'w',zipfile.ZIP_DEFLATED) as bundle:
        for file in files:bundle.write(file,str(file.relative_to(root)))
    index=folder/'index.json';record=json.loads(index.read_text()) if index.exists() else {}
    record[digest]={'archive':destination.name,'archive_sha256':hashlib.sha256(destination.read_bytes()).hexdigest(),
                   'python_files':len(files),
                   'purpose':'Exact source of computed states before postprocessing/provenance fixes. Extract into a separate STEP06 copy to reproduce/resume those original source identities.'}
    index.write_text(json.dumps(record,indent=2));return record[digest]


if __name__=='__main__':print(json.dumps(archive(Path(__file__).resolve().parents[1]),indent=2))
