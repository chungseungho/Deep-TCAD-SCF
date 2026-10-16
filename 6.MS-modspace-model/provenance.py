"""Source/config identity excludes generated output and tests by definition."""
import hashlib,json
from pathlib import Path


def fingerprint(c):
    return hashlib.sha256(json.dumps(c.to_dict(),sort_keys=True).encode()).hexdigest()


def source_digest():
    root=Path(__file__).resolve().parent;h=hashlib.sha256()
    for f in sorted(root.rglob('*.py')):
        if any(p in ('output','results','tests','validation','__pycache__','.venv','venv','.git','.pytest_cache') for p in f.relative_to(root).parts): continue
        h.update(str(f.relative_to(root)).encode());h.update(f.read_bytes())
    return h.hexdigest()
