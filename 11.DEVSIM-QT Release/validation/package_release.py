"""Create a complete, hash-manifested source/results archive without caches."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {"__pycache__", ".pytest_cache", ".git", ".venv", "node_modules"}


def package(destination):
    destination = Path(destination).resolve()
    if destination.is_relative_to(ROOT):
        raise ValueError("Archive destination must be outside the source directory")
    files = [p for p in sorted(ROOT.rglob("*")) if p.is_file()
             and not EXCLUDED.intersection(p.relative_to(ROOT).parts)
             and p.name not in {"PACKAGE_MANIFEST.json", "PACKAGING.json"}]
    if any(p.is_symlink() for p in files):
        raise ValueError("Do not bundle unreviewed symlink targets")
    prohibited = [str(p.relative_to(ROOT)) for p in files
                  if p.suffix.lower() in {".f", ".f90", ".f95", ".pem", ".key"}
                  or p.name.lower().endswith((".tar", ".tar.gz", ".tgz"))]
    if prohibited:
        raise ValueError(f"Review excluded private/source artifacts: {prohibited}")
    metadata = {"directory": ROOT.name, "classification": "RELEASE_CANDIDATE",
                "excluded_directory_names": sorted(EXCLUDED),
                "native_libraries_bundled": False,
                "siesta_originals_bundled": False,
                "source_results_logs_checkpoints_included": True}
    (ROOT / "PACKAGING.json").write_text(json.dumps(metadata, indent=2) + "\n")
    files.append(ROOT / "PACKAGING.json")
    manifest = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(files)}
    (ROOT / "PACKAGE_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    files.append(ROOT / "PACKAGE_MANIFEST.json")
    with zipfile.ZipFile(destination, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(files):
            z.write(p, Path(ROOT.name) / p.relative_to(ROOT))
    with zipfile.ZipFile(destination) as z:
        if z.testzip() is not None:
            raise ValueError("Archive CRC verification failed")
        for name, sha in manifest.items():
            if hashlib.sha256(z.read(ROOT.name + "/" + name)).hexdigest() != sha:
                raise ValueError("Archive manifest mismatch: " + name)
    return {"archive": str(destination), "files": len(files),
            "bytes": destination.stat().st_size,
            "sha256": hashlib.sha256(destination.read_bytes()).hexdigest()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("destination")
    args = parser.parse_args()
    print(json.dumps(package(args.destination), indent=2))
