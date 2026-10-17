#!/usr/bin/env python3
"""Audit saved evidence and package identity; this is not physics certification.

Run without options to print an audit without changing the delivered files.
--baseline-root optionally checks the original restored STEP07 payload.
--write-report is intended for the pre-package audit, before manifest creation.
"""
import argparse
import csv
import hashlib
import json
import math
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from provenance import source_digest


def read(path):
    return json.loads((ROOT / path).read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def audit(baseline_root=None):
    checks = []
    current = source_digest()
    snapshots = {}

    def check(name, function):
        try:
            details = function()
            checks.append({'name': name, 'pass': True, 'details': details})
        except Exception as error:
            checks.append({'name': name, 'pass': False,
                           'error': f'{type(error).__name__}: {error}'})

    def sources():
        files = [ROOT / 'baseline/step07_production_c0fbe6dc.zip',
                 *sorted((ROOT / 'results/source_snapshots').glob('*.zip'))]
        for path in files:
            digest = hashlib.sha256()
            with zipfile.ZipFile(path) as archive:
                require(archive.testzip() is None, str(path))
                names = sorted(archive.namelist())
                require(all(n.endswith('.py') for n in names), 'Unexpected snapshot file')
                for name in names:
                    digest.update(name.encode())
                    digest.update(archive.read(name))
                value = digest.hexdigest()
                snapshots[value] = str(path.relative_to(ROOT))
                if path.parent.name == 'source_snapshots':
                    require(path.stem == value, 'Snapshot filename/digest mismatch')
                if value == current:
                    for name in names:
                        require((ROOT / name).read_bytes() == archive.read(name), name)
        require(current in snapshots, 'Current production snapshot missing')
        tests = read('results/tests_batched/summary.json')
        require(tests['source_sha256'] == current, 'Tests refer to different core')
        require(read('STATUS.json')['source_sha256'] == current, 'Status/core mismatch')
        return {'current_core_sha256': current, 'snapshots': snapshots,
                'digest_scope': 'production Python; validation/tests covered by package manifest'}

    def junit():
        summary = read('results/tests_batched/summary.json')
        total = 0
        for group in summary['groups']:
            path = ROOT / 'results/tests_batched' / (group['group'] + '.xml')
            suites = ET.parse(path).getroot().findall('.//testsuite')
            counts = {key: sum(int(s.get(key, '0')) for s in suites)
                      for key in ('tests', 'failures', 'errors', 'skipped')}
            require(counts['tests'] == group['tests'], 'JUnit test count mismatch')
            require(not any(counts[k] for k in ('failures', 'errors', 'skipped')),
                    'Incomplete or failed tests')
            require(group['exit_code'] == 0 and counts['tests'] > 0, 'Test completion')
            total += counts['tests']
        require(total == 64 == summary['total_tests'] and summary['pass'], 'Total tests')
        return {'completed_passed': total, 'skipped': 0}

    def states():
        records = []
        for path in sorted((ROOT / 'results').glob('*/states/*.json')):
            state = json.loads(path.read_text())
            setting = json.loads((path.parents[1] / 'settings.json').read_text())
            config = setting['config']
            fingerprint = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
            require(fingerprint == setting['configuration_fingerprint']
                    == state['configuration_fingerprint'], f'Config identity: {path}')
            require(state['source_sha256'] in snapshots, f'Unknown core: {path}')
            require(state['source_sha256'] == setting['source_sha256'], f'Core identity: {path}')
            require(state['quantitative_validation_pass'] is False, 'Unexpected certification')
            with np.load(path.with_suffix('.npz'), allow_pickle=False) as data:
                for name in ('phi_v', 'electron_density_cm3', 'poisson_density_cm3'):
                    require(np.isfinite(data[name]).all(), f'Nonfinite {name}: {path}')
                require(np.min(data['electron_density_cm3']) >= -1e-7, f'Negative density: {path}')
                require(str(data['source_sha256']) == state['source_sha256'], 'NPZ core mismatch')
                require(str(data['configuration_fingerprint']) == fingerprint, 'NPZ config mismatch')
                if state['numerical_status'] == 'SCF_CONVERGED':
                    require(bool(data['converged']), 'JSON/NPZ convergence mismatch')
                    require(state['phi_residual_v'] < config['phi_tol'], 'Potential residual')
                    require(state['density_relative_residual'] < config['charge_tol'], 'Density residual')
                    require(state['gauss_relative_error'] < 1e-7, 'Gauss residual')
            records.append({'path': str(path.relative_to(ROOT)),
                            'numerical_status': state['numerical_status'],
                            'physical_certification': False})
        require(len(records) == 9, 'Unexpected recorded state count')
        return {'state_count': len(records), 'states': records,
                'scope': 'saved numeric fields and identities, not a fresh SCF solve'}

    def failures():
        count = 0
        for folder in ('baseline_ms', 'baseline_negf', 'repaired_ms',
                       'repaired_negf', 'repaired_negf_vg04', 'final_ms',
                       'final_negf_equilibrium'):
            for path in (ROOT / 'results' / folder).glob('idvds.csv'):
                for row in csv.DictReader(path.open()):
                    if row['numerical_status'] != 'SCF_CONVERGED':
                        require(math.isnan(float(row['id_a'])), f'Failed bias has current: {path}')
                        count += 1
        require(count > 0, 'No rejected-bias evidence')
        for folder in ('repaired_negf', 'repaired_negf_vg04'):
            mark = read(f'results/{folder}/INTERRUPTED.json')
            require(mark['status'] == 'INTERRUPTED_NOT_A_COMPLETED_SWEEP', folder)
        oom = read('results/gaa_modes_high/RESOURCE_FAILURE.json')
        require(oom['exit_code'] == -9, 'OOM record')
        incomplete = read('results/gaa_modes_high/summary.json')
        require(incomplete['complete'] is False and len(incomplete['cases']) == 2,
                'OOM run marked completed')
        return {'rejected_bias_rows_with_nan_current': count,
                'interrupted_runs': 2, 'oom_runs': 1}

    def experiments():
        expected = {'uniform': 80, 'constriction': 45, 'gaa_modes': 14,
                    'gaa_modes_batched': 6, 'matched_mu_modes': 14,
                    'origin_sensitivity': 21}
        for folder, count in expected.items():
            data = read(f'results/{folder}/summary.json')
            require(len(data['cases']) == count, f'Case count: {folder}')
            require(data.get('complete', True), f'Incomplete: {folder}')
            require(data['source_sha256'] in snapshots, f'Unknown source: {folder}')
            execution = read(f'results/{folder}/execution.execution.json')
            require(execution['exit_code'] == 0, f'Execution completion: {folder}')
        require(read('results/constriction/summary.json')['full_basis_identity_pass'], 'Full basis')
        chain = read('results/open_chain/summary.json')
        require(chain['pass'] and len(chain['open_chain']) == 15, 'Open chain reference')
        require(len(chain['finite_bath_contact_sweep']) == 16, 'Finite bath cases')
        require(read('results/open_chain/execution.execution.json')['exit_code'] == 0, 'Chain completion')
        origin = read('results/origin_sensitivity/summary.json')
        require(origin['unknown_occupations_applied_to_density'] is False, 'Artificial filling')
        require(read('results/state_audit/summary.json')['active_ambiguous_state_records'] == 186,
                'State audit count')
        return {**expected, 'open_chain': 15, 'finite_bath': 16,
                'scope': 'completion evidence; no GAA plateau assertion'}

    def certification():
        status = read('STATUS.json')
        for flag in ('quantitative_validation_pass', 'ready_for_training',
                     'full_task_physics_goal_complete', 'siesta_latest_attachment_read'):
            require(status[flag] is False, flag)
        require(status['biased_gaa_certified_points'] == 0, 'Biased GAA status')
        plan = read('hpc_plan/plan.json')
        require(not plan['executed'] and plan['jobs_submitted'] == 0, 'HPC status')
        require(len(plan['cases']) == 17 and all(not c['executed'] for c in plan['cases']),
                'HPC plan entries')
        return {'quantitative_validation_pass': False, 'ready_for_training': False,
                'hpc_configurations_generated': 17, 'hpc_jobs_submitted': 0}

    def original():
        inventory = read('baseline/step07_file_inventory.json')['files']
        require(len(inventory) == 917, 'STEP07 inventory size')
        if baseline_root is None:
            return {'inventory_files': 917, 'original_directory_checked': False,
                    'reason': 'optional external STEP07 directory not supplied'}
        for name, expected in inventory:
            path = Path(baseline_root) / name
            require(path.is_file() and sha(path) == expected, f'Original STEP07 changed: {name}')
        return {'inventory_files': 917, 'original_directory_checked': True,
                'all_original_payloads_unchanged': True}

    def payload():
        docs = ['README.md', 'ROOT_CAUSE_REPORT.md', 'MSDFT_SOURCE_COMPARISON.md',
                'MODE_CONVERGENCE.md', 'CONTACT_CONVERGENCE.md', 'MS_NEGF_COMPARISON.md',
                'VALIDATION.md', 'HOW_TO_RUN.md', 'HPC_PLAN.md', 'job.sh']
        require(all((ROOT / p).is_file() for p in docs), 'Required deliverable absent')
        excluded = {'.venv', 'venv', '.git', '__pycache__', '.pytest_cache'}
        files = [p for p in ROOT.rglob('*') if p.is_file()
                 and not excluded.intersection(p.relative_to(ROOT).parts)]
        require(not any(p.suffix.lower() in ('.f', '.f90', '.f95', '.gz', '.tar')
                        for p in files), 'Unexpected private-source/archive payload')
        figures = list((ROOT / 'results/figures').glob('*.png'))
        require(len(figures) == 7, 'Figure count')
        for path in figures:
            require(path.read_bytes()[:8] == b'\x89PNG\r\n\x1a\n', 'Invalid PNG signature')
        manifest = ROOT / 'PACKAGE_MANIFEST.json'
        if manifest.exists():
            entries = json.loads(manifest.read_text())['files']
            require({x['path'] for x in entries} == {
                str(p.relative_to(ROOT)) for p in files if p != manifest}, 'Manifest file list')
            for entry in entries:
                path = ROOT / entry['path']
                require(path.stat().st_size == entry['size_bytes'] and sha(path) == entry['sha256'],
                        f'Manifest hash: {entry["path"]}')
            return {'figures': 7, 'private_source_files_in_payload': 0,
                    'manifest_checked': True, 'manifest_files': len(entries)}
        return {'figures': 7, 'private_source_files_in_payload': 0,
                'manifest_checked': False, 'reason': 'pre-package audit; rerun after unpacking'}

    for name, function in [('source_snapshots', sources), ('completed_tests', junit),
                           ('saved_states', states), ('honest_failure_status', failures),
                           ('experiment_completion', experiments),
                           ('certification_and_hpc_status', certification),
                           ('original_step07_integrity', original), ('package_payload', payload)]:
        check(name, function)
    return {'scope': 'STEP08 artifact integrity and recorded evidence, not physical certification',
            'source_sha256': current, 'pass': all(c['pass'] for c in checks), 'checks': checks}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-root', type=Path)
    parser.add_argument('--write-report', type=Path)
    args = parser.parse_args()
    result = audit(args.baseline_root)
    text = json.dumps(result, indent=2, ensure_ascii=False) + '\n'
    if args.write_report:
        args.write_report.parent.mkdir(parents=True, exist_ok=True)
        args.write_report.write_text(text)
    print(text)
    sys.exit(0 if result['pass'] else 1)
