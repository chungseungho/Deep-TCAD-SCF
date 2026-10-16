#!/usr/bin/env python3
"""Merge disjoint, identically configured SCF sweep partitions without fitting.

Only completed bias rows are merged. Accepted rows must have their actual
saved state. This is file organization after all partition solvers stop.
"""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


def read_rows(path):
    return list(csv.DictReader(path.open())) if path.exists() else []


def write_rows(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def merge(target, parts, require_complete=False):
    target = Path(target).resolve()
    sources = [target] + [Path(part).resolve() for part in parts]
    if len(set(sources)) != len(sources):
        raise ValueError('Partition paths must be distinct')
    reference = json.loads((target/'settings.json').read_text())
    records, allrows, histories, seen = [], [], [], set()
    for directory in sources:
        settings = json.loads((directory/'settings.json').read_text())
        for key in ['config', 'fingerprint', 'source_sha256']:
            if settings[key] != reference[key]:
                raise ValueError(f'Incompatible {key}: {directory}')
        table = directory/'idvds.csv'
        rows = read_rows(table)
        for row in rows:
            if row['key'] in seen:
                raise ValueError('Overlapping bias partitions: '+row['key'])
            seen.add(row['key'])
            row['execution_partition']=directory.name
            if row['status'] == 'SCF_CONVERGED_PROTOTYPE':
                for suffix in ['.json', '.npz']:
                    if not (directory/'states'/(row['key']+suffix)).is_file():
                        raise ValueError('Missing accepted state: '+row['key']+suffix)
            allrows.append(row)
        histories += [{**row,'execution_partition':directory.name} for row in read_rows(directory/'scf_history.csv')]
        records.append({'directory':str(directory.relative_to(target.parent)),
                        'completed_biases':len(rows),
                        'bias_table_sha256':hashlib.sha256(table.read_bytes()).hexdigest(),
                        'argv':settings.get('argv', []),
                        'source_sha256':settings['source_sha256'],
                        'configuration_fingerprint':settings['fingerprint']})
    expected = {f'vg{g:.6f}_vd{i*.05:.6f}' for g in [.4, .8, 1.2] for i in range(25)}
    if seen - expected:
        raise ValueError('Unexpected bias keys: '+str(seen-expected))
    if require_complete and seen != expected:
        raise ValueError('Incomplete sweep: '+str(expected-seen))
    for directory in sources[1:]:
        completed_keys={row['key'] for row in read_rows(directory/'idvds.csv')}
        for subdir in ['states', 'checkpoints']:
            for source in sorted((directory/subdir).glob('*')):
                if not source.is_file():
                    continue
                if source.stem.removesuffix('_spectra') not in completed_keys:
                    continue
                destination = target/subdir/source.name
                if destination.exists():
                    raise ValueError('Existing state would be overwritten: '+str(destination))
                destination.parent.mkdir(exist_ok=True)
                shutil.copy2(source, destination)
    allrows.sort(key=lambda row:(float(row['vg_v']),float(row['vd_v'])))
    write_rows(target/'idvds.csv', allrows)
    write_rows(target/'scf_history.csv', histories)
    report = {'summary_kind':'disjoint_executed_sweep_partitions',
              'requested_points':75, 'finished_points':len(allrows),
              'converged_points':sum(row['status']=='SCF_CONVERGED_PROTOTYPE' for row in allrows),
              'failed_points':sum(row['status']!='SCF_CONVERGED_PROTOTYPE' for row in allrows),
              'missing_bias_keys':sorted(expected-seen),
              'quantitative_device_validation':False,
              'partitions':records,
              'note':'Same Hamiltonian, material, grid, contact and SCF tolerances; different numerical continuation seeds. No interpolation or current replacement.'}
    (target/'partition_execution.json').write_text(json.dumps(report,indent=2))
    (target/'run_summary.json').write_text(json.dumps(report,indent=2))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target',required=True)
    parser.add_argument('--parts',nargs='+',required=True)
    parser.add_argument('--require-complete',action='store_true')
    args = parser.parse_args()
    print(json.dumps(merge(args.target,args.parts,args.require_complete),indent=2))
