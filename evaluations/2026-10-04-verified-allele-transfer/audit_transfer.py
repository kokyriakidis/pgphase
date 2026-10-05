#!/usr/bin/env python3
"""Audit exact per-read allele preservation; no truth or competitor inputs."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def load_matrix(path):
    variants, metadata, calls, primary, gauges = {}, {}, {}, {}, {}
    with path.open() as stream:
        for line in stream:
            fields = line.rstrip('\n').split('\t')
            if fields[0] == 'VAR':
                variants[int(fields[1])] = {
                    'key': (int(fields[2]), fields[3], int(fields[6]), fields[7]),
                    'flags': int(fields[4]),
                }
            elif fields[0] == '#META':
                metadata[int(fields[1])] = dict(item.split('=') for item in fields[2:])
            elif fields[0] == 'READ':
                gauges[fields[1]] = tuple(fields[6:8])
            elif fields[0] == 'OBS':
                key = (fields[1], variants[int(fields[2])]['key'])
                primary[key], calls[key] = int(fields[3]), int(fields[5])
    for index, variant in variants.items():
        variant['msa'] = metadata.get(index, {}).get('msa_verified') == '1'
    return {value['key']: value for value in variants.values()}, calls, primary, gauges


def audit_overlay(before, after):
    variants, old_calls, old_primary, old_gauges = before
    _, new_calls, new_primary, new_gauges = after
    counts, examples = Counter(), []
    for (name, key), allele in old_calls.items():
        variant = variants[key]
        if allele < 0 or not (variant['msa'] or variant['flags'] in (4, 8)):
            continue
        current = new_calls.get((name, key), -1)
        if current == allele:
            counts['retained'] += 1
            continue
        reason = 'lost' if current < 0 else 'changed'
        counts[reason] += 1
        counts[('msa_' if variant['msa'] else 'clean_') + reason] += 1
        if len(examples) < 25:
            examples.append({
                'read': name, 'key': key, 'before': allele, 'after': current,
                'primary_before': old_primary.get((name, key), -1),
                'primary_after': new_primary.get((name, key), -1), 'msa': variant['msa'],
            })
    return {
        'counts': dict(counts), 'examples': examples,
        'changed_read_gauges': sum(old_gauges.get(name) != new_gauges.get(name)
                                   for name in old_gauges.keys() | new_gauges.keys()),
    }


def audit_source(folder):
    expected, actual, census = defaultdict(set), {}, Counter()
    traces = sorted(folder.glob('matrix.chunk0.transfer*.tsv'))
    for path in traces:
        with path.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                census[row['status']] += 1
                if row['status'] != 'mapped':
                    continue
                # The C++ trace resolves canonical sequence keys to destination
                # indices, including graph-walk aliases. Comparing raw graph
                # keys with BAM sequence keys would manufacture missing calls.
                key = (path.name, row['qname'], int(row['destination_index']))
                expected[key].add(int(row['source_allele']))
                actual[key] = int(row['bam_allele'])
    counts, examples = Counter(), []
    for key, alleles in expected.items():
        if len(alleles) > 1:
            counts['overlapping_source_conflict'] += 1
            continue
        allele = next(iter(alleles))
        current = actual[key]
        if current == allele:
            counts['retained'] += 1
            continue
        reason = 'lost' if current < 0 else 'changed'
        counts[reason] += 1
        if len(examples) < 25:
            examples.append({
                'stage': key[0], 'read': key[1], 'destination_index': key[2],
                'source': allele, 'destination': current,
            })
    return {'available': bool(traces), 'source_census': dict(census),
            'counts': dict(counts), 'examples': examples}


def audit_admission(folder):
    counts = Counter()
    for path in sorted(folder.glob('*.msa-admission.tsv')):
        with path.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                kind = 'fixed_msa' if row['update_counts'] == '1' else 'physical'
                counts[kind + ':' + row['status']] += 1
    return dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True,
                        help='Parent of completed owning-chunk dump directories')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--strict', action='store_true',
                        help='Require source traces and fail on admitted call loss or overwrite')
    args = parser.parse_args()
    results = []
    for folder in sorted(args.root.iterdir(), key=lambda path: path.name):
        if not folder.is_dir():
            continue
        before = folder / 'matrix.chunk0.bam-overlay-input.tsv'
        after = folder / 'matrix.chunk0.bam-overlay-output.tsv'
        if not (before.exists() and after.exists()):
            continue
        results.append({
            'window': folder.name, 'source_to_input': audit_source(folder),
            'msa_admission': audit_admission(folder),
            'late_whole_bam_overlay': audit_overlay(load_matrix(before), load_matrix(after)),
        })
    if not results:
        parser.error(f'no completed BAM overlay dumps found under {args.root}')
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    for result in results:
        print(result['window'], result['source_to_input']['counts'],
              result['late_whole_bam_overlay']['counts'])
    if args.strict:
        for result in results:
            if not result['source_to_input']['available']:
                parser.error(f"source transfer trace missing for window {result['window']}")
            if result['source_to_input']['source_census'].get('read_not_transferred', 0):
                raise SystemExit(1)
            for stage in ('source_to_input', 'late_whole_bam_overlay'):
                counts = result[stage]['counts']
                if counts.get('lost', 0) or counts.get('changed', 0):
                    raise SystemExit(1)


if __name__ == '__main__':
    main()
