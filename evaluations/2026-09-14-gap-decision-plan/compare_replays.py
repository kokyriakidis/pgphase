#!/usr/bin/env python3
"""Compare regional replay assignments without treating file headers as evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics

import pysam


def signature(path, bam=False):
    if bam:
        with pysam.AlignmentFile(path, 'rb') as stream:
            rows = [(r.query_name, r.reference_id, r.reference_start, r.flag,
                     r.cigarstring or '', r.get_tag('HP') if r.has_tag('HP') else 0,
                     r.get_tag('PS') if r.has_tag('PS') else -1)
                    for r in stream.fetch(until_eof=True)]
    else:
        rows = [line.rstrip() for line in path.read_text().splitlines()
                if line and not line.startswith('#')]
    return len(rows), hashlib.sha256(json.dumps(sorted(rows)).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--left', type=Path, required=True)
    parser.add_argument('--right', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    a = json.loads((args.left / 'manifest.json').read_text())
    b = json.loads((args.right / 'manifest.json').read_text())
    if a['targets'] != b['targets'] or a['binary_sha256'] != b['binary_sha256']:
        raise SystemExit('Comparison requires identical targets and binary hashes')
    checks, times = [], [[], []]
    for chrom, left, right in a['targets']:
        name = f"{chrom.replace('#', '_')}_{left}_{right}"
        results = [json.loads((root / name / 'result.json').read_text())
                   for root in (args.left, args.right)]
        for arm in ('baseline', 'graph_bam'):
            paths = [root / name / arm for root in (args.left, args.right)]
            read_a, read_b = [signature(p / 'phased.bam', True) for p in paths]
            vcf_a, vcf_b = [signature(p / 'native.vcf') for p in paths]
            checks.append(dict(target=name, arm=arm, reads=read_a[0],
                               same_reads=read_a == read_b,
                               same_variants=vcf_a == vcf_b,
                               same_endpoint_status=results[0][arm + '_status'] == results[1][arm + '_status'],
                               endpoint_status=results[1][arm + '_status']))
            for index in (0, 1):
                times[index].append(results[index][arm + '_seconds'])
    passed = all(r['same_reads'] and r['same_variants'] and r['same_endpoint_status']
                 for r in checks)
    report = dict(passed=passed, binary_sha256=a['binary_sha256'],
                  left=str(args.left), right=str(args.right),
                  left_workers=[a['args']['jobs'], a['args']['threads']],
                  right_workers=[b['args']['jobs'], b['args']['threads']],
                  median_seconds_per_arm=[statistics.median(t) for t in times],
                  checks=checks,
                  limitation='Regional replay only; does not prove chromosome-context order independence or truth accuracy.')
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(passed=passed, comparisons=len(checks),
                          median_seconds_per_arm=report['median_seconds_per_arm'])))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
