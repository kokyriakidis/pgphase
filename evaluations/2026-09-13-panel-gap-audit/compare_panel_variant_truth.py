#!/usr/bin/env python3
"""Compare local panel VCFs to truth with WhatsHap; counts overlap across windows."""
import argparse
import csv
import json
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--label', required=True)
parser.add_argument('--baseline-label', required=True)
parser.add_argument('--whatshap', default='whatshap')
parser.add_argument('--workers', type=int, default=4)
args = parser.parse_args()
report = Path(__file__).resolve().parent / args.label
output = Path('/tmp') / ('pgphase-' + args.label)
data = Path.home() / 'Downloads/pgphase-eval-data/results/chr12-18-20-comparison'
manifest = json.loads((report / 'manifest.json').read_text())
metrics = ('covered_variants', 'all_assessed_pairs', 'all_switches',
           'all_switchflips', 'blockwise_hamming', 'blockwise_hamming_rate')


def compare(case):
    region = case['region'].replace('CHM13#0#', '')
    chrom, interval = region.split(':')
    beg, end = map(int, interval.split('-'))
    directory = output / case['name'] / 'variant_truth'
    directory.mkdir(exist_ok=True)
    truth = directory / 'truth.vcf'
    with pysam.VariantFile(str(data / chrom / 'truth.vcf.gz')) as source:
        with pysam.VariantFile(str(truth), 'w', header=source.header) as dest:
            for record in source.fetch(chrom, beg - 1, end):
                dest.write(record)
    result = {'region': case['name']}
    for prefix, label in [('before', args.baseline_label), ('after', args.label)]:
        path = directory / (prefix + '.tsv')
        vcf = Path('/tmp') / ('pgphase-' + label) / case['name'] / 'recovery2/shared.vcf'
        command = [args.whatshap, 'compare', '--ignore-sample-name',
                   '--tsv-pairwise', str(path), str(truth), str(vcf)]
        (directory / (prefix + '.command.txt')).write_text(shlex.join(command) + '\n')
        with (directory / (prefix + '.log')).open('w') as log:
            subprocess.run(command, stdout=log, stderr=log, check=True)
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        if len(rows) != 1:
            raise ValueError(f'Expected one chromosome comparison for {case["name"]}: {rows}')
        for metric in metrics:
            result[prefix + '_' + metric] = rows[0][metric]
    return result


with ThreadPoolExecutor(max_workers=args.workers) as pool:
    results = list(pool.map(compare, manifest))
with (report / 'variant_truth_comparison.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(results)
print(json.dumps({'cases': len(results), 'hamming_count_increases': [
    row for row in results if int(row['after_blockwise_hamming']) > int(row['before_blockwise_hamming'])
]}, indent=2))
