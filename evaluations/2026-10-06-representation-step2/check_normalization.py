#!/usr/bin/env python3
"""Compare the production identity helper with bcftools and frozen source records."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import subprocess
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--vcf', type=Path, default=Path('test_data/tmp_representation_step1/full-current/phased.vcf'))
parser.add_argument('--reference', type=Path, default=Path('test_data/chm13v2.0.chr20.renamed.fa'))
parser.add_argument('--bcftools', default='/home/kokyriakidis/micromamba/envs/bench-phasers/bin/bcftools')
parser.add_argument('--out', type=Path, default=Path('evaluations/2026-10-06-representation-step2/normalization-checks.json'))
args = parser.parse_args()
original = [line.split('\t') for line in args.vcf.read_text().splitlines() if not line.startswith('#')]
normalized = subprocess.run([args.bcftools, 'norm', '-f', str(args.reference), '--old-rec-tag', 'ORIGINAL', '-Ov', str(args.vcf)], text=True, capture_output=True, check=True)
rows = [line.split('\t') for line in normalized.stdout.splitlines() if not line.startswith('#')]
source_by_key = {(row[0], row[1], row[3], row[4]): i for i, row in enumerate(original)}
assert len(source_by_key) == len(original), 'frozen input must have unique raw records'
expected = [None] * len(original)
changed = 0
for row in rows:
    info = dict(field.split('=', 1) for field in row[7].split(';') if '=' in field)
    chrom, pos, ref, alt = info['ORIGINAL'].split('|')[:4] if 'ORIGINAL' in info else (row[0], row[1], row[3], row[4])
    index = source_by_key[chrom, pos, ref, alt]
    changed += 'ORIGINAL' in info
    # bcftools supplies the leftmost anchored representation. Independently
    # remove its anchor/context; do not left-align with the helper being tested.
    pos, ref, alt = int(row[1]), row[3].upper(), row[4].upper()
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    pos += prefix
    ref, alt = ref[prefix:], alt[prefix:]
    kind = 8 if len(ref) == len(alt) == 1 else 1 if len(alt) > len(ref) else 2
    expected[index] = f'{pos if kind == 8 else pos - 1}\t{kind}\t{len(ref)}\t{alt or "."}'
assert all(value is not None for value in expected)
input_text = ''.join(f'{row[1]}\t{row[3]}\t{row[4]}\n' for row in original)
start = time.monotonic()
actual = subprocess.run(['./test_allele_identity', '--normalize', str(args.reference)], input=input_text, text=True, capture_output=True, check=True).stdout.splitlines()
seconds = time.monotonic() - start
assert len(actual) == len(expected)
differences = [dict(index=i, original=original[i][:5], expected=want, actual=got) for i, (want, got) in enumerate(zip(expected, actual)) if want != got]
groups = defaultdict(list)
for row, identity in zip(original, actual):
    sample = dict(zip(row[8].split(':'), row[9].split(':')))
    groups[identity].append(dict(pos=int(row[1]), ref=row[3], alt=row[4], gt=sample.get('GT'), ps=sample.get('PS'), ad=sample.get('AD')))
counts = Counter(actual)
duplicate_groups = sum(count > 1 for count in counts.values())
result = dict(records=len(actual), bcftools_realigned_records=changed, differences=len(differences), examples=differences[:20], normalized_duplicate_groups=duplicate_groups, normalized_extra_rows=sum(count - 1 for count in counts.values()), duplicate_gt_ps_disagreements=sum(len({(row['gt'], row['ps']) for row in group}) > 1 for group in groups.values() if len(group) > 1), duplicate_ad_disagreements=sum(len({row['ad'] for row in group}) > 1 for group in groups.values() if len(group) > 1), duplicate_groups=[dict(identity=identity, originals=records) for identity, records in groups.items() if len(records) > 1], helper_seconds=seconds, source_sha256=hashlib.sha256(args.vcf.read_bytes()).hexdigest(), helper_sha256=hashlib.sha256(Path('test_allele_identity').read_bytes()).hexdigest(), oracle_version=subprocess.check_output([args.bcftools, '--version'], text=True).splitlines()[0], oracle_stderr=normalized.stderr.strip())
args.out.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
assert not differences
assert duplicate_groups == sum(count > 1 for count in Counter(expected).values())
