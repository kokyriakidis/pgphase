#!/usr/bin/env python3
"""Audit unchanged variants/tags and measure new assignments against parents."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True)
parser.add_argument('--after', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
truth = {f[0]: f[1] == 'PATERNAL' for line in open('test_data/derived/chr20_truth_hap.tsv')
         if len(f := line.rstrip().split('\t')) == 2}

def reads(folder):
    tags, votes = {}, defaultdict(Counter)
    with pysam.AlignmentFile(str(folder / 'phased.bam'), check_sq=False) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            hp = read.get_tag('HP') if read.has_tag('HP') else 0
            ps = read.get_tag('PS') if read.has_tag('PS') else 0
            tags[read.query_name] = (hp, ps)
            if hp in (1, 2) and ps > 0 and read.query_name in truth:
                votes[ps][(hp == 1) != truth[read.query_name]] += 1
    orientation = {ps: c.most_common()[0][0] for ps, c in votes.items()}
    status = {name: (((hp == 1) != truth[name]) == orientation[ps]
                     if hp in (1, 2) and ps in orientation else None)
              for name, (hp, ps) in tags.items() if name in truth}
    return tags, status

before, old_status = reads(args.before)
after, new_status = reads(args.after)
def variants(folder):
    path = folder / 'phased.vcf'
    if not path.exists():
        path = folder / 'native.vcf'
    return [line for line in path.read_text().splitlines() if not line.startswith('#')]

old_variants, new_variants = variants(args.before), variants(args.after)
assert old_variants == new_variants, 'Variant records or phase blocks changed'
assert before.keys() == after.keys(), 'Primary output read set changed'
changed = [name for name in before if before[name] != after[name]]
assert all(before[name][0] == 0 and after[name][0] in (1, 2) and
           0 < after[name][1] < 1_000_000_000 for name in changed), 'Existing assignment changed'
assert all(new_status[n] == s for n, s in old_status.items() if s is not None), 'Existing parental correctness changed'
result = {'variant_records_unchanged': len(old_variants), 'primary_reads_unchanged': len(before),
          'existing_assignments_preserved': sum(tag[0] in (1, 2) for tag in before.values()),
          'new_assignments': len(changed),
          'truth_before': dict(Counter(str(s) for s in old_status.values())),
          'truth_after': dict(Counter(str(s) for s in new_status.values())),
          'changed_reads': [{'read': name, 'before': before[name], 'after': after[name],
                            'correct': new_status.get(name)} for name in changed]}
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'changed_reads'}, indent=2))
