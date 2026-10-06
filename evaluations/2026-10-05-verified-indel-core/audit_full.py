#!/usr/bin/env python3
"""Verify that core materialization preserves variants and parental assignments."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True)
parser.add_argument('--after', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--end', type=int, help='Audit complete chunks before this endpoint')
args = parser.parse_args()
truth = {f[0]: f[1] == 'PATERNAL' for line in open('test_data/derived/chr20_truth_hap.tsv')
         if len(f := line.rstrip().split('\t')) == 2}
eligible = None
if args.end:
    with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as bam:
        eligible = {read.query_name for read in bam.fetch('CHM13#0#chr20', 0, args.end)
                    if not read.is_secondary and not read.is_supplementary
                    and read.reference_end <= args.end}


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
    if eligible is not None:
        tags = {name: tag for name, tag in tags.items() if name in eligible}
        status = {name: value for name, value in status.items() if name in eligible}
    return tags, status


def variants(folder):
    return [line for line in (folder / 'phased.vcf').read_text().splitlines()
            if not line.startswith('#') and
            (args.end is None or int(line.split('\t')[1]) <= args.end)]


old_variants, new_variants = variants(args.before), variants(args.after)
assert old_variants == new_variants, 'Variant records or phase blocks changed'
before, old_status = reads(args.before)
after, new_status = reads(args.after)
assert before.keys() == after.keys(), 'Primary read set changed'
changed = [name for name in before if before[name] != after[name]]
assert all(before[name][1] >= 1_000_000_000 and 0 < after[name][1] < 1_000_000_000
           for name in changed), 'An existing core or unphased assignment changed'
assert all(before[name][0] == after[name][0] for name in changed), 'An HP tag changed'
assert {n for n, s in old_status.items() if s is None} == {
    n for n, s in new_status.items() if s is None}, 'Abstentions changed'
assert sum(s is True for s in new_status.values()) >= sum(
    s is True for s in old_status.values()), 'Overall parental correctness fell'
result = {'variant_records_unchanged': len(old_variants),
          'primary_reads_unchanged': len(before),
          'truth_before': dict(Counter(str(s) for s in old_status.values())),
          'truth_after': dict(Counter(str(s) for s in new_status.values())),
          'truth_transitions': {str(k): v for k, v in Counter(
              (old_status[n], new_status[n]) for n in old_status).items()},
          'rescues_materialized': len(changed),
          'changed_reads': [{'read': name, 'before': before[name], 'after': after[name],
                             'correct': new_status.get(name)} for name in changed]}
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'changed_reads'}, indent=2))
