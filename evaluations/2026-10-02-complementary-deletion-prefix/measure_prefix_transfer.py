#!/usr/bin/env python3
"""Evaluate the BAM-prefix transfer; parental truth is evaluation-only input."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True, help='Full chr20 output directory')
parser.add_argument('--after', type=Path, required=True, help='Full chr20 output directory')
parser.add_argument('--matrix', type=Path, required=True, help='Owning 64 Mb recovery-final matrix')
parser.add_argument('--bam', default='test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
parser.add_argument('--truth', default='test_data/derived/chr20_truth_hap.tsv')
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
truth = {f[0]: f[1] == 'PATERNAL' for line in open(args.truth)
         if len(f := line.rstrip().split('\t')) == 2}
left, right = 64128828, 64134226
names, spanning = set(), set()
with pysam.AlignmentFile(args.bam) as bam:
    for read in bam.fetch('CHM13#0#chr20', left - 1, right):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        names.add(read.query_name)
        if not (read.is_qcfail or read.is_duplicate) and \
                30 <= read.mapping_quality < 255 and \
                read.reference_start <= left - 1 and read.reference_end >= right:
            spanning.add(read.query_name)


def measure(directory):
    tags, votes, local_votes = {}, defaultdict(Counter), defaultdict(Counter)
    with pysam.AlignmentFile(str(directory / 'phased.bam'), check_sq=False) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary or \
                    not read.has_tag('HP') or not read.has_tag('PS'):
                continue
            hp, ps, name = read.get_tag('HP'), read.get_tag('PS'), read.query_name
            if hp not in (1, 2) or ps <= 0 or name not in truth:
                continue
            bit = (hp == 1) != truth[name]
            tags[name] = (ps, bit)
            votes[ps][bit] += 1
            if name in names:
                local_votes[ps][bit] += 1
    majority = {ps: counts[True] > counts[False] for ps, counts in votes.items()}
    correctness = {name: bit == majority[ps] for name, (ps, bit) in tags.items()}
    n = sum(sum(v.values()) for v in local_votes.values())
    c = sum(max(v.values()) for v in local_votes.values())
    dominant = max(local_votes.values(), key=lambda v: sum(v.values()))
    return correctness, {'phased': n, 'correct': c, 'discordant': n - c,
                         'dominant_correct': max(dominant.values()),
                         'truth_scorable_overlaps': len(names & truth.keys()),
                         'read_phase_sets': len(local_votes)}


before, local_before = measure(args.before)
after, local_after = measure(args.after)
transitions = Counter((before.get(n), after.get(n))
                      for n in before.keys() | after.keys())
variants, observations = {}, defaultdict(dict)
for line in args.matrix.open():
    fields = line.rstrip('\n').split('\t')
    if fields[0] == 'VAR':
        variants[int(fields[1])] = fields
    elif fields[0] == 'OBS':
        observations[fields[1]][int(fields[2])] = int(fields[5])
snp = [i for i, v in variants.items() if int(v[2]) == left and v[3] == 'X']
pair = sorted((i for i, v in variants.items() if int(v[2]) == right + 1 and v[3] == 'D'),
              key=lambda i: int(variants[i][6]))
assert len(snp) == 1 and len(pair) == 2
exclusive = Counter()
ambiguous = Counter()
for name in spanning:
    calls = observations.get(name, {})
    triple = tuple(calls.get(i, -1) for i in [snp[0], *pair])
    if triple[0] not in (0, 1):
        continue
    if triple[1:] in ((1, 0), (0, 1)):
        exclusive[triple] += 1
    else:
        ambiguous[triple] += 1
result = {'gap': [left, right], 'mapq30_spanning_reads': len(spanning),
          'exclusive_alt_pairs': {'/'.join(map(str, k)): v for k, v in exclusive.items()},
          'abstaining_pairs': {'/'.join(map(str, k)): v for k, v in ambiguous.items()},
          'local_before': local_before, 'local_after': local_after,
          'full_chr20_truth_transitions':
              {f'{a}->{b}': n for (a, b), n in transitions.items()}}
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
