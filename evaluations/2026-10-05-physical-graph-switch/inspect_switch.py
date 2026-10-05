#!/usr/bin/env python3
"""Inspect primary BAM evidence for the graph block's internal SNP switch."""
import argparse
from collections import Counter
import json
import math
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bam', default='test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
parser.add_argument('--vcf', default='test_data/tmp_gap_fix55/full-final/phased.vcf')
parser.add_argument('--output', required=True)
args = parser.parse_args()
positions = (62716326, 62718395, 62722021, 62722416)
rows = {}
with pysam.VariantFile(args.vcf) as vcf:
    for record in vcf:
        if record.pos not in positions:
            continue
        sample = next(iter(record.samples.values()))
        assert sample.phased and len(record.ref) == len(record.alts[0]) == 1
        rows[record.pos] = (record.ref, record.alts[0], sample['GT'][0])
assert set(rows) == set(positions)
result = []
with pysam.AlignmentFile(args.bam) as bam:
    for left, right in zip(positions, positions[1:]):
        votes, allele_counts, seen = Counter(), Counter(), set()
        log_odds = 0.0
        for read in bam.fetch('CHM13#0#chr20', left - 1, right):
            if read.flag & (4 | 256 | 2048 | 512 | 1024) or not 30 <= read.mapping_quality < 255:
                continue
            if read.query_name in seen:
                continue
            aligned = {p: q for q, p in read.get_aligned_pairs()
                       if p is not None and q is not None}
            qleft, qright = aligned.get(left - 1), aligned.get(right - 1)
            if qleft is None or qright is None:
                continue
            qualities = (read.query_qualities[qleft], read.query_qualities[qright])
            if any(not 30 <= q < 255 for q in qualities):
                continue
            bases = (read.query_sequence[qleft], read.query_sequence[qright])
            if any(base not in rows[pos][:2] for base, pos in zip(bases, (left, right))):
                continue
            seen.add(read.query_name)
            haps = [(base == rows[pos][1]) == (rows[pos][2] == 1)
                    for base, pos in zip(bases, (left, right))]
            votes['same' if haps[0] == haps[1] else 'cross'] += 1
            allele_counts['hap1' if haps[0] else 'hap2'] += 1
            error = sum(10 ** (-q / 10) for q in qualities) + 2 * 10 ** (-read.mapping_quality / 10)
            log_odds += (1 if haps[0] != haps[1] else -1) * math.log((1 - error) / error)
        result.append({'left': left, 'right': right, 'votes': dict(votes),
                       'represented_haps': dict(allele_counts),
                       'log_odds_for_reversal': log_odds,
                       'two_sided_unanimous_count_p': math.ldexp(2.0, -sum(votes.values()))})
with open(args.output, 'w') as out:
    json.dump(result, out, indent=2)
    out.write('\n')
print(json.dumps(result, indent=2))
