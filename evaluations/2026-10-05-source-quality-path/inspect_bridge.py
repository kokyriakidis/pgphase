#!/usr/bin/env python3
"""Record primary physical SNP calls and source assignments at the new gap."""
import argparse
import json
from pathlib import Path
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--bam', default='test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
args = parser.parse_args()
source_reads = {}
for line in args.source.read_text().splitlines():
    fields = line.split('\t')
    if fields[0] == 'READ':
        source_reads[fields[1]] = {'hp': int(fields[6]), 'ps': int(fields[7])}
left, right = 13752640, 13773452
pairs = []
with pysam.AlignmentFile(args.bam) as bam:
    for read in bam.fetch('CHM13#0#chr20', left - 1, right):
        if read.flag & (4 | 256 | 2048 | 1024 | 512) or read.mapping_quality < 30 or read.mapping_quality == 255:
            continue
        calls = {}
        for query_pos, ref_pos in read.get_aligned_pairs():
            if query_pos is not None and ref_pos in (left - 1, right - 1):
                calls[ref_pos + 1] = (read.query_sequence[query_pos], int(read.query_qualities[query_pos]))
        if set(calls) != {left, right}:
            continue
        if calls[left][0] not in 'CT' or calls[right][0] not in 'GC':
            continue
        if any(q < 20 or q == 255 for _, q in calls.values()):
            continue
        pairs.append({'qname': read.query_name, 'mapq': read.mapping_quality,
                      'calls': calls, 'source': source_reads.get(read.query_name),
                      'error_bound': sum(10 ** (-q / 10) for _, q in calls.values()) +
                                     2 * 10 ** (-read.mapping_quality / 10)})
args.output.write_text(json.dumps({'left': left, 'right': right, 'primary_pairs': pairs}, indent=2) + '\n')
assert len(pairs) == 2
assert all(p['source'] == {'hp': 1, 'ps': 13751618} for p in pairs)
assert all(p['calls'][left][0] == 'C' and p['calls'][right][0] == 'G' for p in pairs)
