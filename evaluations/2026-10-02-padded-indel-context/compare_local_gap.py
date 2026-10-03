#!/usr/bin/env python3
"""Score reads overlapping a gap without importing truth into phasing."""
import argparse
from collections import Counter, defaultdict
import json
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bam', required=True, help='Input BAM defining overlap read names')
parser.add_argument('--truth', required=True)
parser.add_argument('--left', type=int, required=True)
parser.add_argument('--right', type=int, required=True)
parser.add_argument('--output-bam', action='append', required=True)
args = parser.parse_args()
truth = {f[0]: f[1] == 'PATERNAL' for line in open(args.truth)
         if len(f := line.rstrip().split('\t')) == 2}
with pysam.AlignmentFile(args.bam) as bam:
    names = {read.query_name for read in bam.fetch('CHM13#0#chr20', args.left - 1, args.right)
             if not read.is_secondary and not read.is_supplementary}
results = []
for path in args.output_bam:
    votes = defaultdict(Counter)
    seen = set()
    with pysam.AlignmentFile(path, check_sq=False) as bam:
        # Competitor BAMs are indexed; pgphase tag-only output has no SQ/index.
        if bam.references and bam.has_index():
            contig = 'chr20' if 'chr20' in bam.references else 'CHM13#0#chr20'
            records = bam.fetch(contig, args.left - 1, args.right)
        else:
            records = iter(bam)
        for read in records:
            name = read.query_name
            if name not in names or name not in truth or name in seen:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if not read.has_tag('HP') or not read.has_tag('PS'):
                continue
            hp, ps = read.get_tag('HP'), read.get_tag('PS')
            if hp not in (1, 2) or ps <= 0:
                continue
            seen.add(name)
            votes[ps][(hp == 1) != truth[name]] += 1
    phased = sum(sum(v.values()) for v in votes.values())
    correct = sum(max(v.values()) for v in votes.values())
    results.append(dict(path=path, phased=phased, correct=correct,
                        discordant=phased - correct, blocks=len(votes),
                        dominant_correct=max((max(v.values()) for v in votes.values()), default=0),
                        overlap_reads=len(names)))
print(json.dumps(results, indent=2))
