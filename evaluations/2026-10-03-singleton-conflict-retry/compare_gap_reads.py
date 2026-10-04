#!/usr/bin/env python3
"""Score gap-overlapping reads using each complete output block's truth gauge.

Evaluation only: truth and competitor output are never inputs to phasing.
"""
import argparse
from collections import Counter, defaultdict
import json
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input-bam', required=True)
parser.add_argument('--truth', required=True)
parser.add_argument('--left', type=int, required=True)
parser.add_argument('--right', type=int, required=True)
parser.add_argument('--output-bam', action='append', required=True)
args = parser.parse_args()
truth = {f[0]: f[1] == 'PATERNAL' for line in open(args.truth)
         if len(f := line.rstrip().split('\t')) == 2}
with pysam.AlignmentFile(args.input_bam) as bam:
    names = {read.query_name for read in
             bam.fetch('CHM13#0#chr20', args.left - 1, args.right)
             if not read.is_secondary and not read.is_supplementary}
results = []
for path in args.output_bam:
    tags, votes = {}, defaultdict(Counter)
    with pysam.AlignmentFile(path, check_sq=False) as bam:
        for read in bam:
            name = read.query_name
            if read.is_secondary or read.is_supplementary or name not in truth:
                continue
            hp = read.get_tag('HP') if read.has_tag('HP') else 0
            ps = read.get_tag('PS') if read.has_tag('PS') else 0
            if hp in (1, 2) and ps > 0:
                tags[name] = (ps, (hp == 1) != truth[name])
    for ps, orientation in tags.values():
        votes[ps][orientation] += 1
    gauges = {ps: count[True] > count[False] for ps, count in votes.items()}
    local_correct, local_phased = Counter(), Counter()
    for name in names & tags.keys():
        ps, orientation = tags[name]
        local_phased[ps] += 1
        local_correct[ps] += orientation == gauges[ps]
    phased, correct = sum(local_phased.values()), sum(local_correct.values())
    results.append(dict(path=path, physical_overlaps=len(names),
                        truth_overlaps=len(names & truth.keys()), phased=phased,
                        correct=correct, discordant=phased - correct,
                        dominant_correct=max(local_correct.values(), default=0),
                        blocks=len(local_phased)))
print(json.dumps(results, indent=2))
