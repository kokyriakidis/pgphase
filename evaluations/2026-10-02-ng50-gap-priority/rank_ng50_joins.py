#!/usr/bin/env python3
"""Rank audited gap nominations by modeled whole-phase-set NG50 gain.

Evaluation only. A simulated union does not certify an allele connection or
its parental orientation. NG50 uses half the reference contig length; N50 uses
half the sum of phase-set spans, including any overlapping extents.
"""
import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path


def load_blocks(path, chrom):
    positions = defaultdict(list)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if fields[0] != chrom:
                continue
            sample = dict(zip(fields[8].split(':'), fields[9].split(':')))
            gt, ps = sample.get('GT', ''), sample.get('PS', '.')
            alleles = gt.split('|')
            if len(alleles) != 2 or '.' in alleles or alleles[0] == alleles[1]:
                continue
            if ps == '.' or int(ps) <= 0:
                continue
            positions[ps].append(int(fields[1]))
    return {ps: (min(pos), max(pos)) for ps, pos in positions.items()}


def nx(spans, denominator):
    accumulated = 0
    for span in sorted(spans, reverse=True):
        accumulated += span
        if 2 * accumulated >= denominator:
            return span
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', type=Path, required=True)
    parser.add_argument('--fai', type=Path, required=True)
    parser.add_argument('--chrom', required=True)
    parser.add_argument('--targets', type=Path, required=True,
                        help='JSON nominations from audit_competitor_gaps.py')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    lengths = {f[0]: int(f[1]) for line in args.fai.read_text().splitlines()
               if len(f := line.split('\t')) >= 2}
    length = lengths[args.chrom]
    blocks = load_blocks(args.vcf, args.chrom)
    spans = [end - beg + 1 for beg, end in blocks.values()]
    baseline = nx(spans, length)
    nominations = defaultdict(list)
    for row in json.loads(args.targets.read_text()):
        nominations[(row['left'], row['right'])].append(row)
    ranked, ambiguous = [], []
    for (left, right), evidence in nominations.items():
        left_blocks = [ps for ps, (_, end) in blocks.items() if end == left]
        right_blocks = [ps for ps, (beg, _) in blocks.items() if beg == right]
        if len(left_blocks) != 1 or len(right_blocks) != 1:
            ambiguous.append(dict(left=left, right=right))
            continue
        left_ps, right_ps = left_blocks[0], right_blocks[0]
        merged = max(blocks[left_ps][1], blocks[right_ps][1]) - \
            min(blocks[left_ps][0], blocks[right_ps][0]) + 1
        joined_spans = [end - beg + 1 for ps, (beg, end) in blocks.items()
                        if ps not in (left_ps, right_ps)] + [merged]
        projected = nx(joined_spans, length)
        ranked.append(dict(left=left, right=right, gap_bp=right - left,
                           left_ps=left_ps, right_ps=right_ps,
                           left_extent=blocks[left_ps], right_extent=blocks[right_ps],
                           merged_span_bp=merged, projected_ng50_bp=projected,
                           ng50_gain_bp=projected - baseline, evidence=evidence))
    ranked.sort(key=lambda row: (-row['ng50_gain_bp'], -row['merged_span_bp'],
                                 row['left'], row['right']))
    report = dict(vcf=str(args.vcf),
                  vcf_sha256=hashlib.sha256(args.vcf.read_bytes()).hexdigest(),
                  chrom=args.chrom, genome_length_bp=length, vcf_blocks=len(blocks),
                  baseline_ng50_bp=baseline, baseline_n50_bp=nx(spans, sum(spans)),
                  ranked=ranked, ambiguous_nominations=ambiguous)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(f"NG50 {baseline:,} bp; {len(ranked)} distinct audited nominations")
    for row in ranked[:10]:
        print(f"{row['left']}-{row['right']}: "
              f"NG50 +{row['ng50_gain_bp']:,} bp; "
              f"merged span {row['merged_span_bp']:,} bp")


if __name__ == '__main__':
    main()
