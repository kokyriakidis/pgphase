#!/usr/bin/env python3
"""Compare phasing outputs exactly and score read tags by parental truth."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import json
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True)
parser.add_argument('--after', type=Path, required=True)
parser.add_argument('--truth', default='test_data/derived/chr20_truth_hap.tsv')
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
truth = {f[0]: f[1] == 'PATERNAL' for line in open(args.truth)
         if len(f := line.rstrip().split('\t')) == 2}


def measure(directory):
    tags, rows, votes, blocks = {}, {}, defaultdict(Counter), defaultdict(list)
    with pysam.AlignmentFile(str(directory / 'phased.bam'), check_sq=False) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            hp = read.get_tag('HP') if read.has_tag('HP') else 0
            ps = read.get_tag('PS') if read.has_tag('PS') else 0
            tags[read.query_name] = (hp, ps)
            if hp in (1, 2) and ps > 0 and read.query_name in truth:
                votes[ps][(hp == 1) != truth[read.query_name]] += 1
    for line in (directory / 'phased.vcf').open():
        if line.startswith('#'):
            continue
        f = line.rstrip().split('\t')
        row = dict(zip(f[8].split(':'), f[9].split(':')))
        key = (int(f[1]), f[3], f[4])
        rows[key] = tuple(f)
        gt, ps = row['GT'], row.get('PS')
        if '|' in gt and len(set(gt.split('|'))) == 2 and ps not in (None, '.', '0'):
            blocks[ps].append(key[0])
    lengths = sorted((max(p) - min(p) + 1 for p in blocks.values()), reverse=True)
    acc, n50 = 0, 0
    for length in lengths:
        acc += length
        if acc >= sum(lengths) / 2:
            n50 = length
            break
    scored = sum(sum(v.values()) for v in votes.values())
    correct = sum(max(v.values()) for v in votes.values())
    return tags, rows, {'output_reads': len(tags), 'scored': scored, 'correct': correct,
                        'discordant': scored - correct, 'concordance': correct / scored,
                        'read_phase_sets': len(votes), 'vcf_keys': len(rows),
                        'vcf_blocks': len(blocks), 'span_n50': n50}


bt, bv, before = measure(args.before)
at, av, after = measure(args.after)
result = {'before': before, 'after': after,
          'changed_read_tags': sum(bt.get(k) != at.get(k) for k in bt.keys() | at.keys()),
          'changed_vcf_rows': sum(bv.get(k) != av.get(k) for k in bv.keys() | av.keys()),
          'lost_keys': len(bv.keys() - av.keys()), 'gained_keys': len(av.keys() - bv.keys())}
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
