#!/usr/bin/env python3
"""Evaluation only: compare VCF extent gaps and parental read orientation.

The 25--30 Mb exclusion avoids centromeric cases. A spanning competitor PS
must have at least 95% purity on ten gap-overlapping reads and on five reads
in each disjoint 10-kb flank, with the same parental majority on both flanks.
Extent coverage is a nomination, not proof of an exact allele connection.
"""

from pathlib import Path
from collections import defaultdict, Counter
import pysam, json
import argparse
parser = argparse.ArgumentParser(description='Audit competitor-spanned pgphase VCF extent gaps against disjoint parental flanks.')
parser.add_argument('--pg-vcf', type=Path, required=True)
parser.add_argument('--pg-bam', type=Path, required=True)
parser.add_argument('--competitor-root', type=Path, required=True)
parser.add_argument('--bam', default='test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
parser.add_argument('--truth', default='test_data/derived/chr20_truth_hap.tsv')
parser.add_argument('--expectations', default='src/test_gap_windows_expect.tsv')
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
truth = {f[0]: f[1] == 'PATERNAL' for x in open(args.truth) if len((f := x.rstrip().split('\t'))) == 2}
root = args.competitor_root

def blocks(path):
    d = defaultdict(list)
    with pysam.VariantFile(str(path)) as v:
        for x in v:
            s = list(x.samples.values())[0]
            gt = s.get('GT')
            ps = s.get('PS')
            if s.phased and ps and gt and (len(set(gt)) > 1):
                d[ps].append(x.pos)
    return sorted(((min(x), max(x), ps) for ps, x in d.items()))

def tags(path):
    d = {}
    with pysam.AlignmentFile(str(path), check_sq=False) as b:
        for x in b:
            if x.is_secondary or x.is_supplementary:
                continue
            if x.has_tag('HP') and x.has_tag('PS') and (x.get_tag('HP') in [1, 2]):
                d[x.query_name] = (x.get_tag('HP'), x.get_tag('PS'))
    return d
pgblocks = blocks(args.pg_vcf)
hi = {arm: (blocks(root / arm / 'phased.vcf.gz'), tags(root / arm / 'phased.bam')) for arm in ['hiphase_dv', 'hiphase_pg_calls']}
pg = tags(args.pg_bam)
panel = set()
for x in open(args.expectations):
    f = x.rstrip().split('\t')
    if len(f) > 2 and f[0] == 'graph':
        panel.add(f[1])
gaps = []
prev = pgblocks[0]
for cur in pgblocks[1:]:
    if cur[0] > prev[1]:
        gaps.append((prev[1], cur[0]))
    if cur[1] > prev[1]:
        prev = cur
results = []
with pysam.AlignmentFile(args.bam) as b:
    for a, z in gaps:
        if a < 30000000 and z > 25000000:
            continue
        known = {x.query_name: (x.reference_start + 1, x.reference_end) for x in b.fetch('CHM13#0#chr20', max(0, a - 10001), z + 10000) if not x.is_secondary and (not x.is_supplementary) and (x.query_name in truth)}
        local = {n for n, (l, h) in known.items() if l <= z and h >= a}
        for arm, (bs, ts) in hi.items():
            for l, h, ps in bs:
                if not (l <= a and h >= z):
                    continue
                votes = [Counter(), Counter()]
                lv = Counter()
                for n, (beg, end) in known.items():
                    if n not in ts or ts[n][1] != ps:
                        continue
                    hp = ts[n][0]
                    bit = (hp == 1) != truth[n]
                    if n in local:
                        lv[bit] += 1
                    if end < a and end >= a - 10000:
                        votes[0][bit] += 1
                    if beg > z and beg <= z + 10000:
                        votes[1][bit] += 1
                if any((sum(x.values()) < 5 or max(x.values()) / sum(x.values()) < 0.95 for x in votes)):
                    continue
                if votes[0].most_common()[0][0] != votes[1].most_common()[0][0]:
                    continue
                if sum(lv.values()) < 10 or max(lv.values()) / sum(lv.values()) < 0.95:
                    continue
                pgb = defaultdict(Counter)
                for n in local:
                    if n in pg:
                        pgb[pg[n][1]][(pg[n][0] == 1) != truth[n]] += 1
                results.append(dict(left=a, right=z, length=z - a, arm=arm, hiphase_phased=sum(lv.values()), hiphase_correct=max(lv.values()), pg_phased=sum((sum(x.values()) for x in pgb.values())), pg_correct=sum((max(x.values()) for x in pgb.values())), flanks=[dict(x) for x in votes], tracked=f'{a}-{z}' in panel))
print(json.dumps(results, indent=2))
args.output.write_text(json.dumps(results, indent=2) + '\n')
print('GAPS', len(gaps), 'TARGETS', len({(x['left'], x['right']) for x in results}), 'UNTRACKED', [(x['left'], x['right'], x['arm']) for x in results if not x['tracked']])
