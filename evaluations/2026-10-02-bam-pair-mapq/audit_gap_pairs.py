#!/usr/bin/env python3
"""Evaluate callable boundary pairs; truth is never an input to phasing."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import json
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--runs-root', type=Path, required=True)
parser.add_argument('--prefix', default='baseline')
parser.add_argument('--bam', default='test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
parser.add_argument('--truth', default='test_data/derived/chr20_truth_hap.tsv')
parser.add_argument('--competitor-root', type=Path, default=Path('/tmp/pgphase-hiphase-comparison-2026-10-01'))
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
truth = {f[0]: f[1] == 'PATERNAL' for line in open(args.truth)
         if len(f := line.rstrip().split('\t')) == 2}


tag_cache = {}


def score(path, names):
    # Competitor BAMs contain full read sequences; decode each file once.
    if path not in tag_cache:
        tags = {}
        with pysam.AlignmentFile(str(path), check_sq=False) as bam:
            for read in bam:
                if read.is_secondary or read.is_supplementary or \
                        read.query_name not in truth or not read.has_tag('HP') or \
                        not read.has_tag('PS'):
                    continue
                hp, ps = read.get_tag('HP'), read.get_tag('PS')
                if hp in (1, 2) and ps > 0:
                    tags[read.query_name] = (hp, ps)
        tag_cache[path] = tags
    counts = defaultdict(Counter)
    for name in names:
        if name in tag_cache[path]:
            hp, ps = tag_cache[path][name]
            counts[ps][(hp == 1) != truth[name]] += 1
    phased = sum(sum(x.values()) for x in counts.values())
    correct = sum(max(x.values()) for x in counts.values())
    return {'phased': phased, 'correct': correct, 'discordant': phased - correct}


results = []
for chunk, left, right in [(4, 4866153, 4874129), (24, 24121713, 24131707),
                           (34, 34094604, 34102867), (41, 41879449, 41880908),
                           (61, 61738239, 61747506), (64, 64128828, 64134226)]:
    directory = args.runs_root / f'{args.prefix}-{chunk}'
    variants, obs = {}, {}
    for line in (directory / 'matrix.chunk0.recovery-final.tsv').open():
        f = line.rstrip('\n').split('\t')
        if f[0] == 'VAR':
            variants[int(f[1])] = f
        elif f[0] == 'OBS':
            obs.setdefault(f[1], {})[int(f[2])] = tuple(map(int, f[3:6]))
    names, callable_names = set(), set()
    with pysam.AlignmentFile(args.bam) as bam:
        for read in bam.fetch('CHM13#0#chr20', left - 1, right):
            if read.flag & (4 | 256 | 2048 | 512 | 1024):
                continue
            names.add(read.query_name)
            if 30 <= read.mapping_quality < 255 and read.reference_start <= left - 1 \
                    and read.reference_end >= right:
                callable_names.add(read.query_name)
    # Internal indel coordinates follow the VCF anchor; SNP coordinates do not.
    def at_anchor(row, anchor):
        return int(row[2]) - (row[3] in ('I', 'D')) == anchor and \
            int(row[8]) > 0 and row[9] != row[10]
    left_rows = [i for i, row in variants.items() if at_anchor(row, left)]
    right_rows = [i for i, row in variants.items() if at_anchor(row, right)]
    pairs = []
    for li in left_rows:
        for ri in right_rows:
            counts = Counter((z[li][2], z[ri][2]) for name, z in obs.items()
                             if name in callable_names and li in z and ri in z
                             and z[li][2] in (0, 1) and z[ri][2] in (0, 1))
            pairs.append({'left_row': variants[li][2:8], 'right_row': variants[ri][2:8],
                          'bam_allele_pairs': {f'{a}/{b}': n for (a, b), n in sorted(counts.items())}})
    result = {'chunk': chunk, 'left': left, 'right': right,
              'mapq30_spanning_primary_reads': len(callable_names), 'pairs': pairs,
              'pgphase': score(directory / 'phased.bam', names)}
    for arm in ('hiphase_dv', 'hiphase_pg_calls'):
        result[arm] = score(args.competitor_root / arm / 'phased.bam', names)
    if chunk == 64:
        counts = Counter((z[left_rows[0]][2], z[right_rows[0]][2], z[right_rows[1]][2])
                         for name, z in obs.items() if name in callable_names
                         and all(i in z for i in left_rows + right_rows)
                         and z[left_rows[0]][2] in (0, 1)
                         and (z[right_rows[0]][2], z[right_rows[1]][2]) in ((0, 1), (1, 0)))
        result['exclusive_deletion_alt_pairs'] = {'/'.join(map(str, k)): n for k, n in sorted(counts.items())}
    results.append(result)
args.output.write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps([{k: v for k, v in result.items() if k != 'pairs'} for result in results], indent=2))
