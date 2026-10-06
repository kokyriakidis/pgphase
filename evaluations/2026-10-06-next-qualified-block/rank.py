#!/usr/bin/env python3
"""Rank current HiPhase spans and audit the next seam in the next split block."""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PG = ROOT/'test_data/tmp_gap_fix81/final/0'
HI = ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv'
truth = {f[0]: f[1] == 'PATERNAL'
         for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}


def blocks(path):
    positions = defaultdict(list)
    with pysam.VariantFile(str(path)) as vcf:
        for row in vcf:
            sample = next(iter(row.samples.values()))
            gt, ps = sample.get('GT'), sample.get('PS')
            if sample.phased and gt and None not in gt and len(set(gt)) > 1 and ps and ps > 0:
                positions[ps].append(row.pos)
    return sorted([{'ps': ps, 'left': min(pos), 'right': max(pos),
                    'length': max(pos)-min(pos)+1, 'rows': len(pos)}
                   for ps, pos in positions.items() if len(pos) > 1],
                  key=lambda b: b['length'], reverse=True)


pg_blocks, hi_blocks = blocks(PG/'phased.vcf'), blocks(HI/'phased.vcf.gz')
ranking = []
for rank, block in enumerate(hi_blocks, 1):
    overlap = sorted([b for b in pg_blocks if b['left'] <= block['right'] and b['right'] >= block['left']],
                     key=lambda b: b['left'])
    ranking.append({**block, 'hiphase_rank': rank,
                    'exact_span_covered': any(b['left'] <= block['left'] and b['right'] >= block['right'] for b in overlap),
                    'pgphase_blocks': overlap})
(OUT/'ranked-blocks.json').write_text(json.dumps(ranking, indent=2)+'\n')
