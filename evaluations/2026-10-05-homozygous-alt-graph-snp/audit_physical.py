#!/usr/bin/env python3
"""Reproduce the original-molecule homozygous ALT certificate, without truth."""
from collections import Counter
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
pos = 36620864
calls = []
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    for read in bam.fetch('CHM13#0#chr20', pos-1, pos):
        if read.flag & (4 | 256 | 2048 | 512 | 1024) or read.mapping_quality < 30 or read.mapping_quality == 255:
            continue
        rp, qp = read.reference_start, 0
        for op, length in read.cigartuples:
            if op in (0, 2, 3, 7, 8):
                if rp <= pos-1 < rp+length:
                    if op == 2:
                        calls.append({'read': read.query_name, 'base': 'DEL'})
                    elif op in (0, 7, 8):
                        qi = qp + pos-1-rp
                        quality = read.query_qualities[qi]
                        if quality >= 20:
                            calls.append({'read': read.query_name, 'base': read.query_sequence[qi], 'quality': quality})
                rp += length
            if op in (0, 1, 4, 7, 8):
                qp += length
counts = Counter(c['base'] for c in calls)
assert counts == {'A': 54}, counts
# All catalog rows are a conservative upper bound on the tested chunk SNPs.
with pysam.VariantFile(str(ROOT/'test_data/chr20.sites.striped.vcf.gz')) as vcf:
    family_upper_bound = sum(1 for _ in vcf)
p = 2.0**-counts['A'] * family_upper_bound
assert p <= .01
result = {'position': pos, 'ref': 'G', 'alt': 'A', 'counts': dict(counts),
          'catalog_family_upper_bound': family_upper_bound, 'conservative_adjusted_p': p,
          'mapq_floor': 30, 'base_quality_floor': 20, 'calls': calls}
(OUT/'physical-evidence.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: v for k, v in result.items() if k != 'calls'}))
