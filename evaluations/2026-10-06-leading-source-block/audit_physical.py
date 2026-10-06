#!/usr/bin/env python3
"""Compare selected repeat SNPs with original primary CIGAR bases and HiPhase variants."""
from collections import Counter
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
report = {}
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    for pos in (7280346, 11813622, 2282198, 42290009, 47003854, 50431256, 62940433, 63734075, 63793665):
        counts = {10: Counter(), 20: Counter()}
        for read in bam.fetch('CHM13#0#chr20', pos-1, pos):
            if read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.mapping_quality < 30 or read.mapping_quality == 255:
                continue
            rp, qp = read.reference_start+1, 0
            for op, length in read.cigartuples:
                if op in (0, 7, 8) and rp <= pos < rp+length:
                    qi = qp+pos-rp
                    for min_bq in counts:
                        if read.query_qualities[qi] >= min_bq:
                            counts[min_bq][read.query_sequence[qi]] += 1
                if op == 2 and rp <= pos < rp+length:
                    for min_bq in counts:
                        counts[min_bq]['DEL'] += 1
                if op in (0, 2, 3, 7, 8): rp += length
                if op in (0, 1, 4, 7, 8): qp += length
        report[str(pos)] = {str(bq): dict(count) for bq, count in counts.items()}
assert report['11813622']['10'] == {'T': 61, 'C': 8}
with pysam.VariantFile(str(ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.vcf.gz')) as vcf:
    report['hiphase_input_variants'] = [str(row).rstrip() for row in vcf.fetch('chr20', 11813610, 11813680)]
(OUT/'physical-evidence.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
