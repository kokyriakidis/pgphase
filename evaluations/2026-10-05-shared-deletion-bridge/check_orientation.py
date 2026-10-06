#!/usr/bin/env python3
"""Check parental orientation on disjoint original primary flanks."""
import json
from pathlib import Path
from collections import Counter
import pysam

left, right = 62408056, 62432427
truth = {f[0]: f[1] == 'PATERNAL' for line in open('test_data/derived/chr20_truth_hap.tsv')
         if len(f := line.rstrip().split('\t')) == 2}
original = pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
flanks = [set(), set()]
for read in original.fetch('CHM13#0#chr20', left - 50000, right + 50000):
    if read.is_secondary or read.is_supplementary or read.query_name not in truth:
        continue
    if read.reference_end < left:
        flanks[0].add(read.query_name)
    if read.reference_start >= right:
        flanks[1].add(read.query_name)
result = []
for name, vcf_path, bam_path, contig in [
    ('pgphase', 'test_data/tmp_gap_fix67/certified/graph/w62408056-62432427/native.vcf',
     'test_data/tmp_gap_fix67/certified/graph/w62408056-62432427/phased.bam', 'CHM13#0#chr20'),
    ('hiphase', 'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.vcf.gz',
     'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.bam', 'chr20')]:
    records = {}
    with pysam.VariantFile(vcf_path) as vcf:
        for row in vcf:
            if row.contig == contig and row.pos in (left, right):
                sample = next(iter(row.samples.values()))
                records[row.pos] = (sample['GT'], sample['PS'])
    assert len(records) == 2 and records[left] == records[right], records
    ps = records[left][1]
    votes = [Counter(), Counter()]
    with pysam.AlignmentFile(bam_path, check_sq=False) as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary or not read.has_tag('HP') or not read.has_tag('PS'):
                continue
            if read.get_tag('PS') != ps or read.get_tag('HP') not in (1, 2):
                continue
            for side in range(2):
                if read.query_name in flanks[side]:
                    votes[side][(read.get_tag('HP') == 1) != truth[read.query_name]] += 1
    orientation = votes[0].most_common(1)[0][0]
    assert all(sum(v.values()) >= 5 and v[orientation] / sum(v.values()) >= .9 for v in votes), votes
    result.append({'tool': name, 'phase_set': ps, 'boundary_genotypes': records[left][0],
                   'disjoint_flanks': [{'scored': sum(v.values()), 'correct': v[orientation],
                                       'concordance': v[orientation] / sum(v.values())} for v in votes],
                   'hap1_parent': 'MATERNAL' if orientation else 'PATERNAL'})
Path('evaluations/2026-10-05-shared-deletion-bridge/orientation.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
