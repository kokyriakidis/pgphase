#!/usr/bin/env python3
"""Prepare same-read HiPhase controls; alternate alleles are evaluation-only."""
import argparse
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--replays', type=Path, default=ROOT/'test_data/tmp_joint_evidence_investigation')
args = parser.parse_args()
directory = args.replays/'hiphase5'
directory.mkdir(exist_ok=True)
left, right = 5250000, 5400000
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as source, \
        pysam.AlignmentFile(str(directory/'input.bam'), 'wb', header=source.header) as destination:
    for read in source.fetch('CHM13#0#chr20', left, right):
        destination.write(read)
pysam.index(str(directory/'input.bam'))
with pysam.VariantFile(str(args.replays/'5/phased.vcf')) as source:
    header = source.header.copy()
    records = [row.copy() for row in source if left <= row.start < right]
with pysam.VariantFile(str(ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.vcf.gz')) as source:
    alternate = [row.copy() for row in source.fetch('chr20', left, right)]

def converted(row):
    sample = next(iter(row.samples.values()))
    result = header.new_record(contig='CHM13#0#chr20', start=row.start, alleles=row.alleles)
    result.samples[0]['GT'] = sample['GT']
    result.samples[0].phased = False
    result.samples[0]['PS'] = None
    return result

insertion = next(row for row in alternate if row.pos == 5315591)
contrast = next(row for row in alternate if row.pos == 5339363)
for mode in ['input', 'dv-input', 'add-marker-input', 'diploid-input', 'both-input']:
    selected = alternate if mode == 'dv-input' else records
    if mode in ('diploid-input', 'both-input'):
        selected = [row for row in selected if row.pos != 5339368] + [contrast]
    if mode in ('add-marker-input', 'both-input'):
        selected = selected + [insertion]
    selected = sorted(selected, key=lambda row: row.pos)
    path = directory/(mode+'.vcf.gz')
    with pysam.VariantFile(str(path), 'wz', header=header) as destination:
        for row in selected:
            if mode != 'dv-input' and row.contig == 'CHM13#0#chr20':
                row.translate(header)
                for sample in row.samples.values():
                    sample.phased = False
                    sample['PS'] = None
                destination.write(row)
            else:
                destination.write(converted(row))
    pysam.tabix_index(str(path), preset='vcf', force=True)
