#!/usr/bin/env python3
"""Inventory unresolved local targets without using competitor calls as rescue input."""
import argparse
import csv
from pathlib import Path

import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--label', default='chr20_multi_final')
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
report = Path(__file__).resolve().parent
root = report.parents[1]
data = Path.home() / 'Downloads/pgphase-eval-data'
fields = ['region', 'target_status', 'left', 'right', 'gap_length',
          'primary_mapq30_spanning', 'intermediate_dv_het_snps',
          'intermediate_dv_het_indels', 'evidence_class']
with pysam.AlignmentFile(str(root / 'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam, \
     pysam.VariantFile(str(data / 'shared_calls/chr20/deepvariant.vcf.gz')) as calls, \
     (report / args.label / 'results.tsv').open() as results, args.output.open('w') as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    for row in csv.DictReader(results, delimiter='\t'):
        if row['recovery2_target'] == 'joined':
            continue
        left, right = int(row['left']), int(row['right'])
        spanning = {r.query_name for r in bam.fetch('CHM13#0#chr20', left - 1, left)
                    if not r.is_secondary and not r.is_supplementary and
                    r.mapping_quality >= 30 and r.reference_end >= right}
        snps = indels = 0
        for record in calls.fetch('chr20', left, right - 1):
            gt = next(iter(record.samples.values())).get('GT')
            if not gt or None in gt or len(set(gt)) < 2:
                continue
            if len(record.ref) == 1 and all(len(allele) == 1 for allele in record.alts):
                snps += 1
            else:
                indels += 1
        if row['recovery2_target'] == 'endpoint_unphased':
            reason = 'endpoint_allele_or_read_tag_audit'
        elif not spanning:
            reason = 'intermediate_evidence_required'
        elif len(spanning) == 1:
            reason = 'single_read_bridge_evidence'
        else:
            reason = 'spanning_read_observation_or_link_conflict'
        writer.writerow(dict(zip(fields, [row['region'], row['recovery2_target'], left, right,
                                         right - left, len(spanning), snps, indels, reason])))
