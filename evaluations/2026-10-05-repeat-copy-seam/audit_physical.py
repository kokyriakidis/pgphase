#!/usr/bin/env python3
"""Reproduce original-sequence contradictions; truth is diagnostic only."""
from collections import Counter
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
MATRIX = ROOT/'test_data/tmp_gap_fix73/baseline/25/matrix.chunk0.flags4236.tsv'
observations = {25853217: {}, 25855631: {}}
indices = {}
skipped = set()
for line in MATRIX.open():
    fields = line.rstrip().split('\t')
    if fields[0] == 'VAR' and int(fields[2]) in (25853217, 25855630):
        indices[int(fields[1])] = 25853217 if int(fields[2]) == 25853217 else 25855631
    elif fields[0] == 'READ' and fields[6] == '1':
        skipped.add(fields[1])
    elif fields[0] == 'OBS' and int(fields[2]) in indices and fields[3] in ('0', '1'):
        observations[indices[int(fields[2])]][fields[1]] = int(fields[3])
for cohort in observations.values():
    for qname in skipped:
        cohort.pop(qname, None)
assert observations[25853217].keys() == observations[25855631].keys()
truth = {f[0]: f[1] for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}


def call(read, position):
    ref, query = read.reference_start + 1, 0
    for op, length in read.cigartuples:
        if op in (0, 7, 8) and ref <= position < ref + length:
            offset = query + position - ref
            return read.query_sequence[offset], read.query_qualities[offset]
        if op == 2 and ref <= position < ref + length:
            return 'DEL', 0
        if op in (0, 7, 8, 2, 3):
            ref += length
        if op in (0, 7, 8, 1, 4):
            query += length
    return 'UNKNOWN', 0


records = {pos: [] for pos in observations}
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    seen = set()
    for read in bam.fetch('CHM13#0#chr20', 25853216, 25855632):
        if read.flag & (4 | 256 | 2048 | 1024 | 512) or read.mapping_quality in (0, 255):
            continue
        if read.query_name in seen:
            continue
        seen.add(read.query_name)
        for pos, cohort in observations.items():
            base, quality = call(read, pos)
            ref, alt = ('G', 'A') if pos == 25853217 else ('T', 'C')
            assert not (base == alt and 30 <= quality < 255 and read.mapping_quality >= 30), read.query_name
            if read.query_name not in cohort:
                continue
            assert not (base == alt and 30 <= quality < 255), read.query_name
            records[pos].append({'qname': read.query_name, 'graph_allele': cohort[read.query_name],
                                 'physical_base': base, 'base_quality': quality,
                                 'mapq': read.mapping_quality,
                                 'adjacent_calls': [call(read, pos-1), call(read, pos+1)],
                                 'diagnostic_parent': truth.get(read.query_name)})

certificate = {}
for pos, rows in records.items():
    ref = 'G' if pos == 25853217 else 'T'
    reference_class = [r for r in rows if r['graph_allele'] == 0 and r['physical_base'] == ref and
                       30 <= r['base_quality'] < 255 and r['mapq'] >= 30]
    masked_alt = [r for r in rows if r['graph_allele'] == 1 and r['physical_base'] == 'DEL' and r['mapq'] >= 30 and
                  all(base == 'T' and 30 <= quality < 255 for base, quality in r['adjacent_calls'])]
    alternate_refs = [r for r in rows if r['graph_allele'] == 1 and r['physical_base'] == ref and
                      30 <= r['base_quality'] < 255 and
                      10**(-r['mapq']/10) + 10**(-r['base_quality']/10) < 0.5]
    probability = 1.0
    for row in alternate_refs:
        probability *= 10**(-row['mapq']/10) + 10**(-row['base_quality']/10)
    certificate[pos] = {'graph_molecules': len(observations[pos]),
                        'graph_allele_counts': dict(Counter(observations[pos].values())),
                        'reference_class_verified': len(reference_class),
                        'alternate_class_deleted': len(masked_alt),
                        'alternate_class_physical_ref': len(alternate_refs),
                        'wrong_alternate_bound': probability,
                        'diagnostic_parent_counts': dict(Counter(r['diagnostic_parent'] for r in rows)),
                        'reads': rows}
assert certificate[25855631]['reference_class_verified'] >= 1
assert certificate[25855631]['alternate_class_deleted'] >= 2
assert certificate[25853217]['reference_class_verified'] >= 1
assert certificate[25853217]['alternate_class_physical_ref'] >= 2
assert certificate[25853217]['wrong_alternate_bound'] <= 0.001
(OUT/'physical-certificate.json').write_text(json.dumps(certificate, indent=2)+'\n')
print(json.dumps({pos: {k: v for k, v in data.items() if k != 'reads'} for pos, data in certificate.items()}, indent=2))
