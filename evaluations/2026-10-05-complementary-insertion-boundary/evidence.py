#!/usr/bin/env python3
"""Audit the two exact insertion ALT molecules and replay-only candidate state."""
import argparse
import json
import math
from pathlib import Path
import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--replay', type=Path, required=True)
parser.add_argument('--bam', default='test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
rows = []
with pysam.AlignmentFile(args.bam) as bam:
    for read in bam.fetch('CHM13#0#chr20', 8977828, 8998973):
        if (read.is_secondary or read.is_supplementary or read.is_duplicate or
                read.is_qcfail or read.is_unmapped or read.mapping_quality < 30 or
                read.mapping_quality == 255):
            continue
        pairs = {r: q for q, r in read.get_aligned_pairs() if r is not None}
        q = pairs.get(8977828)
        if q is None or read.query_qualities[q] < 30:
            continue
        base, baseq = read.query_sequence[q], read.query_qualities[q]
        if base not in ('T', 'A'):
            continue
        refpos, querypos = read.reference_start, 0
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                refpos += length
                querypos += length
            elif op == 1:
                allele = read.query_sequence[querypos:querypos + length]
                if refpos == 8998972 and allele in ('AC', 'ACAC'):
                    quality = min(read.query_qualities[querypos - 1:querypos + length + 1])
                    if quality >= 30:
                        rows.append(dict(qname=read.query_name, snp=base, snp_quality=baseq,
                                         insertion=allele, insertion_quality=quality,
                                         mapq=read.mapping_quality))
                querypos += length
            elif op in (2, 3):
                refpos += length
            elif op == 4:
                querypos += length
assert len({r['qname'] for r in rows}) == len(rows) == 2
assert {(r['snp'], r['insertion']) for r in rows} == {('T', 'ACAC'), ('A', 'AC')}
log_odds = sum(math.log((1 - p) / p) for r in rows
               for p in [10 ** (-r['snp_quality'] / 10) +
                         10 ** (-r['insertion_quality'] / 10) +
                         2 * 10 ** (-r['mapq'] / 10)])
assert log_odds >= math.log(999)
variants = [line.rstrip().split('\t') for line in
            (args.replay / 'matrix.chunk0.recovery-final.tsv').open()
            if line.startswith('VAR\t')]
insertions = [v for v in variants if v[2] == '8998973' and v[7] in ('AC', 'ACAC')]
assert len(insertions) == 2 and insertions[0][8] == insertions[1][8]
assert {v[9] for v in insertions} == {'0', '1'}
assert all(v[11:13] == ['1', '1'] for v in insertions)
shadow = [v for v in variants if v[2] == '9048832' and v[7] == '>115281574>115281604']
assert len(shadow) == 5 and [int(v[8]) > 0 for v in shadow] == [False, False, True, False, False]
result = dict(physical_molecules=rows, log_odds=log_odds,
              wrong_parity_probability=1 / (1 + math.exp(log_odds)),
              complementary_source_rows=insertions, phased_shadow_rows=shadow)
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
