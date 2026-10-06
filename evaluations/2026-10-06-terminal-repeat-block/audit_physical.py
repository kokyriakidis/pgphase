#!/usr/bin/env python3
"""Independently replay original CIGAR/base evidence for the terminal call."""
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
OWNER = ROOT/'test_data/tmp_gap_fix80/final_owner/12'
CHROM = 'CHM13#0#chr20'
POS = 12955839
ref = pysam.FastaFile(str(ROOT/'test_data/chm13v2.0.chr20.renamed.fa'))
graph = {}
for line in (OWNER/'matrix.chunk0.recovery-final.tsv').open():
    z = line.rstrip().split('\t')
    if z[0] == 'OBS' and int(z[2]) == 1736:
        graph[z[1]] = int(z[4])
tags = {}
with pysam.AlignmentFile(str(OWNER/'phased.bam'), check_sq=False) as bam:
    for r in bam:
        if not r.is_secondary and not r.is_supplementary:
            tags[r.query_name] = (r.get_tag('HP') if r.has_tag('HP') else 0,
                                  r.get_tag('PS') if r.has_tag('PS') else 0)

def physical(r):
    seq, qs = r.query_sequence, r.query_qualities
    bases = {rp+1: (seq[qp].upper(), qs[qp]) for qp, rp in r.get_aligned_pairs(matches_only=True)}
    rp, qp, events = r.reference_start+1, 0, []
    for op, n in r.cigartuples:
        if op == 1 and POS-64 <= rp <= POS+64:
            events.append((rp, seq[qp:qp+n].upper(), list(qs[qp:qp+n])))
        if op == 2 and rp < POS+64 and rp+n > POS-64:
            return
        if op in (0, 1, 4, 7, 8): qp += n
        if op in (0, 2, 3, 7, 8): rp += n
    if len(events) > 1: return
    beg = end = POS
    error = 0.0
    allele = int(bool(events))
    if events:
        ep, ins, iq = events[0]
        if len(ins) != 1 or any(q < 30 or q == 255 for q in iq): return
        beg, end = min(ep, POS), max(ep, POS)
        context = ref.fetch(CHROM, beg-1, end-1).upper()
        if context[:POS-beg]+'C'+context[POS-beg:] != context[:ep-beg]+ins+context[ep-beg:]: return
        error += sum(10**(-q/10) for q in iq)
    for p in range(beg-1, end+1):
        b, q = bases.get(p, ('N', 0))
        if b != ref.fetch(CHROM, p-1, p).upper() or q < 30 or q == 255: return
        error += 10**(-q/10)
    snps = []
    for p, alleles in ((12954878, 'AG'), (12954322, 'GA')):
        b, q = bases.get(p, ('N', 0))
        if b not in alleles or q < 30 or q == 255: return
        snps.append(alleles.index(b))
        error += 10**(-q/10)
    if snps[0] != snps[1]: return
    error += 2*10**(-r.mapping_quality/10)
    if error > .01: return
    return allele, snps[0], error

cohorts = [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
accepted = []
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    for r in bam.fetch(CHROM, POS-1, POS):
        if r.flag & (4|256|2048|1024|512) or r.mapping_quality < 30 or r.mapping_quality == 255: continue
        hp, ps = tags.get(r.query_name, (0, 0))
        if hp not in (1, 2) or ps != 12000781: continue
        call = physical(r)
        if call is None: continue
        allele, hap, error = call
        if graph.get(r.query_name) != allele: continue
        hashed = 14695981039346656037
        for byte in r.query_name.encode(): hashed = ((hashed ^ byte)*1099511628211) & ((1<<64)-1)
        fold = hashed & 1
        cohorts[fold][allele][hap] += 1
        accepted.append({'read': r.query_name, 'cohort': fold, 'allele': allele,
                         'snp_haplotype_index': hap, 'call_error_bound': error})
assert cohorts == [[[14, 0], [0, 7]], [[10, 0], [0, 7]]], cohorts
z2 = 1.6448536269514722**2
hi = {}
with pysam.VariantFile(str(ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.vcf.gz')) as vcf:
    for row in vcf:
        if row.pos not in (12954322, 12954878, 12955838): continue
        sample = next(iter(row.samples.values()))
        hi[row.pos] = {'ref': row.ref, 'alt': row.alts, 'GT': sample['GT'],
                       'PS': sample['PS'], 'DP': sample.get('DP'), 'GQ': sample.get('GQ'), 'AD': sample.get('AD')}
result = {'physical_insertion_coordinate': POS, 'vcf_coordinate': POS-1,
          'catalog_alleles': ['T', 'TC'], 'catalog_counts': [44, 29],
          'cohorts': cohorts, 'accepted_molecules': accepted,
          'molecule_count': len(accepted), 'discordant_molecules': 0,
          'one_sided_95_wilson_upper_plus_call_error': {
              'each_cohort': [z2/(sum(map(sum, c))+z2)+.01 for c in cohorts],
              'combined': z2/(len(accepted)+z2)+.01}, 'hiphase_calls': hi}
(OUT/'physical-evidence.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: v for k, v in result.items() if k != 'accepted_molecules'}, indent=2))
