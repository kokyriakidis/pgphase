#!/usr/bin/env python3
"""Reconstruct the tandem-repeat calibration from original CIGAR and bases."""
import json
import math
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
CHROM = 'CHM13#0#chr20'
POS, FLANK = 882278, 16
reference = pysam.FastaFile(str(ROOT/'test_data/chm13v2.0.chr20.renamed.fa'))


def edit_distance(a, b):
    previous = list(range(len(b)+1))
    for i, x in enumerate(a, 1):
        current = [i]
        for j, y in enumerate(b, 1):
            current.append(min(previous[j]+1, current[j-1]+1, previous[j-1]+(x != y)))
        previous = current
    return previous[-1]


def physical(read):
    beg, end = POS-FLANK, POS+FLANK
    rp, qp, qb, qe, net, error = read.reference_start+1, 0, None, None, 0, 0.0
    events = []
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            if rp <= beg < rp+length: qb = qp+beg-rp
            if rp <= end < rp+length: qe = qp+end-rp+1
        if op == 1 and beg <= rp <= end:
            qualities = read.query_qualities[qp:qp+length]
            if any(q < 10 or q == 255 for q in qualities): return
            net += length
            error += sum(10**(-q/10) for q in qualities)
            events.append([rp, read.query_sequence[qp:qp+length]])
        if op == 2 and rp < end and rp+length > beg:
            net -= min(rp+length, end)-max(rp, beg)
        if op == 3 and rp < end and rp+length > beg: return
        if op in (0, 1, 4, 7, 8): qp += length
        if op in (0, 2, 3, 7, 8): rp += length
    if net < 0 or net > 50 or qb is None or qe is None: return
    qualities = [read.query_qualities[qb], read.query_qualities[qe-1]]
    if any(q < 20 or q == 255 for q in qualities): return
    ref = reference.fetch(CHROM, beg-1, end).upper()
    query = read.query_sequence[qb:qe].upper()
    distances = [edit_distance(query, ref[:FLANK]+alt+ref[FLANK:]) for alt in ('TC', 'TCTC')]
    if distances[0] == distances[1]: return
    allele = int(distances[1] < distances[0])
    if net != 3 and allele != int(net > 3): return
    if distances[allele] > abs(net-(2 if allele == 0 else 4))+1: return
    error += sum(10**(-q/10) for q in qualities)
    return {'class': allele, 'net_length': net, 'distances': distances,
            'events': events, 'call_error': error}


calibration, parity, odds = [[0, 0], [0, 0]], [0, 0], 0.0
accepted = []
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    seen = set()
    for read in bam.fetch(CHROM, 863405, 890263):
        if read.flag & (4|256|2048|1024|512) or read.mapping_quality < 30 or read.mapping_quality == 255:
            continue
        if read.query_name in seen: continue
        seen.add(read.query_name)
        call = physical(read)
        if call is None: continue
        bases = {rp+1: (read.query_sequence[qp].upper(), read.query_qualities[qp])
                 for qp, rp in read.get_aligned_pairs(matches_only=True)}
        error = call['call_error'] + 2*10**(-read.mapping_quality/10)
        repeat_hap = 1-call['class']
        if read.reference_start+1 > 863406:
            observed = [bases.get(p, ('N', 0)) for p in (890261, 890262)]
            sequence = ''.join(b for b, q in observed)
            quality = min(q for b, q in observed)
            if sequence not in ('TT', 'AC') or quality < 20 or any(q == 255 for b, q in observed): continue
            error += 2*10**(-quality/10)
            if error > .01: continue
            gauge_hap = int(sequence == 'TT')
            calibration[repeat_hap][gauge_hap] += 1
            kind = 'calibration'
        else:
            base, quality = bases.get(863406, ('N', 0))
            if base not in 'CA' or quality < 20 or quality == 255: continue
            error += 10**(-quality/10)
            if error > .01: continue
            gauge_hap = int(base == 'A')
            flipped = repeat_hap != gauge_hap
            parity[flipped] += 1
            odds += (1 if flipped else -1)*math.log((1-error)/error)
            kind = 'bridge'
        accepted.append({'read': read.query_name, 'kind': kind, 'repeat_haplotype_index': repeat_hap,
                         'snp_haplotype_index': gauge_hap, 'error_bound': error, **call})
assert calibration == [[10, 0], [1, 12]], calibration
assert parity == [0, 1], parity
assert abs(odds-7.23801) < .00001, odds
n, bad, z = sum(map(sum, calibration)), 1, 1.6448536269514722
p = bad/n
wilson = (p+z*z/(2*n)+z*math.sqrt(p*(1-p)/n+z*z/(4*n*n)))/(1+z*z/n)
result = {'physical_insertion_coordinate': POS, 'vcf_coordinate': POS-1,
          'alleles': ['TC', 'TCTC'], 'calibration': calibration, 'bridge_parity_counts': parity,
          'bridge_log_odds': odds, 'parity_error': 1/(1+math.exp(abs(odds))),
          'joint_error_bound': wilson+.01+1/(1+math.exp(abs(odds))),
          'calibration_molecules': n, 'accepted': accepted}
assert result['joint_error_bound'] <= .2
(OUT/'physical-evidence.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k != 'accepted'}, indent=2))
