#!/usr/bin/env python3
"""Reproduce the physical repeat certificate without consulting parental truth."""
from collections import Counter
import json
import math
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
POS = 24121714
RIGHT = 24142287
LEFT = [(24103779, 'C', 'T'), (24105188, 'A', 'G')]
ref = pysam.FastaFile(str(ROOT/'test_data/chm13v2.0.chr20.renamed.fa'))
bam = pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'))


def distance(a, b):
    row = list(range(len(b)+1))
    for i, x in enumerate(a, 1):
        nxt = [i]
        for j, y in enumerate(b, 1):
            nxt.append(min(nxt[-1]+1, row[j]+1, row[j-1]+(x != y)))
        row = nxt
    return row[-1]


def bases(read):
    return {p+1: (read.query_sequence[q], read.query_qualities[q], q)
            for q, p in read.get_aligned_pairs() if q is not None and p is not None}


def snp(call, position, reference, alternate):
    base, quality, _ = call.get(position, ('N', 0, -1))
    return (0 if base == reference else 1 if base == alternate else -1), quality


def eligible(read):
    return not read.flag & (4 | 256 | 512 | 1024 | 2048) and 30 <= read.mapping_quality < 255


def repeat(read, call):
    beg, end = POS-16, POS+16
    if beg not in call or end not in call:
        return None
    qb, qe = call[beg][2], call[end][2]+1
    if min(call[beg][1], call[end][1]) < 20 or max(call[beg][1], call[end][1]) == 255:
        return None
    rp, qp, observed, error = read.reference_start+1, 0, 0, 0.0
    for op, length in read.cigartuples:
        if op == 1 and beg <= rp <= end:
            observed += length
            qualities = read.query_qualities[qp:qp+length]
            if min(qualities) < 10 or max(qualities) == 255:
                return None
            error += sum(10**(-q/10) for q in qualities)
        if op == 2 and rp < end and rp+length > beg:
            observed -= min(rp+length, end)-max(rp, beg)
        if op == 3 and rp < end and rp+length > beg:
            return None
        if op in (0, 1, 4, 7, 8): qp += length
        if op in (0, 2, 3, 7, 8): rp += length
    if not 0 < observed <= 64:
        return None
    reference = ref.fetch('CHM13#0#chr20', beg-1, end).upper()
    query = read.query_sequence[qb:qe]
    scores = [distance(query, reference[:16]+'T'*n+reference[16:]) for n in (4, 8)]
    if scores[0] == scores[1]:
        return None
    allele = int(scores[1] < scores[0])
    length_distances = [abs(observed-n) for n in (4, 8)]
    if length_distances[0] != length_distances[1] and allele != int(length_distances[1] < length_distances[0]):
        return None
    if scores[allele] > length_distances[allele]+1:
        return None
    error += 10**(-call[beg][1]/10)+10**(-call[end][1]/10)+2*10**(-read.mapping_quality/10)
    return 1-allele, observed, error, scores


gauge = [[0, 0], [0, 0]]
parity = [[0, 0], [0, 0]]
errors = [0.0, 0.0]
calibration, bridge, seen = [], [], set()
for read in bam.fetch('CHM13#0#chr20', POS-2, RIGHT):
    if not eligible(read) or read.query_name in seen:
        continue
    seen.add(read.query_name)
    call = bases(read)
    observed = repeat(read, call)
    if observed is None:
        continue
    hap, length, error, scores = observed
    rc, rq = snp(call, RIGHT, 'A', 'G')
    record = {'qname': read.query_name, 'length': length, 'class_hap_index': hap, 'scores': scores}
    if rc >= 0 and 20 <= rq < 255:
        error += 10**(-rq/10)
        if error < 0.5:
            flip = hap != rc
            parity[hap][int(flip)] += 1
            errors[hap] = max(errors[hap], error)
            bridge.append({**record, 'right_allele': rc, 'flip': flip, 'error': error})
        continue
    if read.reference_end >= RIGHT:
        continue
    votes = []
    for position, reference, alternate in LEFT:
        sc, sq = snp(call, position, reference, alternate)
        if sc >= 0 and 20 <= sq < 255:
            votes.append(sc)
            error += 10**(-sq/10)
    if votes and error <= 0.01:
        snp_hap = votes[-1] if len(set(votes)) == 1 else 1-hap
        gauge[hap][snp_hap] += 1
        calibration.append({**record, 'snp_hap_index': snp_hap, 'error': error})
assert not ({r['qname'] for r in calibration} & {r['qname'] for r in bridge})


def wilson(discordant, total):
    z = 1.6448536269514722
    rate = discordant/total
    return (rate+z*z/(2*total)+z*math.sqrt(rate*(1-rate)/total+z*z/(4*total*total)))/(1+z*z/total)


bounds = [wilson(gauge[h][1-h], sum(gauge[h]))+0.01+errors[h] for h in (0, 1)]
assert gauge == [[4, 0], [0, 5]], gauge
assert parity == [[0, 1], [0, 1]], parity
assert math.prod(bounds) <= 0.20
same = gauge[0][0]+gauge[1][1]
cross = gauge[0][1]+gauge[1][0]
association = 2*sum(math.comb(same+cross, i) for i in range(same, same+cross+1))/2**(same+cross)
assert association <= 0.01
edge, edge_seen = Counter(), set()
for read in bam.fetch('CHM13#0#chr20', LEFT[0][0]-1, LEFT[1][0]):
    if not eligible(read) or read.query_name in edge_seen:
        continue
    edge_seen.add(read.query_name)
    call = bases(read)
    calls = [snp(call, *site) for site in LEFT]
    if all(c >= 0 and 30 <= q < 255 for c, q in calls):
        edge['agree' if calls[0][0] == calls[1][0] else 'contrary'] += 1
assert edge['agree'] >= 2 and edge['contrary'] == 0
report = {'gauge': gauge, 'parity_by_class': parity, 'class_max_call_errors': errors,
          'class_augmented_wilson_bounds': bounds, 'joint_wrong_parity_bound': math.prod(bounds),
          'diploid_association_p': association, 'recovered_snp_physical_edge': dict(edge),
          'calibration_reads': calibration, 'bridge_reads': bridge, 'parental_truth_used': False}
(OUT/'physical-certificate.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
