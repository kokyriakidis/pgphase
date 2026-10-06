#!/usr/bin/env python3
"""Reconstruct the two certificates from original BAM sequence and quality, without parental truth."""
import json
import math
from pathlib import Path
import pysam
ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
CHROM = 'CHM13#0#chr20'
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
def base(read, pos):
    rp, qp = read.reference_start+1, 0
    for op, n in read.cigartuples:
        if op in (0,7,8) and rp <= pos < rp+n:
            qi = qp+pos-rp
            return read.query_sequence[qi], read.query_qualities[qi], qi
        if op in (0,1,4,7,8): qp += n
        if op in (0,2,3,7,8): rp += n
    return '', 255, -1
def repeat(read, beg, end, pos, lengths, deletion=False):
    _, bq, qb = base(read, beg)
    _, eq, qe = base(read, end)
    if qb < 0 or qe <= qb or min(bq,eq) < 20 or max(bq,eq) == 255: return None
    rp, qp, net, error = read.reference_start+1, 0, 0, 10**(-bq/10)+10**(-eq/10)
    for op, n in read.cigartuples:
        if op == 1 and beg <= rp <= end:
            net += n
            if not deletion:
                qualities = read.query_qualities[qp:qp+n]
                if min(qualities) < 10 or max(qualities) == 255: return None
                error += sum(10**(-q/10) for q in qualities)
        if op == 2 and rp < end and rp+n > beg: net -= min(rp+n,end)-max(rp,beg)
        if op == 3 and rp < end and rp+n > beg: return None
        if op in (0,1,4,7,8): qp += n
        if op in (0,2,3,7,8): rp += n
    reference = ref.fetch(CHROM,beg-1,end).upper()
    expected = [reference[:pos-beg]+reference[pos-beg+n:] if deletion else
                reference[:pos-beg]+'TG'*(n//2)+reference[pos-beg:] for n in lengths]
    distances = [distance(read.query_sequence[qb:qe+1], sequence) for sequence in expected]
    if distances[0] == distances[1]: return None
    allele = int(distances[1] < distances[0])
    observed = -net if deletion else net
    nearest = [abs(observed-n) for n in lengths]
    if nearest[0] == nearest[1] or allele != int(nearest[1] < nearest[0]): return None
    if deletion and nearest[allele] > 2: return None
    if not deletion and net < 0 and nearest[allele] > 4: return None
    if distances[allele] > nearest[allele]+1: return None
    return allele, error, observed, distances
reports = {}
for name, params in {'deletion':(15015554,15039543,15019240,15019317,15019295,[0,6]),
                     'insertion':(15087221,15115387,15101247,15101332,15101263,[2,8])}.items():
    left, right, beg, end, pos, lengths = params
    gauge, parity, records = [[0,0],[0,0]], [0,0], []
    log_odds = 0.0
    for read in bam.fetch(CHROM,left-1,right):
        if read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.is_unmapped or not 30 <= read.mapping_quality < 255: continue
        call = repeat(read,beg,end,pos,lengths,name=='deletion')
        if call is None: continue
        allele, error, net, distances = call
        error += 2*10**(-read.mapping_quality/10)
        calibration = read.reference_end < right if name=='deletion' else read.reference_start+1 > left
        snp = left if (name=='deletion' and calibration) or (name=='insertion' and not calibration) else right
        sequence, quality, _ = base(read,snp)
        error += 10**(-quality/10)
        if quality < 20 or quality == 255 or error > .01: continue
        # Baseline source gauges: deletion HP1=REF (upstream C); insertion HP1=short (downstream G).
        if name=='deletion': hap = 0 if sequence == ('C' if calibration else 'A') else 1 if sequence == ('G' if calibration else 'G') else -1
        else:
            hap1_base, hap2_base = ('G','A') if calibration else ('A','G')
            hap = 0 if sequence == hap1_base else 1 if sequence == hap2_base else -1
        if hap < 0: continue
        marker_hap = allele
        if calibration: gauge[marker_hap][hap] += 1
        else:
            reverse = marker_hap != hap
            parity[reverse] += 1
            log_odds += (1 if reverse else -1)*math.log((1-error)/error)
        records.append({'qname':read.query_name,'role':'calibration' if calibration else 'bridge',
                        'allele':allele,'snp_hap':hap,'net_length':net,'sequence_distances':distances,'error':error})
    reports[name] = {'gauge':gauge,'parity':parity,'log_odds':log_odds,'molecules':records}
assert reports['deletion']['gauge'] == [[22,0],[0,27]], reports['deletion']
assert reports['deletion']['parity'] == [2,0]
assert reports['insertion']['gauge'] == [[8,0],[0,6]], reports['insertion']
assert reports['insertion']['parity'] == [0,12]
# Independently verify the retained source I1 gauge in the joined owner.
def source_insertion(read, pos):
    symbol = ref.fetch(CHROM,pos-1,pos).upper()
    beg, end = pos, pos
    while ref.fetch(CHROM,beg-2,beg-1).upper() == symbol: beg -= 1
    while ref.fetch(CHROM,end-1,end).upper() == symbol: end += 1
    assert end-beg < 64
    rp, qp, net, events, eq = read.reference_start+1, 0, 0, 0, 255
    flanks, pure = [0,0], True
    for op, length in read.cigartuples:
        if op == 1 and beg <= rp <= end:
            net += length; events += 1
            pure &= all(x == symbol for x in read.query_sequence[qp:qp+length])
            eq = min(eq,*read.query_qualities[qp:qp+length])
        if op == 2 and rp < end and rp+length > beg:
            if rp < beg or rp+length > end or qp == 0 or qp >= read.query_length: return None
            net -= length; events += 1
            eq = min(eq,read.query_qualities[qp-1],read.query_qualities[qp])
        if op == 3 and rp < end+16 and rp+length > beg-16: return None
        if op in (0,7,8):
            for at in range(max(rp,beg),min(rp+length,end)):
                pure &= read.query_sequence[qp+at-rp] == symbol
            for side,(left,right) in enumerate(((beg-16,beg),(end,end+16))):
                for at in range(max(rp,left),min(rp+length,right)):
                    qi = qp+at-rp
                    if read.query_qualities[qi] != 255 and read.query_sequence[qi] == ref.fetch(CHROM,at-1,at).upper():
                        flanks[side] = max(flanks[side],read.query_qualities[qi])
        if op in (0,1,4,7,8): qp += length
        if op in (0,2,3,7,8): rp += length
    allele = int(net > 0)
    if not pure or events > 1 or abs(net-allele) > 1 or min(flanks) < 20: return None
    quality = min(flanks) if events == 0 else eq
    error = sum(10**(-q/10) for q in flanks)+10**(-quality/10)+2*10**(-read.mapping_quality/10)
    if quality < 20 or quality == 255 or error > .01: return None
    return allele,quality,error,net
cohorts = [[[0,0],[0,0]],[[0,0],[0,0]]]
records = []
for read in bam.fetch(CHROM,15109300,15109301):
    if read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.is_unmapped or not 30 <= read.mapping_quality < 255: continue
    call = source_insertion(read,15109301)
    if call is None: continue
    allele,quality,error,net = call
    symbol,sq,_ = base(read,15115387)
    error += 10**(-sq/10)
    if quality < 30 or sq < 30 or sq == 255 or symbol not in ('G','A') or error > .01: continue
    # In the joined owner hap1 is paternal: source ALT and downstream SNP ALT A.
    hap = int(symbol != 'A')
    hashed = 14695981039346656037
    for byte in read.query_name.encode(): hashed = ((hashed ^ byte)*1099511628211) & ((1<<64)-1)
    fold = hashed & 1
    cohorts[fold][allele][hap] += 1
    records.append({'qname':read.query_name,'cohort':fold,'allele':allele,'snp_hap':hap,'net_length':net,'error':error})
assert cohorts == [[[1,8],[8,0]],[[0,8],[8,0]]],cohorts
n,wrong,z = len(records),1,1.6448536269514722
rate = wrong/n
bound = (rate+z*z/(2*n)+z*math.sqrt(rate*(1-rate)/n+z*z/(4*n*n)))/(1+z*z/n)+.01
assert bound <= .20
reports['retained_source_insertion'] = {'position':15109301,'cohorts':cohorts,'hap1_allele':1,
    'discordant':wrong,'total':n,'joint_error_bound':bound,'molecules':records}
(OUT/'physical-certificates.json').write_text(json.dumps(reports,indent=2)+'\n')
print(json.dumps({name:{k:v for k,v in report.items() if k!='molecules'} for name,report in reports.items()},indent=2))
