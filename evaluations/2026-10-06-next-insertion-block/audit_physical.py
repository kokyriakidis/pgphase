#!/usr/bin/env python3
"""Reconstruct the mixed insertion/deletion calibration from original CIGAR and bases."""
import json
import math
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
CHROM = 'CHM13#0#chr20'
POS, FLANK = 55309790, 16
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
    if abs(net) > 64 or qb is None or qe is None: return
    qualities = [read.query_qualities[qb], read.query_qualities[qe-1]]
    if any(q < 20 or q == 255 for q in qualities): return
    ref = reference.fetch(CHROM, beg-1, end).upper()
    query = read.query_sequence[qb:qe].upper()
    distances = [edit_distance(query, expected) for expected in (ref[:FLANK]+'T'+ref[FLANK:], ref[:FLANK]+ref[FLANK+2:])]
    if distances[0] == distances[1]: return
    allele = int(distances[1] < distances[0])
    if allele != int(abs(net+2) < abs(net-1)): return
    if distances[allele] > abs(net-(1 if allele == 0 else -2))+1: return
    error += sum(10**(-q/10) for q in qualities)
    return {'class': allele, 'net_length': net, 'distances': distances,
            'events': events, 'call_error': error}


# Read graph channels and core membership, independently of parental truth.
observations={}
for line in (ROOT/'test_data/tmp_gap_fix83/baseline/55/matrix.chunk0.recovery-final.tsv').open():
 f=line.rstrip('\n').split('\t')
 if f[0]=='OBS' and f[2] in ('530','531'):observations.setdefault(f[1],{})[f[2]]=int(f[4])
core=set()
for line in (ROOT/'test_data/tmp_gap_fix83/baseline/55/matrix.chunk0.recovery-final.tsv').open():
 f=line.rstrip('\n').split('\t')
 if f[0]=='READ' and f[5]=='0' and f[6] in ('1','2') and f[7]=='55260093':core.add(f[1])
cohorts=[[[0,0],[0,0]],[[0,0],[0,0]]];rows=[]
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
 for read in bam.fetch(CHROM,55309797,55309798):
  if read.flag&(4|256|2048|1024|512) or read.mapping_quality<30 or read.mapping_quality==255 or read.query_name not in core:continue
  call=physical(read)
  if call is None:continue
  graph=observations.get(read.query_name,{})
  if '530' not in graph or '531' not in graph:continue
  a,b=graph['530'],graph['531']
  if (a==1)==(b==1) or (a!=1 if call['class']==0 else b!=1):continue
  bases={rp+1:(read.query_sequence[qp].upper(),read.query_qualities[qp]) for qp,rp in read.get_aligned_pairs(matches_only=True)}
  a,aq=bases.get(55309475,('N',0));b,bq=bases.get(55308854,('N',0))
  if a not in 'GA' or b not in 'CT' or aq<30 or bq<30 or aq==255 or bq==255 or (a=='A')!=(b=='T'):continue
  err=call['call_error']+10**(-aq/10)+10**(-bq/10)+2*10**(-read.mapping_quality/10)
  if err>.01:continue
  h=14695981039346656037
  for byte in read.query_name.encode():h=((h^byte)*1099511628211)&((1<<64)-1)
  cohorts[h&1][call['class']][int(a=='A')]+=1
  rows.append({'read':read.query_name,'physical':call,'error':err,'cohort':h&1})
assert cohorts==[[[7,0],[0,10]],[[9,0],[0,8]]],cohorts
background=reference.fetch(CHROM,55309789,55309798).upper()
assert 'T'+background[:8]==background[:8]+'T'
assert background[2:8]==background[:6]
def wilson(n):
 z=1.6448536269514722
 return z*z/(n+z*z)
result={'source_physical_position':55309790,'graph_physical_insertion_position':55309798,'graph_physical_deletion_position':55309796,
 'source_alleles':['INS_T','DEL_TT'],'same_complete_edited_sequences':True,'cohorts':cohorts,
 'cohort_error_bounds':[wilson(17)+.01]*2,'combined_error_bound':wilson(34)+.01,'rows':rows}
assert all(x<=.2 for x in result['cohort_error_bounds']) and result['combined_error_bound']<=.1
(OUT/'physical-evidence.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
