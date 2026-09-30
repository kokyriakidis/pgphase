#!/usr/bin/env python3
"""Audit exact-key HiPhase joins of pgphase chromosome phase-set seams."""
from collections import defaultdict
from pathlib import Path
import argparse
import bisect
import gzip
import pysam

ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pg-vcf', type=Path,
                    default=Path('/tmp/pgphase-repeat-pair-guarded-full/phased.vcf'))
parser.add_argument('--panel', type=Path,
                    default=Path(__file__).with_name('panel_before.tsv'))
parser.add_argument('--out', type=Path,
                    default=Path(__file__).with_name('audit.tsv'))
args=parser.parse_args()
PG=args.pg_vcf
HI=Path('/tmp/hiphase-on-pgphase-final-chr20/phased.vcf')
HINPUT=Path('/tmp/hiphase-on-pgphase-final-chr20/input.vcf.gz')
BAM=Path('/tmp/hiphase-on-pgphase-final-chr20/phased.bam')
TRUTH=ROOT/'test_data/derived/chr20_truth_hap.tsv'
PANEL=args.panel
OUT=args.out

def open_vcf(path):
    return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)
def read_vcf(path):
    rows=[]; blocks=defaultdict(list); phased={}
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith('#'): continue
            f=line.rstrip().split('\t'); key=(int(f[1]),f[3],f[4])
            rows.append(key)
            fmt=dict(zip(f[8].split(':'),f[9].split(':')))
            gt=fmt.get('GT','');ps=fmt.get('PS','.')
            if '|' in gt and len(set(gt.split('|')))==2 and ps not in ('.','0',''):
                blocks[ps].append(key);phased[key]=ps
    return sorted(rows),blocks,phased

pgrows,pgblocks,pgkey=read_vcf(PG)
hirows,hiblocks,hikey=read_vcf(HI)
inputrows,_,_=read_vcf(HINPUT)
truth={}
with TRUTH.open() as fh:
    for line in fh:
        f=line.rstrip().split('\t')
        if len(f)==2 and f[1] in ('MATERNAL','PATERNAL'):
            truth[f[0]]=f[1][0]
panel=set()
for line in PANEL.open():
    if line and line[0].isdigit():
        f=line.split('\t');panel.add((int(f[0]),int(f[1])))
spans=[]
for ps,keys in pgblocks.items():
    low=min(keys,key=lambda k:k[0]);high=max(keys,key=lambda k:k[0])
    spans.append((low[0],high[0],ps,low,high))
spans.sort()
head=['left','right','gap_bp','left_ref','left_alt','right_ref','right_alt','hiphase_ps','local_reads','local_correct','left_reads','left_correct','right_reads','right_correct','same_orientation','same_calls','panel']
lines=[]
with pysam.AlignmentFile(BAM) as bam:
    for a,b in zip(spans,spans[1:]):
        if a[1]>=b[0] or (a[1]>=26000000 and b[0]<=29500000):
            continue
        left=a[4];right=b[3];hi_ps=hikey.get(left)
        if not hi_ps or hi_ps!=hikey.get(right):continue
        low=left[0]-50000; high=right[0]+50000
        def keys_between(rows):
            i=bisect.bisect_left(rows,(low,'',''));j=bisect.bisect_right(rows,(high,'~','~'))
            return rows[i:j]
        same_calls=keys_between(pgrows)==keys_between(inputrows)
        votes=[[0,0] for _ in range(3)]
        seen=set()
        for read in bam.fetch('CHM13#0#chr20',max(0,left[0]-10000),right[0]+10000):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.query_name in seen:
                continue
            seen.add(read.query_name)
            if not read.has_tag('HP') or not read.has_tag('PS') or str(read.get_tag('PS'))!=hi_ps:
                continue
            parent=truth.get(read.query_name)
            if parent is None:continue
            agree=(read.get_tag('HP')==1)==(parent=='M')
            ix=0 if agree else 1
            beg=read.reference_start+1;end=read.reference_end
            if end>=left[0]-10000 and beg<=right[0]+10000:votes[0][ix]+=1
            if end<=left[0] and end>=left[0]-10000:votes[1][ix]+=1
            if beg>=right[0] and beg<=right[0]+10000:votes[2][ix]+=1
        n=[sum(v) for v in votes];correct=[max(v) for v in votes]
        orientation=(votes[1][0]>votes[1][1])==(votes[2][0]>votes[2][1]) if n[1] and n[2] else False
        lines.append([left[0],right[0],right[0]-left[0],left[1],left[2],right[1],right[2],hi_ps,n[0],correct[0],n[1],correct[1],n[2],correct[2],int(orientation),int(same_calls),int((left[0],right[0]) in panel)])
with OUT.open('w') as fh:
    fh.write('\t'.join(head)+'\n')
    for row in lines:fh.write('\t'.join(map(str,row))+'\n')
qualified=[r for r in lines if r[8]>=20 and r[9]/r[8]>=.98 and r[10]>=5 and r[11]/r[10]>=.95 and r[12]>=5 and r[13]/r[12]>=.95 and r[14] and r[15] and not r[16]]
print('exact joins',len(lines),'same calls',sum(r[15] for r in lines),'qualified new',len(qualified),'output',OUT)
for r in qualified:
    print(r[0],r[1],r[2],r[8:17])
