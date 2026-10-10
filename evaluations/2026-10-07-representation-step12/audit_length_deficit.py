#!/usr/bin/env python3
"""Check why the 65 Mb repeat context cannot be repaired by parent alignment alone."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import pysam

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
FOLDER=ROOT/'test_data/tmp_representation_step12/accepted-final/matrix65-final'

def rows(kind):
    with (FOLDER/f'matrix.chunk0.{kind}.tsv').open() as stream:
        return [r for r in csv.DictReader(stream,delimiter='\t') if r['physical']=='3']
context=rows('composition-contexts')[0]
beg,end=int(context['beg']),int(context['end'])
cohort={r['read']:r for r in rows('physical-cohort')}
maximum=max(len(r['sequence']) for r in rows('path-alleles'))
residual_path=ROOT/'test_data/tmp_representation_step11/cached/matrix65-final.residuals.tsv'
previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step11/manifest.json').read_text())
record=next(r for r in previous['files'] if r['path']==str(residual_path.relative_to(ROOT)))
assert hashlib.sha256(residual_path.read_bytes()).hexdigest()==record['sha256']
with residual_path.open() as stream:
    residuals={r['read']:int(r['path_min']) for r in csv.DictReader(stream,delimiter='\t') if r['physical']=='3'}
assert residuals.keys()==cohort.keys()
verified=[]
original={}
with pysam.AlignmentFile(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as bam:
    for read in bam.fetch('CHM13#0#chr20',beg-1,end):
        if read.is_secondary or read.is_supplementary or read.query_name not in cohort:continue
        assert read.query_name not in original
        aligned={r:q for q,r in read.get_aligned_pairs(matches_only=True) if beg-1<=r<end}
        lo,hi=aligned[beg-1],aligned[end-1]+1
        row=cohort[read.query_name]
        assert lo==int(row['query_beg']) and hi==int(row['query_end'])
        assert read.query_sequence[lo:hi]==row['sequence']
        inserted=deleted=0
        cursor=read.reference_start
        events=[]
        for op,length in read.cigartuples:
            if op==1 and beg<=cursor<=end-1:
                inserted+=length;events.append(dict(boundary=cursor,length=length))
            if op==2:
                deleted+=max(0,min(cursor+length,end)-max(cursor,beg-1))
            if op in (0,2,3,7,8):cursor+=length
        assert len(row['sequence'])==end-beg+1+inserted-deleted
        original[read.query_name]=(read.reference_start,read.reference_end,read.cigarstring,hashlib.sha256(read.query_sequence.encode()).hexdigest())
        bound=max(0,len(row['sequence'])-maximum)
        assert residuals[read.query_name]>=bound
        verified.append(dict(read=read.query_name,query_length=len(row['sequence']),internal_inserted_bases=inserted,
                             internal_deleted_bases=deleted,insertions=events,length_lower_bound=bound,path_residual=residuals[read.query_name]))
assert original.keys()==cohort.keys()
seen=set()
with pysam.AlignmentFile(ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.bam') as bam:
    for read in bam:
        if read.is_secondary or read.is_supplementary or read.query_name not in original:continue
        assert read.query_name not in seen
        assert original[read.query_name]==(read.reference_start,read.reference_end,read.cigarstring,hashlib.sha256(read.query_sequence.encode()).hexdigest())
        seen.add(read.query_name)
assert seen==original.keys()
result=dict(physical=3,beg=beg,end=end,molecules=len(cohort),reference_length=end-beg+1,
            maximum_parent_length=max(len(r['sequence']) for r in rows('physical-alleles')),maximum_path_length=maximum,
            query_lengths=dict(sorted(Counter(len(r['sequence']) for r in cohort.values()).items())),
            queries_longer_than_catalog=sum(v['length_lower_bound']>0 for v in verified),
            unavoidable_length_residual=sum(v['length_lower_bound'] for v in verified),path_residual=sum(residuals.values()),
            original_diploid_cost=int(rows('physical-genotypes')[0]['cost']),same_hiphase_input_alignments=True,
            raw_query_slices_and_cigar_length_balance_verified=True,reads=verified)
(OUT/'length-deficit-checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='reads'},indent=2))
