#!/usr/bin/env python3
"""Measure full original-BAM coverage omitted by contrast-observed saved cohorts."""
import csv
import json
from pathlib import Path

import pysam

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
results=[]
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
 for name in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
  folder=ROOT/'test_data/tmp_representation_step6/accepted'/name
  contexts=list(csv.DictReader((folder/'matrix.chunk0.joint-contexts.tsv').open(),delimiter='\t'))
  scored={}
  for row in csv.DictReader((folder/'matrix.chunk0.joint-sequences.tsv').open(),delimiter='\t'):
   if row['status']=='scored':scored.setdefault(row['locus'],set()).add(row['read'])
  loci=[]
  for c in contexts:
   if c['status']!='valid':continue
   lo,hi=int(c['beg'])-1,int(c['end'])-1
   covered={}
   for read in bam.fetch('CHM13#0#chr20',lo,hi+1):
    # Same primary/flag/recovery-MAPQ floor as the owning whole-BAM solve.
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.mapping_quality<1:
     continue
    if read.reference_start>lo or read.reference_end<=hi:continue
    ref=read.reference_start
    skipped=False
    for op,length in read.cigartuples:
     if op==3 and ref<=hi and ref+length>lo:skipped=True;break
     if op in (0,2,3,7,8):ref+=length
    if skipped:continue
    anchors={r:q for q,r in read.get_aligned_pairs() if r in (lo,hi) and q is not None}
    if set(anchors)!={lo,hi}:continue
    start,end=anchors[lo],anchors[hi]+1
    sequence=read.query_sequence[start:end].upper()
    if not sequence or len(sequence)>4096 or set(sequence)-set('ACGT'):continue
    if read.query_name in covered:
     covered[read.query_name]=False
    else:covered[read.query_name]=True
   eligible={n for n,unique in covered.items() if unique}
   old=scored.get(c['locus'],set())
   assert old<=eligible,(name,c['locus'],old-eligible)
   loci.append(dict(locus=int(c['locus']),contrast_observed_molecules=len(old),fully_covered_primary_molecules=len(eligible),omitted_molecules=len(eligible-old),omitted_names=sorted(eligible-old)))
  report=dict(owner=name,contexts=len(loci),contexts_with_omitted_molecules=sum(r['omitted_molecules']>0 for r in loci),omitted_molecule_locus_records=sum(r['omitted_molecules'] for r in loci),loci=loci)
  results.append(report)
  print(json.dumps({k:v for k,v in report.items() if k!='loci'}),flush=True)
(OUT/'cohort-coverage-checks.json').write_text(json.dumps(dict(interpretation='Current diagnostic fits condition on the original contrast-observed cohort; no complete sample genotype or phasing certificate claimed',owners=results),indent=2)+'\n')
