#!/usr/bin/env python3
"""Orient HiPhase globally and verify identical original alignments for cached molecule states."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pysam

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
names=set()
for name in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
 p=ROOT/'test_data/tmp_representation_step8/accepted'/name/'matrix.chunk0.joint-cohort.tsv'
 names.update(r['read'] for r in csv.DictReader(p.open(),delimiter='\t') if r['status']=='scored')
def alignment(r):
 return (r.reference_name.removeprefix('CHM13#0#'),r.reference_start,r.reference_end,r.cigarstring,
         hashlib.sha256(r.query_sequence.encode()).hexdigest(),
         hashlib.sha256(bytes(r.query_qualities) if r.query_qualities is not None else b'').hexdigest(),
         r.mapping_quality,r.flag)
original={}
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
 for lo,hi in [(4000000,5000000),(5000000,6000000),(65000000,66000000)]:
  for r in bam.fetch('CHM13#0#chr20',lo,hi):
   if not r.is_secondary and not r.is_supplementary and r.query_name in names:
    a=alignment(r)
    assert r.query_name not in original or original[r.query_name]==a
    original[r.query_name]=a
assert set(original)==names
hiphase=ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.bam'
verified=set()
with pysam.AlignmentFile(str(hiphase)) as bam, pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as source:
 assert tuple(n.removeprefix('CHM13#0#') for n in source.references)==bam.references
 assert source.lengths==bam.lengths
 for r in bam:
  if not r.is_secondary and not r.is_supplementary and r.query_name in names:
   assert r.query_name not in verified and alignment(r)==original[r.query_name],r.query_name
   verified.add(r.query_name)
assert verified==names
spec=importlib.util.spec_from_file_location('read_audit',ROOT/'evaluations/2026-10-06-shared-insertion-source/read_audit.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
truth={f[0]:f[1]=='PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
       if len(f:=line.rstrip().split('\t'))==2 and f[1] in ('PATERNAL','MATERNAL')}
tags,status=audit.assignments(hiphase,truth)
assert names<=tags.keys() and names<=status.keys()
work=ROOT/'test_data/tmp_representation_step8'; cache=work/'hiphase-molecule-tags.json'
cache.write_text(json.dumps([{n:tags[n] for n in names},{n:status[n] for n in names}]))
report=dict(verified_identical_primary_alignments=len(verified),status_cache=str(cache.relative_to(ROOT)),global_primary_reads=len(tags),global_orientation=True,contig_alias='CHM13#0#chr20 -> chr20',sequence_qualities_mapq_flags_verified=True,hiphase=str(hiphase.relative_to(ROOT)),helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/'hiphase-molecule-checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
