#!/usr/bin/env python3
"""Verify pair enumeration and every valid sequence with independent trim placements."""
from collections import Counter,defaultdict
import csv
import hashlib
import itertools
import json
from pathlib import Path
import time

import pysam

ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step10';BASE=ROOT/'test_data/tmp_representation_step9/accepted'
previous=json.loads((ROOT/'evaluations/2026-10-07-representation-step9/manifest.json').read_text())
known={r['path']:r['sha256'] for r in previous['files']}

def rows(path):
 with path.open() as h:return list(csv.DictReader(h,delimiter='\t'))

def placements(pos,ref,alt):
 # Enumerate all maximal shared-padding removals, without left alignment.
 if ref==alt:return [None]
 prefix=0
 while prefix<min(len(ref),len(alt)) and ref[prefix]==alt[prefix]:prefix+=1
 suffix=0
 while suffix<min(len(ref),len(alt)) and ref[-1-suffix]==alt[-1-suffix]:suffix+=1
 removed=min(len(ref),len(alt),prefix+suffix)
 choices=[]
 for left in range(max(0,removed-suffix),min(prefix,removed)+1):
  right=removed-left
  choices.append((pos+left,len(ref)-removed,alt[left:len(alt)-right if right else len(alt)]))
 return choices

def render(reference,beg,edits):
 events=sorted(set(e for e in edits if e is not None));cursor=beg;chunks=[];insertions=set()
 for pos,length,alt in events:
  if pos<cursor or length==0 and pos in insertions:return None
  chunks.extend([reference[cursor-beg:pos-beg],alt]);cursor=pos+length
  if length==0:insertions.add(pos)
 return ''.join(chunks)+reference[cursor-beg:]

reports=[];start=time.monotonic()
with pysam.FastaFile(str(ROOT/'test_data/chm13v2.0.chr20.renamed.fa')) as fasta:
 for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
  folder=WORK/'accepted'/owner;before=BASE/owner
  for pattern in ['matrix.chunk0.joint-*.tsv','matrix.chunk0.physical-*.tsv']:
   for p in before.glob(pattern):
    assert hashlib.sha256(p.read_bytes()).hexdigest()==known[str(p.relative_to(ROOT))],p
    assert p.read_bytes()==(folder/p.name).read_bytes(),(owner,p.name,'baseline changed')
  contexts={int(r['physical']):r for r in rows(folder/'matrix.chunk0.composition-contexts.tsv')}
  physical=rows(folder/'matrix.chunk0.physical-contexts.tsv')
  assert len(contexts)==len(physical)
  for c in physical:
   ctx=contexts[int(c['physical'])];assert all(ctx[k]==v for k,v in c.items())
   assert ctx['reference']==fasta.fetch('CHM13#0#chr20',int(ctx['beg'])-1,int(ctx['end'])).upper()
  alleles=defaultdict(list)
  for r in rows(folder/'matrix.chunk0.physical-alleles.tsv'):
   assert int(r['allele'])==len(alleles[int(r['physical'])]);alleles[int(r['physical'])].append(r['sequence'])
  parent_candidates=defaultdict(set);loci=defaultdict(set)
  for r in rows(folder/'matrix.chunk0.joint-parents.tsv'):loci[int(r['locus'])].add(int(r['candidate']))
  for r in rows(folder/'matrix.chunk0.physical-members.tsv'):parent_candidates[int(r['physical'])].update(loci[int(r['locus'])])
  catalog=rows(folder/'matrix.chunk0.site-catalog.tsv')
  assert len({(r['candidate'],r['alt_index']) for r in catalog})==len(catalog)
  assert {r['candidate'] for r in catalog}=={r['candidate'] for r in rows(folder/'matrix.chunk0.allele-contrasts.tsv')}
  by_candidate=defaultdict(list)
  for r in catalog:by_candidate[int(r['candidate'])].append(r)
  for r in rows(folder/'matrix.chunk0.joint-parents.tsv'):
   source=by_candidate[int(r['candidate'])]
   assert all(s['pos']==r['pos'] and s['ref']==r['ref'] for s in source)
   assert [s['alt'] for s in source]==r['alts'].split(',')
  expected_sites=[];expected_pairs=[]
  for id_,ctx in contexts.items():
   beg,end=int(ctx['beg']),int(ctx['end'])
   for r in catalog:
    pos=int(r['pos'])
    if int(r['candidate']) in parent_candidates[id_] or not r['ref'] or pos<=beg or pos+len(r['ref'])-1>=end:continue
    expected_sites.append(dict(physical=str(id_),**r))
    for parent in range(len(alleles[id_])):expected_pairs.append((id_,parent,int(r['candidate']),int(r['alt_index'])))
  assert rows(folder/'matrix.chunk0.composition-sites.tsv')==expected_sites
  compositions=rows(folder/'matrix.chunk0.compositions.tsv')
  assert [(int(r['physical']),int(r['parent_allele']),int(r['candidate']),int(r['alt_index'])) for r in compositions]==expected_pairs
  sites={(int(r['physical']),int(r['candidate']),int(r['alt_index'])):r for r in expected_sites}
  full={id_:set(v) for id_,v in alleles.items()};statuses=Counter();checked=0
  for r in compositions:
   id_=int(r['physical']);ctx=contexts[id_];beg=int(ctx['beg']);ref=ctx['reference']
   parent=alleles[id_][int(r['parent_allele'])];site=sites[id_,int(r['candidate']),int(r['alt_index'])]
   pos=int(site['pos']);raw_ref,raw_alt=site['ref'].upper(),site['alt'].upper();statuses[r['status']]+=1
   if r['status']!='valid':assert r['sequence']=='.';continue
   assert set(raw_ref+raw_alt)<=set('ACGT') and ref[pos-beg:pos-beg+len(raw_ref)]==raw_ref
   parent_choices=placements(beg,ref,parent);neighbor_choices=placements(pos,raw_ref,raw_alt)
   # For valid pairs the conservative interval guard makes all padding placements
   # interchangeable; independently check both extreme anchored placements.
   choices=[list(dict.fromkeys([parent_choices[0],parent_choices[-1]])),list(dict.fromkeys([neighbor_choices[0],neighbor_choices[-1]]))]
   expected={render(ref,beg,pair) for pair in itertools.product(*choices)}
   assert expected=={r['sequence']},(owner,r,expected)
   assert len(r['sequence'])<=4096
   full[id_].add(r['sequence']);checked+=1
  hypotheses=defaultdict(list)
  for r in rows(folder/'matrix.chunk0.composed-alleles.tsv'):
   id_=int(r['physical']);assert int(r['allele'])==len(hypotheses[id_]);hypotheses[id_].append(r['sequence'])
  assert {id_:sorted(v) for id_,v in full.items()}==dict(hypotheses)
  for kind in ['composition-sites','compositions','composed-alleles']:
   assert (folder/f'matrix.chunk0.{kind}.tsv').read_bytes()==(WORK/'cached'/f'{owner}.{kind}.tsv').read_bytes(),(owner,kind,'cached mismatch')
  report=dict(owner=owner,contexts=len(contexts),catalog_alternatives=len(catalog),neighbor_alternatives=len(expected_sites),composition_attempts=len(compositions),statuses=dict(statuses),independently_verified_valid_sequences=checked,
   old_hypotheses=sum(map(len,alleles.values())),composed_hypotheses=sum(map(len,hypotheses.values())),new_unique_sequences=sum(len(full[id_]-set(v)) for id_,v in alleles.items()),legacy_and_physical_state_identical=True)
  reports.append(report);print(json.dumps(report),flush=True)
for record in previous['input_metadata']:
 stat=(ROOT/record['path']).stat();assert stat.st_size==record['size'] and stat.st_mtime_ns==record['mtime_ns']
(OUT/'composition-checks.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,owners=reports),indent=2)+'\n')
