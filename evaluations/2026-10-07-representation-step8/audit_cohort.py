#!/usr/bin/env python3
"""Independently census full parent coverage and verify saved raw BAM slices."""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
WORK = ROOT / 'test_data/tmp_representation_step8/accepted'

def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

reports = []
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    for owner, start, end in [('matrix4-final',4000000,5000000), ('matrix5-final',5000000,6000000),
                              ('matrix65-final',65000000,66000000), ('whole-final',4000000,5000000)]:
        folder = WORK/owner
        loaded = [r for r in bam.fetch('CHM13#0#chr20',start,end)
                  if not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail)
                  and r.mapping_quality >= 1]
        multiplicity = Counter(r.query_name for r in loaded)
        raw = {r.query_name:r for r in loaded if multiplicity[r.query_name] == 1}
        saved = {(r['locus'],r['read']):r for r in rows(folder/'matrix.chunk0.joint-cohort.tsv')}
        assert len(saved) == len(rows(folder/'matrix.chunk0.joint-cohort.tsv'))
        original = {(r['locus'],r['read']):r for r in rows(folder/'matrix.chunk0.joint-sequences.tsv') if r['status']=='scored'}
        loci = []
        expected_keys = set()
        for c in rows(folder/'matrix.chunk0.joint-contexts.tsv'):
            if c['status'] != 'valid':
                assert not any(k[0]==c['locus'] for k in saved)
                continue
            lo, hi = int(c['beg'])-1, int(c['end'])-1
            expected = {}
            for name, r in raw.items():
                if r.reference_start > lo or r.reference_end <= hi:
                    continue
                ref = r.reference_start
                skipped = False
                for op, length in r.cigartuples:
                    if op == 3 and ref <= hi and ref+length > lo:
                        skipped = True
                        break
                    if op in (0,2,3,7,8):ref += length
                if skipped:continue
                anchors = {ref:q for q,ref in r.get_aligned_pairs() if ref in (lo,hi) and q is not None}
                if set(anchors) != {lo,hi}:continue
                beg, end_ = anchors[lo], anchors[hi]+1
                sequence = r.query_sequence[beg:end_].upper()
                if not sequence or len(sequence)>4096 or set(sequence)-set('ACGT'):continue
                key = c['locus'],name
                expected_keys.add(key)
                assert key in saved, (owner,key,'omitted complete primary')
                s = saved[key]
                assert s['status']=='scored' and int(s['query_beg'])==beg and int(s['query_end'])==end_
                assert s['sequence']==sequence and int(s['mapq'])==r.mapping_quality
                assert bytes.fromhex(s['qualities_hex'])==bytes(r.query_qualities[beg:end_])
                if key in original:assert s==original[key], (owner,key,'source provenance changed')
                else:assert s['source_allele']=='.', (owner,key,'invented source allele')
                expected[name]=s
            old = {n for l,n in original if l==c['locus']}
            assert old <= expected.keys()
            loci.append(dict(locus=int(c['locus']),old_molecules=len(old),full_molecules=len(expected),added_molecules=len(expected.keys()-old)))
        assert set(saved)==expected_keys, (owner,'ineligible extra molecules',set(saved)-expected_keys)
        costs = rows(folder/'matrix.chunk0.joint-costs.tsv')
        assert len(costs)==len(saved) and {(r['locus'],r['read']) for r in costs}==expected_keys
        hypotheses = defaultdict(list)
        for r in rows(folder/'matrix.chunk0.joint-alleles.tsv'):
            hypotheses[r['locus']].append(r['sequence'])
        physical = defaultdict(list)
        for c in rows(folder/'matrix.chunk0.joint-contexts.tsv'):
            if c['status']=='valid':physical[(c['beg'],c['end'],tuple(hypotheses[c['locus']]))].append(c['locus'])
        fits = {r['locus']:{k:v for k,v in r.items() if k!='locus'} for r in rows(folder/'matrix.chunk0.joint-genotypes.tsv')}
        shared_contexts = []
        for members in physical.values():
            if len(members)<2:continue
            shared_contexts.append(members)
            vectors = [{r['read']:r['costs'] for r in costs if r['locus']==l} for l in members]
            assert all(v==vectors[0] for v in vectors), (owner,members,'source-selected contrast changed complete costs')
            assert all(fits[l]==fits[members[0]] for l in members), (owner,members,'source-selected contrast changed genotype')
        report = dict(owner=owner, shared_physical_context_groups=shared_contexts,contexts=len(loci),molecule_locus_records=len(saved),
                      added_records=len(saved)-len(original),affected_contexts=sum(r['added_molecules']>0 for r in loci),
                      raw_sequences_qualities_bounds_mapq_verified=True,old_provenance_identical=True,loci=loci)
        reports.append(report)
        print(json.dumps({k:v for k,v in report.items() if k!='loci'}),flush=True)
assert [r['added_records'] for r in reports] == [82,72,97,0]
(OUT/'cohort-checks.json').write_text(json.dumps(dict(owners=reports,added_records=sum(r['added_records'] for r in reports)),indent=2)+'\n')
