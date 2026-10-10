#!/usr/bin/env python3
"""Independently rebuild physical hypotheses and score original saved BAM slices."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import time

import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--work', type=Path, default=ROOT/'test_data/tmp_representation_step6/accepted')
args = parser.parse_args()

def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def key(text):
    if text == 'REF':
        return (0,)
    pos, kind, ref_len, alt = text.split(':', 3)
    return (1, int(pos), int(kind), int(ref_len), alt)

def distance(a, b):
    row = list(range(len(b) + 1))
    for i, x in enumerate(a):
        diagonal, row[0] = row[0], i + 1
        for j, y in enumerate(b):
            old = row[j+1]
            row[j+1] = min(row[j] + 1, row[j+1] + 1, diagonal + (x != y))
            diagonal = old
    return row[-1]

reports = []
started = time.monotonic()
with pysam.FastaFile(str(ROOT/'test_data/chm13v2.0.chr20.renamed.fa')) as fasta, \
     pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    for name in ['matrix4-final', 'matrix5-final', 'matrix65-final', 'whole-final']:
        folder = args.work/name
        groups = defaultdict(list)
        contrasts = {r['candidate']: r for r in rows(folder/'matrix.chunk0.allele-contrasts.tsv')}
        for r in contrasts.values():
            if r['allele0'] != '.':
                groups[(key(r['allele0']), key(r['allele1']))].append(r['candidate'])
        pairs = sorted(k for k, members in groups.items() if len(members) > 1)
        parents = defaultdict(list)
        for r in rows(folder/'matrix.chunk0.joint-parents.tsv'):
            parents[int(r['locus'])].append(r)
        contexts = rows(folder/'matrix.chunk0.joint-contexts.tsv')
        sequences = rows(folder/'matrix.chunk0.joint-sequences.tsv')
        molecules = rows(folder/'matrix.chunk0.joint-molecules.tsv')
        assert {(r['locus'], r['read'], r['source_allele']) for r in sequences} == \
               {(r['locus'], r['read'], r['allele']) for r in molecules}
        assert len(contexts) == len(pairs) and len(sequences) == len(molecules)
        models = {}
        for r in contexts:
            locus = int(r['locus'])
            assert {p['candidate'] for p in parents[locus]} == set(groups[pairs[locus]])
            if r['status'] != 'valid':
                continue
            lo, hi = int(r['beg']), int(r['end'])
            reference = fasta.fetch('CHM13#0#chr20', lo-1, hi).upper()
            assert reference == r['reference']
            selected = []
            begins = [int(p['pos']) for p in parents[locus]]
            ends = [int(p['pos']) + len(p['ref']) - 1 for p in parents[locus]]
            for edit in pairs[locus]:
                if edit == (0,):
                    selected.append(reference)
                    continue
                _, pos, kind, length, alt = edit
                physical = pos if kind == 8 else pos + 1
                begins.append(physical)
                ends.append(physical + max(1, length) - 1)
                offset = physical - lo
                selected.append(reference[:offset] + alt + reference[offset + length:])
            assert lo == min(begins) - 16 and hi == max(ends) + 16
            assert selected == [r['allele0'], r['allele1']]
            alternatives = {reference}
            for p in parents[locus]:
                offset = int(p['pos']) - lo
                assert reference[offset:offset+len(p['ref'])] == p['ref'].upper()
                for alt in p['alts'].split(','):
                    alternatives.add(reference[:offset] + alt.upper() + reference[offset+len(p['ref']):])
            alternatives.difference_update(selected)
            assert alternatives == set(r['other_alleles'].split(',')) - {''}
            models[r['locus']] = (selected, alternatives)
        scored = [r for r in sequences if r['status'] == 'scored']
        names = {r['read'] for r in scored}
        begin = min(int(r['beg']) for r in contexts if r['status'] == 'valid')
        end = max(int(r['end']) for r in contexts if r['status'] == 'valid')
        alignments = {}
        for alignment in bam.fetch('CHM13#0#chr20', begin-1, end):
            if alignment.is_secondary or alignment.is_supplementary or alignment.query_name not in names:
                continue
            assert alignment.query_name not in alignments, 'ambiguous original primary read'
            alignments[alignment.query_name] = alignment
        assert names == set(alignments)
        anchors = defaultdict(set)
        for r in scored:
            c = contexts[int(r['locus'])]
            anchors[r['read']].update([int(c['beg'])-1, int(c['end'])-1])
        mappings = {}
        for name_, a in alignments.items():
            mappings[name_] = {r: q for q, r in a.get_aligned_pairs()
                               if q is not None and r in anchors[name_]}
        decisions = Counter()
        for r in scored:
            a = alignments[r['read']]
            c = contexts[int(r['locus'])]
            lo, hi = int(c['beg']), int(c['end'])
            qs, qe = int(r['query_beg']), int(r['query_end'])
            assert mappings[r['read']][lo-1] == qs and mappings[r['read']][hi-1] + 1 == qe
            assert a.query_sequence[qs:qe].upper() == r['sequence']
            quals = bytes(a.query_qualities[qs:qe]) if a.query_qualities is not None else bytes([255])*(qe-qs)
            assert quals.hex() == r['qualities_hex'] and a.mapping_quality == int(r['mapq'])
            selected, others = models[r['locus']]
            costs = [distance(r['sequence'], allele) for allele in selected]
            other_cost = min((distance(r['sequence'], allele) for allele in others), default=-1)
            winner = 0 if costs[0] < costs[1] else 1
            nearest = winner if costs[winner] < costs[1-winner] and (other_cost < 0 or costs[winner] < other_cost) else -1
            assert costs == [int(r['distance0']), int(r['distance1'])]
            assert other_cost == int(r['other_distance']) and nearest == int(r['nearest'])
            if r['source_allele'] == '-2':
                if nearest >= 0:
                    decisions['selected_'+str(nearest)] += 1
                elif other_cost >= 0 and other_cost < min(costs):
                    decisions['other_allele_closer'] += 1
                else:
                    decisions['tie'] += 1
        report = dict(owner=name, contexts=len(contexts), molecules=len(sequences),
                      states=dict(Counter(r['status'] for r in sequences)),
                      source_conflicts=sum(r['source_allele']=='-2' for r in sequences),
                      conflict_rankings=dict(decisions),
                      original_sequence_quality_cigar_verified=len(scored),
                      independent_scalar_dp_verified=len(scored))
        reports.append(report)
        print(json.dumps(report), flush=True)
(OUT/'sequence-checks.json').write_text(json.dumps(dict(seconds=time.monotonic()-started, owners=reports), indent=2)+'\n')
