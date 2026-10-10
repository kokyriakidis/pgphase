#!/usr/bin/env python3
"""Certify exported matches by counting all optimal alignments independently."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
WORK = ROOT / 'test_data/tmp_representation_step12'
BASE = ROOT / 'test_data/tmp_representation_step11/accepted-final'

def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

def costs_and_counts(a, b):
    costs = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    ways = [[1] * (len(b) + 1) for _ in range(len(a) + 1)]
    costs[0] = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        costs[i][0] = i
        for j in range(1, len(b) + 1):
            edges = [(costs[i-1][j]+1, ways[i-1][j]),
                     (costs[i][j-1]+1, ways[i][j-1]),
                     (costs[i-1][j-1]+(a[i-1]!=b[j-1]), ways[i-1][j-1])]
            best = min(c for c, _ in edges)
            costs[i][j] = best
            ways[i][j] = sum(n for c, n in edges if c == best)
    return costs, ways

def oracle(ref, alt):
    costs, ways = costs_and_counts(ref, alt)
    backwards, suffixes = costs_and_counts(ref[::-1], alt[::-1])
    optimum, total = costs[-1][-1], ways[-1][-1]
    matches = []
    for i, base in enumerate(ref):
        for j, alt_base in enumerate(alt):
            if base == alt_base and costs[i][j] + backwards[len(ref)-i-1][len(alt)-j-1] == optimum:
                through = ways[i][j] * suffixes[len(ref)-i-1][len(alt)-j-1]
                if through == total:
                    matches.append((i, j))
    spans = []
    for i, j in matches:
        if spans and spans[-1][0] + spans[-1][2] == i and spans[-1][1] + spans[-1][2] == j:
            spans[-1][2] += 1
        else:
            spans.append([i, j, 1])
    return optimum, spans, total

reports = []
started = time.monotonic()
for owner in ['matrix4-final', 'matrix5-final', 'matrix65-final', 'whole-final']:
    folder = WORK / 'accepted-final' / owner
    prefix = WORK / 'cached' / owner
    subprocess.run([str(ROOT/'test_allele_context'), '--maps', str(folder), str(prefix)], check=True)
    previous = list((BASE/owner).glob('matrix.chunk0.*.tsv'))
    for path in previous:
        assert path.read_bytes() == (folder/path.name).read_bytes(), (owner, path.name)
    for kind in ['composition-sites','compositions','composed-alleles','path-status','sequence-paths','path-alleles','parent-maps','matched-subpaths']:
        assert (folder/f'matrix.chunk0.{kind}.tsv').read_bytes() == Path(str(prefix)+f'.{kind}.tsv').read_bytes(), (owner, kind)
    contexts = {r['physical']:r for r in table(folder/'matrix.chunk0.composition-contexts.tsv')}
    alleles = {(r['physical'],r['allele']):r['sequence'] for r in table(folder/'matrix.chunk0.physical-alleles.tsv')}
    maps = {(r['physical'],r['parent_allele']):r for r in table(folder/'matrix.chunk0.parent-maps.tsv')}
    spans = defaultdict(list)
    for r in table(folder/'matrix.chunk0.matched-subpaths.tsv'):
        spans[(r['physical'],r['parent_allele'])].append([int(r[k]) for k in ['ref_beg','alt_beg','length']])
    assert maps.keys() == alleles.keys()
    ambiguous, matched, bases = 0, 0, 0
    for key, alt in alleles.items():
        ref = contexts[key[0]]['reference']
        distance, expected, total = oracle(ref, alt)
        assert maps[key]['status'] == 'complete'
        assert int(maps[key]['distance']) == distance and spans[key] == expected, (owner, key)
        length = sum(s[2] for s in expected)
        assert int(maps[key]['matched_bases']) == length and int(maps[key]['subpaths']) == len(expected)
        ambiguous += total > 1
        matched += length
        bases += len(ref)
    sites = {(r['physical'],r['candidate'],r['alt_index']):r for r in table(folder/'matrix.chunk0.composition-sites.tsv')}
    nested = Counter()
    examples = []
    for r in table(folder/'matrix.chunk0.compositions.tsv'):
        if r['status'] != 'overlap': continue
        key = (r['physical'],r['parent_allele'])
        site = sites[(r['physical'],r['candidate'],r['alt_index'])]
        context = contexts[key[0]]
        offset = int(site['pos']) - int(context['beg'])
        certified = any(i <= offset and offset + len(site['ref']) <= i+n for i,j,n in spans[key])
        nested['overlap_attempts'] += 1
        nested['whole_raw_ref_certified' if certified else 'unresolved_raw_ref'] += 1
        if certified and len(examples) < 5:
            examples.append(dict(physical=key[0],parent=key[1],candidate=site['candidate'],pos=site['pos'],ref=site['ref'],alt_index=site['alt_index']))
    reports.append(dict(owner=owner,contexts=len(contexts),parents=len(alleles),complete=len(maps),multiple_optimal_alignments=ambiguous,
        reference_bases=bases,certified_matched_bases=matched,matched_subpaths=sum(map(len,spans.values())),
        previous_tables_byte_identical=True,raw_replay_byte_identical=True,overlaps=dict(nested),examples=examples))
result = dict(oracle='all-optimal path counts with arbitrary precision, matched-edge count equals total',
              wall_seconds=time.monotonic()-started,owners=reports)
(OUT/'map-checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
