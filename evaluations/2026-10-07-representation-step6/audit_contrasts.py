#!/usr/bin/env python3
"""Validate saved canonical molecule evidence against original observation channels."""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
reports = []
def allele_key(value):
    if value == 'REF':
        return (0,)
    pos, kind, length, sequence = value.split(':', 3)
    return (1, int(pos), int(kind), int(length), sequence)

for name in ['matrix4-final', 'matrix5-final', 'matrix65-final', 'whole-final']:
    work = ROOT / 'test_data/tmp_representation_step6/accepted' / name
    contrast_rows = list(csv.DictReader((work / 'matrix.chunk0.allele-contrasts.tsv').open(), delimiter='\t'))
    explicit = {int(row['candidate']): row for row in contrast_rows if row['local0'] != '.'}
    loci = defaultdict(list)
    for ci, row in explicit.items():
        loci[allele_key(row['allele0']), allele_key(row['allele1'])].append(ci)
        assert row['non_selected_alt_class'] == '0'
        assert row['local0'] != row['local1']
    assert all(row['local0'] == '.' for row in contrast_rows if row['non_selected_alt_class'] == '1')
    if name == 'whole-final':
        assert sum(row['allele0'] != 'REF' and row['allele1'] != 'REF' for row in explicit.values()) == 3
        assert sum(row['non_selected_alt_class'] == '1' for row in contrast_rows) == 99
    duplicate_groups = [loci[key] for key in sorted(loci) if len(loci[key]) > 1]
    original = {}
    for line in (work / 'matrix.chunk0.bam-overlay-output.tsv').open():
        if line.startswith('OBS\t'):
            _, read, ci, primary, graph, bam = line.rstrip().split('\t')
            original[read, int(ci)] = tuple(map(int, [primary, graph, bam]))
    rows = list(csv.DictReader((work / 'matrix.chunk0.joint-molecules.tsv').open(), delimiter='\t'))
    assert len({(r['locus'], r['read']) for r in rows}) == len(rows)
    reasons = Counter()
    conflicting_sources = Counter()
    for row in rows:
        members = duplicate_groups[int(row['locus'])]
        calls = [tuple(map(int, value.split(':'))) for value in row['observations'].split(',') if value]
        assert len(calls) == len(set(calls)) and calls
        represented = set()
        for ci, call, qi in calls:
            assert ci in members
            physical = original[row['read'], ci]
            if call == -2:
                assert -2 in physical
            else:
                assert call in (0, 1)
                assert int(explicit[ci]['local' + str(call)]) in physical
            represented.add((ci, call))
        for ci in members:
            physical = original.get((row['read'], ci))
            if physical is None:
                continue
            local = [int(explicit[ci]['local0']), int(explicit[ci]['local1'])]
            for call in physical[1:]:
                if call in local:
                    assert (ci, local.index(call)) in represented
                elif call == -2:
                    assert (ci, -2) in represented
            if physical[1] in local and physical[2] in local and physical[1] != physical[2]:
                reasons['opposing_graph_bam_on_same_candidate'] += 1
        unique = {call for _, call, _ in calls}
        expected = -2 if -2 in unique or {0, 1} <= unique else next(iter(unique))
        assert int(row['allele']) == expected
        if expected == -2:
            signature = tuple(sorted({(explicit[ci]['bam_injected'], explicit[ci]['category'], call)
                                      for ci, call, qi in calls}))
            conflicting_sources[signature] += 1
            assert int(row['query_index']) == 0
            if -2 in unique:
                reasons['contains_existing_conflict'] += 1
            else:
                reasons['opposing_retained_calls'] += 1
        else:
            indices = {qi for _, _, qi in calls}
            assert int(row['query_index']) == (next(iter(indices)) if len(indices) == 1 else 0)
    report = dict(owner=name, candidates=len(contrast_rows), explicit_pairs=len(explicit),
        selected_alt_alt_pairs=sum(row['allele0'] != 'REF' and row['allele1'] != 'REF' for row in explicit.values()),
        collapsed_other_alt_classes=sum(row['non_selected_alt_class'] == '1' for row in contrast_rows),
        duplicate_complete_pairs=len(duplicate_groups), joint_molecule_loci=len(rows),
        molecule_calls=dict(Counter(row['allele'] for row in rows)),
        unique_conflicting_reads=len({row['read'] for row in rows if row['allele'] == '-2'}),
        conflict_reasons=dict(reasons),
        conflicting_source_patterns=[dict(sources=[dict(bam_injected=bam, category=category, canonical_allele=call)
            for bam, category, call in signature], molecule_loci=count)
            for signature, count in conflicting_sources.items()], all_calls_trace_to_original_channels=True,
        all_retained_graph_bam_pair_calls_represented=True,
        independent_reduction_and_query_coordinates_match=True)
    reports.append(report)
(OUT / 'contrast-checks.json').write_text(json.dumps(reports, indent=2) + '\n')
print(json.dumps(reports, indent=2))
