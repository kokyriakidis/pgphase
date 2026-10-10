#!/usr/bin/env python3
"""Measure retained evidence and diagnostic SNP votes; never emit phasing."""
import argparse
import importlib.util
import json
from collections import Counter, defaultdict
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--replays', type=Path, default=ROOT/'test_data/tmp_joint_evidence_investigation')
args = parser.parse_args()
spec = importlib.util.spec_from_file_location('read_audit', ROOT/'evaluations/2026-10-06-shared-insertion-source/read_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth = {f[0]: f[1] == 'PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
hi_tags, hi_status = json.loads((ROOT/'test_data/tmp_gap_fix78/hi-tags.json').read_text())
results = []
for chunk, left, right in [(4, 4766928, 4792960), (5, 5309406, 5345085), (6, 6513891, 6516221)]:
    directory = args.replays/str(chunk)
    tags, status = audit.assignments(directory/'phased.bam', truth)
    votes = defaultdict(Counter)
    for q, (hp, ps) in tags.items():
        if q in truth and hp in (1, 2) and ps > 0:
            votes[ps][(hp == 1) != truth[q]] += 1
    orientation = {ps: v[True] >= v[False] for ps, v in votes.items()}
    with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
        names = {r.query_name for r in bam.fetch('CHM13#0#chr20', left-1, right)
                 if not r.is_secondary and not r.is_supplementary and r.query_name in truth}
    core = Counter(tags[q][1] for q in names if q in tags and tags[q][0] in (1, 2) and 0 < tags[q][1] < 1_000_000_000).most_common(1)[0][0]
    fixed = {}
    with pysam.VariantFile(str(directory/'phased.vcf')) as vcf:
        for row in vcf:
            sample = next(iter(row.samples.values()))
            if len(row.ref) == 1 and row.alts and len(row.alts) == 1 and len(row.alts[0]) == 1 and sample.phased and sample.get('PS') == core and set(sample.get('GT', ())) == {0, 1}:
                fixed[row.pos] = sample['GT'][0]
    variants, reads, observations = {}, {}, defaultdict(dict)
    matrix = directory/'matrix.chunk0.bam-overlay-output.tsv'
    for line in matrix.open():
        f = line.rstrip().split('\t')
        if f[0] == 'VAR':
            variants[int(f[1])] = {'pos': int(f[2]), 'type': f[3], 'flags': int(f[4]), 'ps': int(f[8]), 'h1': int(f[9]), 'h2': int(f[10])}
        elif f[0] == 'READ':
            reads[f[1]] = {'skipped': int(f[5]) != 0}
        elif f[0] == 'OBS':
            observations[f[1]][int(f[2])] = tuple(map(int, f[3:6]))
    anchors = {i: fixed[v['pos']] for i, v in variants.items() if v['type'] == 'X' and v['pos'] in fixed}
    targets = {q for q in names if hi_status.get(q) == 'correct' and (status.get(q) != 'correct' or tags.get(q, (0, 0))[1] != core)}
    modes = {}
    for mode in ['effective', 'graph', 'bam', 'agreement_union']:
        def call(channels):
            effective, graph, bam = channels
            if mode != 'agreement_union':
                return {'effective': effective, 'graph': graph, 'bam': bam}[mode]
            if graph in (0, 1) and bam in (0, 1) and graph != bam:
                return -1
            return graph if graph in (0, 1) else bam
        profiles = {q: {i: a for i, channels in p.items() if (a := call(channels)) in (0, 1)}
                    for q, p in observations.items() if not reads[q]['skipped']}
        seen = defaultdict(set)
        for p in profiles.values():
            for i, a in p.items():
                seen[i].add(a)
        # Connectivity is an optimistic availability diagnostic, not a join
        # certificate: repeat rows and correlated calls can supply false edges.
        nodes = {i for i, alleles in seen.items() if alleles == {0, 1} and not
                 (variants[i]['ps'] > 0 and variants[i]['h1'] == variants[i]['h2'] and variants[i]['h1'] >= 0)} | set(anchors)
        parent = {i: i for i in nodes}
        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for p in profiles.values():
            known = [i for i in p if i in nodes]
            for i in known[1:]:
                parent[find(i)] = find(known[0])
        core_roots = {find(i) for i in anchors}
        detail, measures = [], Counter()
        for q in sorted(targets):
            p = profiles.get(q, {})
            counts = Counter(a == anchors[i] for i, a in p.items() if i in anchors)
            connected = any(i in nodes and find(i) in core_roots for i in p)
            measures['targets'] += 1
            measures['missing_read_profile'] += q not in reads
            measures['no_callable_observations'] += not p
            measures['potentially_connected'] += connected
            predicted = 0 if counts[True] == counts[False] else 1 if counts[True] > counts[False] else 2
            correct = bool(predicted and ((predicted == 1) != truth[q]) == orientation[core])
            if predicted:
                measures['direct_snp_vote_correct' if correct else 'direct_snp_vote_wrong'] += 1
                measures['conflicting_direct_snp_votes'] += bool(counts[True] and counts[False])
            else:
                measures['direct_snp_vote_abstains'] += 1
            detail.append({'read': q, 'pgphase_status': status.get(q, 'unphased'), 'tags': tags.get(q),
                           'callable_sites': len(p), 'core_snp_votes': dict(counts), 'potentially_connected': connected,
                           'diagnostic_hap': predicted, 'diagnostic_correct': correct})
        modes[mode] = {'counts': dict(measures), 'reads': detail}
    results.append({'window': f'{left}-{right}', 'owning_region': [chunk*1_000_000+1, (chunk+1)*1_000_000],
                    'matrix': str(matrix), 'dominant_core': core, 'scorable': len(names),
                    'pgphase_counts': dict(Counter(status.get(q, 'unphased') for q in names)),
                    'pgphase_core_correct': sum(status.get(q) == 'correct' and tags[q][1] == core for q in names),
                    'hiphase_correct': sum(hi_status.get(q) == 'correct' for q in names),
                    'modes': modes})
(OUT/'matrix-audit.json').write_text(json.dumps(results, indent=2)+'\n')
for row in results:
    print(row['window'], row['pgphase_counts'], 'core', row['pgphase_core_correct'], 'hiphase', row['hiphase_correct'])
    for mode, measurements in row['modes'].items():
        print(mode, measurements['counts'])
