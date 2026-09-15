#!/usr/bin/env python3
"""Replay frozen gap proposals. Scores are diagnostics, not calibrated probabilities."""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import time

csv.field_size_limit(16 * 1024 * 1024)
FEATURE_SCHEMA = 2  # Bump when observation grouping or feature extraction changes.
VARIANT_SNP = 8  # VariantType::Snp follows BAM_CDIFF.


def orientation_costs(votes):
    left, right = votes
    ld, lo = left[0] + left[3], left[1] + left[2]
    rd, ro = right[0] + right[3], right[1] + right[2]
    return min(lo + ro, ld + rd), min(lo + rd, ld + ro)


def classify(votes, margin):
    if any(sum(side) == 0 for side in votes):
        return 'INSUFFICIENT'
    same, flip = orientation_costs(votes)
    if abs(same - flip) + 1e-9 < margin:
        return 'CONFLICTING' if same > 0 and flip > 0 else 'INSUFFICIENT'
    return 'SAME' if same < flip else 'FLIP'


def summarize(votes, contributions, margin):
    same, flip = orientation_costs(votes)
    decision = classify(votes, margin)
    reversed_reads, lost_support = [], []
    if decision in ('SAME', 'FLIP'):
        for molecule, side, cell, weight in contributions:
            reduced = [list(v) for v in votes]
            reduced[side][cell] -= weight
            after = classify(reduced, margin)
            if after in ('SAME', 'FLIP') and after != decision:
                reversed_reads.append(molecule)
            elif after not in ('SAME', 'FLIP'):
                lost_support.append(molecule)
    return dict(decision=decision, same_cost=same, flip_cost=flip,
                margin=abs(same - flip), votes=votes,
                reversing_molecules=reversed_reads,
                support_critical_molecules=lost_support)


def resolve_components(edges):
    """Reject inconsistent components; tree edges still require local validation."""
    adjacency = collections.defaultdict(list)
    for left, right, flip in edges:
        adjacency[left].append((right, bool(flip)))
        adjacency[right].append((left, bool(flip)))
    seen, components = set(), []
    for root in sorted(adjacency):
        if root in seen:
            continue
        parity, pending, conflicts = {root: False}, [root], set()
        while pending:
            node = pending.pop()
            seen.add(node)
            for other, flip in sorted(adjacency[node]):
                expected = parity[node] ^ flip
                if other in parity:
                    if parity[other] != expected:
                        conflicts.add((min(node, other), max(node, other)))
                else:
                    parity[other] = expected
                    pending.append(other)
        components.append(dict(root=root, consistent=not conflicts,
                               nodes=sorted(parity),
                               flips={str(n): parity[n] for n in sorted(parity)} if not conflicts else {},
                               conflicts=[list(e) for e in sorted(conflicts)]))
    return components


def direct_observation_support(gap, sites, reads, observations):
    """Diagnostic per-molecule votes from frozen block sites, including HP0 reads.

    Overlap grouping prevents duplicate local events from multiplying votes;
    it does not claim sequence equivalence or calibrate allele errors.
    """
    events, end, event = {}, -1, -1
    for vi, site in sorted(sites.items(), key=lambda item: (item[1]['pos'] - (item[1]['type'] != VARIANT_SNP), item[0])):
        beg = site['pos'] - (site['type'] != VARIANT_SNP)
        stop = site['pos'] + site['length'] if site['type'] != VARIANT_SNP else site['pos']
        if beg > end:
            event += 1
        end = max(end, stop)
        events[vi] = event
    result = {}
    for channel in ('hybrid', 'graph'):
        molecule_events = collections.defaultdict(lambda: collections.defaultdict(set))
        for ri, vi, hybrid, graph, qi, quality in observations:
            read, site = reads[ri], sites[vi]
            if read['skipped'] or (read['flag'] >= 0 and read['flag'] & (4 | 256 | 2048)):
                continue
            side = 0 if site['ps'] == gap[1] else 1 if site['ps'] == gap[2] else -1
            allele = hybrid if channel == 'hybrid' else graph
            if side < 0 or allele < 0 or site['h1'] < 0 or site['h2'] < 0 or site['h1'] == site['h2']:
                continue
            hap = 0 if allele == site['h1'] else 1 if allele == site['h2'] else -1
            if hap >= 0:
                molecule_events[ri][side, events[vi]].add(hap)
        votes = [0, 0]
        bridges, ambiguous = [], []
        for ri, evidence in sorted(molecule_events.items()):
            counts = [[0, 0], [0, 0]]
            conflicts = 0
            for (side, event), haps in evidence.items():
                if len(haps) == 1:
                    counts[side][next(iter(haps))] += 1
                else:
                    conflicts += 1
            if not all(sum(c) for c in counts):
                continue
            # Keep inconsistent events visible. This is a diagnostic, not a
            # license for a majority of correlated calls to override them.
            if conflicts or any(min(c) > 0 for c in counts):
                ambiguous.append(dict(molecule=reads[ri]['molecule'], counts=counts,
                                      conflicting_events=conflicts))
                continue
            flip = (counts[0][1] > 0) != (counts[1][1] > 0)
            votes[int(flip)] += 1
            bridges.append(dict(molecule=reads[ri]['molecule'], flip=flip,
                                original_hp=reads[ri]['hp'], counts=counts))
        result[channel] = dict(same=votes[0], flip=votes[1], bridges=bridges,
                               ambiguous_molecules=ambiguous)
    return result


def prepare_gap(prefix):
    evidence_path = Path(str(prefix) + '.evidence.tsv')
    reads, molecules = {}, set()
    sites, observations = {}, []
    gap = None
    schema_seen = False
    with evidence_path.open() as stream:
        for row in csv.reader(stream, delimiter='\t'):
            if not row:
                continue
            if row[0] == 'SCHEMA':
                if schema_seen or row[1] != '1':
                    raise ValueError(f'Unsupported or duplicate schema in {evidence_path}')
                schema_seen = True
            if row[0] == 'GAP':
                gap = list(map(int, row[1:]))
            elif row[0] == 'SITE':
                sites[int(row[1])] = dict(pos=int(row[2]), type=int(row[3]), length=int(row[4]),
                                         ps=int(row[10]), h1=int(row[11]), h2=int(row[12]))
            elif row[0] == 'OBS':
                observations.append(tuple(map(int, row[1:])))
            elif row[0] == 'READ':
                ri, source, name = int(row[1]), int(row[2]), row[3]
                molecule = (source, name)
                if molecule in molecules or ri in reads:
                    raise ValueError(f'Duplicate molecule or read index in {evidence_path}')
                molecules.add(molecule)
                reads[ri] = dict(molecule=f'{source}:{name}', mapq=int(row[6]),
                                 hp=int(row[8]), ps=int(row[9]), skipped=bool(int(row[7])),
                                 flag=int(row[10]))
    if not schema_seen:
        raise ValueError(f'Missing schema in {evidence_path}')
    if gap is None:
        raise ValueError(f'Missing gap in {evidence_path}')
    raw_votes, contributors = {}, collections.defaultdict(list)
    seen = set()
    with Path(str(prefix) + '.reads.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            view, tier, ri, hp, ps, skipped = (int(row[k]) for k in
                ('VIEW', 'TIER', 'READ', 'HP', 'PS', 'SKIPPED'))
            identity = view, tier, ri
            if identity in seen:
                raise ValueError(f'Duplicate proposal molecule in {prefix}')
            seen.add(identity)
            original = reads[ri]
            if skipped or hp not in (1, 2) or ps < 0 or original['hp'] not in (1, 2):
                continue
            side = 0 if original['ps'] == gap[1] else 1 if original['ps'] == gap[2] else -1
            if side < 0:
                continue
            key = view, tier, ps
            cell = 2 * (original['hp'] - 1) + hp - 1
            raw_votes.setdefault(key, [[0] * 4, [0] * 4])[side][cell] += 1
            contributors[key].append((original['molecule'], side, cell, original['mapq']))
    expected = {}
    legacy = {}
    with Path(str(prefix) + '.votes.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            key = tuple(int(row[k]) for k in ('VIEW', 'TIER', 'PS'))
            if key in expected:
                raise ValueError(f'Duplicate vote matrix in {prefix}')
            expected[key] = [[int(row[side + cell]) for cell in ('11', '12', '21', '22')]
                             for side in ('L', 'R')]
            legacy[key[:2]] = dict(joined=bool(int(row['JOINED'])), flip=bool(int(row['FLIP'])))
    if expected != raw_votes:
        raise ValueError(f'Exported votes differ from molecule replay: {prefix}')
    paths = [Path(str(prefix) + suffix) for suffix in ('.evidence.tsv', '.reads.tsv', '.votes.tsv')]
    manifest = prefix.parent / f'tid{gap[0]}.manifest.tsv'
    if manifest.exists():
        paths.append(manifest)
    return dict(snapshot=str(prefix.parent), gap=gap,
                proposals=[dict(view=k[0], tier=k[1], ps=k[2], molecules=contributors[k])
                           for k in sorted(raw_votes)],
                direct_observations=direct_observation_support(gap, sites, reads, observations),
                legacy_trials=[dict(view=k[0], tier=k[1], **v) for k, v in sorted(legacy.items())],
                input_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def haplotype_margin(votes, flip):
    costs = []
    for left_flip in (False, True):
        right_flip = left_flip ^ flip
        margins = []
        cost = 0
        for v, swapped in zip(votes, (left_flip, right_flip)):
            margins.extend((v[1] - v[0], v[2] - v[3]) if swapped else
                           (v[0] - v[1], v[3] - v[2]))
            cost += v[0] + v[3] if swapped else v[1] + v[2]
        costs.append((cost, -min(margins)))
    return -min(costs)[1]


def score_gap(features, margin=2.0, weighting='uniform', min_haplotype_margin=0.0):
    proposals = []
    for p in features['proposals']:
        votes, contributions = [[0.0] * 4, [0.0] * 4], []
        for name, side, cell, q in p['molecules']:
            weight = 1.0 if weighting == 'uniform' else (0.0 if q == 255 or q < 0 else 1.0 - 10.0 ** (-q / 10.0))
            votes[side][cell] += weight
            contributions.append((name, side, cell, weight))
        summary = summarize(votes, contributions, margin)
        summary['raw_decision'] = summary['decision']
        summary['haplotype_margin'] = haplotype_margin(votes, summary['flip_cost'] < summary['same_cost'])
        if summary['decision'] in ('SAME', 'FLIP') and summary['haplotype_margin'] + 1e-9 < min_haplotype_margin:
            summary['decision'] = 'INSUFFICIENT'
        proposals.append(dict(view=p['view'], tier=p['tier'], ps=p['ps'], **summary))
    directions = {p['decision'] for p in proposals if p['decision'] in ('SAME', 'FLIP')}
    # Alternative tiers and graph/BAM views share molecules: never add their scores.
    decision = ('CONFLICTING' if len(directions) > 1 else next(iter(directions))
                if directions else 'CONFLICTING' if any(p['decision'] == 'CONFLICTING' for p in proposals)
                else 'INSUFFICIENT')
    return dict(features, decision=decision, proposals=proposals)


def replay(prefix, margin=2.0, weighting='uniform'):
    return score_gap(prepare_gap(prefix), margin, weighting)


def load_features(cache, paths):
    data = json.loads(cache.read_text())
    if data.get('schema') != FEATURE_SCHEMA:
        raise ValueError(f'Unsupported feature cache: {cache}')
    features = data['gaps']
    cached_paths = {str(Path(f['snapshot']) / name) for f in features
                    for name in f['input_sha256'] if name.endswith('.evidence.tsv')}
    if cached_paths != {str(p) for p in paths}:
        raise ValueError(f'Feature cache targets changed; use a new cache path: {cache}')
    for f in features:
        for name, digest in f['input_sha256'].items():
            path = Path(f['snapshot']) / name
            if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise ValueError(f'Frozen input changed; use a new feature cache: {path}')
    return features


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--feature-cache', type=Path, help='Reuse verified scoring features across parameter trials')
    parser.add_argument('--min-margin', type=float, default=2.0)
    parser.add_argument('--weighting', choices=('uniform', 'mapq'), default='uniform')
    parser.add_argument('--min-haplotype-margin', type=float, default=2.0,
                        help='Required net support from each original haplotype on both flanks')
    args = parser.parse_args()
    if not math.isfinite(args.min_margin) or args.min_margin <= 0:
        parser.error('--min-margin must be finite and positive')
    if not math.isfinite(args.min_haplotype_margin) or args.min_haplotype_margin < 0:
        parser.error('--min-haplotype-margin must be finite and nonnegative')
    paths = sorted(args.audit.resolve().rglob('*.evidence.tsv'))
    if not paths:
        parser.error('no frozen gap evidence found')
    started = time.monotonic()
    cache_reused = args.feature_cache is not None and args.feature_cache.exists()
    if cache_reused:
        features = load_features(args.feature_cache, paths)
    else:
        features = [prepare_gap(Path(str(p)[:-len('.evidence.tsv')])) for p in paths]
        if args.feature_cache:
            args.feature_cache.parent.mkdir(parents=True, exist_ok=True)
            with args.feature_cache.open('x') as stream:
                json.dump(dict(schema=FEATURE_SCHEMA, gaps=features), stream)
                stream.write('\n')
    prepared_at = time.monotonic()
    results = [score_gap(f, args.min_margin, args.weighting, args.min_haplotype_margin) for f in features]
    graph_edges = collections.defaultdict(list)
    for r in results:
        if r['decision'] in ('SAME', 'FLIP'):
            graph_edges[(r['snapshot'], r['gap'][0])].append(
                (r['gap'][1], r['gap'][2], r['decision'] == 'FLIP'))
    components = [dict(snapshot=k[0], tid=k[1], components=resolve_components(edges))
                  for k, edges in sorted(graph_edges.items())]
    result = dict(schema=1, seconds=time.monotonic() - started,
                  feature_cache_reused=cache_reused, preparation_seconds=prepared_at - started,
                  scoring_seconds=time.monotonic() - prepared_at,
                  tentative_components=components,
                  weighting=args.weighting, min_margin=args.min_margin,
                  min_haplotype_margin=args.min_haplotype_margin, gaps=results,
                  counts=dict(collections.Counter(r['decision'] for r in results)),
                  limitation='Diagnostic proposal consistency only. Shared-source votes are not independent. '
                  'No event-level error calibration, new phasing decisions, or truth-based acceptance is performed.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('gaps', 'tentative_components')}))


if __name__ == '__main__':
    main()
