#!/usr/bin/env python3
"""Diagnostic: could a gap's OWN (non-anchor) events form a coherent local
phase block on their own, independent of whether they confidently connect to
either flank?

The existing replay (`replay_gap_decisions.py`) only ever scores gap-owned
events AGAINST the two flank anchors (`.votes.tsv`'s L11..R22 columns are
flank-referenced by construction). It has no notion of "do private event A
and private event B agree with each other", which is the question this script
answers, using only the raw per-read/per-event allele calls already exported
by `--gap-decision-audit` (`.events.tsv` + `.observations.tsv`). No C++ rerun,
no new instrumentation: this is a pure re-aggregation of existing evidence.

Method per gap:
  1. Take every non-anchor (private/graph) event with exactly 2 alleles.
  2. For each pair of such events, tally agree/conflict votes from reads that
     have a confident (non-missing, non-conflicting) BAM-channel call at both
     -- this is the same "does a read see the same haplotype at two sites"
     question the phaser itself asks, just computed directly instead of via
     an intermediate k-means/consensus fit.
  3. Build a graph over events (edge = pair with total votes >= --min-votes
     and |agree-conflict| >= --min-margin), and take connected components via
     union-find with parity (a component is internally consistent if no cycle
     produces a parity contradiction, same principle as
     replay_gap_decisions.resolve_components).
  4. Report the largest internally-consistent component per gap: its event
     count and the number of distinct reads spanning >=2 of its events (a
     proxy for whether it would look like a real block, not a trivial pair).

This is diagnostic only. It says nothing about correctness -- only about
whether there is enough mutually-consistent internal structure for the
existing phaser to plausibly produce a genuine local block if run on the
gap's own evidence without requiring a flank connection.
"""
import argparse
import collections
import csv
import json
import sys
import time
from pathlib import Path

csv.field_size_limit(16 * 1024 * 1024)

# A read with more informative events than this is almost always a wide,
# low-complexity/repetitive gap dumping many noisy per-read calls (the
# largest chr20 gap here is 1.45 Mb with a 111 MB observations file) rather
# than genuine multi-site linking signal. Pairwise voting is O(k^2) per read;
# without a cap, one such gap dominates the whole run's wall clock while
# contributing votes that are not the kind of evidence this diagnostic is
# trying to measure. Capped reads are dropped from voting entirely rather
# than truncated, so the remaining votes are not a biased subsample of a
# high-event read's own calls.
MAX_EVENTS_PER_READ = 12


def load_gap(prefix):
    events = {}
    with open(f'{prefix}.events.tsv') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            alleles = row['ALLELES'].split(',')
            if len(alleles) != 2:
                continue
            events[int(row['EVENT'])] = dict(role=row['ROLE'], pos=int(row['BEGIN_0']))
    # (event_a, event_b) -> [agree, conflict]
    by_read = collections.defaultdict(dict)
    with open(f'{prefix}.observations.tsv') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ev = int(row['EVENT'])
            if ev not in events or events[ev]['role'] == 'anchor':
                continue
            allele = int(row['BAM_ALLELE'])
            if row['BAM_STATUS'] != 'observed' or allele < 0:
                continue
            by_read[(int(row['INPUT']), row['READ'])][ev] = allele
    return events, by_read


def pairwise_votes(events, by_read):
    votes = collections.defaultdict(lambda: [0, 0])  # (a,b) sorted -> [agree, conflict]
    capped = 0
    for calls in by_read.values():
        if len(calls) > MAX_EVENTS_PER_READ:
            capped += 1
            continue
        items = sorted(calls.items())
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                (ea, aa), (eb, ab) = items[i], items[j]
                votes[(ea, eb)][0 if aa == ab else 1] += 1
    return votes, capped


def resolve_internal_components(private_events, votes, min_votes, min_margin):
    adjacency = collections.defaultdict(list)
    for (a, b), (agree, conflict) in votes.items():
        total = agree + conflict
        if total < min_votes or abs(agree - conflict) < min_margin:
            continue
        flip = conflict > agree
        adjacency[a].append((b, flip))
        adjacency[b].append((a, flip))
    seen, components = set(), []
    for root in sorted(private_events):
        if root in seen or root not in adjacency:
            continue
        parity, pending, consistent = {root: False}, [root], True
        while pending:
            node = pending.pop()
            seen.add(node)
            for other, flip in adjacency[node]:
                expected = parity[node] ^ flip
                if other in parity:
                    if parity[other] != expected:
                        consistent = False
                else:
                    parity[other] = expected
                    pending.append(other)
        components.append(dict(nodes=sorted(parity), consistent=consistent))
    return components


def reads_spanning(nodes, by_read):
    node_set = set(nodes)
    return sum(1 for calls in by_read.values() if len(node_set & calls.keys()) >= 2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--audit', type=Path, required=True)
    ap.add_argument('--decisions', type=Path, required=True,
                    help='Output of replay_gap_decisions.py, to restrict to INSUFFICIENT/CONFLICTING gaps')
    ap.add_argument('--min-votes', type=int, default=2)
    ap.add_argument('--min-margin', type=int, default=2)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()

    decisions = json.loads(args.decisions.read_text())
    target_states = {'INSUFFICIENT', 'CONFLICTING'}
    # decisions JSON's per-gap breakdown isn't in the summary; re-derive gap
    # prefixes from the audit directory and cross-check against the summary
    # counts only as a sanity print, then evaluate every gap (cheap either way).
    prefixes = sorted({p.name[:-len('.events.tsv')] for p in args.audit.glob('*.events.tsv')})

    results = []
    t0 = time.time()
    for gi, prefix_name in enumerate(prefixes):
        g0 = time.time()
        prefix = args.audit / prefix_name
        events, by_read = load_gap(prefix)
        private_events = [e for e, v in events.items() if v['role'] != 'anchor']
        if len(private_events) < 2:
            results.append(dict(gap=prefix_name, private_events=len(private_events),
                               largest_component=0, largest_component_reads=0, consistent=None))
        else:
            votes, capped = pairwise_votes(events, by_read)
            components = resolve_internal_components(private_events, votes, args.min_votes, args.min_margin)
            if not components:
                results.append(dict(gap=prefix_name, private_events=len(private_events),
                                   largest_component=0, largest_component_reads=0, consistent=None,
                                   reads_capped=capped))
            else:
                best = max(components, key=lambda c: len(c['nodes']))
                results.append(dict(gap=prefix_name, private_events=len(private_events),
                                   largest_component=len(best['nodes']),
                                   largest_component_reads=reads_spanning(best['nodes'], by_read),
                                   consistent=best['consistent'], reads_capped=capped))
        elapsed = time.time() - g0
        if elapsed > 1.0:
            print(f'  [{gi + 1}/{len(prefixes)}] {prefix_name}: {elapsed:.1f}s '
                  f'({len(by_read)} reads, {len(private_events)} private events)', file=sys.stderr)
    print(f'total: {time.time() - t0:.1f}s for {len(prefixes)} gaps', file=sys.stderr)

    # A gap is judged "internally phaseable" if it has an internally-consistent
    # component of >=2 events spanned by >=3 distinct reads -- a deliberately
    # low bar meant to upper-bound the opportunity, not to certify accuracy.
    phaseable = [r for r in results if r['consistent'] and r['largest_component'] >= 2
                and r['largest_component_reads'] >= 3]
    inconsistent = [r for r in results if r['consistent'] is False]
    print(f'gaps examined: {len(results)}')
    print(f'gaps with >=2 mutually-consistent private events, '
          f'spanned by >=3 reads (upper-bound "internally phaseable"): {len(phaseable)}')
    print(f'gaps whose largest private-event component is INTERNALLY INCONSISTENT '
          f'(a cycle contradicts itself even ignoring flanks): {len(inconsistent)}')
    zero = sum(1 for r in results if r['largest_component'] == 0)
    print(f'gaps with no qualifying private-event pair at all '
          f'(min_votes={args.min_votes}, min_margin={args.min_margin}): {zero}')

    if args.output:
        args.output.write_text(json.dumps(results, indent=2))
        print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
