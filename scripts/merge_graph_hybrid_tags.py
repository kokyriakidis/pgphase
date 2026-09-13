#!/usr/bin/env python3
"""Preserve graph HP/PS tags and add only graph-supported hybrid assignments."""

import argparse
from collections import Counter, defaultdict

import pysam


def phase(record):
    if not record.has_tag("HP") or not record.has_tag("PS"):
        return None
    hp = record.get_tag("HP")
    ps = record.get_tag("PS")
    return (hp, ps) if hp in (1, 2) and ps >= 0 else None


def clear_phase(record):
    if record.has_tag("HP"):
        record.set_tag("HP", None)
    if record.has_tag("PS"):
        record.set_tag("PS", None)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--graph-bam", required=True)
    ap.add_argument("--hybrid-bam", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-shared-reads", type=int, default=10)
    ap.add_argument("--min-vote-margin", type=int, default=5)
    ap.add_argument("--min-purity", type=float, default=0.9)
    ap.add_argument("--require-both-haplotypes", action="store_true")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    if args.min_shared_reads < 1 or args.min_vote_margin < 0:
        ap.error("read thresholds must be non-negative and support at least one read")
    if not 0.5 <= args.min_purity <= 1.0:
        ap.error("--min-purity must be between 0.5 and 1.0")
    if args.threads < 1:
        ap.error("--threads must be at least 1")

    graph_phase = {}
    with pysam.AlignmentFile(args.graph_bam, "rb", check_sq=False,
                             threads=args.threads) as graph:
        for record in graph.fetch(until_eof=True):
            assignment = phase(record)
            if assignment is not None:
                graph_phase[record.query_name] = assignment

    votes = defaultdict(Counter)
    shared_haps = defaultdict(lambda: defaultdict(set))
    with pysam.AlignmentFile(args.hybrid_bam, threads=args.threads) as hybrid:
        for record in hybrid.fetch(until_eof=True):
            hybrid_assignment = phase(record)
            graph_assignment = graph_phase.get(record.query_name)
            if hybrid_assignment is None or graph_assignment is None:
                continue
            hybrid_hp, hybrid_ps = hybrid_assignment
            graph_hp, graph_ps = graph_assignment
            parity = 0 if hybrid_hp == graph_hp else 1
            votes[hybrid_ps][(graph_ps, parity)] += 1
            shared_haps[hybrid_ps][(graph_ps, parity)].add(hybrid_hp)

    accepted = {}
    rejected = Counter()
    for hybrid_ps, block_votes in votes.items():
        ranked = block_votes.most_common()
        (graph_ps, parity), support = ranked[0]
        runner_up = ranked[1][1] if len(ranked) > 1 else 0
        total = sum(block_votes.values())
        if support < args.min_shared_reads:
            rejected["support"] += 1
        elif support - runner_up < args.min_vote_margin:
            rejected["margin"] += 1
        elif support / total < args.min_purity:
            rejected["purity"] += 1
        elif (args.require_both_haplotypes and
              len(shared_haps[hybrid_ps][(graph_ps, parity)]) < 2):
            rejected["one_haplotype"] += 1
        else:
            accepted[hybrid_ps] = (graph_ps, parity)

    graph_locked = 0
    hybrid_added = 0
    unphased = 0
    with pysam.AlignmentFile(args.hybrid_bam, threads=args.threads) as hybrid:
        with pysam.AlignmentFile(args.output, "wb", header=hybrid.header,
                                 threads=args.threads) as out:
            for record in hybrid.fetch(until_eof=True):
                graph_assignment = graph_phase.get(record.query_name)
                if graph_assignment is not None:
                    hp, ps = graph_assignment
                    record.set_tag("HP", hp, value_type="i")
                    record.set_tag("PS", ps, value_type="i")
                    graph_locked += 1
                else:
                    hybrid_assignment = phase(record)
                    anchor = (accepted.get(hybrid_assignment[1])
                              if hybrid_assignment is not None else None)
                    if anchor is None:
                        clear_phase(record)
                        unphased += 1
                    else:
                        hp, _ = hybrid_assignment
                        graph_ps, parity = anchor
                        record.set_tag("HP", 3 - hp if parity else hp,
                                       value_type="i")
                        record.set_tag("PS", graph_ps, value_type="i")
                        hybrid_added += 1
                out.write(record)

    pysam.index(args.output)
    print(f"graph_locked={graph_locked} hybrid_added={hybrid_added} "
          f"unphased={unphased} accepted_blocks={len(accepted)} "
          f"rejected_blocks={sum(rejected.values())} "
          f"rejected_reasons={dict(rejected)}")


if __name__ == "__main__":
    main()
