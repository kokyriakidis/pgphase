#!/usr/bin/env python3
"""Preserve graph HP/PS tags and add only graph-supported hybrid assignments."""

import argparse
from collections import Counter, defaultdict

import pysam


class ParityUnionFind:
    def __init__(self):
        self.parent = {}
        self.parity = {}

    def find(self, item):
        if item not in self.parent:
            self.parent[item] = item
            self.parity[item] = 0
        if self.parent[item] != item:
            parent, parent_parity = self.find(self.parent[item])
            self.parity[item] ^= parent_parity
            self.parent[item] = parent
        return self.parent[item], self.parity[item]

    def union(self, left, right, relation):
        left_root, left_parity = self.find(left)
        right_root, right_parity = self.find(right)
        if left_root == right_root:
            return (left_parity ^ right_parity) == relation
        self.parent[right_root] = left_root
        self.parity[right_root] = left_parity ^ right_parity ^ relation
        return True


def merge_graph_links(accepted, max_distance):
    merged_graph = ParityUnionFind()
    conflicts = 0
    merge_edges = 0
    distance_rejected = 0
    report = []
    for hybrid_ps, links in accepted.items():
        anchor_ps, anchor_parity, _ = links[0]
        for graph_ps, parity, support in links[1:]:
            distance = abs(anchor_ps - graph_ps)
            if max_distance > 0 and distance > max_distance:
                distance_rejected += 1
                report.append((hybrid_ps, anchor_ps, graph_ps,
                               anchor_parity ^ parity, support,
                               distance, False, "distance"))
                continue
            merge_edges += 1
            relation = anchor_parity ^ parity
            merged = merged_graph.union(anchor_ps, graph_ps, relation)
            if not merged:
                conflicts += 1
            report.append((hybrid_ps, anchor_ps, graph_ps, relation,
                           support, distance, merged,
                           "accepted" if merged else "parity_conflict"))
    return merged_graph, merge_edges, conflicts, distance_rejected, report


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
    ap.add_argument("--min-output-phase-set-reads", type=int, default=0)
    ap.add_argument(
        "--merge-graph-phase-sets", action="store_true",
        help=("Merge graph phase sets when one accepted hybrid block has "
              "independently strong links to both."),
    )
    ap.add_argument(
        "--max-graph-bridge-distance", type=int, default=0,
        help=("Maximum distance between graph phase-set anchors for a merge; "
              "0 disables the distance gate."),
    )
    ap.add_argument("--graph-bridge-report")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    if (args.min_shared_reads < 1 or args.min_vote_margin < 0 or
            args.min_output_phase_set_reads < 0 or
            args.max_graph_bridge_distance < 0):
        ap.error("thresholds must be non-negative and shared support at least one read")
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
    hybrid_only_reads = defaultdict(set)
    with pysam.AlignmentFile(args.hybrid_bam, threads=args.threads) as hybrid:
        for record in hybrid.fetch(until_eof=True):
            hybrid_assignment = phase(record)
            graph_assignment = graph_phase.get(record.query_name)
            if hybrid_assignment is None:
                continue
            if graph_assignment is None:
                hybrid_only_reads[hybrid_assignment[1]].add(record.query_name)
                continue
            hybrid_hp, hybrid_ps = hybrid_assignment
            graph_hp, graph_ps = graph_assignment
            parity = 0 if hybrid_hp == graph_hp else 1
            votes[hybrid_ps][(graph_ps, parity)] += 1
            shared_haps[hybrid_ps][(graph_ps, parity)].add(hybrid_hp)

    accepted = {}
    rejected = Counter()
    for hybrid_ps, block_votes in votes.items():
        if not args.merge_graph_phase_sets:
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
                accepted[hybrid_ps] = [(graph_ps, parity, support)]
            continue

        by_graph_ps = defaultdict(Counter)
        for (graph_ps, parity), support in block_votes.items():
            by_graph_ps[graph_ps][parity] += support
        accepted_links = []
        for graph_ps, parity_votes in by_graph_ps.items():
            parity, support = parity_votes.most_common(1)[0]
            runner_up = parity_votes[1 - parity]
            total = sum(parity_votes.values())
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
                accepted_links.append((graph_ps, parity, support))
        if accepted_links:
            accepted_links.sort(key=lambda link: link[2], reverse=True)
            accepted[hybrid_ps] = accepted_links

    if args.merge_graph_phase_sets:
        (merged_graph, merge_edges, merge_conflicts, distance_rejected,
         bridge_report) = merge_graph_links(
             accepted, args.max_graph_bridge_distance)
    else:
        merged_graph = ParityUnionFind()
        merge_edges = 0
        merge_conflicts = 0
        distance_rejected = 0
        bridge_report = []

    if args.graph_bridge_report:
        with open(args.graph_bridge_report, "w") as report:
            report.write("HYBRID_PS\tANCHOR_GRAPH_PS\tGRAPH_PS\tRELATION\t"
                         "SUPPORT\tDISTANCE\tACCEPTED\tREASON\n")
            for row in bridge_report:
                report.write("\t".join(map(str, row)) + "\n")

    def merged_assignment(hp, graph_ps):
        root, parity = merged_graph.find(graph_ps)
        return (3 - hp if parity else hp), root

    phase_set_reads = Counter()
    for hp, graph_ps in graph_phase.values():
        _, root = merged_assignment(hp, graph_ps)
        phase_set_reads[root] += 1
    for hybrid_ps, links in accepted.items():
        graph_ps, _, _ = links[0]
        _, root = merged_assignment(1, graph_ps)
        phase_set_reads[root] += len(hybrid_only_reads[hybrid_ps])

    def output_phase_set(ps):
        return (ps if phase_set_reads[ps] >= args.min_output_phase_set_reads
                else None)

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
                    hp, ps = merged_assignment(hp, ps)
                    ps = output_phase_set(ps)
                    if ps is None:
                        clear_phase(record)
                        unphased += 1
                    else:
                        record.set_tag("HP", hp, value_type="i")
                        record.set_tag("PS", ps, value_type="i")
                        graph_locked += 1
                else:
                    hybrid_assignment = phase(record)
                    links = (accepted.get(hybrid_assignment[1])
                             if hybrid_assignment is not None else None)
                    if not links:
                        clear_phase(record)
                        unphased += 1
                    else:
                        hp, _ = hybrid_assignment
                        graph_ps, parity, _ = links[0]
                        hp = 3 - hp if parity else hp
                        hp, graph_ps = merged_assignment(hp, graph_ps)
                        graph_ps = output_phase_set(graph_ps)
                        if graph_ps is None:
                            clear_phase(record)
                            unphased += 1
                        else:
                            record.set_tag("HP", hp, value_type="i")
                            record.set_tag("PS", graph_ps, value_type="i")
                            hybrid_added += 1
                out.write(record)

    pysam.index(args.output)
    print(f"graph_locked={graph_locked} hybrid_added={hybrid_added} "
          f"unphased={unphased} accepted_blocks={len(accepted)} "
          f"graph_merge_edges={merge_edges} merge_conflicts={merge_conflicts} "
          f"distance_rejected={distance_rejected} "
          f"rejected_blocks={sum(rejected.values())} "
          f"rejected_reasons={dict(rejected)}")


if __name__ == "__main__":
    main()
