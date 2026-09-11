#!/usr/bin/env python3
"""Tests for the AT-synthesis logic in pantree_to_pgphase_catalog.py.

pantree itself is not required: the tests stub the two graph methods the walk
builder calls.  What is actually being checked is the boundary contract pgphase
enforces in graph_site_validation_skip_reason() -- every allele walk has >=2
handles and all alleles share their first and last handle -- across the three
orientation cases of PangenomeGraph.walk_with_variants, whose returned walk
excludes `first`, and excludes `last` except on a + to - strand flip.

Run:  python3 scripts/test_pantree_to_pgphase_catalog.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pantree_to_pgphase_catalog import (  # noqa: E402
    allele_walk_nodes,
    build_allele_traversals,
)
from compare_catalogs import parse_walk  # noqa: E402


def node_recover(node_id):
    """Minimal stand-in for pantree.utils.node_recover."""
    node, direction = node_id.rsplit("_", 1)
    return (">" if direction == "+" else "<") + node


class StubGraph:
    """Serves fixed interiors for (first, last, variant_edges) lookups."""

    def __init__(self, interiors, branch_point):
        self._interiors = interiors
        self.edges = {}
        self._branch_point = branch_point

    def direction(self, node):
        return 1 if node.endswith("_+") else -1

    def walk_with_variants(self, first, last, variant_edges):
        key = (first, last, tuple(variant_edges))
        return list(self._interiors[key])

    def set_edge(self, edge):
        self.edges[edge] = {"branch_point": self._branch_point}


FAILURES = []


def check(condition, label):
    if condition:
        print(f"  ok    {label}")
    else:
        print(f"  FAIL  {label}")
        FAILURES.append(label)


def test_flank_reattachment_forward():
    """+ to +: walk_with_variants drops both ends, so both are re-added."""
    graph = StubGraph({("1_+", "5_+", ()): ["2_+", "3_+"]}, "1_+")
    nodes = allele_walk_nodes(graph, "1_+", "5_+", [])
    check(nodes == ["1_+", "2_+", "3_+", "5_+"], "forward walk re-adds both flanks")


def test_flank_reattachment_strand_flip():
    """+ to -: the '-' last node is already present; do not duplicate it."""
    graph = StubGraph({("1_+", "5_-", ()): ["2_+", "5_-"]}, "1_+")
    nodes = allele_walk_nodes(graph, "1_+", "5_-", [])
    check(nodes == ["1_+", "2_+", "5_-"], "strand-flip walk does not duplicate last")


def test_flank_reattachment_reverse():
    """- to -: both ends dropped again, both re-added."""
    graph = StubGraph({("4_-", "1_-", ()): ["3_-", "2_-"]}, "1_-")
    nodes = allele_walk_nodes(graph, "4_-", "1_-", [])
    check(nodes == ["4_-", "3_-", "2_-", "1_-"], "reverse walk re-adds both flanks")


def test_snp_traversals_share_boundaries():
    """A SNP bubble: ref through 2, alt through 4, shared flanks 1 and 3.

    The tree holds 1->2->3 and 1->4, so the variant edge is the one closing the
    bubble, (4,3); its branch point is 1.
    """
    edge = ("4_+", "3_+")
    graph = StubGraph(
        {
            ("1_+", "3_+", ()): ["2_+"],
            ("1_+", "3_+", (edge,)): ["4_+"],
        },
        branch_point="1_+",
    )
    graph.set_edge(edge)
    result = build_allele_traversals(graph, edge, node_recover)
    check(result is not None, "SNP bubble yields traversals")
    if result is None:
        return
    ref_at, alt_at = result
    check(ref_at == ">1>2>3", f"ref AT is >1>2>3 (got {ref_at})")
    check(alt_at == ">1>4>3", f"alt AT is >1>4>3 (got {alt_at})")

    ref_walk, alt_walk = parse_walk(ref_at), parse_walk(alt_at)
    check(ref_walk is not None and alt_walk is not None, "both AT strings parse")
    check(len(ref_walk) >= 2 and len(alt_walk) >= 2,
          "both walks have >=2 handles (allele_walk_without_boundaries)")
    check(ref_walk[0] == alt_walk[0] and ref_walk[-1] == alt_walk[-1],
          "alleles share first and last handle (incompatible_allele_boundaries)")
    check(ref_walk != alt_walk, "alleles are distinct (too_few_unique_alleles)")


def test_deletion_traversals():
    """A deletion: alt skips the deleted node entirely, flanks still shared."""
    edge = ("1_+", "3_+")
    graph = StubGraph(
        {
            ("1_+", "3_+", ()): ["2_+"],
            ("1_+", "3_+", (edge,)): [],
        },
        branch_point="1_+",
    )
    graph.set_edge(edge)
    result = build_allele_traversals(graph, edge, node_recover)
    check(result is not None, "deletion yields traversals")
    if result is None:
        return
    ref_at, alt_at = result
    check(ref_at == ">1>2>3", f"deletion ref AT is >1>2>3 (got {ref_at})")
    check(alt_at == ">1>3", f"deletion alt AT is >1>3 (got {alt_at})")
    check(len(parse_walk(alt_at)) == 2, "deletion alt walk still has 2 handles")


def test_degenerate_traversal_rejected():
    """Identical ref and alt walks are unusable as a site; reject them here."""
    edge = ("1_+", "3_+")
    graph = StubGraph(
        {
            ("1_+", "3_+", ()): ["2_+"],
            ("1_+", "3_+", (edge,)): ["2_+"],
        },
        branch_point="1_+",
    )
    graph.set_edge(edge)
    check(build_allele_traversals(graph, edge, node_recover) is None,
          "degenerate ref==alt traversal rejected")


def test_reverse_strand_variant_edge():
    """direction(u) == -1: the walk runs from u back to the branch point."""
    edge = ("4_-", "2_-")
    graph = StubGraph(
        {
            ("4_-", "1_-", ()): ["3_-"],
            ("4_-", "1_-", (edge,)): ["2_-"],
        },
        branch_point="1_-",
    )
    graph.set_edge(edge)
    result = build_allele_traversals(graph, edge, node_recover)
    check(result is not None, "reverse-strand edge yields traversals")
    if result is None:
        return
    ref_at, alt_at = result
    check(ref_at == "<4<3<1", f"reverse ref AT is <4<3<1 (got {ref_at})")
    check(alt_at == "<4<2<1", f"reverse alt AT is <4<2<1 (got {alt_at})")
    ref_walk, alt_walk = parse_walk(ref_at), parse_walk(alt_at)
    check(ref_walk[0] == alt_walk[0] and ref_walk[-1] == alt_walk[-1],
          "reverse-strand alleles share flanks")


def main():
    print("pantree_to_pgphase_catalog AT synthesis")
    test_flank_reattachment_forward()
    test_flank_reattachment_strand_flip()
    test_flank_reattachment_reverse()
    test_snp_traversals_share_boundaries()
    test_deletion_traversals()
    test_degenerate_traversal_rejected()
    test_reverse_strand_variant_edge()
    if FAILURES:
        print(f"\n{len(FAILURES)} failure(s)")
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
