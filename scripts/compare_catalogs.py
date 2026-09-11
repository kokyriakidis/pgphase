#!/usr/bin/env python3
"""Differential comparison of two pgphase site catalogs, stratified by region.

Built to answer one question: does a pantree reference-tree catalog contribute
phasing anchors that the `vg deconstruct` superbubble catalog does not, and are
those anchors concentrated in the segmental duplications where pgphase's residual
switch errors live (CHECKPOINT.md, "haplotype ambiguity in segmental duplications")?

Comparing the two catalogs by POS/REF/ALT does not work: deconstruct draws allele
boundaries at snarl endpoints while pantree draws them at the reference-tree branch
point, so the same underlying variation gets different coordinates, different
flanking context, and different allele strings.  The one representation-independent
unit both catalogs agree on is the **variant edge**: the graph edge an alternate
allele traverses that its reference allele does not.  That is pantree's native
definition, and it is recoverable from any `AT` field by diffing the alt walk
against the ref walk.  This script decomposes both catalogs to canonical variant
edges and compares those sets.

Usage:
    python3 compare_catalogs.py \
        --a chr20.deconstruct.sites.vcf.gz --a-label deconstruct \
        --b chr20.pantree.sites.vcf.gz     --b-label pantree \
        --bed data/annotations/hg002v1.1.segdups.bed --bed-label segdup \
        --contig chr20 \
        --out-prefix out/chr20_catalog_cmp

Outputs:
    <prefix>.summary.txt    human-readable counts and per-Mb densities
    <prefix>.strata.tsv     machine-readable counts by region stratum and class
    <prefix>.b_only.tsv     every variant edge unique to catalog B
    <prefix>.windows.tsv    per-window edge counts (with --window)

Only needs the Python standard library.
"""

import argparse
import bisect
import gzip
import os
import sys
from collections import defaultdict


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalize_chrom(chrom):
    """Pangenome path name -> plain contig, matching GraphSite::ref_contig."""
    return chrom.rsplit("#", 1)[-1]


# --- walk parsing ----------------------------------------------------------

def parse_walk(text):
    """'>1>2<3' -> [('1', False), ('2', False), ('3', True)]; None if unparseable."""
    if not text or text in (".", "*"):
        return None
    steps = []
    i = 0
    n = len(text)
    while i < n:
        marker = text[i]
        if marker not in ("<", ">"):
            return None
        j = i + 1
        while j < n and text[j] not in ("<", ">"):
            j += 1
        node = text[i + 1:j]
        if not node:
            return None
        steps.append((node, marker == "<"))
        i = j
    return steps or None


def canonical_edge(step_a, step_b):
    """Orientation-independent identifier for the edge step_a -> step_b.

    A walk and its reverse complement traverse the same edges, so canonicalize
    by taking the lexicographically smaller of the two renderings.
    """
    def render(step):
        node, rev = step
        return ("<" if rev else ">") + node

    forward = render(step_a) + render(step_b)
    backward = render((step_b[0], not step_b[1])) + render((step_a[0], not step_a[1]))
    return min(forward, backward)


def walk_edges(steps):
    return {canonical_edge(steps[i], steps[i + 1]) for i in range(len(steps) - 1)}


def classify_allele(ref, alt):
    """SNP / INS / DEL / MNP / OTHER from the VCF allele pair."""
    if ref in (".", "", "N") or alt in (".", "", "*"):
        return "NONREF"
    if len(ref) == len(alt):
        return "SNP" if len(ref) == 1 else "MNP"
    if alt.startswith("<"):
        return "OTHER"
    return "INS" if len(alt) > len(ref) else "DEL"


# --- region index ----------------------------------------------------------

class RegionIndex:
    """Sorted-interval membership test over a BED, per contig."""

    def __init__(self):
        self.starts = defaultdict(list)
        self.ends = defaultdict(list)

    @classmethod
    def from_bed(cls, path):
        index = cls()
        raw = defaultdict(list)
        with open_maybe_gz(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("track", "browser", "#")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    raw[fields[0]].append((int(fields[1]), int(fields[2])))
                except ValueError:
                    continue
        # Merge so overlapping BED entries do not distort the covered-bp total.
        for contig, spans in raw.items():
            spans.sort()
            merged = []
            for start, end in spans:
                if merged and start <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            index.starts[contig] = [s for s, _ in merged]
            index.ends[contig] = [e for _, e in merged]
        return index

    def contains(self, contig, pos_1based):
        starts = self.starts.get(contig)
        if not starts:
            return False
        pos = pos_1based - 1
        i = bisect.bisect_right(starts, pos) - 1
        return i >= 0 and pos < self.ends[contig][i]

    def covered_bp(self, contig):
        starts = self.starts.get(contig)
        if not starts:
            return 0
        return sum(e - s for s, e in zip(starts, self.ends[contig]))


# --- catalog loading -------------------------------------------------------

class EdgeRecord:
    __slots__ = ("pos", "site_id", "klass", "af", "ref", "alt")

    def __init__(self, pos, site_id, klass, af, ref, alt):
        self.pos = pos
        self.site_id = site_id
        self.klass = klass
        self.af = af
        self.ref = ref
        self.alt = alt


def parse_info(info):
    out = {}
    for item in info.split(";"):
        if not item:
            continue
        key, _, value = item.partition("=")
        out[key] = value
    return out


def load_catalog(path, contig_filter):
    """Decompose a sites VCF into {canonical variant edge -> EdgeRecord}.

    A site's alternate allele contributes every edge on its walk that the
    reference allele's walk does not contain.  Sites without a usable AT field
    are counted separately and reported -- they cannot be compared, and pgphase
    cannot genotype them either.
    """
    edges = {}
    stats = defaultdict(int)
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            contig = normalize_chrom(fields[0])
            if contig_filter and contig != contig_filter:
                continue
            stats["sites"] += 1
            pos = int(fields[1])
            site_id = fields[2]
            ref = fields[3]
            alts = fields[4].split(",")
            info = parse_info(fields[7])

            at_raw = info.get("AT") or info.get("UT") or info.get("LVAT")
            if not at_raw:
                stats["sites_without_at"] += 1
                continue
            walk_texts = at_raw.split(",")
            if len(walk_texts) != len(alts) + 1:
                stats["sites_at_count_mismatch"] += 1
                continue
            walks = [parse_walk(text) for text in walk_texts]
            if walks[0] is None:
                stats["sites_unparseable_ref_walk"] += 1
                continue

            afs = info.get("AF", "").split(",")
            ref_edges = walk_edges(walks[0])
            for allele_i, alt in enumerate(alts, start=1):
                alt_walk = walks[allele_i]
                if alt_walk is None:
                    stats["alleles_unparseable_walk"] += 1
                    continue
                novel = walk_edges(alt_walk) - ref_edges
                if not novel:
                    stats["alleles_without_novel_edge"] += 1
                    continue
                try:
                    af = float(afs[allele_i - 1])
                except (IndexError, ValueError):
                    af = None
                klass = info.get("VT") or classify_allele(ref, alt)
                stats["alleles"] += 1
                for edge in novel:
                    stats["edges_seen"] += 1
                    # First site to claim an edge owns it.  Overlapping snarls in
                    # a deconstruct catalog routinely re-derive the same edge.
                    if edge not in edges:
                        edges[edge] = EdgeRecord(pos, site_id, klass, af, ref, alt)
    stats["edges_unique"] = len(edges)
    return edges, stats


# --- reporting -------------------------------------------------------------

def stratum_of(record, contig, region, bed_label):
    if region is None:
        return "all"
    return bed_label if region.contains(contig, record.pos) else f"non_{bed_label}"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--a", required=True, help="Catalog A sites VCF (baseline)")
    parser.add_argument("--b", required=True, help="Catalog B sites VCF (candidate)")
    parser.add_argument("--a-label", default="A")
    parser.add_argument("--b-label", default="B")
    parser.add_argument("--bed", help="BED of regions to stratify by (e.g. segdups)")
    parser.add_argument("--bed-label", default="segdup")
    parser.add_argument("--contig", help="Restrict to this contig (e.g. chr20)")
    parser.add_argument("--contig-length", type=int,
                        help="Contig length in bp; enables per-Mb densities")
    parser.add_argument("--window", type=int,
                        help="Emit per-window edge counts at this window size (bp)")
    parser.add_argument("--out-prefix", required=True, help="Output file prefix")
    args = parser.parse_args()

    out_dir = os.path.dirname(os.path.abspath(args.out_prefix))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    region = RegionIndex.from_bed(args.bed) if args.bed else None
    bed_label = args.bed_label

    print(f"[compare] loading {args.a_label}: {args.a}", file=sys.stderr)
    a_edges, a_stats = load_catalog(args.a, args.contig)
    print(f"[compare] loading {args.b_label}: {args.b}", file=sys.stderr)
    b_edges, b_stats = load_catalog(args.b, args.contig)

    a_keys = set(a_edges)
    b_keys = set(b_edges)
    shared = a_keys & b_keys
    a_only = a_keys - b_keys
    b_only = b_keys - a_keys

    strata = ["all"] if region is None else [bed_label, f"non_{bed_label}"]
    classes_seen = set()
    # counts[stratum][class][bucket] -> int
    counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    def tally(edge_set, edges, bucket):
        for edge in edge_set:
            record = edges[edge]
            contig = args.contig or "?"
            stratum = stratum_of(record, contig, region, bed_label)
            classes_seen.add(record.klass)
            counts[stratum][record.klass][bucket] += 1
            counts[stratum]["ALL"][bucket] += 1

    tally(shared, a_edges, "shared")
    tally(a_only, a_edges, "a_only")
    tally(b_only, b_edges, "b_only")

    # --- strata.tsv
    with open(f"{args.out_prefix}.strata.tsv", "w") as fh:
        fh.write("stratum\tclass\tshared\t%s_only\t%s_only\t%s_total\t%s_total\n"
                 % (args.a_label, args.b_label, args.a_label, args.b_label))
        for stratum in strata:
            for klass in ["ALL"] + sorted(classes_seen):
                row = counts[stratum][klass]
                a_total = row["shared"] + row["a_only"]
                b_total = row["shared"] + row["b_only"]
                if a_total == 0 and b_total == 0:
                    continue
                fh.write(f"{stratum}\t{klass}\t{row['shared']}\t{row['a_only']}\t"
                         f"{row['b_only']}\t{a_total}\t{b_total}\n")

    # --- b_only.tsv
    with open(f"{args.out_prefix}.b_only.tsv", "w") as fh:
        fh.write("edge\tpos\tsite_id\tclass\taf\tstratum\tref\talt\n")
        for edge in sorted(b_only, key=lambda e: (b_edges[e].pos, e)):
            record = b_edges[edge]
            stratum = stratum_of(record, args.contig or "?", region, bed_label)
            af = "." if record.af is None else f"{record.af:.4g}"
            fh.write(f"{edge}\t{record.pos}\t{record.site_id}\t{record.klass}\t{af}\t"
                     f"{stratum}\t{record.ref[:40]}\t{record.alt[:40]}\n")

    # --- windows.tsv
    if args.window:
        window_counts = defaultdict(lambda: defaultdict(int))
        for bucket, edge_set, edges in (("a_only", a_only, a_edges),
                                        ("b_only", b_only, b_edges),
                                        ("shared", shared, a_edges)):
            for edge in edge_set:
                window = (edges[edge].pos - 1) // args.window
                window_counts[window][bucket] += 1
        with open(f"{args.out_prefix}.windows.tsv", "w") as fh:
            fh.write(f"window_start\twindow_end\tshared\t{args.a_label}_only\t{args.b_label}_only\n")
            for window in sorted(window_counts):
                row = window_counts[window]
                fh.write(f"{window * args.window}\t{(window + 1) * args.window}\t"
                         f"{row['shared']}\t{row['a_only']}\t{row['b_only']}\n")

    # --- summary.txt
    lines = []
    lines.append(f"Catalog comparison: {args.a_label} (A) vs {args.b_label} (B)")
    lines.append(f"  A: {args.a}")
    lines.append(f"  B: {args.b}")
    if args.contig:
        lines.append(f"  contig: {args.contig}")
    if args.bed:
        lines.append(f"  stratifying BED: {args.bed} (label '{bed_label}')")
    lines.append("")
    lines.append("Parse stats (sites -> comparable variant edges)")
    header = f"  {'metric':<32}{args.a_label:>14}{args.b_label:>14}"
    lines.append(header)
    for key in ["sites", "sites_without_at", "sites_at_count_mismatch",
                "sites_unparseable_ref_walk", "alleles", "alleles_unparseable_walk",
                "alleles_without_novel_edge", "edges_unique"]:
        lines.append(f"  {key:<32}{a_stats.get(key, 0):>14}{b_stats.get(key, 0):>14}")
    lines.append("")
    lines.append("Variant-edge overlap")
    lines.append(f"  shared                          {len(shared):>14}")
    lines.append(f"  {args.a_label + ' only':<32}{len(a_only):>14}")
    lines.append(f"  {args.b_label + ' only':<32}{len(b_only):>14}")
    if a_keys:
        lines.append(f"  B adds {len(b_only) / len(a_keys) * 100:.2f}% on top of {args.a_label}")
    lines.append("")

    for stratum in strata:
        lines.append(f"[{stratum}]")
        lines.append(f"  {'class':<10}{'shared':>12}{args.a_label + '_only':>20}"
                     f"{args.b_label + '_only':>20}{'B gain':>10}")
        for klass in ["ALL"] + sorted(classes_seen):
            row = counts[stratum][klass]
            a_total = row["shared"] + row["a_only"]
            b_total = row["shared"] + row["b_only"]
            if a_total == 0 and b_total == 0:
                continue
            gain = f"{(row['b_only'] / a_total * 100):.1f}%" if a_total else "n/a"
            lines.append(f"  {klass:<10}{row['shared']:>12}{row['a_only']:>20}"
                         f"{row['b_only']:>20}{gain:>10}")
        if region is not None and args.contig:
            covered = region.covered_bp(args.contig)
            if stratum == bed_label:
                span = covered
            elif args.contig_length:
                span = max(args.contig_length - covered, 0)
            else:
                span = 0
            if span:
                row = counts[stratum]["ALL"]
                per_mb = (row["shared"] + row["b_only"]) / (span / 1e6)
                add_mb = row["b_only"] / (span / 1e6)
                lines.append(f"  span {span / 1e6:.2f} Mb -> {per_mb:.0f} B edges/Mb "
                             f"({add_mb:.0f}/Mb added over {args.a_label})")
        lines.append("")

    lines.append("Interpretation guide")
    if region is None:
        lines.append("  Re-run with --bed pointing at a segmental-duplication BED. Without the")
        lines.append("  stratification this only says how many edges differ, not whether the")
        lines.append("  added ones land where pgphase's switch errors actually are.")
    else:
        lines.append(f"  The number that matters is '{args.b_label}_only' for class SNP inside")
        lines.append(f"  '{bed_label}', as anchors/Mb. Those are the paralog-distinguishing sites")
        lines.append("  the superbubble catalog cannot represent. A large gain there is the only")
        lines.append("  result that justifies the preprocessing cost; a gain concentrated outside")
        lines.append("  the BED, or in INS/DEL/REP classes, is not.")

    summary = "\n".join(lines) + "\n"
    with open(f"{args.out_prefix}.summary.txt", "w") as fh:
        fh.write(summary)
    print(summary)


if __name__ == "__main__":
    main()
