#!/usr/bin/env python3
"""Extract linear heterozygous sites private to gaps between graph phase blocks.

The graph-phased VCF defines the authoritative blocks. The output contains only
biallelic heterozygous records from the linear VCF that are absent by exact
POS/REF/ALT match from the graph catalog and lie strictly between merged graph
blocks. Input phasing and PS are cleared so downstream gap phasing cannot inherit
the linear caller's phase. Optionally, a surjected BAM can restrict output to
private-site chains with read overlaps connecting both neighboring graph blocks.
"""

import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import pysam


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_graph_alleles(path, contig):
    alleles = set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[0].split("#")[-1] != contig:
                continue
            for alt in fields[4].split(","):
                alleles.add((int(fields[1]), fields[3], alt))
    return alleles


def graph_phase_blocks(path, contig):
    by_ps = defaultdict(list)
    with pysam.VariantFile(path) as vcf:
        for rec in vcf.fetch(contig):
            if len(rec.samples) != 1:
                raise ValueError("graph VCF must contain exactly one sample")
            sample = rec.samples[0]
            gt = sample.get("GT")
            ps = sample.get("PS")
            if not sample.phased or ps is None or gt is None \
                    or len(gt) != 2 or gt[0] == gt[1]:
                continue
            by_ps[int(ps)].append((rec.start, rec.stop))

    blocks = []
    for positions in by_ps.values():
        blocks.append((min(p[0] for p in positions),
                       max(p[1] for p in positions)))
    blocks.sort()

    merged = []
    for start, end in blocks:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def block_gaps(blocks, min_gap):
    return [(blocks[i][1], blocks[i + 1][0])
            for i in range(len(blocks) - 1)
            if blocks[i + 1][0] - blocks[i][1] >= min_gap]


def in_gap(pos, gaps, gap_i):
    while gap_i < len(gaps) and pos >= gaps[gap_i][1]:
        gap_i += 1
    return gap_i < len(gaps) and gaps[gap_i][0] <= pos < gaps[gap_i][1], gap_i


def bam_contig_name(bam, contig):
    if contig in bam.references:
        return contig
    matches = [name for name in bam.references
               if name.split("#")[-1] == contig]
    if len(matches) != 1:
        raise ValueError(
            f"cannot uniquely map VCF contig {contig!r} to BAM references")
    return matches[0]


def bridge_path(bam, bam_contig, gap, records, min_reads, min_mapq):
    """Return records on the strongest supported left-to-right overlap path."""
    if not records:
        return []

    # Positions are 0-based. Put the anchors just inside the neighboring graph
    # blocks; they represent the evidence that a private chain must reach.
    nodes = [max(0, gap[0] - 1)] + [rec.start for rec in records] + [gap[1]]
    edge_support = defaultdict(int)
    for read in bam.fetch(bam_contig, nodes[0], nodes[-1] + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary \
                or read.mapping_quality < min_mapq:
            continue
        covered = [i for i, pos in enumerate(nodes)
                   if read.reference_start <= pos < read.reference_end]
        for oi, left in enumerate(covered):
            for right in covered[oi + 1:]:
                edge_support[(left, right)] += 1

    last = len(nodes) - 1
    # Maximize the weakest edge first, then prefer fewer private anchors. A
    # direct block-to-block edge does not establish that a private site helps.
    scores = [None] * len(nodes)
    previous = [-1] * len(nodes)
    scores[0] = (10**9, 0, 0)
    for right in range(1, len(nodes)):
        for left in range(right):
            support = edge_support.get((left, right), 0)
            if support < min_reads or scores[left] is None:
                continue
            if left == 0 and right == last:
                continue
            score = (min(scores[left][0], support),
                     scores[left][1] - 1,
                     scores[left][2] + support)
            if scores[right] is None or score > scores[right]:
                scores[right] = score
                previous[right] = left

    if scores[last] is None:
        return []
    path = []
    node = last
    while previous[node] >= 0:
        if node != last:
            path.append(records[node - 1])
        node = previous[node]
    path.reverse()
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--graph-phased-vcf", required=True)
    ap.add_argument("--graph-sites", required=True)
    ap.add_argument("--linear-vcf", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--gaps-bed")
    ap.add_argument("--bam", help="surjected BAM/CRAM for bridge-only filtering")
    ap.add_argument("--min-bridge-reads", type=int, default=2)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-gq", type=int, default=10)
    ap.add_argument("--clean-snps-only", action="store_true",
                    help="accept only CLEAN single-nucleotide records")
    ap.add_argument("--exclude-graph-positions", action="store_true",
                    help="reject records at any position present in the graph catalog")
    ap.add_argument("--min-vaf", type=float, default=0.0,
                    help="minimum alternate allele fraction")
    ap.add_argument("--max-vaf", type=float, default=1.0,
                    help="maximum alternate allele fraction")
    ap.add_argument("--min-gap", type=int, default=1)
    args = ap.parse_args()
    if args.min_bridge_reads < 1:
        ap.error("--min-bridge-reads must be at least 1")
    if args.min_mapq < 0:
        ap.error("--min-mapq must be non-negative")
    if not 0.0 <= args.min_vaf <= args.max_vaf <= 1.0:
        ap.error("require 0 <= --min-vaf <= --max-vaf <= 1")

    blocks = graph_phase_blocks(args.graph_phased_vcf, args.contig)
    gaps = block_gaps(blocks, args.min_gap)
    graph_alleles = read_graph_alleles(args.graph_sites, args.contig)
    graph_positions = {key[0] for key in graph_alleles}

    if args.gaps_bed:
        with Path(args.gaps_bed).open("w") as out:
            for start, end in gaps:
                out.write(f"{args.contig}\t{start}\t{end}\n")

    selected = 0
    rejected_graph = 0
    rejected_outside = 0
    gap_i = 0
    gap_records = [[] for _ in gaps]
    with pysam.VariantFile(args.linear_vcf) as vin:
        if len(vin.header.samples) != 1:
            ap.error("linear VCF must contain exactly one sample")
        header = vin.header.copy()
        if "PS" not in header.formats:
            header.formats.add("PS", 1, "Integer", "Phase set identifier")
        for rec in vin.fetch(args.contig):
            sample = rec.samples[0]
            gt = sample.get("GT")
            if gt is None or len(gt) != 2 or set(gt) != {0, 1} \
                    or len(rec.alts or ()) != 1:
                continue
            if rec.filter.keys() not in ([], ["PASS"]):
                continue
            gq = sample.get("GQ")
            if gq is None or gq < args.min_gq:
                continue
            if args.clean_snps_only and (len(rec.ref) != 1
                                         or len(rec.alts[0]) != 1
                                         or "CLEAN" not in rec.info):
                continue
            if args.min_vaf > 0.0 or args.max_vaf < 1.0:
                ad = sample.get("AD")
                if ad is None or len(ad) < 2 or ad[0] + ad[1] == 0:
                    continue
                vaf = ad[1] / (ad[0] + ad[1])
                if vaf < args.min_vaf or vaf > args.max_vaf:
                    continue
            key = (rec.pos, rec.ref, rec.alts[0])
            if key in graph_alleles:
                rejected_graph += 1
                continue
            if args.exclude_graph_positions and rec.pos in graph_positions:
                rejected_graph += 1
                continue
            inside, gap_i = in_gap(rec.start, gaps, gap_i)
            if not inside:
                rejected_outside += 1
                continue
            sample["GT"] = tuple(gt)
            sample.phased = False
            sample["PS"] = None
            gap_records[gap_i].append(rec.copy())

        bridge_gaps = 0
        if args.bam:
            with pysam.AlignmentFile(args.bam) as bam:
                bam_contig = bam_contig_name(bam, args.contig)
                filtered = []
                for gap, records in zip(gaps, gap_records):
                    path = bridge_path(bam, bam_contig, gap, records,
                                       args.min_bridge_reads, args.min_mapq)
                    if path:
                        bridge_gaps += 1
                    filtered.append(path)
                gap_records = filtered

        with pysam.VariantFile(args.output, "w", header=header) as out:
            for records in gap_records:
                for rec in records:
                    out.write(rec)
                    selected += 1

    print(f"merged_graph_blocks={len(blocks)} gaps={len(gaps)} private_gap_sites={selected} "
          f"bridge_gaps={bridge_gaps if args.bam else 'NA'} "
          f"graph_matches={rejected_graph} outside_gaps={rejected_outside}")


if __name__ == "__main__":
    main()
