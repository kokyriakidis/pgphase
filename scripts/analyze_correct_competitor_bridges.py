#!/usr/bin/env python3
"""Find accurate competitor blocks that bridge pgphase graph block breaks."""

import argparse
import bisect
import csv
import gzip
from collections import Counter
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def norm_chrom(chrom):
    return chrom.split("#")[-1]


def read_blocks(path, chrom):
    with path.open() as fh:
        rows = []
        for row in csv.DictReader(fh, delimiter="\t"):
            if norm_chrom(row["chromosome"]) != chrom:
                continue
            rows.append({
                "ps": row["phase_set"],
                "start": int(row["from"]),
                "end": int(row["to"]),
                "variants": int(row["variants"]),
            })
    return sorted(rows, key=lambda row: (row["start"], row["end"]))


def graph_gaps(blocks):
    gaps = []
    furthest = blocks[0]
    for block in blocks[1:]:
        if block["start"] > furthest["end"]:
            gaps.append({
                "left_ps": furthest["ps"],
                "right_ps": block["ps"],
                "start": furthest["end"],
                "end": block["start"],
            })
        if block["end"] > furthest["end"]:
            furthest = block
    return gaps


def read_phase_accuracy(path, chrom):
    values = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if norm_chrom(row["chrom"]) == chrom:
                values[row["phase_set"]] = {
                    "accuracy": float(row["accuracy"]),
                    "reads": int(row["n_reads"]),
                    "discordant": int(row["discordant"]),
                }
    return values


def read_switch_intervals(path, chrom):
    intervals = []
    with path.open() as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and norm_chrom(fields[0]) == chrom:
                intervals.append((int(fields[1]), int(fields[2])))
    return sorted(intervals)


def overlaps(intervals, start, end):
    return any(left < end and right > start for left, right in intervals)


def is_het(sample):
    gt = sample.split(":", 1)[0].replace("|", "/").split("/")
    return len(gt) == 2 and "." not in gt and gt[0] != gt[1]


def vcf_positions(path, chrom, phased_only=False):
    positions = set()
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if norm_chrom(fields[0]) != chrom or len(fields) < 10:
                continue
            if not is_het(fields[9]):
                continue
            if phased_only and "|" not in fields[9].split(":", 1)[0]:
                continue
            positions.add(int(fields[1]))
    return sorted(positions)


def catalog_positions(path, chrom):
    positions = set()
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 2)
            if norm_chrom(fields[0]) == chrom:
                positions.add(int(fields[1]))
    return sorted(positions)


def tsv_positions(path, chrom=None, phased_only=False):
    positions = set()
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row_chrom = row.get("CHROM")
            if chrom and row_chrom and norm_chrom(row_chrom) != chrom:
                continue
            if phased_only and row.get("PHASE_SET") in (None, "", "-1", "0"):
                continue
            positions.add(int(row["POS"]))
    return sorted(positions)


def tsv_label_index(path, label_field, chrom=None):
    rows = []
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row_chrom = row.get("CHROM")
            if chrom and row_chrom and norm_chrom(row_chrom) != chrom:
                continue
            rows.append((int(row["POS"]), row[label_field]))
    rows.sort()
    return [row[0] for row in rows], rows


def labels_between(index, start, end):
    positions, rows = index
    left = bisect.bisect_right(positions, start)
    right = bisect.bisect_left(positions, end)
    return Counter(label for _, label in rows[left:right])


def count_between(positions, start, end):
    return bisect.bisect_left(positions, end) - bisect.bisect_right(positions, start)


def classify(counts):
    if counts["graph_phased_sites"]:
        return "block_stitching"
    if counts["clean_het_snps"] or counts["clean_het_indels"]:
        return "clean_candidate_unphased"
    if counts["repeat_het_indels"]:
        return "repeat_indels_excluded"
    if counts["graph_candidates"]:
        return "other_candidate_unphased"
    if counts["catalog_sites"]:
        return "catalog_site_not_candidate"
    return "catalog_absent"


def parse_tool(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError("tool must be LABEL=DIR")
    label, directory = value.split("=", 1)
    return label, Path(directory)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--root", type=Path, required=True,
                        help="chromosome result directory")
    parser.add_argument("--truth-vcf", type=Path, required=True)
    parser.add_argument("--shared-vcf", type=Path, required=True)
    parser.add_argument("--catalog-vcf", type=Path, required=True)
    parser.add_argument("--competitor", action="append", type=parse_tool,
                        required=True, metavar="LABEL=DIR")
    parser.add_argument("--min-accuracy", type=float, default=0.99)
    parser.add_argument("--min-reads", type=int, default=50)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    chrom = args.chromosome
    graph_blocks = read_blocks(args.root / "graph.blocks.tsv", chrom)
    gaps = graph_gaps(graph_blocks)
    position_sets = {
        "truth_hets": vcf_positions(args.truth_vcf, chrom),
        "shared_hets": vcf_positions(args.shared_vcf, chrom),
        "catalog_sites": catalog_positions(args.catalog_vcf, chrom),
        "graph_candidates": tsv_positions(args.root / "graph/candidates.tsv", chrom),
        "graph_phased_sites": tsv_positions(
            args.root / "graph/phase_sites.tsv", phased_only=True),
        "private_sites": vcf_positions(args.root / "private_gap_sites.vcf", chrom),
    }
    candidate_index = tsv_label_index(
        args.root / "graph/candidates.tsv", "CATEGORY", chrom)
    filtered_index = tsv_label_index(
        args.root / "graph/filtered_sites.tsv", "REASON", chrom)
    rows = []
    for tool, directory in args.competitor:
        blocks = read_blocks(args.root / f"{tool}.blocks.tsv", chrom)
        accuracy = read_phase_accuracy(
            args.root / "eval" / f"{tool}_reads/per_phase_set.tsv", chrom)
        switches = read_switch_intervals(args.root / f"{tool}.switch_errors.bed", chrom)
        phased_positions = vcf_positions(directory / "phased.vcf.gz", chrom, True)
        for block in blocks:
            block_accuracy = accuracy.get(block["ps"])
            if block_accuracy is None \
                    or block_accuracy["accuracy"] < args.min_accuracy \
                    or block_accuracy["reads"] < args.min_reads:
                continue
            for gap in gaps:
                if block["start"] > gap["start"] or block["end"] < gap["end"]:
                    continue
                if overlaps(switches, gap["start"], gap["end"]):
                    continue
                counts = {name: count_between(values, gap["start"], gap["end"])
                          for name, values in position_sets.items()}
                counts["competitor_phased_sites"] = count_between(
                    phased_positions, gap["start"], gap["end"])
                candidate_labels = labels_between(
                    candidate_index, gap["start"], gap["end"])
                filtered_labels = labels_between(
                    filtered_index, gap["start"], gap["end"])
                counts.update({
                    "clean_het_snps": candidate_labels["CLEAN_HET_SNP"],
                    "clean_het_indels": candidate_labels["CLEAN_HET_INDEL"],
                    "repeat_het_indels": candidate_labels["REP_HET_INDEL"],
                    "filtered_ref_only": filtered_labels["ref_only"],
                    "filtered_high_af": filtered_labels["high_af"],
                    "filtered_low_af": filtered_labels["low_af"],
                    "filtered_low_depth": filtered_labels["low_depth"],
                    "filtered_no_reads": filtered_labels["no_reads_in_chunk"],
                })
                rows.append({
                    "chromosome": chrom,
                    "tool": tool,
                    "competitor_ps": block["ps"],
                    "competitor_accuracy": block_accuracy["accuracy"],
                    "competitor_reads": block_accuracy["reads"],
                    "competitor_discordant_reads": block_accuracy["discordant"],
                    "gap_start": gap["start"],
                    "gap_end": gap["end"],
                    "gap_bp": gap["end"] - gap["start"],
                    "left_graph_ps": gap["left_ps"],
                    "right_graph_ps": gap["right_ps"],
                    **counts,
                    "primary_reason": classify(counts),
                })
    rows.sort(key=lambda row: (row["chromosome"], row["tool"],
                               -row["gap_bp"], row["gap_start"]))
    fields = list(rows[0]) if rows else []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
        if fields:
            writer.writeheader()
            writer.writerows(rows)

    summary = Counter((row["tool"], row["primary_reason"]) for row in rows)
    print(f"{chrom}: {len(rows)} correct competitor bridges over {len(gaps)} graph gaps")
    for (tool, reason), count in sorted(summary.items()):
        print(f"  {tool}\t{reason}\t{count}")


if __name__ == "__main__":
    main()
