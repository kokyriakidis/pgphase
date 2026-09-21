#!/usr/bin/env python3
"""Check whether standalone BAM spans the graph-only gaps shorter than 10 kb."""

import argparse
import csv
import gzip
import re
import subprocess
from pathlib import Path


REFERENCE_CIGAR_OPS = frozenset("MDN=X")


def read_phase_blocks(path, contig):
    """Return PS bounds and phased heterozygous positions from a VCF."""
    blocks = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or fields[0] != contig:
                continue
            format_keys = fields[8].split(":")
            sample = fields[9].split(":")
            if "PS" not in format_keys or not sample or "|" not in sample[0]:
                continue
            alleles = sample[0].split("|")
            if len(alleles) != 2 or alleles[0] == alleles[1]:
                continue
            ps_index = format_keys.index("PS")
            if ps_index >= len(sample) or sample[ps_index] in ("", ".", "0"):
                continue
            ps = sample[ps_index]
            pos = int(fields[1])
            block = blocks.setdefault(ps, {"beg": pos, "end": pos, "sites": []})
            block["beg"] = min(block["beg"], pos)
            block["end"] = max(block["end"], pos)
            block["sites"].append(pos)
    return blocks


def read_truth_status(path):
    with gzip.open(path, "rt") as stream:
        return {
            row["read_name"]: row["status"].lower()
            for row in csv.DictReader(stream, delimiter="\t")
        }


def reference_end(pos, cigar):
    consumed = sum(
        int(length)
        for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if operation in REFERENCE_CIGAR_OPS
    )
    return pos + consumed - 1


def alignment_tag(fields, name):
    prefix = name + ":i:"
    for tag in fields[11:]:
        if tag.startswith(prefix):
            return tag[len(prefix):]
    return None


def spanning_reads(bam, contig, left, right, phase_sets):
    """Return primary HP-tagged reads physically covering both gap flanks."""
    result = subprocess.run(
        ["samtools", "view", str(bam), f"{contig}:{left}-{right}"],
        check=True,
        capture_output=True,
        text=True,
    )
    reads = set()
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        flag = int(fields[1])
        if flag & (4 | 256 | 2048):
            continue
        pos = int(fields[3])
        if pos > left or reference_end(pos, fields[5]) < right:
            continue
        hp = alignment_tag(fields, "HP")
        ps = alignment_tag(fields, "PS")
        if hp in ("1", "2") and ps in phase_sets:
            reads.add(fields[0])
    return reads


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gaps", type=Path, required=True)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--per-read", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--vcf-contig", default="chr20")
    parser.add_argument("--bam-contig", default="chr20")
    args = parser.parse_args()

    blocks = read_phase_blocks(args.vcf, args.vcf_contig)
    truth_status = read_truth_status(args.per_read)
    with args.gaps.open() as stream:
        gaps = list(csv.DictReader(stream, delimiter="\t"))

    rows = []
    for gap in gaps:
        left = int(gap["graph_left"])
        right = int(gap["graph_right"])
        covering = {
            ps: block
            for ps, block in blocks.items()
            if block["beg"] <= left and block["end"] >= right
        }
        reads = spanning_reads(args.bam, args.bam_contig, left, right,
                               frozenset(covering)) if covering else set()
        scored = [
            truth_status[name]
            for name in reads
            if truth_status.get(name) in ("concordant", "discordant")
        ]
        concordant = sum(status == "concordant" for status in scored)
        discordant = sum(status == "discordant" for status in scored)
        local_consistent = max(concordant, discordant)
        chosen_ps = min(covering, key=lambda ps: covering[ps]["end"] -
                        covering[ps]["beg"]) if covering else "."
        chosen = covering.get(chosen_ps)
        rows.append({
            "gap_beg": gap["gap_beg"],
            "gap_end": gap["gap_end"],
            "gap_bp": gap["gap_bp"],
            "graph_left": left,
            "graph_right": right,
            "bam_spans": int(bool(covering)),
            "bam_ps": chosen_ps,
            "bam_ps_beg": chosen["beg"] if chosen else 0,
            "bam_ps_end": chosen["end"] if chosen else 0,
            "bam_het_sites": len(chosen["sites"]) if chosen else 0,
            "spanning_reads": len(reads),
            "scored_reads": len(scored),
            "concordant_reads": concordant,
            "discordant_reads": discordant,
            # Haplotype numbers are arbitrary within a local block. The first
            # value retains the frozen whole-block orientation; local purity
            # flips it when needed and measures whether the gap itself is
            # crossed consistently.
            "truth_concordance": (f"{concordant / len(scored):.4f}"
                                  if scored else "."),
            "local_phase_purity": (f"{local_consistent / len(scored):.4f}"
                                   if scored else "."),
            "recovery_spans": gap["recovery_spans"],
        })

    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    closed = [row for row in rows if row["bam_spans"]]
    scored = [row for row in closed if row["scored_reads"]]
    print(f"BAM spans {len(closed)}/{len(rows)} short graph gaps")
    print(f"spanned bases: {sum(int(row['gap_bp']) for row in closed)}/"
          f"{sum(int(row['gap_bp']) for row in rows)}")
    print(f"spans with scored physical crossing reads: {len(scored)}/{len(closed)}")
    concordant = sum(row["concordant_reads"] for row in scored)
    consistent = sum(max(row["concordant_reads"], row["discordant_reads"])
                     for row in scored)
    total = sum(row["scored_reads"] for row in scored)
    print(f"whole-block-oriented crossing-read concordance: {concordant}/{total} "
          f"({100 * concordant / total:.2f}%)" if total else
          "crossing-read concordance: no scored reads")
    if total:
        print(f"orientation-independent local phase purity: {consistent}/{total} "
              f"({100 * consistent / total:.2f}%)")


if __name__ == "__main__":
    main()
