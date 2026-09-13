#!/usr/bin/env python3
"""Aggregate chromosome comparison outputs into stable TSV files."""

import csv
import json
import os
import re
import sys
from statistics import median
from pathlib import Path


TOOLS = ("graph", "bam", "hybrid", "graph_lock", "whatshap", "whatshap_opt", "hiphase", "longphase")


def parse_elapsed(value):
    parts = value.split(":")
    seconds = float(parts[-1])
    if len(parts) >= 2:
        seconds += int(parts[-2]) * 60
    if len(parts) == 3:
        seconds += int(parts[0]) * 3600
    return round(seconds, 2)


def parse_resources(path):
    values = {"phasing_elapsed_seconds": "NA", "max_rss_kb": "NA"}
    if not path.exists():
        return values
    text = path.read_text()
    elapsed = re.search(r"Elapsed \(wall clock\) time.*\):\s*([0-9:.]+)$",
                        text, re.MULTILINE)
    rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    if elapsed:
        values["phasing_elapsed_seconds"] = parse_elapsed(elapsed.group(1))
    if rss:
        values["max_rss_kb"] = int(rss.group(1))
    return values


def parse_compare(path):
    text = path.read_text()
    values = {"assessed_pairs": "NA", "switches": "NA", "flips": "NA",
              "switchflips": "NA", "hamming": "NA"}
    patterns = {
        "assessed_pairs": r"ALL INTERSECTION BLOCKS:.*?phased pairs of variants assessed:\s*(\d+)",
        "switches": r"ALL INTERSECTION BLOCKS:.*?switch errors:\s*(\d+)",
        "hamming": r"ALL INTERSECTION BLOCKS:.*?Block-wise Hamming distance:\s*(\d+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
        if match:
            values[key] = match.group(1)
    decomposition = re.search(
        r"ALL INTERSECTION BLOCKS:.*?switch/flip decomposition:\s*(\d+)/(\d+)",
        text, re.DOTALL | re.IGNORECASE)
    if decomposition:
        switches, flips = map(int, decomposition.groups())
        values["flips"] = str(flips)
        values["switchflips"] = str(switches + flips)
    return values


def parse_stats(path):
    values = {"phased_variants": "NA", "blocks": "NA", "n50_bp": "NA"}
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    row = next((r for r in rows if r.get("chromosome") == "ALL"), rows[-1] if rows else {})
    aliases = {
        "phased_variants": ("phased", "variants_phased"),
        "blocks": ("blocks", "phase_blocks"),
        "n50_bp": ("block_n50", "N50", "n50"),
    }
    for output_key, candidates in aliases.items():
        for candidate in candidates:
            if candidate in row:
                values[output_key] = row[candidate]
                break
    return values


def parse_gap_summary(path):
    scalar = {}
    sections = {}
    section = None
    for line in path.read_text().splitlines():
        if not line:
            section = None
            continue
        fields = line.split("\t")
        if len(fields) == 1:
            section = fields[0]
        elif len(fields) == 2 and fields[1].isdigit():
            if section:
                sections.setdefault(section, {})[fields[0]] = int(fields[1])
            else:
                scalar[fields[0]] = int(fields[1])
    return scalar, sections


def main():
    if len(sys.argv) < 3:
        raise SystemExit("usage: summarize.py OUT_ROOT CHROM [CHROM ...]")
    root = Path(sys.argv[1])
    output_dir = Path(os.environ.get("SUMMARY_OUT_DIR", root))
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for chrom in sys.argv[2:]:
        for tool in TOOLS:
            prefix = root / chrom / tool
            compare = parse_compare(root / chrom / f"{tool}.compare.txt")
            stats = parse_stats(root / chrom / f"{tool}.stats.tsv")
            ngc = json.loads((root / chrom / f"{tool}.ngc50.json").read_text())
            read = json.loads((root / chrom / "eval" / f"{tool}_reads" / "summary.json").read_text())
            rows.append({
                "chromosome": chrom,
                "tool": tool,
                **compare,
                **stats,
                "chromosome_ngc50_bp": ngc.get("ngc50_bp", "NA"),
                "phased_reads": read.get("total_phased_reads", "NA"),
                "evaluated_reads": read.get("total_reads_evaluated", "NA"),
                "discordant_reads": read.get("discordant_reads", "NA"),
                "read_hamming_rate": read.get("hamming_error_rate", "NA"),
                "read_switchflips": read.get("switchflip_errors", "NA"),
                "read_block_n50_bp": read.get("phase_block_n50_bp", "NA"),
                **parse_resources(prefix / "resources.txt"),
                "vcf": str(prefix / ("shared.vcf.gz" if tool in {"graph", "bam", "hybrid", "graph_lock"} else "phased.vcf.gz")),
            })
    fields = list(rows[0]) if rows else []
    with (output_dir / "results.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    pooled = []
    summed = ("assessed_pairs", "switches", "flips", "switchflips", "hamming",
              "phased_variants", "blocks", "phased_reads", "evaluated_reads",
              "discordant_reads", "read_switchflips")
    for tool in TOOLS:
        tool_rows = [row for row in rows if row["tool"] == tool]
        aggregate = {"tool": tool, "chromosomes": len(tool_rows)}
        for field in summed:
            aggregate[field] = sum(int(row[field]) for row in tool_rows)
        aggregate["variant_hamming_rate"] = (
            aggregate["hamming"] / aggregate["assessed_pairs"]
            if aggregate["assessed_pairs"] else "NA")
        aggregate["read_hamming_rate"] = (
            aggregate["discordant_reads"] / aggregate["evaluated_reads"]
            if aggregate["evaluated_reads"] else "NA")
        aggregate["median_chromosome_ngc50_bp"] = int(median(
            int(row["chromosome_ngc50_bp"]) for row in tool_rows))
        aggregate["median_read_block_n50_bp"] = int(median(
            int(row["read_block_n50_bp"]) for row in tool_rows))
        numeric_times = [float(row["phasing_elapsed_seconds"]) for row in tool_rows
                         if row["phasing_elapsed_seconds"] != "NA"]
        numeric_rss = [int(row["max_rss_kb"]) for row in tool_rows
                       if row["max_rss_kb"] != "NA"]
        aggregate["total_phasing_elapsed_seconds"] = (
            round(sum(numeric_times), 2) if numeric_times else "NA")
        aggregate["max_rss_kb"] = max(numeric_rss) if numeric_rss else "NA"
        pooled.append(aggregate)
    pooled_fields = list(pooled[0]) if pooled else []
    with (output_dir / "pooled.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, delimiter="\t", fieldnames=pooled_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(pooled)

    gap_rows = []
    for chrom in sys.argv[2:]:
        scalar, sections = parse_gap_summary(root / chrom / "gap_site_loss.summary.txt")
        reasons = sections.get("primary_reason", {})
        recovered_sites = sections.get("recovered_by_phased_site", {})
        recovered_blocks = sections.get("recovered_inside_phase_block", {})
        gap_rows.append({
            "chromosome": chrom,
            "truth_hets": scalar.get("truth_hets", "NA"),
            "truth_hets_outside_graph_blocks": scalar.get("truth_hets_outside_blocks", "NA"),
            "competitor_covered_gap_hets": scalar.get(
                "truth_hets_outside_pgphase_inside_competitor_blocks", "NA"),
            "no_catalog_record": reasons.get("no_catalog_record", 0),
            "ref_only": reasons.get("ref_only", 0),
            "no_reads_in_chunk": reasons.get("no_reads_in_chunk", 0),
            "catalog_retained_unphased": reasons.get("catalog_site_retained_unphased", 0),
            "hybrid_recovered_phased_sites": recovered_sites.get("hybrid", 0),
            "graph_lock_recovered_phased_sites": recovered_sites.get("graph_lock", 0),
            "hybrid_recovered_inside_blocks": recovered_blocks.get("hybrid", 0),
            "graph_lock_recovered_inside_blocks": recovered_blocks.get("graph_lock", 0),
        })
    gap_fields = list(gap_rows[0]) if gap_rows else []
    with (output_dir / "gap_summary.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, delimiter="\t", fieldnames=gap_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(gap_rows)


if __name__ == "__main__":
    main()
