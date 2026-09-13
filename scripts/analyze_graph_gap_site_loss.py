#!/usr/bin/env python3
"""Explain truth hets in pgphase graph phase-block gaps.

This is a lightweight, dependency-free counterpart to the larger benchmark
notebooks: it asks, for truth heterozygous sites that fall outside pgphase phase
blocks, whether the graph catalog has a nearby snarl record and whether the
current graph run emitted or filtered a nearby candidate.  When competitor VCFs
are supplied, the target is narrowed either to truth hets inside competitor
merged phase blocks or to truth hets near a phased competitor variant. Recovery
VCFs can then measure how many of those target sites another pgphase mode calls,
phases, or covers with a phase block.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
from collections import Counter
from pathlib import Path


def parent_site_id(site_id: str) -> str:
    if ":" not in site_id:
        return site_id
    parent, suffix = site_id.rsplit(":", 1)
    return parent if suffix.isdigit() else site_id


def open_maybe_gzip(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def norm_chrom(chrom: str) -> str:
    return chrom.split("#")[-1]


def is_het_gt(gt: str) -> bool:
    gt = gt.split(":", 1)[0].replace("|", "/")
    alleles = [a for a in gt.split("/") if a != "."]
    return len(alleles) == 2 and alleles[0] != alleles[1]


def gt_alt_indices(gt: str) -> set[int]:
    gt = gt.split(":", 1)[0].replace("|", "/")
    out: set[int] = set()
    for allele in gt.split("/"):
        if allele.isdigit() and int(allele) > 0:
            out.add(int(allele))
    return out


def parse_truth_hets(path: Path, contig: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if norm_chrom(fields[0]) != contig:
                continue
            if len(fields) < 10 or not is_het_gt(fields[9]):
                continue
            alts = fields[4].split(",")
            for ai in sorted(gt_alt_indices(fields[9])):
                if ai - 1 >= len(alts):
                    continue
                rows.append({
                    "chrom": contig,
                    "pos": int(fields[1]),
                    "ref": fields[3],
                    "alt": alts[ai - 1],
                })
    rows.sort(key=lambda r: int(r["pos"]))
    return rows


def parse_phase_blocks(path: Path, contig: str) -> list[tuple[int, int]]:
    by_ps: dict[int, list[int]] = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if norm_chrom(fields[0]) != contig or len(fields) < 10:
                continue
            fmt = fields[8].split(":")
            sample = fields[9].split(":")
            data = dict(zip(fmt, sample))
            gt = data.get("GT", "")
            if "|" not in gt or not is_het_gt(gt):
                continue
            ps_s = data.get("PS")
            if not ps_s or ps_s == ".":
                continue
            try:
                ps = int(ps_s)
                pos = int(fields[1])
            except ValueError:
                continue
            by_ps.setdefault(ps, []).append(pos)

    blocks = [(min(v), max(v) + 1) for v in by_ps.values() if v]
    blocks.sort()
    merged: list[tuple[int, int]] = []
    for beg, end in blocks:
        if not merged or beg > merged[-1][1]:
            merged.append((beg, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def parse_phased_positions(path: Path, contig: str) -> list[int]:
    positions: set[int] = set()
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if norm_chrom(fields[0]) != contig or len(fields) < 10:
                continue
            fmt = fields[8].split(":")
            sample = fields[9].split(":")
            data = dict(zip(fmt, sample))
            gt = data.get("GT", "")
            if "|" not in gt or not is_het_gt(gt):
                continue
            ps_s = data.get("PS")
            if not ps_s or ps_s == ".":
                continue
            try:
                positions.add(int(fields[1]))
            except ValueError:
                continue
    return sorted(positions)


def parse_het_positions(path: Path, contig: str) -> list[int]:
    positions: set[int] = set()
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if norm_chrom(fields[0]) != contig or len(fields) < 10:
                continue
            fmt = fields[8].split(":")
            sample = fields[9].split(":")
            data = dict(zip(fmt, sample))
            if is_het_gt(data.get("GT", "")):
                positions.add(int(fields[1]))
    return sorted(positions)


def in_any_block(pos: int, blocks: list[tuple[int, int]]) -> bool:
    starts = [b for b, _ in blocks]
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and blocks[i][0] <= pos < blocks[i][1]


def near_any_position(pos: int, positions: list[int], window: int) -> bool:
    i = bisect.bisect_left(positions, pos - window)
    return i < len(positions) and positions[i] <= pos + window


def positions_in_complement(
    rows: list[dict[str, object]], blocks: list[tuple[int, int]]
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        pos = int(row["pos"])
        if not in_any_block(pos, blocks):
            out.append(row)
    return out


def positions_in_any_competitor_block(
    rows: list[dict[str, object]], competitor_blocks: dict[str, list[tuple[int, int]]]
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        pos = int(row["pos"])
        labels = [
            label for label, blocks in competitor_blocks.items()
            if in_any_block(pos, blocks)
        ]
        if labels:
            row = dict(row)
            row["competitors_covering"] = ",".join(labels)
            out.append(row)
    return out


def positions_near_any_competitor_site(
    rows: list[dict[str, object]],
    competitor_positions: dict[str, list[int]],
    window: int,
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        pos = int(row["pos"])
        labels = [
            label for label, positions in competitor_positions.items()
            if near_any_position(pos, positions, window)
        ]
        if labels:
            row = dict(row)
            row["competitors_covering"] = ",".join(labels)
            out.append(row)
    return out


def parse_catalog(path: Path, contig: str) -> tuple[list[int], dict[int, list[tuple[str, int]]]]:
    by_pos: dict[int, list[tuple[str, int]]] = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if norm_chrom(fields[0]) != contig:
                continue
            pos = int(fields[1])
            n_alt = 0 if fields[4] == "." else len(fields[4].split(","))
            by_pos.setdefault(pos, []).append((fields[2], n_alt))
    return sorted(by_pos), by_pos


def parse_filtered(path: Path) -> tuple[list[int], dict[int, Counter[str]]]:
    by_pos: dict[int, Counter[str]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pos = int(row["POS"])
            by_pos.setdefault(pos, Counter())[row["REASON"]] += 1
    return sorted(by_pos), by_pos


def parse_filtered_by_site(path: Path) -> dict[str, Counter[str]]:
    by_site: dict[str, Counter[str]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid = parent_site_id(row["SITE_ID"])
            by_site.setdefault(sid, Counter())[row["REASON"]] += 1
    return by_site


def parse_candidates(path: Path) -> tuple[list[int], dict[int, Counter[str]]]:
    by_pos: dict[int, Counter[str]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pos = int(row["POS"])
            category = row["CATEGORY"]
            phased = "phased" if row["PHASE_SET"] not in ("", "-1", "0") else "unphased"
            by_pos.setdefault(pos, Counter())[f"{category}:{phased}"] += 1
    return sorted(by_pos), by_pos


def parse_phase_sites(path: Path) -> dict[str, Counter[str]]:
    by_site: dict[str, Counter[str]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid = parent_site_id(row["SITE_ID"])
            phased = "phased" if row["PHASE_SET"] not in ("", "-1", "0") else "unphased"
            by_site.setdefault(sid, Counter())[phased] += 1
    return by_site


def nearby_values(positions: list[int], by_pos: dict[int, object], pos: int, window: int):
    left = bisect.bisect_left(positions, pos - window)
    right = bisect.bisect_right(positions, pos + window)
    return [(p, by_pos[p]) for p in positions[left:right]]


def best_catalog_class(records: list[tuple[int, list[tuple[str, int]]]]) -> str:
    if not records:
        return "no_catalog_record"
    max_alt = max(n_alt for _, site_records in records for _, n_alt in site_records)
    if max_alt >= 8:
        return "catalog_8plus_alt"
    if max_alt >= 4:
        return "catalog_4to7_alt"
    if max_alt >= 2:
        return "catalog_2to3_alt"
    return "catalog_biallelic"


def top_counter(counter: Counter[str]) -> str:
    if not counter:
        return ""
    return counter.most_common(1)[0][0]


def top_informative_reason(counter: Counter[str]) -> str:
    for reason in (
        "high_af",
        "low_af",
        "low_depth",
        "ref_only",
        "multiallelic_unsupported",
        "no_reads_in_chunk",
    ):
        if counter[reason] > 0:
            return reason
    return top_counter(counter)


def classify_site(
    catalog_hits: list[tuple[int, list[tuple[str, int]]]],
    filtered_hits: list[tuple[int, Counter[str]]],
    candidate_hits: list[tuple[int, Counter[str]]],
) -> str:
    if not catalog_hits:
        return "no_catalog_record"
    if candidate_hits:
        merged: Counter[str] = Counter()
        for _, counts in candidate_hits:
            merged.update(counts)
        phased = sum(v for k, v in merged.items() if k.endswith(":phased"))
        return "candidate_nearby_phased" if phased else "candidate_nearby_unphased"
    if filtered_hits:
        merged: Counter[str] = Counter()
        for _, counts in filtered_hits:
            merged.update(counts)
        return top_counter(merged)
    return "catalog_nearby_no_candidate_or_filter_row"


def classify_site_by_id(
    catalog_hits: list[tuple[int, list[tuple[str, int]]]],
    filtered_by_site: dict[str, Counter[str]],
    phase_by_site: dict[str, Counter[str]],
) -> str:
    if not catalog_hits:
        return "no_catalog_record"
    site_ids = {
        parent_site_id(site_id)
        for _, site_records in catalog_hits
        for site_id, _ in site_records
    }
    retained: Counter[str] = Counter()
    for site_id in site_ids:
        retained.update(phase_by_site.get(site_id, Counter()))
    if retained:
        return "catalog_site_retained_phased" if retained["phased"] else "catalog_site_retained_unphased"

    filtered: Counter[str] = Counter()
    for site_id in site_ids:
        filtered.update(filtered_by_site.get(site_id, Counter()))
    if filtered:
        return top_informative_reason(filtered)
    return "catalog_nearby_no_site_id_accounting_row"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth-vcf", type=Path, required=True)
    ap.add_argument("--phased-vcf", type=Path, required=True)
    ap.add_argument("--catalog-vcf", type=Path, required=True)
    ap.add_argument("--candidates-tsv", type=Path, required=True)
    ap.add_argument("--filtered-tsv", type=Path, required=True)
    ap.add_argument(
        "--phase-sites-tsv",
        type=Path,
        help="Optional retained graph-site TSV from collect-graph-variation --phase-sites-out.",
    )
    ap.add_argument("--contig", default="chr20")
    ap.add_argument("--window", type=int, default=25)
    ap.add_argument(
        "--competitor-vcf",
        action="append",
        default=[],
        metavar="LABEL=PATH",
        help="Optional phased competitor VCF. May be repeated.",
    )
    ap.add_argument(
        "--competitor-target",
        choices=("block", "site"),
        default="block",
        help="With competitor VCFs, use phase blocks or nearby phased sites as the target.",
    )
    ap.add_argument(
        "--competitor-site-window",
        type=int,
        default=25,
        help="Window for --competitor-target site.",
    )
    ap.add_argument(
        "--recovery-vcf",
        action="append",
        default=[],
        metavar="LABEL=PATH",
        help=("Optional phased VCF to test for recovery of target sites. "
              "May be repeated."),
    )
    ap.add_argument(
        "--recovery-site-window",
        type=int,
        default=25,
        help="Window for matching phased sites in --recovery-vcf outputs.",
    )
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    truth = parse_truth_hets(args.truth_vcf, args.contig)
    blocks = parse_phase_blocks(args.phased_vcf, args.contig)
    gap_truth = positions_in_complement(truth, blocks)
    competitor_blocks: dict[str, list[tuple[int, int]]] = {}
    competitor_positions: dict[str, list[int]] = {}
    for spec in args.competitor_vcf:
        if "=" not in spec:
            raise SystemExit(f"--competitor-vcf must be LABEL=PATH, got {spec!r}")
        label, path_s = spec.split("=", 1)
        path = Path(path_s)
        competitor_blocks[label] = parse_phase_blocks(path, args.contig)
        competitor_positions[label] = parse_phased_positions(path, args.contig)
    recovery_blocks: dict[str, list[tuple[int, int]]] = {}
    recovery_het_positions: dict[str, list[int]] = {}
    recovery_positions: dict[str, list[int]] = {}
    for spec in args.recovery_vcf:
        if "=" not in spec:
            raise SystemExit(f"--recovery-vcf must be LABEL=PATH, got {spec!r}")
        label, path_s = spec.split("=", 1)
        path = Path(path_s)
        recovery_blocks[label] = parse_phase_blocks(path, args.contig)
        recovery_het_positions[label] = parse_het_positions(path, args.contig)
        recovery_positions[label] = parse_phased_positions(path, args.contig)
    if competitor_blocks:
        if args.competitor_target == "block":
            target_truth = positions_in_any_competitor_block(gap_truth, competitor_blocks)
        else:
            target_truth = positions_near_any_competitor_site(
                gap_truth, competitor_positions, args.competitor_site_window
            )
    else:
        target_truth = gap_truth
    cat_pos, cat_by_pos = parse_catalog(args.catalog_vcf, args.contig)
    filt_pos, filt_by_pos = parse_filtered(args.filtered_tsv)
    filt_by_site = parse_filtered_by_site(args.filtered_tsv) if args.phase_sites_tsv else {}
    cand_pos, cand_by_pos = parse_candidates(args.candidates_tsv)
    phase_by_site = parse_phase_sites(args.phase_sites_tsv) if args.phase_sites_tsv else {}

    primary = Counter()
    catalog_class = Counter()
    rows: list[dict[str, object]] = []
    competitor_cover = Counter()
    recovery_het_cover = Counter()
    recovery_site_cover = Counter()
    recovery_block_cover = Counter()
    for row in target_truth:
        pos = int(row["pos"])
        cat_hits = nearby_values(cat_pos, cat_by_pos, pos, args.window)
        filt_hits = nearby_values(filt_pos, filt_by_pos, pos, args.window)
        cand_hits = nearby_values(cand_pos, cand_by_pos, pos, args.window)
        why = classify_site_by_id(cat_hits, filt_by_site, phase_by_site) \
            if args.phase_sites_tsv else classify_site(cat_hits, filt_hits, cand_hits)
        cclass = best_catalog_class(cat_hits)
        primary[why] += 1
        catalog_class[cclass] += 1
        if "competitors_covering" in row:
            for label in str(row["competitors_covering"]).split(","):
                competitor_cover[label] += 1
        recovered_as_het = [
            label for label, positions in recovery_het_positions.items()
            if near_any_position(pos, positions, args.recovery_site_window)
        ]
        recovered_by_site = [
            label for label, positions in recovery_positions.items()
            if near_any_position(pos, positions, args.recovery_site_window)
        ]
        recovered_by_block = [
            label for label, recovery_blocks_for_label in recovery_blocks.items()
            if in_any_block(pos, recovery_blocks_for_label)
        ]
        for label in recovered_as_het:
            recovery_het_cover[label] += 1
        for label in recovered_by_site:
            recovery_site_cover[label] += 1
        for label in recovered_by_block:
            recovery_block_cover[label] += 1
        if args.out:
            rows.append({
                **row,
                "primary_reason": why,
                "catalog_class": cclass,
                "near_catalog_records": len(cat_hits),
                "near_filtered_positions": len(filt_hits),
                "near_candidate_positions": len(cand_hits),
                "recovered_as_het": ",".join(recovered_as_het),
                "recovered_by_phased_site": ",".join(recovered_by_site),
                "recovered_inside_phase_block": ",".join(recovered_by_block),
            })

    print(f"truth_hets\t{len(truth)}")
    print(f"phase_blocks_merged\t{len(blocks)}")
    print(f"truth_hets_outside_blocks\t{len(gap_truth)}")
    for label, comp_blocks in competitor_blocks.items():
        print(f"competitor_blocks_merged:{label}\t{len(comp_blocks)}")
    for label, comp_positions in competitor_positions.items():
        print(f"competitor_phased_positions:{label}\t{len(comp_positions)}")
    for label, positions in recovery_positions.items():
        print(f"recovery_phased_positions:{label}\t{len(positions)}")
    if competitor_blocks:
        if args.competitor_target == "block":
            print(f"truth_hets_outside_pgphase_inside_competitor_blocks\t{len(target_truth)}")
        else:
            print(f"truth_hets_outside_pgphase_near_competitor_sites\t{len(target_truth)}")
        print("\ncompetitor_cover")
        for key, val in competitor_cover.most_common():
            print(f"{key}\t{val}")
    else:
        print(f"target_truth_hets\t{len(target_truth)}")
    if recovery_positions:
        print("\nrecovered_as_het")
        for label in recovery_het_positions:
            print(f"{label}\t{recovery_het_cover[label]}")
        print("\nrecovered_by_phased_site")
        for label in recovery_positions:
            print(f"{label}\t{recovery_site_cover[label]}")
        print("\nrecovered_inside_phase_block")
        for label in recovery_blocks:
            print(f"{label}\t{recovery_block_cover[label]}")
    print("\nprimary_reason")
    for key, val in primary.most_common():
        print(f"{key}\t{val}")
    print("\ncatalog_class")
    for key, val in catalog_class.most_common():
        print(f"{key}\t{val}")

    if args.out:
        with args.out.open("w", newline="") as fh:
            fieldnames = [
                "chrom", "pos", "ref", "alt", "primary_reason", "catalog_class",
                "competitors_covering",
                "near_catalog_records", "near_filtered_positions",
                "near_candidate_positions", "recovered_as_het",
                "recovered_by_phased_site",
                "recovered_inside_phase_block",
            ]
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote\t{args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
