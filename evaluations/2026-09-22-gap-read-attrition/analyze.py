#!/usr/bin/env python3
"""Measure which read evidence reaches pgphase and HiPhase at missed chr20 gaps.

Coordinates in the target table are one based.  A qname bridges a gap when one
accepted alignment reaches its left edge and one accepted alignment reaches its
right edge; those can be the same alignment or separate primary/supplementary
records.  This matches HiPhase's decision to collapse records by qname while
also exposing the evidence pgphase loses by accepting primary records only.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pysam


def load_targets(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        rows = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(rows, delimiter="\t"))


def load_truth(path: Path) -> dict[str, int]:
    truth: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            qname, label = line.rstrip().split("\t")[:2]
            truth[qname] = 1 if label == "MATERNAL" else 2
    return truth


def load_tags(path: Path) -> dict[str, tuple[int, int]]:
    tags: dict[str, tuple[int, int]] = {}
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        for record in bam.fetch(until_eof=True):
            if not record.has_tag("HP") or not record.has_tag("PS"):
                continue
            hp = int(record.get_tag("HP"))
            ps = int(record.get_tag("PS"))
            if hp in (1, 2) and ps > 0:
                tags.setdefault(record.query_name, (hp, ps))
    return tags


def truth_purity(qnames: set[str], tags: dict[str, tuple[int, int]],
                 truth: dict[str, int]) -> tuple[int, float]:
    # HP labels are arbitrary independently in every phase set. Orient each PS
    # to truth before summing, matching the read evaluator used by the tests.
    by_phase_set: dict[int, list[int]] = defaultdict(lambda: [0, 0])
    for qname in qnames:
        if qname not in tags or qname not in truth:
            continue
        hp, phase_set = tags[qname]
        by_phase_set[phase_set][hp != truth[qname]] += 1
    scored = sum(same + swapped for same, swapped in by_phase_set.values())
    correct = sum(max(same, swapped) for same, swapped in by_phase_set.values())
    return scored, (correct / scored if scored else 0.0)


def accepted(record: pysam.AlignedSegment, min_mapq: int,
             allow_supplementary: bool) -> bool:
    return not (
        record.is_unmapped
        or record.is_secondary
        or record.is_qcfail
        or record.is_duplicate
        or (record.is_supplementary and not allow_supplementary)
        or record.mapping_quality < min_mapq
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--raw-bam", type=Path, required=True)
    parser.add_argument("--graph-bam", type=Path, required=True)
    parser.add_argument("--bam-bam", type=Path, required=True)
    parser.add_argument("--hiphase-bam", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    targets = load_targets(args.targets)
    truth = load_truth(args.truth)
    graph_tags = load_tags(args.graph_bam)
    bam_tags = load_tags(args.bam_bam)
    hiphase_tags = load_tags(args.hiphase_bam)

    # Per target and qname, retain whether accepted records reach each flank.
    pg_reach = [defaultdict(lambda: [False, False]) for _ in targets]
    hi_reach = [defaultdict(lambda: [False, False]) for _ in targets]
    primary_span = [set() for _ in targets]
    mapq_bins = [defaultdict(int) for _ in targets]
    supplementary_records = [0 for _ in targets]

    with pysam.AlignmentFile(args.raw_bam, "rb") as bam:
        for record in bam.fetch(until_eof=True):
            if record.reference_name is None or "chr20" not in record.reference_name:
                continue
            start = record.reference_start + 1
            end = record.reference_end or start
            for index, row in enumerate(targets):
                left = int(row["gap_left"])
                right = int(row["gap_right"])
                if end < left or start > right:
                    continue
                if not record.is_secondary and not record.is_supplementary and not record.is_unmapped:
                    if start <= left and end >= right:
                        primary_span[index].add(record.query_name)
                        mq = record.mapping_quality
                        bucket = "0" if mq == 0 else "1_4" if mq < 5 else "5_29" if mq < 30 else "ge30"
                        mapq_bins[index][bucket] += 1
                if record.is_supplementary:
                    supplementary_records[index] += 1
                if accepted(record, 1, False):
                    reach = pg_reach[index][record.query_name]
                    reach[0] |= start <= left <= end
                    reach[1] |= start <= right <= end
                if accepted(record, 5, True):
                    reach = hi_reach[index][record.query_name]
                    reach[0] |= start <= left <= end
                    reach[1] |= start <= right <= end

    columns = [
        "gap_left", "gap_right", "gap_bp", "support_class",
        "primary_span", "mapq0", "mapq1_4", "mapq5_29", "mapq_ge30",
        "pg_eligible_bridge", "hi_eligible_bridge", "hi_supplementary_gain",
        "graph_tagged", "bam_tagged", "graph_and_bam_tagged", "hiphase_tagged",
        "graph_scored", "graph_purity", "bam_scored", "bam_purity",
        "hiphase_scored", "hiphase_purity", "supplementary_records",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for index, row in enumerate(targets):
            pg_bridge = {q for q, reach in pg_reach[index].items() if all(reach)}
            hi_bridge = {q for q, reach in hi_reach[index].items() if all(reach)}
            graph_q = pg_bridge & graph_tags.keys()
            bam_q = pg_bridge & bam_tags.keys()
            both_q = graph_q & bam_q
            hi_q = hi_bridge & hiphase_tags.keys()
            graph_scored, graph_purity = truth_purity(graph_q, graph_tags, truth)
            bam_scored, bam_purity = truth_purity(bam_q, bam_tags, truth)
            hi_scored, hi_purity = truth_purity(hi_q, hiphase_tags, truth)
            writer.writerow({
                "gap_left": row["gap_left"],
                "gap_right": row["gap_right"],
                "gap_bp": row["gap_bp"],
                "support_class": row["support_class"],
                "primary_span": len(primary_span[index]),
                "mapq0": mapq_bins[index]["0"],
                "mapq1_4": mapq_bins[index]["1_4"],
                "mapq5_29": mapq_bins[index]["5_29"],
                "mapq_ge30": mapq_bins[index]["ge30"],
                "pg_eligible_bridge": len(pg_bridge),
                "hi_eligible_bridge": len(hi_bridge),
                "hi_supplementary_gain": len(hi_bridge - pg_bridge),
                "graph_tagged": len(graph_q),
                "bam_tagged": len(bam_q),
                "graph_and_bam_tagged": len(both_q),
                "hiphase_tagged": len(hi_q),
                "graph_scored": graph_scored,
                "graph_purity": f"{graph_purity:.4f}",
                "bam_scored": bam_scored,
                "bam_purity": f"{bam_purity:.4f}",
                "hiphase_scored": hi_scored,
                "hiphase_purity": f"{hi_purity:.4f}",
                "supplementary_records": supplementary_records[index],
            })


if __name__ == "__main__":
    main()
