#!/usr/bin/env python3
"""Summarize the fixed chr20 gap experiment, including actual allele observations."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import pysam


GAP_BEG = 15019294
GAP_END = 15130077
PRIVATE_SNPS = (15039543, 15039634, 15048452, 15055707)


def table(path):
    with path.open() as source:
        return list(csv.DictReader(source, delimiter="\t"))


def assignments(path):
    with pysam.AlignmentFile(path) as source:
        return {r.query_name: (r.get_tag("HP"), r.get_tag("PS"))
                for r in source if r.has_tag("HP") and r.has_tag("PS")
                and r.get_tag("HP") in (1, 2)}


def matrix_links(path, arm):
    variants = {}
    observations = defaultdict(dict)
    with path.open() as source:
        for line in source:
            fields = line.rstrip().split("\t")
            if fields[0] == "VAR":
                variants[int(fields[1])] = int(fields[2])
            elif fields[0] == "OBS" and int(fields[3]) >= 0:
                observations[int(fields[2])][fields[1]] = int(fields[3])
    rows = []
    for left_pos, right_pos in ((15019256, 15039543), (15095642, 15109301)):
        for left in (i for i, pos in variants.items() if pos == left_pos):
            for right in (i for i, pos in variants.items() if pos == right_pos):
                shared = observations[left].keys() & observations[right].keys()
                pairs = Counter((observations[left][r], observations[right][r]) for r in shared)
                rows.append([arm, left_pos, right_pos, len(shared),
                             pairs[0, 0], pairs[0, 1], pairs[1, 0], pairs[1, 1]])
    return rows


def write_table(path, header, rows):
    with path.open("w") as dest:
        writer = csv.writer(dest, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    args = parser.parse_args()
    out = args.out
    frozen = args.data_root / "results/chr12-18-20-comparison/chr20"
    arms = ("graph", "hiphase", "native", "clean", "msa", "msa_margin1",
            "stitched_msa", "stitched_msa_margin1")
    summaries = []
    for arm in arms:
        data = json.loads((out / f"{arm}.eval/summary.json").read_text())
        summaries.append([arm, data["total_reads_evaluated"], data["discordant_reads"],
                          data["total_phase_sets"], data["switchflip_errors"]])
    write_table(out / "results.tsv", ["arm", "evaluated_reads", "discordant_reads",
                                      "phase_sets", "read_switchflips"], summaries)
    candidates = []
    counts = {}
    for arm in ("graph", "native", "clean", "msa", "msa_margin1"):
        path = frozen / "graph/candidates.tsv" if arm == "graph" else out / f"{arm}.tsv"
        rows = [r for r in table(path) if GAP_BEG <= int(r["POS"]) <= GAP_END]
        counts[arm] = dict(Counter(r["CATEGORY"] for r in rows))
        for r in rows:
            candidates.append([arm] + [r[k] for k in
                              ("POS", "TYPE", "REF", "ALT", "DP", "AF", "CATEGORY", "PHASE_SET")])
    write_table(out / "gap_candidates.tsv",
                ["source", "pos", "type", "ref", "alt", "depth", "af", "category", "ps"], candidates)
    catalog = args.repo / "test_data/chr20.sites.striped.vcf.gz"
    with pysam.TabixFile(str(catalog)) as source:
        contig = next(c for c in source.contigs if c.split("#")[-1] == "chr20")
        records = [r.split("\t") for r in source.fetch(contig, GAP_BEG - 1, GAP_END)]
    catalog_positions = Counter(int(r[1]) for r in records)
    graph = assignments(out / "graph.bam")
    stitched = assignments(out / "stitched_msa_margin1.bam")
    transforms = defaultdict(set)
    missing = []
    for name, (hp, ps) in graph.items():
        if name not in stitched:
            missing.append(name)
            continue
        new_hp, new_ps = stitched[name]
        transforms[ps].add((new_ps, int(hp != new_hp)))
    preserved = not missing and all(len(values) == 1 for values in transforms.values())
    audit = {"gap_1based": [GAP_BEG, GAP_END], "catalog_overlapping_records": len(records),
             "catalog_exact_positions_at_private_snps":
                 {pos: catalog_positions[pos] for pos in PRIVATE_SNPS},
             "candidate_categories": counts,
             "graph_reads": len(graph), "graph_reads_missing_after_stitch": len(missing),
             "graph_blocks_preserved_up_to_uniform_orientation": preserved,
             "graph_block_transforms": {ps: sorted(values) for ps, values in transforms.items()}}
    (out / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    links = []
    for arm in ("msa", "msa_margin1"):
        links += matrix_links(out / f"{arm}.matrix.chunk0.flags5004.tsv", arm)
    write_table(out / "boundary_observations.tsv",
                ["arm", "left_pos", "right_pos", "shared_observed_reads", "00", "01", "10", "11"], links)
    if not preserved:
        raise SystemExit("stitched output lost or internally reoriented a graph block")
    print((out / "results.tsv").read_text(), end="")


if __name__ == "__main__":
    main()
