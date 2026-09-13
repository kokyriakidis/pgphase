#!/usr/bin/env python3
"""Audit the fixed chr20:47.67-47.76 Mb second-pass experiment."""

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from summarize_chr20_gap_rephasing import assignments, table, write_table

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--data-root", type=Path, required=True)
args = parser.parse_args()
out = args.out
repo = Path(__file__).resolve().parents[2]
frozen = args.data_root / "results/chr12-18-20-comparison/chr20"
rows = []
for arm in ("graph", "hiphase", "native", "clean", "msa", "msa_margin1",
            "stitched_msa", "stitched_msa_margin1"):
    d = json.loads((out / f"{arm}.eval/summary.json").read_text())
    rows.append([arm] + [d[k] for k in ("total_reads_evaluated", "discordant_reads",
                                       "total_phase_sets", "switchflip_errors")])
write_table(out / "results.tsv", ["arm", "evaluated_reads", "discordant_reads",
                                  "phase_sets", "read_switchflips"], rows)
counts = {}
rows = []
for arm in ("graph", "native", "clean", "msa", "msa_margin1"):
    path = frozen / "graph/candidates.tsv" if arm == "graph" else out / f"{arm}.tsv"
    candidates = [r for r in table(path) if 47671540 <= int(r["POS"]) <= 47762233]
    counts[arm] = dict(Counter(r["CATEGORY"] for r in candidates))
    for r in candidates:
        rows.append([arm] + [r[k] for k in ("POS", "TYPE", "REF", "ALT", "DP", "AF",
                                             "CATEGORY", "PHASE_SET")])
write_table(out / "gap_candidates.tsv",
            ["source", "pos", "type", "ref", "alt", "depth", "af", "category", "ps"], rows)
rows = []
for arm in ("msa", "msa_margin1"):
    variants = {}
    observations = defaultdict(dict)
    with (out / f"{arm}.matrix.chunk0.flags5004.tsv").open() as source:
        for line in source:
            f = line.rstrip().split("\t")
            if f[0] == "VAR":
                variants[int(f[1])] = int(f[2])
            elif f[0] == "OBS" and int(f[3]) >= 0:
                observations[int(f[2])][f[1]] = int(f[3])
    for left, right in ((47694120, 47713869), (47751481, 47755891), (47746345, 47755891)):
        for a in (i for i, p in variants.items() if p == left):
            for b in (i for i, p in variants.items() if p == right):
                shared = observations[a].keys() & observations[b].keys()
                pairs = Counter((observations[a][r], observations[b][r]) for r in shared)
                rows.append([arm, left, right, len(shared),
                             pairs[0, 0], pairs[0, 1], pairs[1, 0], pairs[1, 1]])
write_table(out / "boundary_observations.tsv",
            ["arm", "left_pos", "right_pos", "shared_observed_reads", "00", "01", "10", "11"], rows)
with pysam.TabixFile(str(repo / "test_data/chr20.sites.striped.vcf.gz")) as source:
    chrom = next(c for c in source.contigs if c.split("#")[-1] == "chr20")
    records = [r.split("\t") for r in source.fetch(chrom, 47671539, 47762233)]
positions = Counter(int(r[1]) for r in records)
graph = assignments(out / "graph.bam")
stitched = assignments(out / "stitched_msa_margin1.bam")
transforms = defaultdict(set)
for name, (hp, ps) in graph.items():
    if name not in stitched:
        raise SystemExit(f"Missing original graph read: {name}")
    new_hp, new_ps = stitched[name]
    transforms[ps].add((new_ps, int(hp != new_hp)))
assert all(len(values) == 1 for values in transforms.values()), "inconsistent graph-block orientation"
assert next(iter(transforms[47666163]))[0] == next(iter(transforms[47762233]))[0], "target gap remains open"
audit = {"gap_1based": [47671540, 47762233], "catalog_overlapping_records": len(records),
         "clean_snp_catalog_position_counts": {p: positions[p] for p in
              (47671540, 47689418, 47713869, 47726592, 47738673, 47762233)},
         "candidate_categories": counts, "original_graph_reads_preserved": len(graph),
         "graph_block_transforms": {ps: sorted(values) for ps, values in transforms.items()},
         "target_graph_phase_sets_joined": True}
(out / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
print((out / "results.tsv").read_text(), end="")
