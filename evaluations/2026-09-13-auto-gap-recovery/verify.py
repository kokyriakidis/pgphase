#!/usr/bin/env python3
"""Verify cumulative tiers and preservation of the initial block orientations."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import pysam


def read_phase(path):
    with pysam.AlignmentFile(path) as source:
        return {r.query_name: (r.get_tag("HP"), r.get_tag("PS")) for r in source
                if r.has_tag("HP") and r.has_tag("PS") and r.get_tag("HP") in (1, 2)}


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
results = []
invariants = {}
for name in ("15m", "47m", "15m_split"):
    case = args.out / name
    initial = read_phase(case / "clean.bam")
    for arm in ("clean", "auto24", "auto1"):
        summary = json.loads((case / f"{arm}.eval/summary.json").read_text())
        results.append([name, arm] + [summary[key] for key in
                       ("total_reads_evaluated", "discordant_reads", "total_phase_sets", "switchflip_errors")])
        if arm == "clean":
            continue
        assert summary["discordant_reads"] == 0 and summary["switchflip_errors"] == 0, (name, arm, "truth regression")
        assert summary["total_phase_sets"] == 1, (name, arm, "regional gap not closed")
        after = read_phase(case / f"{arm}.bam")
        transforms = defaultdict(set)
        for read, (hp, ps) in initial.items():
            assert read in after, (name, arm, "lost initial read", read)
            new_hp, new_ps = after[read]
            transforms[ps].add((new_ps, int(hp != new_hp)))
        assert all(len(values) == 1 for values in transforms.values()), (name, arm, "mixed original block")
        with (case / f"{arm}.tiers.tsv").open() as source:
            rows = list(csv.DictReader(source, delimiter="\t"))
        assert rows, (name, arm, "no automatically detected gaps")
        by_gap = defaultdict(list)
        for row in rows:
            by_gap[row["CHROM"], row["GAP_LEFT"], row["GAP_RIGHT"]].append(row)
        for gap, tiers in by_gap.items():
            assert [int(row["TIER"]) for row in tiers] == list(range(1, len(tiers) + 1)), gap
            assert all(row["STATUS"] != "joined" for row in tiers[:-1]), gap
            snps = 0
            for row in tiers:
                tier = int(row["TIER"])
                assert tier == 3 or int(row["MSA_HET_INDELS"]) == 0, (gap, "early indel")
                assert int(row["MSA_HET_SNPS"]) >= snps, (gap, "lost SNP tier")
                snps = int(row["MSA_HET_SNPS"])
        invariants[f"{name}/{arm}"] = {
            "original_reads_preserved": len(initial),
            "block_transforms": {ps: sorted(values) for ps, values in transforms.items()},
            "gaps_attempted": len(by_gap),
            "gaps_joined": sum(rows[-1]["STATUS"] == "joined" for rows in by_gap.values()),
        }
with (args.out / "results.tsv").open("w") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(["region", "arm", "evaluated_reads", "discordant_reads", "phase_sets", "read_switchflips"])
    writer.writerows(results)
(args.out / "invariants.json").write_text(json.dumps(invariants, indent=2) + "\n")
print((args.out / "results.tsv").read_text(), end="")
