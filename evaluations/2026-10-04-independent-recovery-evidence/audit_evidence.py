#!/usr/bin/env python3
"""Check independent allele preservation and working conflict persistence.

Only diagnostic matrices are read; this audit needs no truth or competitor.
Counts include each solve/stage, so they are not unique chromosome calls.
"""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def rows(path):
    with path.open() as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def matrix(path):
    sites, calls = {}, {}
    with path.open() as stream:
        for line in stream:
            f = line.rstrip("\n").split("\t")
            if f[0] == "VAR":
                sites[f[1]] = tuple(f[i] for i in (2, 3, 6, 7))
            elif f[0] == "OBS":
                calls[(f[1], sites[f[2]])] = int(f[5])
    return calls


def source_key(row, allele):
    return tuple(row[name] for name in
                 ("solve", "pos", "type", "ref_len", "alt", "qname")) + (int(row[allele]),)


def audit(folder):
    counts, errors, census = Counter(), [], Counter()
    for trace in sorted(folder.glob("*.transfer*.tsv")):
        suffix = ".bam-source-evidence-retry.tsv" if trace.name.endswith("transfer-retry.tsv") \
                 else ".bam-source-evidence.tsv"
        stem = trace.name.split(".transfer")[0]
        snapshot = folder / (stem + suffix)
        if not snapshot.exists():
            errors.append({"missing_snapshot": snapshot.name})
            continue
        independent = set()
        for row in rows(snapshot):
            if row["kind"] == "call":
                independent.add(source_key(row, "allele"))
            else:
                counts["saved_" + row["kind"]] += 1
        counts["independent_calls"] += len(independent)
        expected, actual = defaultdict(set), {}
        for row in rows(trace):
            census[row["status"]] += 1
            if row["phase_anchor"] == "1":
                counts["checked_source_calls"] += 1
                if source_key(row, "source_allele") not in independent:
                    counts["lost_independent_calls"] += 1
                    if len(errors) < 20:
                        errors.append({"missing_source_call": source_key(row, "source_allele")})
            if row["status"] != "mapped":
                continue
            key = row["qname"], row["destination_index"]
            expected[key].add(int(row["source_allele"]))
            actual[key] = int(row["bam_allele"])
        for key, alleles in expected.items():
            value = actual[key]
            if value == -2:
                counts["working_abstentions"] += 1
            elif value in alleles:
                counts["working_retained"] += 1
                counts["owner_resolved_disagreement"] += len(alleles) > 1
            else:
                counts["invalid_working_calls"] += 1
                if len(errors) < 20:
                    errors.append({"invalid_working": key, "sources": sorted(alleles),
                                   "actual": value})
            counts["source_disagreements"] += len(alleles) > 1
    for path in sorted(folder.glob("*.msa-admission.tsv")):
        for row in rows(path):
            kind = "msa" if row["update_counts"] == "1" else "physical"
            census[kind + ":" + row["status"]] += 1
            if kind == "msa" and row["status"] != "eligible":
                counts["rejected_verified_recalls"] += 1
    for before in sorted(folder.glob("*.bam-overlay-input.tsv")):
        after = before.with_name(before.name.replace("bam-overlay-input", "bam-overlay-output"))
        old, new = matrix(before), matrix(after)
        for key, value in old.items():
            if value < 0 and value != -2:
                continue
            counts["overlay_conflicts" if value == -2 else "overlay_calls"] += 1
            if new.get(key, -1) != value:
                counts["overwritten_overlay_calls"] += 1
                if len(errors) < 20:
                    errors.append({"overwrite": key, "before": value, "after": new.get(key, -1)})
    return {"directory": str(folder), "counts": dict(counts),
            "census": dict(census), "errors": errors}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()
    folders = sorted({path.parent for path in args.root.rglob("*.transfer*.tsv")})
    if not folders:
        parser.error(f"no transfer diagnostics under {args.root}")
    results = [audit(folder) for folder in folders]
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    for result in results:
        print(result["directory"], result["counts"])
    if args.strict and any(result["errors"] or
                           result["counts"].get("rejected_verified_recalls", 0)
                           for result in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
