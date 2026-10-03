#!/usr/bin/env python3
"""Score fresh chr20 runs on a common BAM population and parental truth map."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[2]
SHARED = Path.home() / "Downloads/pgphase-eval-data/shared_calls/chr20"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("run_root", type=Path)
args = parser.parse_args()
truth = {}
for line in (ROOT / "test_data/derived/chr20_truth_hap.tsv").read_text().splitlines():
    name, parent = line.split("\t")
    if parent in ("PATERNAL", "MATERNAL"):
        truth[name] = parent == "PATERNAL"


def input_population(path):
    """Hash SAM records with only contig prefixes normalized, retaining tags."""
    names, primary, flags = set(), set(), Counter()
    digest = hashlib.sha256()
    with pysam.AlignmentFile(path) as bam:
        for read in bam:
            fields = read.to_string().split("\t")
            for column in (2, 6):
                fields[column] = fields[column].split("#")[-1]
            digest.update(("\t".join(fields) + "\n").encode())
            flags[read.flag] += 1
            if read.is_unmapped or read.reference_name.split("#")[-1] != "chr20":
                continue
            names.add(read.query_name)
            if not read.is_secondary and not read.is_supplementary:
                primary.add(read.query_name)
    return names, primary, dict(flags), digest.hexdigest()


def read_tags(path):
    """Prefer a primary alignment; no alignment is counted as another read."""
    selected, records = {}, Counter()
    with pysam.AlignmentFile(path, check_sq=False) as bam:
        tag_only = not bam.references
        for read in bam:
            rank = 2 if read.is_secondary else 1 if read.is_supplementary else 0
            records[rank] += 1
            if read.is_unmapped and not tag_only:
                continue
            hp = read.get_tag("HP") if read.has_tag("HP") else 0
            ps = read.get_tag("PS") if read.has_tag("PS") else 0
            chrom = "chr20" if tag_only else read.reference_name.split("#")[-1]
            tag = (chrom, int(ps), int(hp))
            old = selected.get(read.query_name)
            if old is None or rank < old[0]:
                selected[read.query_name] = (rank, tag)
            elif rank == old[0]:
                assert tag == old[1], f"Conflicting equally ranked tags: {read.query_name}"
    return {name: tag for name, (_, tag) in selected.items()}, dict(records)


def evaluate(tags, population):
    votes, phased, scored = defaultdict(lambda: [0, 0]), set(), set()
    for name in population:
        chrom, ps, hp = tags.get(name, ("chr20", 0, 0))
        if hp not in (1, 2) or ps <= 0:
            continue
        phased.add(name)
        if name in truth:
            votes[(chrom, ps)][int((hp == 1) != truth[name])] += 1
            scored.add(name)
    correct = sum(max(v) for v in votes.values())
    errors = len(scored) - correct
    flips = {block: v[1] > v[0] for block, v in votes.items()}
    correct_names = {
        name for name in scored
        if ((tags[name][2] == 1) != truth[name]) == flips[tags[name][:2]]
    }
    assert len(correct_names) == correct
    return dict(input_reads=len(population), phased_reads=len(phased),
                phased_percent=100 * len(phased) / len(population),
                unphased_reads=len(population) - len(phased),
                truth_scored=len(scored), truth_unavailable=len(phased - scored),
                truth_correct=correct, truth_errors=errors,
                accuracy_percent=100 * correct / len(scored) if scored else 0,
                correct_yield_percent=100 * correct / len(population),
                read_phase_sets=len(votes)), phased, correct_names


def vcf_metrics(path):
    blocks, calls = defaultdict(list), {}
    hets = phased = 0
    with pysam.VariantFile(path) as vcf:
        for record in vcf:
            sample = record.samples[0]
            gt = sample.get("GT")
            key = (record.chrom.split("#")[-1], record.pos, record.ref, record.alts)
            assert key not in calls, f"Duplicate exact VCF key: {key}"
            calls[key] = tuple(sorted(gt, key=lambda x: -1 if x is None else x))
            if gt is None or len(gt) != 2 or None in gt or gt[0] == gt[1]:
                continue
            hets += 1
            ps = sample.get("PS")
            if sample.phased and ps is not None and ps > 0:
                phased += 1
                blocks[(key[0], ps)].append(record.pos)
    spans = sorted((max(p) - min(p) + 1 for p in blocks.values()), reverse=True)
    half, running, n50 = sum(spans) / 2, 0, 0
    for span in spans:
        running += span
        if running >= half:
            n50 = span
            break
    return dict(vcf_records=len(calls), hets=hets, phased_hets=phased,
                vcf_blocks=len(blocks), span_n50_bp=n50,
                largest_block_bp=max(spans, default=0)), calls


def resources(path):
    fields = {}
    for line in path.read_text().splitlines():
        if ": " in line:
            name, value = line.strip().rsplit(": ", 1)
            fields[name] = value
    elapsed = fields["Elapsed (wall clock) time (h:mm:ss or m:ss)"]
    seconds = 0.0
    for part in elapsed.split(":"):
        seconds = seconds * 60 + float(part)
    assert fields["Exit status"] == "0", path
    return dict(wall_seconds=seconds,
                user_seconds=float(fields["User time (seconds)"]),
                system_seconds=float(fields["System time (seconds)"]),
                max_rss_kib=int(fields["Maximum resident set size (kbytes)"]))


annotated = ROOT / "test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
population, primary, flags, digest = input_population(annotated)
linear_population, linear_primary, linear_flags, linear_digest = input_population(
    SHARED / "HG002.chr20.normalized.bam")
assert (population, primary, flags, digest) == (
    linear_population, linear_primary, linear_flags, linear_digest), "BAM records differ"
with pysam.FastaFile(str(ROOT / "test_data/chm13v2.0.chr20.renamed.fa")) as a, \
     pysam.FastaFile(str(SHARED / "chm13v2.0.chr20.normalized.fa")) as b:
    assert a.fetch("CHM13#0#chr20") == b.fetch("chr20"), "References differ"

result = dict(input_unique_reads=len(population), primary_unique_reads=len(primary),
              normalized_sam_sha256=digest, bam_records_identical=True,
              reference_sequences_identical=True, input_flag_counts=flags, tools={})
tags_by_tool, phased_by_tool, correct_by_tool, calls_by_tool = {}, {}, {}, {}
for tool in ("pgphase", "hiphase_dv", "hiphase_pg_calls"):
    directory = args.run_root / tool
    tags, counts = read_tags(directory / "phased.bam")
    metrics, phased_names, correct_names = evaluate(tags, population)
    primary_metrics, _, _ = evaluate(tags, primary)
    vcf_path = directory / ("phased.vcf" if tool == "pgphase" else "phased.vcf.gz")
    variants, calls = vcf_metrics(vcf_path)
    metrics.update(variants)
    metrics.update(resources(directory / "resources.txt"))
    metrics["primary_only"] = primary_metrics
    metrics["output_alignment_counts_by_rank"] = counts
    result["tools"][tool] = metrics
    tags_by_tool[tool], phased_by_tool[tool] = tags, phased_names
    correct_by_tool[tool], calls_by_tool[tool] = correct_names, calls

_, transferred_calls = vcf_metrics(args.run_root / "hiphase_pg_calls/input.vcf.gz")
assert calls_by_tool["pgphase"] == transferred_calls == calls_by_tool["hiphase_pg_calls"], \
    "Same-callset arm changed keys or unordered genotypes"
result["same_callset_preserved"] = True
materialized_tags, _ = read_tags(args.run_root / "pgphase/haplotagged.bam")
assert {name: tag for name, tag in tags_by_tool["pgphase"].items()
        if tag[2] in (1, 2) and tag[1] > 0} == {
            name: tag for name, tag in materialized_tags.items()
            if tag[2] in (1, 2) and tag[1] > 0}, "Tag materialization changed assignments"
materialization = resources(args.run_root / "pgphase/materialize_resources.txt")
result["tools"]["pgphase"]["materialization"] = materialization
result["tools"]["pgphase"]["full_bam_wall_seconds"] = (
    result["tools"]["pgphase"]["wall_seconds"] + materialization["wall_seconds"])
result["materialized_tags_preserved"] = True
result["overlap"] = {}
for competitor in ("hiphase_dv", "hiphase_pg_calls"):
    pg, other = phased_by_tool["pgphase"], phased_by_tool[competitor]
    common = pg & other
    result["overlap"][competitor] = dict(
        both_phased=len(common), pgphase_only=len(pg - other),
        competitor_only=len(other - pg),
        pgphase_only_correct=len((pg - other) & correct_by_tool["pgphase"]),
        competitor_only_correct=len((other - pg) & correct_by_tool[competitor]),
        common_full_block_accuracy={
            tool: 100 * len(common & correct_by_tool[tool]) / len(common)
            for tool in ("pgphase", competitor)},
        common_independently_oriented={
            tool: evaluate(tags_by_tool[tool], common)[0]
            for tool in ("pgphase", competitor)})
(args.run_root / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
