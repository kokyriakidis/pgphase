#!/usr/bin/env python3
"""Audit output allele identity and saved internal candidate-key collisions."""
import collections
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

import pysam


ROOT = Path(__file__).resolve().parents[2]
VCF = ROOT / "test_data/tmp_gap_fix91/final_current/0/phased.vcf"
REF = ROOT / "test_data/chm13v2.0.chr20.renamed.fa"
BCFTOOLS = os.environ.get(
    "BCFTOOLS", "/home/kokyriakidis/micromamba/envs/bench-phasers/bin/bcftools")


def record_description(record):
    sample = next(iter(record.samples.values()))
    return dict(pos=record.pos, ref=record.ref, alt=record.alts,
                gt=sample.get("GT"), phased=sample.phased,
                ps=sample.get("PS"), ad=sample.get("AD"))


def identity(record):
    return record.contig, record.pos, record.ref, record.alts


original = list(pysam.VariantFile(str(VCF)))
assert all(len(r.alts) == 1 for r in original)
with tempfile.TemporaryDirectory(prefix="pgphase-representation-") as directory:
    normalized_path = Path(directory) / "normalized.vcf"
    command = [BCFTOOLS, "norm", "-f", str(REF), "-m", "-any",
               "--old-rec-tag", "PGPHASE_ORIGINAL", "-Ov",
               "-o", str(normalized_path), str(VCF)]
    normalization = subprocess.run(command, check=True, capture_output=True, text=True)
    normalized = list(pysam.VariantFile(str(normalized_path)))
assert len(original) == len(normalized)
raw_counts = collections.Counter(identity(r) for r in original)
assert max(raw_counts.values()) == 1
original_by_identity = {identity(r): r for r in original}
normalized_groups = collections.defaultdict(list)
for after in normalized:
    source = after.info.get("PGPHASE_ORIGINAL")
    if source:
        chrom, pos, ref, alt = source.split("|")[:4]
        before = original_by_identity[(chrom, int(pos), ref, tuple(alt.split(",")))]
    else:
        before = original_by_identity[identity(after)]
    assert before.contig == after.contig
    assert record_description(before)["gt"] == record_description(after)["gt"]
    normalized_groups[identity(after)].append(
        dict(original=record_description(before), normalized=record_description(after)))
duplicates = [g for g in normalized_groups.values() if len(g) > 1]
result = dict(vcf=str(VCF.relative_to(ROOT)),
              vcf_sha256=hashlib.sha256(VCF.read_bytes()).hexdigest(),
              bcftools_version=subprocess.run([BCFTOOLS, "--version"], check=True,
                  capture_output=True, text=True).stdout.splitlines()[0],
              normalization_stderr=normalization.stderr.strip(),
              rows=len(original),
              raw_duplicate_groups=sum(n > 1 for n in raw_counts.values()),
              normalized_duplicate_groups=len(duplicates),
              normalized_extra_rows=sum(len(g) - 1 for g in duplicates),
              duplicate_groups=duplicates, chunks=[])
for chunk in (4, 5, 6):
    path = ROOT / f"test_data/tmp_joint_evidence_investigation/{chunk}/matrix.chunk0.bam-overlay-output.tsv"
    groups = collections.defaultdict(list)
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if fields[0] == "VAR":
            groups[(fields[2], fields[3], fields[6], fields[7])].append(fields)
    collisions = [g for g in groups.values() if len(g) > 1]
    result["chunks"].append(dict(chunk=chunk,
        matrix_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        candidates=sum(map(len, groups.values())),
        raw_key_collision_groups=len(collisions),
        extra_rows=sum(len(g) - 1 for g in collisions),
        example=[dict(index=int(a[1]), pos=int(a[2]), type=a[3], ref_len=int(a[6]),
                      alt=a[7], ps=int(a[8]), ref_cov=int(a[14]), alt_cov=int(a[15]))
                 for a in next((g for g in collisions if g[0][2] == "5331128"),
                               collisions[0])]))
print(json.dumps(result, indent=2))
