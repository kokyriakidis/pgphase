#!/usr/bin/env python3
"""Manage immutable competitor baselines and incremental pgphase panel runs."""

import argparse
import csv
import hashlib
import json
import os
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from string import Template


SAMPLE_BYTES = 1024 * 1024


def sampled_sha256(path):
    size = path.stat().st_size
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        offsets = sorted({0, max(0, size // 2 - SAMPLE_BYTES // 2),
                          max(0, size - SAMPLE_BYTES)})
        for offset in offsets:
            fh.seek(offset)
            digest.update(offset.to_bytes(8, "little"))
            digest.update(fh.read(SAMPLE_BYTES))
    return digest.hexdigest()


def full_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * SAMPLE_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(path, deep=False):
    path = Path(path).resolve()
    stat = path.stat()
    value = {
        "path": str(path),
        "size": stat.st_size,
        "sampled_sha256": sampled_sha256(path),
    }
    if deep:
        value["sha256"] = full_sha256(path)
    return value


def expand(value, context):
    previous = None
    current = value
    while current != previous:
        previous = current
        current = Template(current).safe_substitute(context)
    return current


def base_context(manifest):
    repo = Path(__file__).resolve().parents[1]
    context = {
        "HOME": str(Path.home()),
        "REPO": str(repo),
        "DATA_ROOT": os.environ.get(
            "DATA_ROOT", str(Path.home() / "Downloads/pgphase-eval-data")),
        "PHASER_BIN": os.environ.get(
            "PHASER_BIN", str(Path.home() / "micromamba/envs/bench-phasers/bin")),
    }
    context["SHARED_VCF_ROOT"] = os.environ.get(
        "SHARED_VCF_ROOT", f"{context['DATA_ROOT']}/shared_calls")
    context["OUT_ROOT"] = os.environ.get(
        "OUT_ROOT", f"{context['DATA_ROOT']}/results/chr12-18-20-comparison")
    for key, value in manifest.get("roots", {}).items():
        context[key.upper()] = expand(value, context)
    return context


def resolved_competitors(manifest):
    base = base_context(manifest)
    resolved = []
    for chrom, chrom_values in manifest["chromosomes"].items():
        chrom_context = dict(base, chrom=chrom)
        chrom_context.update({key: expand(value, chrom_context)
                              for key, value in chrom_values.items()})
        chrom_dir = f"{base['OUT_ROOT']}/{chrom}"
        for tool, spec in manifest["competitors"].items():
            if chrom not in spec.get("chromosomes", manifest["chromosomes"]):
                continue
            context = dict(chrom_context, tool=tool, chrom_dir=chrom_dir,
                           tool_dir=f"{chrom_dir}/{tool}")
            commands = [expand(command, context) for command in spec["commands"]]
            artifacts = [expand(path, context) for path in spec.get(
                "artifacts", manifest["competitor_artifacts"])]
            inputs = [expand(path, context) for path in spec.get("inputs", (
                "${linear_reference}", "${linear_bam}", "${shared_vcf}",
                "${truth_bam}", "${chrom_dir}/truth.vcf.gz",
            ))]
            resolved.append({
                "chromosome": chrom,
                "tool": tool,
                "version": spec["version"],
                "commands": commands,
                "inputs": inputs,
                "artifacts": artifacts,
            })
    return resolved


def resolved_chromosomes(manifest):
    base = base_context(manifest)
    values = {}
    for chrom, spec in manifest["chromosomes"].items():
        context = dict(base, chrom=chrom)
        context.update({key: expand(value, context) for key, value in spec.items()})
        values[chrom] = context
    return values


def competitors_for_chromosome(manifest, chrom):
    return [tool for tool, spec in manifest["competitors"].items()
            if chrom in spec.get("chromosomes", manifest["chromosomes"])]


def competitor_spec_hash(entries):
    payload = [{key: entry[key] for key in
                ("chromosome", "tool", "version", "commands", "inputs", "artifacts")}
               for entry in entries]
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def load_manifest(path):
    return json.loads(path.read_text())


def freeze(args, manifest):
    entries = resolved_competitors(manifest)
    locked = []
    for entry in entries:
        print(f"freeze {entry['chromosome']} {entry['tool']}")
        missing = [path for path in entry["inputs"] + entry["artifacts"]
                   if not Path(path).is_file()]
        if missing:
            raise FileNotFoundError("missing frozen artifact(s): " + ", ".join(missing))
        locked.append({
            **entry,
            "input_fingerprints": [fingerprint(path, args.deep)
                                   for path in entry["inputs"]],
            "artifact_fingerprints": [fingerprint(path, args.deep)
                                      for path in entry["artifacts"]],
        })
    lock = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "competitor_spec_sha256": competitor_spec_hash(entries),
        "deep_checksums": args.deep,
        "entries": locked,
    }
    args.lock.write_text(json.dumps(lock, indent=2) + "\n")
    print(f"wrote {args.lock}")


def verify_fingerprints(records, deep):
    errors = []
    for expected in records:
        path = Path(expected["path"])
        if not path.is_file():
            errors.append(f"missing: {path}")
            continue
        observed = fingerprint(path, deep and "sha256" in expected)
        for key in ("size", "sampled_sha256"):
            if observed[key] != expected[key]:
                errors.append(f"changed {key}: {path}")
        if deep and "sha256" in expected \
                and observed.get("sha256") != expected["sha256"]:
            errors.append(f"changed sha256: {path}")
    return errors


def verify(args, manifest, quiet=False):
    if not args.lock.is_file():
        raise FileNotFoundError(
            f"competitor lock missing: {args.lock}; run freeze-competitors once")
    lock = json.loads(args.lock.read_text())
    entries = resolved_competitors(manifest)
    errors = []
    if lock.get("competitor_spec_sha256") != competitor_spec_hash(entries):
        errors.append("manifest competitor specification differs from lock")
    for entry in lock.get("entries", []):
        errors.extend(verify_fingerprints(entry["input_fingerprints"], args.deep))
        errors.extend(verify_fingerprints(entry["artifact_fingerprints"], args.deep))
    if errors:
        raise RuntimeError("competitor verification failed:\n  " + "\n  ".join(errors))
    if not quiet:
        mode = "full SHA-256" if args.deep else "size + sampled SHA-256"
        print(f"verified {len(lock['entries'])} frozen competitor runs ({mode})")
    return lock


def show_commands(args, manifest):
    for entry in resolved_competitors(manifest):
        print(f"# {entry['chromosome']} {entry['tool']} {entry['version']}")
        for command in entry["commands"]:
            print(command)


def run_panel(args, manifest):
    lock = verify(args, manifest, quiet=True)
    repo = Path(__file__).resolve().parents[1]
    driver = repo / "evaluations/2026-09-12-chr12-18-20-comparison/commands.sh"
    env = dict(os.environ)
    env["CHROMS"] = " ".join(args.chromosome or manifest["chromosomes"].keys())
    env["RUN_COMPETITORS"] = "0"
    print(f"competitors locked at {lock['created_at']}; running pgphase/evaluation only")
    subprocess.run([str(driver)], cwd=repo, env=env, check=True)


def analyze_bridges(args, manifest):
    verify(args, manifest, quiet=True)
    repo = Path(__file__).resolve().parents[1]
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for chrom, context in resolved_chromosomes(manifest).items():
        root = Path(context["OUT_ROOT"]) / chrom
        output = output_dir / f"{chrom}.correct_bridges.tsv"
        command = [
            sys.executable, str(repo / "scripts/analyze_correct_competitor_bridges.py"),
            "--chromosome", chrom, "--root", str(root),
            "--truth-vcf", str(root / "truth.vcf.gz"),
            "--shared-vcf", context["shared_vcf"],
            "--catalog-vcf", context["catalog_vcf"],
            "--output", str(output),
        ]
        for tool in competitors_for_chromosome(manifest, chrom):
            command.extend(("--competitor", f"{tool}={root / tool}"))
        subprocess.run(command, cwd=repo, check=True)
        with output.open() as fh:
            all_rows.extend(csv.DictReader(fh, delimiter="\t"))

    combined = output_dir / "correct_bridges.tsv"
    fields = list(all_rows[0]) if all_rows else []
    with combined.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
        if fields:
            writer.writeheader()
            writer.writerows(all_rows)

    grouped = {}
    for row in all_rows:
        key = (row["tool"], row["primary_reason"])
        value = grouped.setdefault(key, {
            "tool": key[0], "primary_reason": key[1], "bridges": 0,
            "gap_bp": 0, "truth_hets": 0, "shared_hets": 0,
            "competitor_phased_sites": 0, "graph_phased_sites": 0,
            "private_sites": 0,
        })
        value["bridges"] += 1
        for field in ("gap_bp", "truth_hets", "shared_hets",
                      "competitor_phased_sites", "graph_phased_sites", "private_sites"):
            value[field] += int(row[field])
    summary_rows = [grouped[key] for key in sorted(grouped)]
    summary = output_dir / "correct_bridge_summary.tsv"
    summary_fields = list(summary_rows[0]) if summary_rows else []
    with summary.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, delimiter="\t", fieldnames=summary_fields, lineterminator="\n")
        if summary_fields:
            writer.writeheader()
            writer.writerows(summary_rows)
    print(f"wrote {combined} and {summary}")


def markdown_table(rows, fields, labels=None):
    labels = labels or fields
    lines = ["| " + " | ".join(labels) + " |",
             "|" + "|".join("---" for _ in fields) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(str(row[field]) for field in fields) + " |")
    return "\n".join(lines)


def report(args, manifest):
    lock = verify(args, manifest, quiet=True)
    repo = Path(__file__).resolve().parents[1]
    evaluation_dir = args.output_dir
    context = base_context(manifest)
    env = dict(os.environ, SUMMARY_OUT_DIR=str(evaluation_dir))
    subprocess.run([
        sys.executable, str(evaluation_dir / "summarize.py"), context["OUT_ROOT"],
        *manifest["chromosomes"].keys()], cwd=repo, env=env, check=True)
    analyze_bridges(args, manifest)

    with (evaluation_dir / "pooled.tsv").open() as fh:
        pooled = list(csv.DictReader(fh, delimiter="\t"))
    with (evaluation_dir / "correct_bridge_summary.tsv").open() as fh:
        bridges = list(csv.DictReader(fh, delimiter="\t"))
    with (evaluation_dir / "correct_bridges.tsv").open() as fh:
        bridge_details = list(csv.DictReader(fh, delimiter="\t"))
    for row in pooled:
        row["read_hamming_pct"] = f"{100 * float(row['read_hamming_rate']):.3f}%"
        row["variant_hamming_pct"] = f"{100 * float(row['variant_hamming_rate']):.3f}%"
        row["median_ngc50_kb"] = f"{int(row['median_chromosome_ngc50_bp']) / 1000:.0f}"
    hiphase_bridges = [row for row in bridge_details if row["tool"] == "hiphase"]
    hiphase_reasons = {}
    for row in hiphase_bridges:
        hiphase_reasons[row["primary_reason"]] = \
            hiphase_reasons.get(row["primary_reason"], 0) + 1
    hiphase_top = sorted(
        hiphase_bridges, key=lambda row: int(row["gap_bp"]), reverse=True)[:15]
    for row in hiphase_top:
        row["region"] = f"{row['chromosome']}:{row['gap_start']}-{row['gap_end']}"
        row["accuracy_pct"] = f"{100 * float(row['competitor_accuracy']):.2f}%"
    hiphase_total = len(hiphase_bridges)
    hiphase_repeat = hiphase_reasons.get("repeat_indels_excluded", 0)
    hiphase_stitch = hiphase_reasons.get("block_stitching", 0)
    hiphase_competitor_sites = sum(
        int(row["competitor_phased_sites"]) for row in hiphase_bridges)
    hiphase_graph_sites = sum(int(row["graph_phased_sites"])
                              for row in hiphase_bridges)
    content = [
        "# Generated Phasing Benchmark Report",
        "",
        f"Baseline frozen: `{lock['created_at']}`",
        f"Panel: `{manifest['panel']}`",
        f"Competitor lock: `{lock['competitor_spec_sha256']}`",
        "",
        "## Pooled Results",
        "",
        markdown_table(
            pooled,
            ["tool", "callset", "chromosomes", "assessed_pairs",
             "variant_hamming_pct", "phased_reads", "read_hamming_pct",
             "median_ngc50_kb"],
            ["method", "callset", "chromosomes", "assessed pairs",
             "variant Hamming", "phased reads", "read Hamming",
             "median chr NGC50 (kb)"]),
        "",
        "The frozen LongPhase baseline is the measured `--pb` SNP mode; it does",
        "not include LongPhase's optional `--indels` mode.",
        "LongcallD is a chr20-only native-caller control. Its contiguity is not",
        "directly comparable to methods evaluated on the shared DeepVariant VCF.",
        "",
        "## Correct Competitor Bridges",
        "",
        "Blocks require >=99% read-truth accuracy, >=50 evaluated reads, and no",
        "variant switch-error interval across the pgphase gap.",
        "",
        markdown_table(
            bridges,
            ["tool", "primary_reason", "bridges", "gap_bp",
             "competitor_phased_sites", "graph_phased_sites", "private_sites"],
            ["tool", "pgphase break reason", "bridges", "gap bp",
             "competitor sites", "graph sites", "private sites"]),
        "",
        "## Why HiPhase Is More Contiguous",
        "",
        (f"HiPhase correctly bridges {hiphase_total} graph breaks under the strict "
         f"accuracy rule. Their median width is "
         f"{median(int(row['gap_bp']) for row in hiphase_bridges):,.0f} bp. "
         f"At these breaks HiPhase phases {hiphase_competitor_sites:,} linking "
         f"sites while pgphase phases {hiphase_graph_sites:,}."),
        "",
        (f"The dominant mechanism is repeat-indel exclusion: {hiphase_repeat} "
         f"bridges ({100 * hiphase_repeat / hiphase_total:.1f}%) contain repeat "
         "heterozygous indel candidates but no graph-phased anchor. Another "
         f"{hiphase_stitch} ({100 * hiphase_stitch / hiphase_total:.1f}%) already "
         "contain graph-phased sites and fail at block stitching. The remaining "
         "breaks are clean candidates left unphased or catalog sites that never "
         "become candidates."),
        "",
        "### Largest Correct HiPhase Bridges",
        "",
        markdown_table(
            hiphase_top,
            ["region", "gap_bp", "accuracy_pct", "primary_reason", "shared_hets",
             "competitor_phased_sites", "graph_phased_sites", "repeat_het_indels"],
            ["region", "gap bp", "read accuracy", "reason", "shared hets",
             "HiPhase sites", "graph sites", "repeat indels"]),
        "",
        "## Reproduction",
        "",
        "```bash",
        "python3 scripts/benchmark_panel.py verify-competitors",
        "python3 scripts/benchmark_panel.py run",
        "python3 scripts/benchmark_panel.py report",
        "```",
        "",
        "Normal `run` mode verifies and reuses frozen competitors. It cannot",
        "invoke competitor phasers unless the lower-level driver is explicitly",
        "called with `RUN_COMPETITORS=1`.",
    ]
    (evaluation_dir / "REPORT.md").write_text("\n".join(content) + "\n")
    print(f"wrote {evaluation_dir / 'REPORT.md'}")


def status(args, manifest):
    lock = verify(args, manifest, quiet=True)
    print(f"panel\t{manifest['panel']}")
    print(f"lock_created\t{lock['created_at']}")
    print(f"frozen_runs\t{len(lock['entries'])}")
    by_tool = {}
    for entry in lock["entries"]:
        by_tool.setdefault(entry["tool"], entry["version"])
    for tool, version in by_tool.items():
        print(f"competitor\t{tool}\t{version}\tfrozen")


def parse_args():
    repo = Path(__file__).resolve().parents[1]
    default_dir = repo / "evaluations/2026-09-12-chr12-18-20-comparison"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=default_dir / "panel.json")
    parser.add_argument("--lock", type=Path, default=default_dir / "competitor_lock.json")
    parser.add_argument("--deep", action="store_true",
                        help="create or verify full-file SHA-256 checksums")
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("freeze-competitors")
    sub.add_parser("verify-competitors")
    sub.add_parser("show-competitor-commands")
    run = sub.add_parser("run")
    run.add_argument("--chromosome", action="append",
                     choices=("chr12", "chr18", "chr20"))
    analyze = sub.add_parser("analyze")
    analyze.add_argument("--output-dir", type=Path, default=default_dir)
    report_parser = sub.add_parser("report")
    report_parser.add_argument("--output-dir", type=Path, default=default_dir)
    sub.add_parser("status")
    return parser.parse_args()


def main():
    args = parse_args()
    manifest = load_manifest(args.manifest)
    actions = {
        "freeze-competitors": freeze,
        "verify-competitors": verify,
        "show-competitor-commands": show_commands,
        "run": run_panel,
        "analyze": analyze_bridges,
        "report": report,
        "status": status,
    }
    actions[args.action](args, manifest)


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
