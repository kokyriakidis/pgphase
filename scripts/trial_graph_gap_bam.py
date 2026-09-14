#!/usr/bin/env python3
"""Compare recovery strategies on a fixed panel, reusing per-region evidence."""

import argparse
import concurrent.futures
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time


def run(command, directory, name):
    (directory / f"{name}.command.json").write_text(json.dumps(command, indent=2) + "\n")
    started = time.monotonic()
    with (directory / f"{name}.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    return round(time.monotonic() - started, 3)


def endpoint_blocks(path, left, right):
    blocks = {left: set(), right: set()}
    with path.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            # Native indel POS names the event; audited VCF POS names its
            # preceding anchor, as in VariantKey::sort_pos().
            pos = int(row["POS"]) - (row["TYPE"] != "SNP")
            if (pos in blocks and row["HAP_ALT"] in ("1", "2") and
                    row["HAP_REF"] in ("1", "2") and row["HAP_ALT"] != row["HAP_REF"] and
                    int(row["PHASE_SET"]) >= 0):
                blocks[pos].add(row["PHASE_SET"])
    if len(blocks[left]) != 1 or len(blocks[right]) != 1:
        return "unresolved_endpoint"
    return "joined" if blocks[left] == blocks[right] else "split"


def trial_region(args, target):
    chrom, left, right = target
    key = f"{chrom.replace('#', '_')}_{left}_{right}"
    directory = args.output / "runs" / args.run_name / key
    directory.mkdir(parents=True, exist_ok=False)
    cache = args.output / "cache" / f"{key}.gapev"
    cache.parent.mkdir(parents=True, exist_ok=True)
    region = f"{args.contig_prefix}{chrom}:{max(1, left - args.flank)}-{right + args.flank}"
    common = [str(args.binary), "collect-hybrid-variation", "--ref", str(args.ref),
              "--bam", str(args.bam), "--graph-sites", str(args.graph_sites),
              "--gaf", str(args.gaf), "-r", region, "--chunk-size", "500000",
              "--threads", str(args.threads), "--link-by-alleles", "--block-link-window", "8",
              "--min-read-margin", "2", "--recover-gaps", "--gap-evidence-cache", str(cache)]
    result = {"chrom": chrom, "left": left, "right": right}
    for arm in ("baseline", "graph_bam"):
        out = directory / arm
        out.mkdir()
        command = common + ["-o", str(out / "candidates.tsv"), "--phased-vcf-out",
                            str(out / "native.vcf"), "-b", str(out / "phased.bam"),
                            "--gap-recovery-report", str(out / "tiers.tsv")]
        if arm == "baseline":
            command.append("--no-graph-gap-bam")
        result[f"{arm}_cache_reused"] = cache.exists()
        result[f"{arm}_seconds"] = run(command, out, "phase")
        result[f"{arm}_status"] = endpoint_blocks(out / "candidates.tsv", left, right)
        with (out / "tiers.tsv").open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        selected = [r for r in rows if r["GRAPH_BAM_PASS"] == "1"]
        result[f"{arm}_selected_reads"] = max(
            (int(r["SELECTED_GRAPH_READS"]) for r in selected), default=0)
        if args.truth_bam:
            evaluation = out / "read_eval"
            evaluation.mkdir()
            evaluator = Path(__file__).resolve().parent / "evaluate_phase_accuracy.py"
            run(["python3", str(evaluator), str(out / "phased.bam"), str(args.truth_bam),
                 "0", "0", "5", "", str(evaluation), "samtools", "", "", "", ""], out, "evaluate")
            summary = json.loads((evaluation / "summary.json").read_text())
            for field in ("discordant_reads", "total_reads_evaluated", "hamming_error_rate",
                          "switchflip_errors"):
                result[f"{arm}_{field}"] = summary[field]
    if args.truth_bam:
        result["discordant_delta"] = result["graph_bam_discordant_reads"] - result["baseline_discordant_reads"]
    (directory / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, required=True,
                        help="TSV with chromosome,left,right; also accepts audit_panel.py missed.tsv")
    parser.add_argument("--binary", type=Path, default=Path("./pgphase"))
    for flag in ("ref", "bam", "graph-sites", "gaf", "output"):
        parser.add_argument(f"--{flag}", type=Path, required=True)
    parser.add_argument("--truth-bam", type=Path, help="Evaluation only; never passed to pgphase")
    parser.add_argument("--run-name", required=True, help="New label for this build; caches are shared between runs")
    parser.add_argument("--contig-prefix", default="CHM13#0#")
    parser.add_argument("--flank", type=int, default=50000)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--limit", type=int, default=0, help="Maximum unique gaps; zero selects all")
    args = parser.parse_args()
    if args.jobs < 1 or args.threads < 1 or args.flank < 1 or args.limit < 0:
        parser.error("jobs, threads and flank must be positive; limit must be nonnegative")
    if Path(args.run_name).name != args.run_name or args.run_name in (".", ".."):
        parser.error("run-name must be a single directory name")
    for name in ("binary", "ref", "bam", "graph_sites", "gaf", "output", "targets", "truth_bam"):
        if getattr(args, name) is not None:
            setattr(args, name, getattr(args, name).resolve())
    with args.targets.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    targets = set()
    for row in rows:
        if "kind" in row and (row["kind"] != "split_block" or
                               row["competitor_pair_truth"] != "concordant" or
                               float(row["competitor_accuracy"]) < 0.99):
            continue
        left, right = int(row["left"]), int(row["right"])
        if right > left:
            targets.add((row["chromosome"], left, right))
    targets = sorted(targets)
    if args.limit:
        targets = targets[:args.limit]
    if not targets:
        parser.error("no eligible gaps in target table")
    run_dir = args.output / "runs" / args.run_name
    run_dir.mkdir(parents=True, exist_ok=False)
    manifest = {"binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(),
                "targets": targets, "args": {k: str(v) if isinstance(v, Path) else v
                                             for k, v in vars(args).items()},
                "note": "Local windows can differ from whole-chromosome phasing. Confirm improvements on chr20."}
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    results, failures = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        pending = {executor.submit(trial_region, args, target): target for target in targets}
        for future in concurrent.futures.as_completed(pending):
            try:
                result = future.result()
                results.append(result)
                print(json.dumps(result), flush=True)
            except Exception as error:
                failures.append({"target": pending[future], "error": str(error)})
                print(json.dumps(failures[-1]), flush=True)
    results.sort(key=lambda r: (r["chrom"], r["left"], r["right"]))
    if results:
        with (run_dir / "results.tsv").open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(results[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(results)
    (run_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
