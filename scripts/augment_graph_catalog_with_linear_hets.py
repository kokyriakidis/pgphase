#!/usr/bin/env python3
"""Add high-confidence linear-VCF hets to a graph-site catalog for experiments.

The synthetic traversal strings make the records structurally valid graph-site
rows, but do not correspond to GAF observations. In hybrid mode their allele
support therefore comes only from the surjected BAM and is re-gated by pgphase.
"""

import argparse
import gzip
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_linear_hets(path, contig, min_gq, include_filtered):
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0].split("#")[-1] != contig or len(fields) < 10:
                continue
            if (not include_filtered and fields[6] not in ("PASS", ".")) \
                    or "," in fields[4]:
                continue
            if not fields[3] or not fields[4] or any(
                    base not in "ACGT" for base in fields[3] + fields[4]):
                continue
            sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
            gt = sample.get("GT", "").replace("|", "/").split("/")
            if len(gt) != 2 or set(gt) != {"0", "1"}:
                continue
            try:
                gq = int(sample.get("GQ", "0"))
            except ValueError:
                gq = 0
            if gq < min_gq:
                continue
            rows.append((int(fields[1]), fields[3], fields[4], gq))
    rows.sort()
    return rows


def read_graph_alleles(path):
    alleles = set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            pos = int(fields[1])
            for alt in fields[4].split(","):
                alleles.add((pos, fields[3], alt))
    return alleles


def synthetic_row(chrom, index, site):
    pos, ref, alt, gq = site
    left = f"linear{index}L"
    right = f"linear{index}R"
    middle = f"linear{index}A"
    site_id = f"linear_{pos}_{index}"
    traversals = f">{left}>{right},>{left}>{middle}>{right}"
    end = pos + max(len(ref), 1) - 1
    info = f"AT={traversals};LV=0;RC={chrom};RS={pos};RD={end};LINEAR_GQ={gq}"
    return "\t".join((chrom, str(pos), site_id, ref, alt, "60", "PASS", info))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-sites", required=True)
    parser.add_argument("--linear-vcf", required=True)
    parser.add_argument("--contig", required=True)
    parser.add_argument("--min-gq", type=int, default=20)
    parser.add_argument("--include-filtered", action="store_true")
    parser.add_argument(
        "--include-graph-matches", action="store_true",
        help="Also inject linear alleles already represented exactly in the graph.",
    )
    parser.add_argument(
        "--count-only", action="store_true",
        help="Report the number of selected linear hets without writing a catalog.",
    )
    parser.add_argument("--output")
    args = parser.parse_args()

    linear = parse_linear_hets(
        args.linear_vcf, args.contig, args.min_gq, args.include_filtered)
    if not args.include_graph_matches:
        graph_alleles = read_graph_alleles(args.graph_sites)
        linear = [site for site in linear if site[:3] not in graph_alleles]
    if args.count_only:
        print(f"linear_hets_selected={len(linear)} min_gq={args.min_gq} "
              f"include_filtered={args.include_filtered} "
              f"include_graph_matches={args.include_graph_matches}")
        return
    if not args.output:
        parser.error("--output is required unless --count-only is used")
    linear_i = 0
    graph_chrom = None
    graph_records = 0

    with open_text(args.graph_sites) as graph, Path(args.output).open("w") as out:
        for line in graph:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    out.write("##INFO=<ID=LINEAR_GQ,Number=1,Type=Integer,"
                              "Description=\"Source linear-call genotype quality\">\n")
                out.write(line)
                continue
            fields = line.split("\t", 2)
            pos = int(fields[1])
            if graph_chrom is None:
                graph_chrom = fields[0]
            while linear_i < len(linear) and linear[linear_i][0] <= pos:
                out.write(synthetic_row(graph_chrom, linear_i, linear[linear_i]) + "\n")
                linear_i += 1
            out.write(line)
            graph_records += 1

        if graph_chrom is None:
            raise SystemExit("graph catalog contains no records")
        while linear_i < len(linear):
            out.write(synthetic_row(graph_chrom, linear_i, linear[linear_i]) + "\n")
            linear_i += 1

    print(f"graph_records={graph_records} linear_hets_added={len(linear)} "
          f"min_gq={args.min_gq} include_filtered={args.include_filtered}")


if __name__ == "__main__":
    main()
