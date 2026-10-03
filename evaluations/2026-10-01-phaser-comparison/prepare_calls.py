#!/usr/bin/env python3
"""Remove phase labels from pgphase calls without changing keys or genotypes."""
import argparse
from pathlib import Path

import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("source", type=Path)
parser.add_argument("output", type=Path, help="BGZF-compressed, indexed VCF")
args = parser.parse_args()
plain = args.output.with_suffix("")
with args.source.open() as source, plain.open("w") as output:
    for line in source:
        if line.startswith("#CHROM"):
            fields = line.rstrip("\n").split("\t")
            assert len(fields) == 10, "Expected one sample"
            fields[9] = "HG002"
        elif line.startswith("#"):
            output.write(line)
            continue
        else:
            fields = line.rstrip("\n").split("\t")
            names = fields[8].split(":")
            values = fields[9].split(":")
            values[names.index("GT")] = values[names.index("GT")].replace("|", "/")
            if "PS" in names:
                values[names.index("PS")] = "."
            fields[9] = ":".join(values)
        output.write("\t".join(fields) + "\n")
pysam.tabix_compress(str(plain), str(args.output), force=True)
pysam.tabix_index(str(args.output), preset="vcf", force=True)
