#!/usr/bin/env python3
"""Copy pgphase's tag-only assignments onto the original BAM, then index it."""
import argparse
from pathlib import Path

import pysam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("tags", type=Path)
parser.add_argument("source", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--io-threads", type=int, default=4)
args = parser.parse_args()
tags = {}
with pysam.AlignmentFile(args.tags, check_sq=False) as source:
    for read in source:
        hp = read.get_tag("HP") if read.has_tag("HP") else 0
        ps = read.get_tag("PS") if read.has_tag("PS") else 0
        if hp in (1, 2) and ps > 0:
            assert read.query_name not in tags, "Duplicate assignment"
            tags[read.query_name] = (hp, ps)
with pysam.AlignmentFile(args.source, threads=args.io_threads) as source, \
     pysam.AlignmentFile(args.output, "wb", template=source,
                         threads=args.io_threads) as output:
    for read in source:
        for tag in ("HP", "PS"):
            if read.has_tag(tag):
                read.set_tag(tag, None)
        phase = tags.get(read.query_name)
        if phase:
            read.set_tag("HP", phase[0], value_type="i")
            read.set_tag("PS", phase[1], value_type="i")
        output.write(read)
pysam.index("-@", str(args.io_threads), str(args.output))
