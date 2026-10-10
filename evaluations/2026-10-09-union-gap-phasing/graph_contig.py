#!/usr/bin/env python3
"""How many links of the graph-only phasing does a run keep? (eval only)
Usage: graph_contig.py INV_DIR RUN_DIR [RUN_DIR ...]
A link is a pair of consecutive phased hets in the same graph-only phase set;
it is kept when one block of the run spans both positions."""
import sys
inv = sys.argv[1]; sys.path.insert(0, inv)
import audit_hybrid as A
g, _ = A.load_vcf(f"{inv}/graph/phased.vcf")
g.sort()
links = [(a[0], b[0]) for a, b in zip(g, g[1:]) if a[2] == b[2] and b[0] > a[0]]
for d in sys.argv[2:]:
    u, _ = A.load_vcf(f"{d}/phased.vcf")
    cover = A.covering(A.blocks_of(u))
    broken = [(a, b) for a, b in links if cover(a, b) is None]
    cen = sum(1 for a, b in broken if 26e6 <= a < 32e6)
    print(f"{d.split('/')[-1]}: graph links {len(links)}, broken {len(broken)} (centromeric {cen}, arms {len(broken) - cen})")
