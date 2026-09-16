# Transferring BAM-channel phasing into the graph's gauge

Instead of discovering sites in a gap, phase the gap window with the BAM channel
that already exists (`collect-bam-variation`, with all its classification and MSA
machinery), then carry the result into the graph's gauge. Read identity is an
exact bridge: a read carrying both a BAM haplotype and a graph phase set says how
the two gauges line up.

```sh
./evaluations/2026-09-15-bam-phase-transfer/run.sh          # 8 gap windows, ~2 s each
python3 evaluations/2026-09-15-bam-phase-transfer/transfer.py \
  --transfer-root /tmp/bam-phase-transfer \
  --graph-reads /tmp/graph-only-mapq/q1/phase_reads.tsv \
  --gaps evaluations/2026-09-15-graph-only-baseline/graph_only_gaps.tsv \
  --truth-map /tmp/truth_hap.tsv --output transfer_results.tsv
```

Windows are the gap plus a 50 kb flank, phased at `-q 1` to match pass 1.

## The anchoring is unambiguous

| gap | gap bp | BAM phase sets | anchors left | anchors right | agreement | residual break | reads gained | accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 10,296,487 | 47,108 | 2 | 302 | 288 | 1.000 | **627** | 66 | 98.5% |
| 7,163,303 | 43,347 | 3 | 263 | 107 | 1.000 | **1,101** | 82 | 97.6% |
| 36,332,599 | 48,420 | 2 | 112 | 228 | 1.000 | **1,111** | 164 | 98.8% |
| 16,673,759 | 41,482 | 2 | 251 | 134 | 1.000 | **2,009** | 175 | 97.7% |
| 60,084,884 | 47,398 | 2 | 131 | 247 | 1.000 | 8,599 | 66 | 98.5% |
| 7,073,919 | 44,697 | 4 | 97 | 262 | 1.000 | 25,203 | 9 | 100.0% |
| 13,429,829 | 201,996 | 6 | 129 | 263 | 1.000 | 174,961 | 80 | 97.5% |
| 35,919,404 | 236,915 | 5 | 255 | 143 | 1.000 | 215,498 | 77 | 93.5% |

Every gap anchors on **both** sides at agreement 1.000, with 97-302 anchor reads
per side. Orientation is never in doubt. Across the eight windows, **719 reads
the graph left unphased receive a haplotype, 702 of them correct against the
diplinator truth (97.64%)**.

## But no BAM phase set spans a gap -- the BAM breaks in the same place

Not one of the eight gaps has a single BAM phase set reaching both graph blocks.
The BAM phase sets run *into* the gap from both sides -- 18 to 33 kb in the
mid-size cases -- and then break. Two independent phasing channels, one driven by
the snarl catalog and one by the pileup, break at the same position, which is
strong evidence the break is a property of the data rather than of the graph
channel. It matches the two other measurements pointing the same way: the
linkage-bottleneck scan, and pileup column discovery admitting only 0-4 columns
inside gaps against 30-103 inside blocks.

So transfer cannot stitch these gaps on its own. What it does instead is shrink
them:

| class | gaps | gap span | residual after transfer |
|---|---:|---:|---:|
| 10-50 kb | 6 | 272.5 kb | 38.7 kb (14.2%) |
| 200 kb+ | 2 | 438.9 kb | 390.5 kb (89.0%) |

Four of the six mid-size gaps come down to **0.6-2.0 kb of residual break** from
41-48 kb originally, a 95-99% reduction. The two large gaps barely move, because
their interiors are phased into BAM phase sets that anchor to *neither* graph
block -- phasing islands that would need chaining, which is the same
both-flanks-link-to-different-sets problem seen in hybrid gap recovery.

## What this changes about the plan

Transfer should run first, before any site discovery. It is pure addition -- the
graph blocks keep their own gauge and only previously untagged reads receive
tags -- so the concordant-to-discordant gate cannot be violated by construction,
and it needs no new variant calling at all. Its cost is that transferred reads
carry 2.4% error against the graph's own 0.85%, so they are lower-quality
coverage, not free coverage.

Then the remaining problem is a **1-2 kb residual break**, not a 48 kb gap. That
is where new evidence has to come from, and at that size the expensive options
become cheap: abPOA consensus over the reads crossing 1 kb, or the GAF `cs:Z:`
private-variant channel keyed on node and offset, are entirely tractable on a
1 kb window and were not on a 48 kb one.

Scope note: eight windows, chosen as the gaps with the most confidently mapped
unphased reads. The per-gap recovery rate (9-175 reads) should not be extrapolated
to all 311 gaps without running them.
