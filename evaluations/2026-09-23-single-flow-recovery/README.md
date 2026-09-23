# Single-flow recovery MEC

Date: 2026-09-23

## Problem

The guarded 14.264 Mb fix had two calls to the exact solver. Any failure of the
complete atomic-block solve could invoke a narrower boundary solve, even though
only variable-limit exhaustion justified changing scope. A tie or conflicting
full solve must not be retried on a smaller matrix.

## Retained design

Recovery now distinguishes resource exhaustion from evidence failure. It first
solves the complete read-connected atomic blocks. The boundary interval is
considered only when that exact problem exceeds 20 internal variables and the
sequence-identical candidate gauge is significant at one-sided binomial
`p <= 0.05`. Full reads, both deterministic halves, the source 2x2 read gauge,
and representation evidence retain their existing agreement requirements.

This is the practical correctness contract. Binary MEC is NP hard, and the
chr20 audit found coupled components as large as 167 variables. An unbounded
branch-and-bound run consumed roughly seven CPU cores for more than 2.5 minutes
without finishing. Unsupported or intractable edges therefore remain separate.

## Rejected simplifications

A joint block-only optimizer and a whole-block majority validator were tested.
The latter left 14 heterozygotes unphased, raised VCF blocks from 419 to 449, and
increased formal truth discordance from 6,258 to 6,416 reads. It was removed.
An unbounded exact solve closed the known unsafe 51.27 Mb boundary before the
independent representation guard and did not finish full chr20 in practical
time. It was also removed.

## Full chr20 validation

The retained single-flow implementation produces a byte-identical phased VCF
and identical SAM records (`49dda8278a7b5b3d4eb0903b1af96f01`) to the prior
guarded baseline. It retains 61,644 phased heterozygotes, 419 VCF phase blocks,
482,085 bp VCF N50, and 225,005 tagged reads. The formal baseline remains 6,258
discordant of 224,985 evaluated reads (97.2185%). The 14.264 Mb target remains
closed; 51.27 Mb remains separate.
