# Preserve insertion source certificates and selected stitch parity

## Corrections

The SNP-to-MSA-insertion stitch still looked up its complete-source certificate
using the insertion's current phase-set label. Transfer can change that label:
it neither identifies the original BAM source nor certifies its allele gauge.
The long-insertion source-path fallback now calls the same original-source
guard used by the corrected deletion bridge. This requires unique adoptable
provenance, a complete original path without weak/quality cuts, and another
distinct-coordinate source anchor in the current block with a consistent gauge.
The ordinary graph-path alternative keeps its existing checks.

The insertion stitch also selects a supported nearby clean-SNP relation when
available, ahead of the physical repeat insertion vote. Its left-suffix branch
previously discarded that selection and merged using the insertion log odds.
Both left-path branches now use the selected parity. A different path-validation
route does not justify reversing an already selected allele relation.

These are truth-agnostic corrections. No allele representation, site admission,
mapping/base-quality threshold, realignment, or region-specific production rule
changes. Four additional adapter checks exercise relabeled long insertions,
an incomplete original source masked by a complete destination, original quality
cuts, and inconsistent source/current gauges. The existing 32.234 Mb owning
regression passes all 27 assertions, including the established upstream join
and parental read orientation. These fixtures do not demonstrate a new native
chr20 failure of the discarded-parity branch; that correction follows from its
inconsistent use of the selected relation in the two return paths.

## Matched default chr20 comparison

Accepted baseline: `/tmp/pgphase-gap-next10/source-full`.
Corrected run: `/tmp/pgphase-gap-next11/insertion`, eight threads, 160.11 s.
Binary SHA256: `018e2db0f041e6ad88f0b13a323714d6e356c9e943c95973a00f99040cc7f075`.
Inputs are the same annotated BAM, graph catalog, GAF and reference, with
default recovery settings. Competitors were not rerun.

| Measurement | Before | After |
|---|---:|---:|
| Primary read names | 256,586 | 256,586 |
| Phased/truth-scored reads | 237,101 | 237,101 |
| Correct / discordant | 229,921 / 7,180 | 229,921 / 7,180 |
| Read concordance | 96.971755% | 96.971755% |
| Read phase sets | 684 | 684 |
| VCF records / blocks | 62,798 / 339 | 62,798 / 339 |
| VCF span N50 | 739,888 bp | 739,888 bp |
| Spanned panel windows | 75 / 93 | 75 / 93 |

Every primary read name and HP/PS assignment, and every VCF key, genotype and
phase-set label, matches exactly. No tracked gap closes in this experiment:
fourteen competitor targets remain open. No window expectation, required-site
entry, or panel membership changes.

The unit tests and the 491 predicate, 27 port-parity and 168,696 upstream-parity
assertions pass without new compiler warnings. Behavior documentation is
updated in `docs/IMPLEMENTATION.md`.

All 105 distinct native panel commands run fresh on the final binary with four
workers in 194.81 s. Binary-hash and exact-command verified scoring passes
**4,082 assertions / 53 cases**, covering all 93 panel windows and owning-chunk
orientation controls. The final standalone unit rerun also passes. No existing
expectation or required-site entry was weakened or regenerated.
