# MSA backfill boundary coordinates (2026-10-01)

## Retained fix

`backfill_msa_observations` compared internal indel positions with recovery
windows expressed as inclusive VCF anchors. An insertion or deletion with
internal position 103 is anchored at 102. A window ending at 102 wrongly
excluded that site; a window starting at 103 wrongly included it.

Membership now uses `VariantKey::sort_pos()` for every candidate. CIGAR allele
calls retain the unchanged internal key. Existing observations, including
explicit ambiguity, are preserved. The change requires no new allocation or
pass, no threshold change, and no truth or competitor input.

The synthetic regression fails five assertions against the old code
(`regression-before.txt`) and passes after the fix. It covers exact insertion
and deletion ALT calls at singleton/right boundaries, REF calls, exclusion
outside either boundary, unchanged SNP coordinates, preserved observations,
quality rejection, and the rebuilt read index. The predicate suite grows from
244 to 294 assertions in 25 cases.

## Matched chromosome result

The baseline includes all preceding uncommitted fixes, including the original
MSA SNP dropout repair at 39.848 Mb. It is not commit `74bc786` alone.
Both binaries use the same annotated chr20 BAM, stripped graph catalog,
coordinate GAF, CHM13 reference and default graph recovery with eight threads.

All 256,570 output read HP/PS assignments and all parsed VCF rows are identical:

| Measure | Before and after |
| --- | ---: |
| Truth-scored phased reads | 237,071 |
| Truth-correct assignments | 229,892 |
| Discordant assignments | 7,179 |
| Read concordance | 96.971793% |
| Read phase sets | 690 |
| VCF keys / blocks | 62,368 / 340 |
| Span N50 | 684,798 bp |

No additional tracked gap closes. Fifteen HiPhase-positive targets remain open;
the permanent panel retains its 92 windows and all prior expectations.
The twelve owning-chunk replays covering those targets also retain their
read scores and VCF genotypes/phase labels.

Before output: `/tmp/pgphase-gap-next3/final-full`; corrected output:
`/tmp/pgphase-gap-next4/anchor-full`. `metrics.json` records the matched measures.
The paired full-run command is the same as
`../2026-10-01-msa-snp-dropout/commands.sh BEFORE_BINARY AFTER_BINARY OUTPUT_DIR`.

## Rejected probes and next blockers

At 54,684,758–54,702,072, all eight primary spanning molecules are represented
in recovery. Seven have callable boundary pairs: six agree and one opposes.
Q30 physical sequence checks retain three ALT and one REF observations;
one of the ALT observations opposes the majority connection. The left graph
block also has a weak internal SNP edge at 54,607,517–54,607,533: GAF supplies
seven agreeing pairs on one haplotype and two opposing pairs. Physical Q30
bases give 31 agreeing and six opposing pairs. The independent BAM source
has an 88-site path with no weak cuts. This is a graph/source certification
and observation issue, not absence of physically spanning molecules.

At 528,850–542,052, the original BAM source already has a weak cut at
528,850 before exact-CIGAR backfill. Backfill does not create that cut.
The raw CIGARs carry a preceding 9–10-base deletion or a 34–39-base insertion;
they do not directly identify the MSA two-base deletion at the left boundary.
Its reference observations cannot stand in for the preceding complex allele.

A temporary probe enabled unplaced-read MSA in every source solve. It joined
both targets, but the 0–1 Mb replay gained 23 discordant assignments and lost
three previously phased heterozygous keys. The 54–55 Mb replay gained 22
correct assignments without adding discordance, but lost five heterozygous
keys. This broad behavior is not retained. A focused probe using majority
absence at a unique MSA deletion beside a clean SNP did not admit a new retry
or change either replay. Neither probe changes production admission criteria.
`probe_summary.tsv` records both comparisons.

## Validation

- Build and shared unit tests pass, with no new compiler warnings.
- Predicates: 294 assertions in 25 cases.
- Port parity: 27 assertions in nine cases.
- Upstream parity: 168,696 assertions in seven cases.
- Permanent gap panel: 3,951 assertions in 49 cases over 92 windows.
- `git diff --check` passes. Existing window expectations are unchanged.
