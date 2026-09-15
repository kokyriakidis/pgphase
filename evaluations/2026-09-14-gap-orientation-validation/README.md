# Remaining chr20 gaps: evidence and orientation validation

The previous panel's “unphased endpoint” at 23,480,815 was a reporting error.
The native VCF already contained `GT=1|2` and a phase set; the candidate TSV's
`HAP_ALT=3,HAP_REF=0` representation cannot determine whether two alternate
alleles are phased. The trial checker now reads native VCF GT/PS, including
alternate FORMAT order. This corrects the baseline to eight split target
rows, representing seven distinct gaps (23.48 Mb has two audited endpoints).

## Confirmed local fixes

At 19.37 Mb the bridge read has an exact `TTCC` insertion with qualities
40,35,40,22 and a Q22 boundary. Two sequence-distinct insertion alleles whose
lengths differ by at least two bases now require every base and boundary to
pass Q20, plus mean inserted-base quality Q30. Other insertion/deletion gates
retain Q30. Unknown qualities, any base below Q20, insufficient mean quality,
and mononucleotide alternatives remain excluded. This change alone joins
126/280 chromosome-wide gaps, versus 125, with identical parental read
accuracy: 2,734 discordant / 189,041 evaluated, 355 switches, 414 flips.

At 23.48 Mb a singleton bridge has a graph-confirmed left SNP and BAM-Q40
right SNPs. The former rule required the confidence source to match across
both sides. Confidence is now assessed separately per flank, retaining the
opposing-read veto. After stitching, an MSA site in a phase set with no tagged
reads can attach to a read-supported block only with independently sufficient
net support from both haplotypes. This handles the nearby `1|2` insertion
without assigning or relabeling any reads. Both audited endpoints join in the
regional test, with 206/206 evaluated reads concordant.

## Why a local success was insufficient

The first full-chromosome combination produced 3,480 discordant / 189,042
reads. An already accepted homopolymer edge at 23,421,003–23,449,834 wrongly
flipped the short left block. The new, locally correct link then propagated
that error into 745 previously concordant right-block reads. This combination
was rejected. `mixed_unvalidated.read_summary.json` preserves its result.

A broader experiment also revisited sparse homozygous MSA calls using missing
BAM observations during general gap recovery. It produced 135 joins but
13,580 discordant reads, 1,171 switches and 1,235 flips. That general promotion
was removed; its output is preserved as `broad_backfill_rejected.read_summary.json`.
Raw allele balance alone is insufficient justification for a new phasing site.

The narrowed validation requires a homopolymer-driven proposed orientation to
also connect both original flanks with the same parity after restoring BAM
SNP observations on graph-selected reads. Both trials are read-only, and the
accepted parity is applied by the final phase-set edge composition. Sparse
MSA genotype backfill is confined to this confirmation trial and cannot add
new gap edges independently.

## Still split in the regional panel

| Gap | Evidence limitation observed |
| --- | --- |
| 1.08 Mb | The only observed right-flank SNP on the bridging read has Q10. |
| 12.72 Mb | A sparse MSA deletion call misses reference observations, but the recovered deletion allele occurs on both tagged right haplotypes. |
| 17.61 Mb | The candidate bridge's right SNP has Q17 and lacks graph confirmation. |
| 47.67 Mb | A repeat insertion is present in raw BAM observations, but does not have tagged-read anchors at that locus. |
| 56.00 Mb | Several overlapping repeat deletions and multiple alleles have incompatible native representations; no accepted connecting orientation yet. |

`regional_panel.tsv` records six joined target rows and five splits, with no
unphased endpoints. These are regional results, not a chromosome-wide claim.
Competitor variants and parental truth were used only for diagnosis/evaluation;
phasing inputs remain the graph, GAF, reference and surjected BAM.

## Final chromosome-wide validation

Output: `/tmp/pgphase-hp-confirm-only`. The same six target rows join and five
remain split in the full native VCF (`chr20_panel.json`). The two new distinct
connections are 19,358,995–19,395,544 and 23,460,963–23,481,134.

| Parental read metric | Previous accepted build | Final build |
| --- | ---: | ---: |
| Evaluated reads | 189,041 | 189,042 |
| Discordant reads | 2,734 | 2,427 |
| Switch errors | 355 | 351 |
| Flip errors | 414 | 414 |
| Recovered initial gaps | 125/280 | 123/280 |

All 186,307 previously concordant reads remain concordant; 307 previously
discordant reads become concordant; no reads are lost. One additional read is
concordant. Of the 307 corrected read statuses, 32 belonged to original phase
set 22,625,517 and 275 to phase set 37,482,721. These improvements include
separating erroneous joins; they are not all new read assignments.

The net decrease in recovered gaps is intentional and explicit: two new edges
are accepted, while four homopolymer edges fail the BAM confirmation:
23,421,003–23,449,834; 37,642,719–37,667,751;
57,841,772–57,866,713; 62,408,056–62,432,427.
The last two splits do not improve the measured read errors and remain a
contiguity cost of the present confirmation rule. No previously accepted
edge changes parity. `edge_changes.json` records all changes.

Evidence loading took 1.552 seconds and cached solving 191.773 seconds in this
single run. Build, C++ unit tests, Python endpoint-parser test and scoped
`git diff --check` pass. No new variant-level Hamming or shared-callset NGC50
measurement is claimed. The five listed original gaps are still unresolved;
no all-gaps-fixed claim is made.
