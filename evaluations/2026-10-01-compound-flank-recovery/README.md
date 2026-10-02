# Recover the missing compound flank at chr20:528,850–542,052

## Root causes and retained corrections

The source boundary inside the graph seam is the MSA deletion at internal
542,053 and clean SNP at 545,002. Its 58 paired MAPQ30 calls contain 42 cross
versus 16 same. A strong majority concealed the excess contradictory indel
calls, so the original recovery did not admit an unplaced-read MSA retry.

A general retry on every such majority regresses unrelated windows. The
retained admission requires a compound flank: a phased, MSA-verified deletion
covers a graph flank anchor, and the BAM source has no SNP row at that base.
Here the nine-base deletion beginning at 528,826 covers the child graph SNP
at 528,827. Existing nondecisive retry admission is unchanged. This is an
admission condition, never evidence for choosing a block orientation.

The successful retry must preserve every previously phased source key and a
single consistent clean-SNP gauge within each original source block. Its
selected, cut-free blocks retain all private flank rows in transfer. Dropping
those rows previously lost the deletion and SNP context just outside the
seam, preventing correct child-SNP projection and left-block attachment.
Shared graph rows keep their original allele representation.

One mismatching MSA observation previously vetoed child-SNP projection. The
physical REF/deletion versus MSA allele table now needs a matching majority
in each allele class and Fisher two-sided p <= 0.01, the existing association
threshold. Physical SNP ALT still vetoes projection; physical/graph allele
association, source haplotype association, and the supported next-site path
remain required. The SNP ALT and complex insertion ALT stay opposite.

A partial read beginning at 528,890 also inherited MSA observations at sites
before its alignment starts. Unplaced-read recovery now removes such calls
and recounts allele and total depths in O(observations + sites). SNP spans and
both surviving indel flanks must be physically covered. Ordinary BAM MSA is
unchanged. Synthetic tests cover boundaries, partial deletions, nonzero profile
starts, unknown calls, query coordinates and depth consistency.

A NoisyCandHom row can retain differing internal consensus alleles. It cannot
supply a direct singleton rescue marker: it needs independent primary-read
association through the inferred-site path. MSA genotypes from the new retry
also require existing singleton confidence or a second independent locus.
The observations remain available for block stitching.

## Targeted HiPhase comparison

HiPhase v1.6.0, source commit 74dbc7eaa46d2f0678cfafe07a3a22b9a4a8c121 (local diagnostic clone), ran with its
default realignment on complete primary alignments overlapping 450–610 kb
and the accepted pgphase VCF sites in that interval, with unphased genotypes.
The diagnostic modification only logs read allele vectors.

For the 11 MAPQ60 molecules crossing both deletion boundaries, original BAM
recovery calls left DEL2 REF on all 11, with six right calls missing. The
physical left CIGARs are compound 34–39-base insertions or 7–10-base deletions;
none is a literal isolated two-base deletion. HiPhase calls six left ALT and
five REF using compound allele context. Neither the read nor its coverage was
lost. Graph matrix READ extents are observation extents, not physical BAM
alignment bounds; they must not be used to diagnose BAM coverage loss.

The retained MSA matrix now calls six left ALT and five REF, and recovers all
six missing right calls. Both deletion calls match HiPhase on all 11 physical
crossing molecules; the committed molecule table preserves those calls.

The retained owning-chunk solve connects the six source/graph flank rows
protected by the new regression. Local gap-overlap scoring uses parental truth
only after phasing. HiPhase phases 129 reads, 127 correct / 2 discordant
(98.4496%). The retained pgphase version keeps its 123 phased gap reads,
all 123 correct, and now places them in one phase set. It does not claim to
match HiPhase's six additional assignments. Disjoint flanks have matching
parental orientation, with no large haplotype switch.

## Rejected experiments

- Broad majority retry: owning chunk 0 changes 3,711 phased / 3,698 correct /
  13 errors to 3,718 / 3,682 / 36 and loses three VCF keys. Rejected.
- Complete private-flank transfer for every ordinary source: restores the
  target's left flank but admits unrelated rows and loses an established
  15 Mb VCF key. Restrict complete transfer to the newly admitted compound
  source and established focused retries.
- Unrestricted conflicting-majority retry with coverage/marker corrections:
  chr20 adds 113 phased reads, 82 correct and 31 discordant; concordance falls
  to 96.9606%. The existing 6.578 Mb and 12.269 Mb accuracy gates fail. Rejected;
  no existing floor is weakened.
- Filtering unplaced MSA source labels before stitching or withholding their
  overlay labels breaks the source/graph calibration and regresses the left
  flank. Removed. The final approach preserves source evidence and controls
  retry admission and singleton confidence.

The newly added closure test originally carried provisional counts from a
rejected trial. Its final measured baseline is 123/123, with a strengthened
100% local accuracy requirement and zero-error ceiling. All pre-existing
accuracy floors remain unchanged. The panel's previously open target becomes
a measured positive span and uses its owning chunk to retain the full flank
context.

## Final verification

| Metric | Accepted baseline | Retained change |
|---|---:|---:|
| Phased/scored reads | 237,101 | 237,098 |
| Truth-correct reads | 229,922 | 229,923 |
| Discordant reads | 7,179 | 7,175 |
| Concordance | 96.972176% | 96.973825% |
| Read phase sets | 682 | 680 |
| VCF rows | 62,811 | 62,850 |
| VCF blocks | 338 | 337 |
| Block span N50 | 739,888 bp | 739,888 bp |
| Tracked spans | 76/93 | 77/93 |

No old VCF key is lost. The target is the only newly spanned panel interval;
twelve competitor targets and four split controls remain open. Total phased
coverage decreases by three reads; there is no coverage-gain claim.

Default eight-thread chr20 completes in 179.69 s. The 105 fresh native panel
commands complete in 234.81 s with four workers. The separate scheduler matrix
regression also runs natively with the same binary SHA, making 106 native
commands in total. SHA/command-verified scoring passes 4,157 assertions in 55
cases. The new closure test fails three assertions on the accepted pre-change
binary and passes all 43 afterward. The source-state regression keeps the
later 57.854 Mb trial and rejects its unsupported join.

All standalone units pass; phase predicates pass 697 assertions/32 cases;
port parity passes 27/9; original upstream C parity passes 168,696/7. The build
passes, with the existing unused SIMD helper warning in vendored abPOA. No
established accuracy floor is lowered. Production uses no truth, competitor,
coordinate or read-name exceptions and adds no alignment stage.

Production SHA256: `0dee719443489efb56f5f33efc69c0c825f93904477af1c25d55cde63206f733`.

`metrics.json`, `verification.json`, `crossing-molecule-calls.tsv`, and
`target-variant-rows.tsv` and `tracked-spans.tsv` preserve the measurements,
allele audit and all 93 before/after span decisions.
