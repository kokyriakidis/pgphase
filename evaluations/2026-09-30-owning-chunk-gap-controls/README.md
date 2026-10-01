# Keep gap regressions in the owning chromosome chunk

The full chr20 graph/BAM VCF already places both boundary alleles in one
phase set at 48,929,511–48,950,388 and 56,064,697–56,083,708. Their
window tests had replayed only 100–120 kb around each seam. Those short
replays lacked the established flanking phase-set context and split both
gaps, so the panel marked two production joins as open. The test now replays
the 48–49 Mb and 56–57 Mb owning chunks and asserts exact VCF boundary
alleles, equal phase-set IDs, and equal phased GTs.

| Gap | Owning-chunk span | Local truth-correct separation | Local read concordance |
|---|---:|---:|---:|
| 48,929,511–48,950,388 | yes | 135/156 (86.54%) | 3,844/3,901 (98.54%) |
| 56,064,697–56,083,708 | yes | 148/167 (88.62%) | 3,914/3,988 (98.14%) |

The focused tests each pass 26 assertions. This is a test-fixture correction;
the accepted full chr20 phased VCF and BAM are unchanged.

## Rejected 36.332 Mb trial

At 36,332,599–36,354,890, recovery exposes a right BAM deletion after the
original seam was recorded. The physical stitch's retry list omits new right
deletions, and the new local seam begins at a left insertion at 36,343,992,
leaving the clean SNP at 36,332,599 outside its search range. A trial routed
this deletion seam and expanded its left search to the original target.
It joined the two whole blocks but dropped owning-window read concordance to
3,320/3,769 (88.09%); the prior floor was 95%. The trial was reverted.

The deletion originates in BAM source PS 36,286,778, but a weak-cut run
attachment places it in graph PS 36,549,255, flipping its local allele
orientation. The full standalone BAM solve also puts the left SNP and this
deletion in PS 36,286,778, but that inherited label is not a certificate
across the targeted source weak cut. Only two MAPQ-60 primary reads span the
exact left SNP and deletion; both carry the left SNP ALT, one has a 6-bp
deletion overlapping the 4-bp target, and the other has a distinct 5-bp
deletion 38 bases later. The physical pair test can use at most the latter
as a deletion REF, providing no independent ALT bridge. A whole-block join
therefore risks reversing the independently accurate right block. The
existing nonspan and read-concordance regression correctly rejects it.

## Rejected left-insertion retry

The original 21.130–21.148 Mb recovery target skips a newly exposed left
BAM insertion and the nearer right clean SNP. Retrying that seam reaches the
existing physical insertion-to-SNP checker, but its four spanning primary
reads supply zero exact paired alleles: one has a different insertion
length, one a nearby deletion, and two have no high-quality usable pair.
A broad retry rule closed instead the adjacent 36,611,593–36,620,864
regression-control seam. Full chr20 kept 62,361 variant keys but lost seven
phased reads and six truth-correct reads (236,866/229,131 to
236,859/229,125). The broad rule was reverted. The control's exact
nonspan expectation remains in the panel.

Validation: `make -j4 pgphase test_gap_windows`, `make unit-tests`, and
`make check` pass. Focused owning-chunk tests for 48.929 and 56.064 Mb
pass 26 assertions each; the 36.611 Mb nonspan control passes 20.
The unchanged production source had passed the 45-case full window panel
before these two replay scopes were corrected. The full panel was not
rerun after this test-only change.
