# Revisit a BAM deletion that becomes a recovery seam boundary

At chr20:35,498,368–35,516,845, the original graph recovery target starts
at the earlier graph SNP 35,490,918. BAM recovery then places a verified
`TA>T` deletion at 35,498,368 in PS 35,479,836. The adjacent right clean SNP
is PS 35,516,845. `collect_phase_set_seams` sees this new boundary, but the
physical stitch previously retried newly exposed seams only for a graph MNP
or a long BAM insertion on the right. Consequently its deletion-to-SNP bridge
never evaluated the actual boundary. The stitch now also retries a new seam
inside the original recovery target when its left boundary is an MSA-verified
BAM deletion. Existing read-quality, likelihood, and block-path checks still
decide whether to join.

The owning-window replay shows why this particular boundary still abstains:
three primary MAPQ-60 reads physically span both sites, but only one has a
callable deletion REF under the existing exact-reference and base-quality
check. One of the others has a nearby insertion and the third has base quality
10 on the deletion's left flank. None calls deletion ALT. The single callable
pair has log odds 6.810, below the required 6.907, and the left graph SNP
path is not certified. The BAM source path itself is complete (no weak or
quality cuts), which is a separate fact from a supported graph SNP path.
The regression for this exact window remains split and passes 20 assertions.

The final full chr20 run at `/tmp/pgphase-left-deletion-final-full/` has a
byte-identical phased VCF and phased BAM compared with
`/tmp/pgphase-verified-msa-output-full/` (62,361 VCF keys). The candidate
TSV differs only in `INIT_CAT` on 2,753 rows, because the already accepted
output-category refinement now preserves BAM source metadata. An earlier
route-only full trial before that metadata refinement likewise left VCF
and BAM byte-identical. The seam fix exposes evidence checks without
changing this fixture's phase assignments.

To test whether the two-allele requirement was excessive, a separate trial
allowed two independent reads of one indel allele when their quality-weighted
likelihood passed the existing 0.001 wrong-parity bound. It joined the
3.964 Mb and 60.7 Mb blocks in the wrong parental orientation. Full chr20
truth-correct reads fell from 229,131 to 228,740 at the same 236,866
truth-scored reads; discordant reads rose 7,735 to 8,126. That trial was
reverted. It shows that independent mapping/base error estimates alone are
overconfident at these repeat indels; both allele classes remain required.
