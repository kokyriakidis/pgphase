# Remaining HiPhase-correct chr20 phase-set gaps

The target inventory in `targets.tsv` compares the full pgphase graph+BAM run
`/tmp/pgphase-next-gap-39-full/` with HiPhase 1.6.0 run on pgphase's variant
VCF (`/tmp/hiphase-on-pgphase-final-chr20/`). Both use the same annotated
HG002 BAM. Exact `(POS, REF, ALT)` boundary rows must be present and phased
in one HiPhase PS but different pgphase PSs. Adjacent pgphase PS extents must
be non-overlapping; boundaries in chr20:26.0–29.5 Mb are excluded.

A target needs at least 20 truth-labeled primary HiPhase reads in its PS over
the gap plus 10 kb on each flank, at least 98% local parental purity, at
least 95% purity on each flank, and the same majority parental orientation on
both flanks. This gives **23 targets**, 22 at least 10 kb long and one shorter.
The existing 15-gap larger regression panel is already closed. These 23 are
newly measured boundaries in the broader full-chromosome comparison, including
the 60,033,052–60,048,237 subgap inside an older test window ending at
60,058,235. The comparator VCF was generated before the latest pgphase
changes, so this inventory covers current boundaries whose exact keys exist
in that HiPhase run; it is not an exhaustive claim about unmatched sites.

`window_truth_reads` counts truth-labeled primary input molecules overlapping
the boundary interval. `hiphase_separated_correct` is the largest number one
HiPhase PS places on the correct parent within that interval, after allowing
one global flip of the PS labels. The ratio is the competitor score copied to
the test panel. All 23 targets were added to
`evaluations/2026-09-16-test-panel/panel.tsv`. The window test now keys cases
by both boundaries so the two 60,033,052 starts have distinct outputs and
expectations. On the baseline short-window replay, eight of those 22
already span locally but remain split in the full-chromosome output. Their
production-context transfer or stitch, rather than local allele discovery,
is the remaining issue.

## First repaired target: 64,144,256–64,144,722

The owning 64–65 Mb replay exposed a BAM source PS with a verified weak cut at
64,140,314. The old seam stitch treated that entire source PS as one block and
attached its post-cut deletion to the left graph block. The nearest deletion
loci vote 14–0 for a connection to the right BAM block, but using that vote to
merge the entire left graph block reduced truth accuracy from 3,010/3,013 to
2,586/3,013. The retained fix leaves graph joins intact and restores only
private BAM sites carried past a zero-observation cut to a local source gauge.
The following component-aware transfer joins both boundary deletions to the
right graph block while graph SNP 64,118,182 stays in the left PS. The owning
chunk remains at 3,010/3,013 truth-correct reads with three discordant reads.
The panel now replays this owning chunk and asserts both the boundary join and
the absence of the wrong whole-block join.

## Second repaired target: 514,902–528,827

The full first-megabase graph solve initially connects these SNPs. Recovery of
the next seam moves the right boundary to a BAM source block and exposes a new
gap after the sole targeted solve had finished. The owning-chunk audit showed
three callable BAM sites inside the new gap, all marked outside its recovery
windows. A bounded second pass over newly exposed, nonoverlapping seams now
injects those three rows; exact 514,902 and 528,827 boundary keys share a PS.

The injected 528,728 C>A SNP has an internal noisy-candidate classification
but is MSA- and alignment-verified. High-quality calls correct 20 wrong-parent
HP labels; consistent calls also tag previously unphased reads. Correctly
separated reads rise from 71/119 to 110/119; HiPhase has 113/119. The owning
1-Mb truth discordance falls from 33 to 13 reads, with 3,711 reads scored.

A second-pass stitch initially split the previously sequence-validated
528,827 SNP/528,828 insertion pair. The retained version saves the original
weak-cut BAM source component and moves only its oriented private rows on the
SNP's side of the cut when the SNP gains an upstream PS. The exact pair remains
connected in both the owning chunk and its dedicated regional regression.

The full chr20 run adds only three VCF keys, closes this target as the third of
23, and changes truth-scored read accuracy from 229,058/236,835 (96.7163%)
to 229,090/236,847 (96.7249%).

## Two more repaired targets: 7,047,080 and 48,971,192

A complete BAM source across 48.971 Mb had exact clean graph SNP anchors in
two padded BAM solves. Duplicate-source suppression discarded both anchor
copies, leaving its left flank unattachable only in the owning chunk. Retaining
exact clean shared SNPs as independent anchors for each source, subject to each
source's own path and read-vote checks, closes that gap and the 7.047 Mb gap.
The owning chunks' parental truth counts remain 3,844/3,901 and 3,926/4,091
correct respectively; their read PS counts each drop by one. Both regressions
now run their owning 1-Mb chunks and assert exact boundary VCF PS equality.

The full chr20 output retains all 62,150 variant keys, closes **5 of 23**
tracked targets, and has 229,091/236,848 truth-correct reads, 7,757
discordant (96.7249%). **18 targets remain open.**

## Sixth repaired target: 56,662,188–56,679,959

The 56–57 Mb owning chunk kept the exact `56662188 G>T` and `56679959
C>T` rows in different phase sets. A nearby clean SNP at 56,658,888 was
being used as the left physical anchor, which excluded reads starting between
it and the actual 56,662,188 boundary. The graph candidate at that boundary
is stored as the three-base snarl allele `TGG>TGT` at 56,662,186; VCF
normalization trims its shared context to `G>T` at 56,662,188. Matching only
single-base graph REF/ALT strings missed this exact SNP.

The recovery bridge now normalizes graph alleles before matching to the BAM
candidate and physically calls the boundary base. The BAM MSA candidate has
consensus haplotype allele IDs 1/2 and is not a clean SNP, so its candidate
IDs are not treated as REF/ALT. Other exact clean SNPs establish the left
block's graph/BAM gauge. The first verified right-source insertion at
56,671,110 (`T>TAC` in VCF) provides a physical bridge: 19 MAPQ/BQ >= 30
molecules call both it and the boundary SNP; 18 vote for the winning
orientation and one disagrees. Both source phase-set paths have no weak cuts,
and each graph flank has multiple exact clean SNP gauge anchors. The stitcher
uses its ordinary conflict checks before merging the two graph blocks.

In the owning chunk, truth-scored phased reads change from 3,913/3,987 to
3,914/3,988 correct; discordant reads remain 74. The target's exact VCF
boundary rows now share a PS. The window panel replays this owning chunk and
requires exact PS equality, preserving the closure as a regression.

With the complete 66,210,255 bp chr20 reference, both runs emit the same
62,150 `(CHROM, POS, REF, ALT)` keys. Tracked target closures rise from 5/23
to **6/23**. Truth-scored phased reads rise from 236,848 to 236,849; correct
reads rise from 229,091 to 229,092; discordant reads remain 7,757. Accuracy
is 96.7249% in both runs at four decimals, and read phase-set count falls
from 729 to 728. These runs are `/tmp/pgphase-sharedsnp-final-full/` and
`/tmp/pgphase-snp-insertion-final-full/`.

A first trial accidentally stopped at 66,000,000 bp and appeared to lose 939
scored reads. The missing molecules were almost entirely in the unprocessed
66.000–66.210 Mb tail. The complete-reference rerun above is the valid
comparison.
