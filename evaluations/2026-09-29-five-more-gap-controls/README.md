# Five more noncentromeric HiPhase joins to track

The current full chr20 graph/recovery VCF has 55 exact-boundary splits that
HiPhase joins on the same input variant callset. Five were absent from the
window panel. Each now runs in its owning 1-Mb graph/BAM chunk. The HiPhase
number below is truth-correct reads in its dominant gap block over all
truth-scorable reads overlapping that gap.

| chr20 gap | HiPhase separated | pgphase separated | pgphase span |
|---|---:|---:|---:|
| 1,180,618–1,194,189 | 125/129 | 66/129 | no |
| 10,727,690–10,746,628 | 133/161 | 66/161 | no |
| 19,373,922–19,395,544 | 111/140 | 59/140 | no |
| 21,594,343–21,612,458 | 118/154 | 59/154 | no |
| 36,332,599–36,354,890 | 144/179 | 144/179 | no |

The counts use the same truth-scored overlap denominator and dominant-block
orientation as the window tests. The committed expectation file rounds their
fractions down for regression floors.

At 1.18 Mb, five MAPQ-60 reads place a 3-bp ATC insertion 18 bp after its
VCF position in the same ATC repeat. The left boundary SNP is a verified
BAM-injected noisy candidate, so the current physical stitch does not use it
as a clean SNP. Eleven reads assigned to the established left phase set call
both right insertion alleles with unanimous HP parity. A trial admitted that
local evidence, but the right graph phase set has a distant unsupported SNP
edge at 1,349,327–1,350,788: GAF votes are 19, 0, and 6 for the two agreeing
haplotypes and reversals. Across that edge, 16 clean BAM pairs came from one
haplotype and four disagreed. The whole right block could not be certified,
so the trial was reverted.

At 12.256 Mb, a 100-kb replay joins the gap while its owning 12–13 Mb chunk
remains split. Both source solves have a weak cut at the boundary SNP.
Several spanning reads place a one-base deletion 19–22 bp later in the same
A run, on both SNP allele classes. A short-window join is therefore not a
safe certificate for the full source. This case was already in the panel.

A broad trial that let verified injected noisy SNPs enter every physical
seam closed 36,332,599–36,354,890, but full chr20 truth-correct reads fell
from 229,108 to 228,964 among the same 236,859 phased reads; 988 VCF sample
fields changed. Of the 148 reads whose truth status changed, 146 lost correctness
and two gained it. The loss is concentrated in a false 37.462–37.467 Mb join:
the trial reversed and absorbed the block beginning at 37,466,820 into the
block beginning at 37,425,328. The seam also contains competing 2- and
3-base deletion rows; whether they cause the false support remains unresolved.
Accepting a verified noisy SNP as a general whole-block bridge is unsafe. That broad trial was reverted; the narrow read-label fix below is the
accepted change. It leaves all 62,154 variant keys unchanged.
A separate repeat-classifier trial that stopped treating deletions as
insertion REF changed no VCF fields and lost four truth-correct reads, so it
was also reverted.

The new 19,373,922–19,395,544 control exposed a separate read-label bug.
Its VCF blocks were split, but 13 reads wholly right of the gap retained the
left read PS. Twelve carried the opposite parental orientation. The trace
showed that weak-cut detachment moved their heterozygous BAM sites to a new
component while the reads kept the old PS through a homozygous row. The fix
moves source reads that call the detached component and no other oriented
source or graph component with those sites, preserving the source HP.
It checks all observation channels and reserves IDs used by variant and read
phase sets. The owning chunk
keeps 4,071 truth-scored phased reads and improves from 4,020 to 4,031
correct; the switch assertion now passes.

On full chr20, the accepted fix keeps all 62,154 VCF rows and sample fields
unchanged. Truth-scored phased reads change 236,859 → 236,858, correct reads
229,108 → 229,125, and discordant reads 7,751 → 7,733. Read phase sets
change 697 → 698. The five cases add six in-gap phased heterozygotes and no
VCF spans to the panel baseline. `make unit-tests` and `make check` pass. The complete window panel passes
3,496 assertions in 42 test cases.
