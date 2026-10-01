# Sparse chr20 gap links and why a lower read threshold is unsafe

This truth-only audit compares the current full-chr20 pgphase graph+BAM output
(`/tmp/pgphase-8166-final-full/`) with HiPhase run on pgphase's variant
callset (`/tmp/hiphase-on-pgphase-final-chr20/`). The input to every phasing
run is the same HiFi BAM; parental labels are used only after phasing.

## Physical coverage is not the missing channel

The 15 still-open, noncentromeric panel gaps each have at least one primary
BAM alignment physically spanning the exact endpoints (1-18 in this screen).
Grouping all MAPQ>=5 primary and supplementary alignments by read name found
**zero** additional split-alignment molecules that bridge either endpoint
pair. HiPhase enables supplementary joins by default, but they do not explain
these 15 joins. At 19,373,922-19,395,544 exactly one primary molecule spans
both boundary positions. Its BAM source calls the left complementary insertion
rows as REF and the right clean SNP as ALT. A second molecule calls the
intermediate 19,377,346 insertion and the right SNP, but starts downstream of
the left insertion. The targeted BAM source matrix contains only two
one-sided intermediate-to-SNP pairs and a weak cut.

Both the standalone BAM VCF and the targeted BAM source label the left
insertion and right SNP with one BAM PS. In the targeted source their ALTs
lie on opposite haplotypes. HiPhase's same-callset block places those ALTs
on the same side and scores
195/203 truth-labelled local reads correctly (left flank 39/43, right flank
45/46). Thus the inherited BAM PS itself has the wrong relation across the
weak cut; transferring it whole would create a large switch.

At 10,727,690-10,746,628, a fresh 10-11 Mb BAM solve still separates the
left clean SNP from the complementary right deletion/insertion block. Only
two exact boundary-pair molecules connect them. At 3,963,987-3,981,464,
the separate six-base graph deletion has seven callable pairs to the right
SNP, split 4:3 by parity; its intermediate graph indels also have mixed
links. Treating the graph deletion as the BAM one-base deletion is both a
representation error and unsupported by these observations.

## HiPhase read-count experiment

The local HiPhase source defaults both `--min-spanning-reads` and
`--min-connecting-reads` to one. Re-running the same chromosome and
callset with `--min-spanning-reads 2` leaves the correct 19.373 Mb join
in place. It splits wrong-parental-orientation joins at 11.456 and 57.085
Mb but leaves wrong joins at 50.099, 60.098, and 60.706 Mb. A single count
threshold therefore does not separate correct from incorrect joins.

The current same-callset audit also finds six unpanelled HiPhase joins longer
than 1 kb where truth-scored reads on both flanks are at least 95% pure yet
the flanks have opposite parental HP orientations inside one HiPhase PS:
11,456,273-11,477,908; 17,839,399-17,852,024;
50,099,530-50,115,566; 57,085,410-57,104,654;
60,098,495-60,114,911; and 60,706,792-60,711,400. The 60.706 Mb
case remains joined by HiPhase even in the two-spanning-read run. These
are concrete anti-join controls against accepting a raw majority or
copying HiPhase's default one-read block rule.

A 5-kb-bin audit of truth-labelled pgphase read tags found no adjacent
high-purity bins with opposite parental orientation inside any current
noncentromeric read PS (minimum five reads and 90% purity in each bin).
This is a screen, not a proof that every individual read tag is correct.

No production stitch threshold was changed, and no new gap was closed by
this audit. The next implementation needs a site-representation or
allele-assignment improvement that creates independent, consistent evidence
across one of these weak cuts.

## Rejected exact-ALT homopolymer backfill trial

At 15,351,845-15,367,755, two MAPQ60 primary reads physically span both
flanks and carry the exact one-base deletion in their CIGAR. The selected
MSA source initially had no allele observation for either at that
homopolymer deletion because `backfill_msa_observations` excludes all
homopolymer indels. One read has deletion-flank base qualities 40/10; the
other has 40/40. Both call the next right clean SNP REF, whereas five of
six other crossing reads call that SNP ALT. The ordinary physical stitch
already checks equivalent deletion events at base quality 30, so this
matrix dropout is not the only blocker: at most one of the two ALT
molecules passes that independent physical quality gate.

I trialled admitting quality-checked exact CIGAR ALT calls, while still
excluding homopolymer REF backfill. Both missing ALT observations appeared
in the source and transferred graph matrices. The 15-16 Mb owning-chunk
VCF stayed byte-identical and the target gap remained split. Read tags
gained four truth-scored reads and lost one, changing local correct /
discordant counts from 4,257 / 22 to 4,258 / 24.

One full chr20 trial of that general rule also closed **none of the 15**
tracked open gaps. It changed variant keys 62,361 -> 62,365, changed 273
shared-key genotypes, and raised read phase sets 691 -> 718. Truth-scored
tagged reads rose 236,880 -> 237,492, but truth-discordant reads rose
7,737 -> 7,970 (correct reads 229,143 -> 229,522). The trial was
reverted. Exact ALT backfill may be useful in a narrower evidence-only
calculation, but admitting it across all recovery phasing matrices has
a measurable accuracy cost without a gap-closure benefit.

## Wider clean-SNP flank screen

A fresh owning 21–22 Mb graph+BAM replay with the current binary confirmed
that 21,594,343–21,612,458 remains split. Only one MAPQ>=30 primary molecule
physically spans both endpoints. Its CIGAR has the exact four-base left
insertion, but two shifted deletions of three and eight bases near the right
repeat deletion. The one paired matrix observation is therefore a real
coverage limit, not a transfer dropout. A one-molecule whole-block stitch
would be especially vulnerable to right-allele representation ambiguity.

For accurate open panel gaps, an additional screen compared the four nearest
phased clean SNPs on each side (within 30 kb) against the original BAM at
MAPQ>=30 and base quality>=30. Most candidate pairs had zero molecules
calling both SNP alleles. The 21,159,070–21,172,487 gap had 16 paired calls
between its left SNP and the right block's first clean SNP: nine support the
current same gauge and seven the cross gauge. All 16 call the same right SNP
allele, so the apparent count does not provide diploid support for either
orientation. These data do not justify relaxing the stitch rule.

The diagnostic script is `/tmp/screen_clean_snp_flank_pairs.py`; the replay
artifacts are `/tmp/pgphase-21594-current2/`. No production code or gap
expectation changed in this follow-up.
