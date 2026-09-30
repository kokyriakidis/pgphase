# Additional chr20 gaps joined by HiPhase on the pgphase callset

`audit.py` compares the accepted pgphase graph and recovery VCF
(`/tmp/pgphase-repeat-pair-guarded-full/phased.vcf`) with HiPhase run on the
same VCF calls (`/tmp/hiphase-on-pgphase-final-chr20/`). It considers adjacent,
nonoverlapping pgphase phase-set extents outside chr20:26–29.5 Mb. The boundary
alleles must have identical VCF POS, REF, and ALT in both outputs and HiPhase
must put them in one phase set. `same_calls` compares all input variant keys
within 50 kb of each boundary. HiPhase read purity is scored against
`test_data/derived/chr20_truth_hap.tsv` after allowing either HP orientation.
`panel_before.tsv` preserves the panel state used to identify new cases.
Run `python3 evaluations/2026-09-29-expanded-hiphase-gaps/audit.py` to
reproduce the baseline `audit.tsv` while the referenced chromosome outputs
remain present. `audit_after.tsv` uses the same command with
`--pg-vcf /tmp/pgphase-newgap-both-full/phased.vcf`,
`--panel evaluations/2026-09-16-test-panel/panel.tsv`, and
`--out evaluations/2026-09-29-expanded-hiphase-gaps/audit_after.tsv`.

There are **56 additional exact-boundary HiPhase joins**, 51 with identical
local calls in the HiPhase input. The two strict additions have at least 20
local HiPhase truth reads, at least 98% local purity, at least 95% purity on
both 10-kb flanks, and the same parental orientation on both flanks:

| Boundary gap | HiPhase local purity | HiPhase separated reads | pgphase result |
|---|---:|---:|---|
| 13,784,702–13,800,824 | 218/218 | 143/155 | Joined; 135/155 separated in the owning chunk |
| 64,140,314–64,144,256 | 103/103 | 62/62 | Joined exact VCF boundary alleles; 45/62 separated |

At 13.8 Mb, the two catalog deletion candidates retain their padded VCF
anchor in the graph chunk. The old physical stitch compared this anchor with
the normalized recovery seam and missed both rows. The repaired stitch
normalizes each selected graph REF/ALT, requires three independent clean
primary deletion pairs with both allele classes, and checks the local graph
SNP links. The left source has an unsupported earlier SNP edge, so only its
certified suffix moves into the right block. The earlier SNP stays split.
Owning-chunk read truth remains 2,775/2,954 correct with no local switch;
the dominant gap block improves from 72/155 to 135/155 correctly separated.

At 64.14 Mb, the BAM solve supplies one left deletion and two overlapping,
complementary right deletion rows. The right rows' REF calls cannot identify
the allele, and the standalone BAM source phase set has mixed parental read
labels in this region. Primary reads spanning both boundaries instead call
the left deletion and use their established right-block HP with a nearby
observed right site. Thirteen callable right-block reads unanimously orient
the left row, with both left allele classes represented. Only that row moves;
the earlier graph SNP and all existing read assignments remain unchanged.
The exact VCF boundary rows now share PS 64,905,482 and have opposite ALT
haplotypes. Eight local truth reads remain unphased, so the read-separation
score does not improve yet.

| Full chr20 measure | Accepted baseline | Both additions |
|---|---:|---:|
| Variant keys | 62,154 | 62,154 |
| Truth-scored phased reads | 236,859 | 236,859 |
| Correct / discordant | 229,108 / 7,751 | 229,108 / 7,751 |
| Read phase sets | 699 | 698 |
| Previously tracked joins | 23/23 | 23/23 |

Only three VCF sample fields change: 13,773,452 G>C and 13,784,702 AT>A
join PS 13,844,727, and 64,140,314 AT>A joins PS 64,905,482. The first
two change haplotype orientation together; the last changes only PS. No truth
read gains or loses an HP tag. Re-auditing the new full VCF against HiPhase
still finds 56 exact HiPhase-joined pgphase splits. Each accepted join
exposes an earlier adjacency. At 13,752,640–13,773,452, HiPhase's left flank
is only 19/36 truth-correct while its right flank is 39/39; their parental
orientations disagree. Keeping that upstream graph prefix separate prevents
a switch. At the other new seam, 64,138,752–64,140,314, HiPhase places
53/53 truth-overlapping reads correctly across the 1,562-bp gap.
It is now an explicit open panel case and owning-chunk assertion. The two
repeat deletions give 14 versus 6 sequence-equivalent paired allele votes;
the clean-SNP bridge gives 3 versus 1, and right-block HP versus left
physical deletion gives 13 versus 8. Both blocks have a consistent parental
orientation in the truth audit, but none of these direct observations alone
supports a safe algorithmic join.

Five more same-callset, noncentromeric HiPhase joins with at least 95% local
truth purity and at least 95% on both concordant flanks are now explicit open
gap tests. The panel records their HiPhase correctly separated read count and
the exact current pgphase nonspan. They are 5,256,785–5,263,741 (145/148
local truth, 65/77 separated), 12,256,072–12,269,535 (217/227, 134/151),
21,823,066–21,844,359 (174/179, 130/160), 35,498,368–35,516,845
(169/176, 93/114), and 50,548,245–50,562,066 (172/177, 101/119).
The 12.256 Mb test replays its full 12–13 Mb owning chunk: a 100-kb replay
joined it, while the full chromosome did not.

Eight further same-callset gaps have at least 80% local HiPhase truth purity
and at least 95% purity on both concordant flanks. All remain open and are
now pinned in the panel:

| Gap (Mb) | HiPhase local truth | HiPhase correctly separated reads |
|---|---:|---:|
| 3.964–3.981 | 180/201 | 104/134 |
| 15.352–15.368 | 200/219 | 118/151 |
| 21.130–21.148 | 210/231 | 128/159 |
| 21.159–21.172 | 202/215 | 108/123 |
| 23.792–23.807 | 161/177 | 86/109 |
| 35.329–35.348 | 163/177 | 98/149 |
| 57.854–57.867 | 158/168 | 94/118 |
| 58.366–58.386 | 156/165 | 104/127 |

The next blockers are concrete. At 5.25 Mb, 35 MAPQ-30 reads physically span
the insertions, but the MSA candidates contain only 9 and 24 allele calls.
A repeat-length trial that treated only exact inserted lengths as ALT left a
6-versus-14 boundary vote conflict. It also reopened the previously verified
21.378 Mb join and raised chromosome discordance by seven reads, so it was
reverted. At 12.256 Mb, 16 primary reads span the SNP and deletion, but only
one directly calls the deletion ALT; a single Q40 SNP pair reaches the next
right SNP 21 kb away. At 21.159 Mb, 16 Q30 primary pairs span the left SNP and
the next right SNP; all 16 call the right SNP REF, divided 9/7 between the
two left SNP alleles. MAPQ 5 admits no additional paired allele call. At
50.548 Mb, complementary deletion lengths have mixed links to the right
SNP. Those observations do not justify a whole-block join.

The two new owning-chunk regressions assert the exact alleles, PS relation,
parental orientation, and that unsupported earlier graph prefixes stay split.
The expanded panel pins 14 newly identified open gaps, including the
exposed adjacent seam, and their per-window read-concordance floors. The required-site check now preserves
all candidates at a shared position; an unphased duplicate can no longer
hide a phased allele when the generated site manifest is expanded.
`make -j4`, `make unit-tests`, `make check`, and the full
`make window-tests` suite pass (3,302 assertions in 42 cases).
