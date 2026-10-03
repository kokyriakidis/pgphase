# Certify graph SNP paths on both sides of an MSA deletion bridge

## Defect and correction

At **chr20:47,751,480–47,762,233**, the MSA deletion-to-right-SNP bridge has
seven callable high-MAPQ molecules, three REF and four ALT. Its signed quality
log odds are 33.8683, above the existing 0.001 wrong-parity gate. Both graph
flanks nevertheless contain internal edges observed on only one haplotype in
GAF: 47,671,540–47,689,418 has votes 0/2/0 and
47,764,156–47,781,656 has votes 3/0/0 (hap1/hap2/reversal).

The insertion bridge already checks such an agreeing left graph edge against
independent clean physical BAM SNP pairs. The deletion bridge omitted that
check, and the right path rejected a one-haplotype edge without checking the
available physical evidence. Two distinct MAPQ-60 molecules corroborate each
edge, with Q35–Q40 bases and zero opposition. Their signed log odds are 15.8691
and 16.2861. These are separate molecules/allele pairs from the boundary vote.

Apply the existing physical edge certificate to deletion boundaries and to
the right graph path. The first weak edge is reported only for missing or
one-sided agreeing graph evidence; a dominant GAF reversal remains a veto.
The prefix before that edge and the suffix after it must pass. Two independent
Q30 primary SNP-pair reads must unanimously agree with the existing block
orientation and meet the existing likelihood bound. The nearest clean suffix
SNP separately certifies the deletion's inherited allele gauge. Both boundary
allele classes and all existing allele/quality/ambiguity checks are preserved.

The implementation adds neither a new search algorithm nor a new alignment,
allele representation, truth input, competitor input, coordinate exception or
relaxed evidence threshold. `physical-edge-evidence.json` records the exact
molecules; `source-trace.log` is exploratory instrumentation from the pre-fix
build, absent from production.

## Results and remaining accuracy issue

| Measurement | Before | After |
|---|---:|---:|
| Local gap blocks | 2 | 1 |
| Local phased / correct / discordant reads | 83 / 70 / 13 | 83 / 70 / 13 |
| Correct local reads in dominant block / input overlaps | 38 / 84 | 70 / 84 |
| Owning 47–48 Mb phased / correct / discordant | 3,993 / 3,957 / 36 | unchanged |
| Owning VCF blocks | 2 | 1 |
| Owning span N50 | 749,000 bp | 997,144 bp |
| chr20 phased / correct / discordant | 237,134 / 229,981 / 7,153 | unchanged |
| chr20 read concordance | 96.983562% | unchanged |
| chr20 read phase sets | 669 | 667 |
| chr20 VCF blocks | 335 | 334 |
| chr20 VCF keys | 63,492 | unchanged |
| chr20 span N50 | 756,878 bp | unchanged |

Every previously phased read retains its individual truth classification in
the chromosome run. There are 1,285 changed HP/PS tuples and 261 changed VCF
rows from label/gauge propagation, with zero lost or gained keys. The correctly
joined deletion and first right SNP carry the same ALT orientation; the next
right SNP is opposite. Exact VCF keys and parental flank agreement are gated.

**This fixes continuity, not local accuracy parity.** Native-DV HiPhase phases
76 overlapping reads with 74 correct and two discordant (97.37%); pgphase still
has 70/83 correct (84.34%). The 13 existing local errors remain, and this change
does not reassign them. Their separate read-tagging investigation must preserve
the newly proven block connection. The competitor files are the saved Oct 1
runs, not new runs against a changed callset.

A fresh noncentromeric extent audit falls from 315 to 314 gaps and from 41 to
**40 distinct competitor-supported nominations**, with 50 records (33 native-DV,
17 saved-pgphase-callset). Nominations are not verified exact allele bridges.

## Other investigations

At 62.408–62.432 Mb, no BAM read physically spans the graph boundary SNPs.
HiPhase uses an intervening deletion; pgphase already retains that deletion as
a repeat indel. Its few physical links are inconsistent, and one spanning
molecule has no callable right SNP in its BAM CIGAR. This turn does not promote
the repeat classification or supply an invented allele.

At 7.28 Mb, the nearest graph SNP has only six eligible graph observations,
while high-quality physical BAM bases at the normalized coordinate have 51
REF and two ALT calls. Its first graph edge is tied (3/0/3). The new check
rejects this edge; it does not reinterpret a tied graph edge as independent
one-haplotype confirmation. That representation/calling discrepancy needs a
separate investigation. The 17.502 Mb deletion bridge remains below its
existing likelihood gate (-6.80953 versus magnitude 6.90675); its threshold
is preserved.

## Regression and validation

Extend the committed panel to **98 coordinate cases / 84 required spans**.
The new owning-chunk test checks the unchanged deletion, exact SNP/deletion
ALT orientation, parental flank agreement, the owning read/error floors, and
local read separation. The accepted pre-fix binary fails six of the 27 new
assertions. Existing expectations are unchanged except the intended total and
new case; the test does not claim the local read accuracy matches HiPhase.

The build has no new warnings; all standalone unit tests pass. Phase predicates
pass 822 assertions / 37 cases, and window tests pass **4,849 assertions / 63
cases**. All 105 fresh native requests are checked against the frozen final
binary and complete in 342.85 s; an additional matrix regression runs natively
with the same binary. A fresh eight-thread chr20 run takes 247.68 s during
validation, not a controlled competitor runtime benchmark. The retained trial
binary and final binary have the same SHA256:

`684a95a62fe10ec32df68c8043bf0fbb730794c167ce54573c1ad2d56c546f9c`

`validation.json`, the comparison JSON files, and gate logs retain measurements.
Reproduce the pipeline gates with `make -j8`, `make unit-tests`,
`./test_phase_predicates` and `make window-tests` after generating the truth map.
The comparison scripts take explicit input/output arguments documented by
`--help` and are evaluation-only.
