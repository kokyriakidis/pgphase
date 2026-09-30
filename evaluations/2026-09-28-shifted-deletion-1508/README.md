# chr20 1.508 Mb: equivalent deletion placements plus a weak left path

The tracked chr20:1,508,171–1,523,721 gap remains split. Its left BAM site is
an MSA-verified four-base deletion at 1,508,172 (`TTAA`); the right boundary
is a clean graph `C>T` SNP. Eleven MAPQ-60 primary reads span both sites.
Two carry the exact deletion at 1,508,172 and SNP ALT, three carry a shifted
four-base deletion at 1,508,174 (`AATT`) and SNP ALT, and six carry neither
deletion and SNP REF. Deleting either `TTAA` or `AATT` from the local reference
segment `AGAAGGGTTAATTCTCTCTGA` yields `AGAAGGGTTCTCTCTGA`. The BAM source
retains its original site row at 1,508,172, but a physical exact-position
CIGAR call does not recognize the three equivalent shifted observations.

That representation fix alone is insufficient for a safe whole-block join.
The left graph phase set fails the existing continuous SNP-path check at
1,349,327–1,350,788: graph read votes are 19, 0 and 6 for the two agreeing
haplotypes and reversal, respectively. The right path passes. A boundary
vote must therefore attach only a validated suffix of the left block, after
splitting at this weak cut, or it could extend an internally unsupported
orientation into the right block. The temporary diagnostics used to locate
this cut were removed; no production join was changed for this gap.

Reproduce the source and path inspection with `collect-graph-variation` over
`CHM13#0#chr20:1000001-2000000`, using the chr20 reference, surjected BAM,
striped site VCF and coordinate GAF in `test_data/`, plus
`--phase-matrix-dump`. The eleven spanning read CIGAR patterns can be checked from the primary BAM
alignments crossing the two boundary positions.

## Resolution and chromosome check

The earlier split diagnosis was conservative. The weak 1,349,327–1,350,788
edge has support from only one haplotype, but a direct graph-observation edge
from 1,349,327 to the next SNP at 1,351,500 supports both haplotypes. The
original BAM bases at that bypass pair give 20 `T/G` and 18 `C/A` calls, with
two `T/A` conflicts. The graph path validator now permits one such weak site
to be bypassed only for a deletion bridge that already has decisive physical
allele evidence. The direct graph edge must pass the existing one-sided
binomial `p <= 0.01` test with at least two reads from each haplotype. It
still rejects reversal-dominant cuts and all further weak sites in that path.

The bridge checks CIGAR deletions within 16 bases by comparing the sequence
left after deleting the observed and catalog reference intervals. It does not
alter the BAM site row. Nearby other indels abstain; both reference flanks and
remaining local bases require quality at least 30. SNP pairs take priority.
The 1.508 Mb seam has three accepted ALT and four REF deletion/SNP pairs,
all supporting one phase relation. An initial trial admitted a 3.964 Mb join
from only one read per allele and raised chr20 discordance. Requiring two
independent reads for each indel allele removed that join while retaining the
1.508 Mb closure.

The owning 1–2 Mb chunk now joins `GTTAA>G` and `C>T` in the same phase set
with their ALT alleles on the same haplotype. Its truth score remains 4,034
correct and 11 discordant of 4,045 scored reads. The short panel window
joins with 109/150 truth-scorable reads correctly separated in one block and
zero discordant among 590 scored reads. The target is added to the owning
chunk and panel regressions, including its normalized deletion position
1,508,172 in the required-site list.

Final full chr20: tracked HiPhase-correct closures rise 13/23 to 14/23, read
phase sets fall 713 to 711, and the graph panel span count rises 41 to 42.
All 62,154 VCF keys and genotypes remain identical to the preceding build.
The same 236,847 truth-scored read names retain their individual correctness
outcomes: 229,090 correct and 7,757 discordant (96.7249%). The broad
skip-one-SNP trial had raised discordance by 535; scoping the bypass to the
MSA deletion bridge and requiring two reads per allele removed that loss.

The final `make window-tests` run passes 1,780 assertions across 34 cases;
`make unit-tests` and `make check` pass.
