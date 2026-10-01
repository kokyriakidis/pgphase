# MSA observation dropout at two open chr20 gaps

This audit used the same BAM, graph catalog, and GAF as the committed gap
panel. All phasing trials were truth agnostic. The trials below were reverted;
no production phasing rule changed.

At 50,548,245–50,562,066, both MSA-verified BAM deletion rows are present,
but none of the 12 MAPQ-30 source reads spanning the deletion and right SNP
has a callable deletion allele. The retry detector skipped this case because
the deletion has two complementary rows at one position and its crossing-read
floor is 20. A trial admitting this complementary pair to the existing MSA
retry restored one-base deletion calls on 10 spanning reads and two-base
calls on eight. Their deletion-to-SNP votes remain mixed (6:4 for the
one-base row and 5:3 for the two-base row), so the output stays split. The
owning-chunk VCF keeps the same 1,199 keys. Those extra observations do not
certify a stitch.

At 21,159,070–21,172,487, 18 MAPQ-30 BAM source reads span the left clean
SNP and right deletion; all call the SNP and none calls the deletion. The
source carries both sites under one inherited PS, but no read calls both
alleles. Lowering the retry's crossing-read floor from 20 to 12 did not
trigger MSA: the graph recovery window is 21,105,723–21,172,501, and its
large grouped source contains more than two phase-set runs. The detector
stops before examining this internal dropout. Forcing unplaced-read MSA on
the whole 21–22 Mb recovery group restored 18 deletion calls. Seventeen
called REF, one ALT; the right allele remains one sided. More seriously,
the grouped retry removed seven established VCF rows, including the
21,172,487 `CA>C` boundary deletion, and gained two other rows (283 to 278
keys). It cannot replace the current source solve. A standalone BAM solve
restricted to 21,129,070–21,202,487 retains the exact deletion, confirming
that region context affects the MSA representation. The grouped trial’s restored spanning calls were 17 REF versus one ALT,
so those observations did not independently establish a diploid bridge to
the right block.

A diagnostic screen of saved recovery matrices for 13 still-open,
HiPhase-correct panel gaps found no path between boundary alleles whose
every edge had at least eight paired calls, both allele classes represented
at least twice, at most one opposing parity call, and at least 90% majority
parity. This screen does not prove that no weaker or representation-aware
path exists. For example, at 21,823,066–21,844,359, the two intermediate
graph repeat deletions have mixed links to both flanks; forcing them into a
chain would risk a haplotype switch.

The `21,159,070–21,172,487` regression now replays the owning 21–22 Mb
chunk and requires its left SNP, deletion, and right SNP to remain present
and phased. It does not require a split, so a future evidence-backed join
can pass. The focused case passes 31 assertions. Broad MSA and the lower
admission floor were both reverted; `src/collect_pipeline.cpp` has no diff.
