# Recovery boundary validation

Date: 2026-10-02. The accepted input state is the uncommitted recovery fix
reported in `../2026-10-01-recovery-anchor-and-single-seam/`, binary SHA256
`598eb8f99b1e14d1ce25eb9a8ef0305bea60b9b44089c2a3f41f9a05c95ac676`.
Inputs and competitors are the identical-BAM chr20 comparison in
`../2026-10-01-phaser-comparison/`. Truth and competitor data are evaluation
only. No production fixture coordinates, read names or truth exceptions.

## Accepted correction: inclusive right anchor

Recovery seams and MSA backfill use closed VCF anchor intervals. The scoped
heterozygote repair used `[beg, end)` and therefore excluded its right boundary.
Use `[beg, end]`, including a singleton seam. Indel membership still uses
`VariantKey::sort_pos()`, preserving the unchanged BAM allele representation.
All allele-depth, transfer-verification and repeat-link requirements remain.

The regression fails on the old implementation at the exact right endpoint.
The expanded predicate test covers SNP, insertion and deletion singleton
boundaries, plus rejection immediately outside either side. This does not
change standalone BAM calling when recovery windows are absent.

Full chr20 output is unchanged: every read's HP/PS pair and every VCF data row
match the accepted input. Phased/scored reads remain 237,118; correct 229,956;
errors 7,162; accuracy 96.979563%; read PS 675; VCF blocks 338; block span N50
756,878 bp. Tracked coordinate spans remain 78/93. This is a latent boundary
correctness fix, not another chr20 gap closure.

## Rejected: remove graph anchors with physical REF-only calls

At 21,172,501 the graph `A>C` boundary has 61 callable physical BAM REF bases
and 17 deletion-overlapped reads, with no physical C calls. Source reads on
both haplotypes carry REF. REF-only loci are absent from BAM candidate rows,
so the existing homozygous-candidate check does not see this evidence.

A trial independently checked physical known MQ30/Q30 calls, requiring both
haplotypes in one source PS and two-sided binomial p <= 0.01. Any callable
non-REF base vetoed it; missing/deleted bases were not counted as REF. The
existing expanded-seam coverage guard stayed in place.

The trial removed the false SNP's phase label and moved the true SNP at
21,183,846 into the left block, with correct relative SNP parity. However,
the deletion at 21,172,487 stayed in a separate block. A coordinate extent
reported the historical gap as spanned even though its deletion was not
connected. **Extent coverage is not proof that both boundary alleles joined.**

| Metric | Accepted input | Rejected REF trial |
| --- | ---: | ---: |
| Full chr20 phased reads | 237,118 | 237,119 |
| Full chr20 correct reads | 229,956 | 229,958 |
| Full chr20 errors | 7,162 | 7,161 |
| Read phase sets | 675 | 677 |
| VCF keys | 63,401 | 63,411 |
| 21 Mb owning-chunk correct | 3,110 | 3,115 |
| 21 Mb owning-chunk errors | 211 | 212 |
| 21.179 Mb dominant correct fraction | at least 0.85 | 0.673913 |

The last row fails an existing committed floor. The trial does not preserve
the neighboring block's read assignments. It is rejected despite a small
chromosome aggregate improvement; no floor is weakened. Local reads over the
true 21,159,070–21,183,846 SNP pair score 151/166 (90.96%) in the trial, versus
149/150 (99.33%) for native HiPhase and 151/164 (92.07%) for same-callset
HiPhase. Removing a false anchor alone does not solve the representation and
block-transfer problem.

The temporary owning-chunk test fails on the accepted input and passes on the
trial, but the complete window suite detects the neighboring regression.
Its logs are preserved here; the trial implementation and changed expectation
are removed. Existing phase-label preservation checks remain.

## Rejected: let source MSA rows bypass transfer verification

BAM source candidates have `msa_verified` before graph transfer supplies
`alignment_verified`. Allowing a source MSA row to bypass the latter guard
can enable depth-based heterozygote seeding earlier. This was tested separately
from the REF-only trial, while retaining all other depth and repeat checks.

It phases 34 extra reads in the 19 Mb owning chunk, but only 17 are correct:
4,072/4,031/41 phased/correct/errors becomes 4,106/4,048/58. The tracked
19,373,922–19,395,544 gap stays split and its local 125/94/31 result does not
improve. The 3 Mb chunk loses two correct assignments. Thus MSA verification
alone is insufficient for this admission. The existing verification guard
stays; this is not presented as a fixed source-provenance bug.

A separate homopolymer-deletion ALT backfill trial made no change at 21 Mb and
was also removed. Physically spanning reads there do not carry the proposed
one-base deletion, so excluding them cannot be repaired by manufacturing REF
or ALT observations.

## Validation

See `validation.json` for the final binary, full output identities and test
results. Native replays are keyed by complete command arguments and final
binary SHA256; cached scoring only copies outputs from those fresh runs.
The retained implementation contains the inclusive endpoint correction only,
in addition to the previously accepted uncommitted fixes. No gap expectation
is relaxed or changed in this round.
