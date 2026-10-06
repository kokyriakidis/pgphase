# Mixed repeat alleles close the next qualifying HiPhase gap

The next qualifying open gap is **chr20:57,085,410–57,104,654 (19,244 bp)**,
inside HiPhase's **56,542,357–57,492,578 (950,222 bp)** block, ranked 24th.
The larger remaining open gaps fail HiPhase's 80% correctness rule on all
original primary truth-scorable overlapping reads. `selection.json` and
`screened-seams.json` preserve the complete ranking and screening.

## Cause and repair

The boundary contains complementary MSA- and alignment-verified BAM alleles,
**57,104,654 T>TAAA / TAA>T**: insertion of AAA and deletion of AA at physical
position 57,104,655. Their original counts are DP67, 44 REF / 23 ALT and
45 REF / 22 ALT. The old repeat bridge accepts insertion pairs only. Worse,
reads with one- or two-base A insertions are encoded as REF at both rows;
this misassigns 16 reads that the original sequence can distinguish.

The physical caller now compares both complete insertion/deletion sequences
in 16-base flanks and handles signed repeat lengths. Separate downstream
molecules calibrate the classes against clean SNPs at 57,119,478 G>C and
57,123,860 A>C. Callable substitutions within 20 kb can corroborate the same
molecule; contrary high-quality calls veto it, and each molecule counts once.
Calibration and bridge cohorts are disjoint. The upstream SNP is
57,085,410 A>C. Original primary, nonduplicate MAPQ30 reads require known
qualities, Q10 inserted bases, Q20 flanks and SNP bases, and summed physical
base/mapping error at most 1% for certification.

Independent original-CIGAR reconstruction in `audit_physical.py` measures
calibration matrix **[[10,0],[0,2]]** and one agreeing bridge, with log odds
**−7.564577**, wrong-parity posterior **0.0518%**, and joint error bound
**19.4499%**. The nearest-SNP-only matrix [[10,0],[0,1]] fails the 20% bound;
the second SNP supplies a genuinely additional independently callable molecule.
The association, Wilson confidence, unanimous bridge and 999:1 odds rules
are preserved. Parental truth never enters production.

After the certified union, reads can correct contradictory two-REF labels
through their own physical class; surviving phased observations still veto
conflicts. Five existing read-only rescues enter the core through the excluded
graph deletion at 57,087,240, after an independent diploid original SNP/deletion
certificate. Each rescue retains its haplotype, matches its original physical
allele at MAPQ30/Q20 with summed base/mapping error at most 1%, and has no
conflicting clean SNP. Unsupported rescues retain their original tags.
Candidate alleles, evidence counts and categories remain unchanged.

## Results

| Gap reads, including abstentions | Before | Fixed | HiPhase |
|---|---:|---:|---:|
| Original primary scorable | 155 | 155 | 155 |
| Correct | 127 (81.9%) | **143 (92.3%)** | 139 (89.7%) |
| Unphased | 10 | 10 | 16 |
| Wrong | 18 | 2 | 0 |
| Correct in one connected core | 74 | **140** | 139 |

The correctness comparison includes every original truth-scorable primary
read, including unphased reads. HiPhase's input sequence, CIGAR, qualities,
flags and mapping qualities are checked against the original BAM.
The repaired gap exceeds HiPhase's correct-read and connected-core counts;
two pre-existing wrong assignments remain, while HiPhase abstains more often.

The native 57–58 Mb owner improves 3,860→3,876 correct and 194→178 wrong,
with 354 abstentions unchanged. The stitched 56–58 Mb continuation improves
7,705→7,721 correct and 268→252 wrong. Both preserve every previously correct
assignment and all variant evidence. One-thread and four-thread owner
candidates, VCF records and BAM assignments are identical. Disjoint left/right
parental cohorts verify that the join preserves orientation.

Full chr20 now has **254 blocks**, down from 255. **N50 remains 955,496 bp**;
the largest block remains 3,012,193 bp. Among 256,610 primary reads, correct
assignments increase 230,749→230,765 and wrong assignments decrease
6,674→6,658, with 19,187 abstentions unchanged. All previously phased and
correct assignments, variant evidence and block extents are preserved.
The merged core spans **56,542,357–57,481,589 (939,233 bp)**.
Across the whole HiPhase target, pgphase still trails by **54 correct reads /
120 connected-core correct reads**; the completed main-gap comparison is
separate from that remaining deficit.

## Remaining terminal seam

This fixes the qualifying main gap; it does **not** claim the entire 950 kb
HiPhase span is covered. The separate **57,481,589–57,492,578 (10,989 bp)**
terminal seam remains open. HiPhase scores **82/109 (75.2%)** there, below the
80% rule. The original physical deletion/insertion bridge has five contrary
molecules among 14 high-quality pairs, so it cannot certify that union.

Investigation also found a padded terminal graph SNP 57,481,589 A>G with
68 high-confidence physical REF observations, no ALT, and three deletions.
The current physical validation misses this single retained mixed-indel
branch. A trial retirement improved some reads but lost a previously correct
assignment; that trial is excluded from this change. The terminal representation
and unsupported union need a separate fix with their own preservation checks.

## Verification

The committed panel records measured exact `spans=1`, the 80% contract,
independently measured HiPhase correct/core counts, both boundary alleles and
two downstream SNPs. The new owning-chunk regression checks original allele
counts, disjoint parental orientation and the preceding-chunk continuation.
It passes **645 assertions**; the old binary fails **11** of those assertions.
The focused saved-state check completes in **0.59 seconds**.

Build, all unit tests, 1,551 predicate assertions and HiFi/ONT golden gates
pass with no new warnings. The HiPhase measurement helper independently
measures all **127** committed windows on identical original alignments.
The complete gap suite passes **17,146 assertions in four test cases**, plus
four cache-helper tests, in **38.71 seconds** from saved states. All **232**
previous native labels / **121** independent requests are audited under the
final binary; three labels containing the same owner improve, and all prior
correct/phased assignments, variant evidence and block extents are preserved.
Final binary SHA256:
`b1d8221cffe4f389a1659b759228653bbb7d627211fa4a305f42f4fcab9c54d6`.

Full-chromosome and complete panel audits are recorded in `results.json`,
`full-preservation.json`, `panel-audit.json` and `test-runtime.json`.

Run `bash replay.sh [output-directory]` for the full chromosome and
`make gap-owner-check GAP=57.085` for the focused saved-state check.
Use `LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu` for pipeline/tests and a Python
with pysam for the audit scripts.
