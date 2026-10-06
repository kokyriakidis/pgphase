# Equivalent terminal alleles complete the next qualifying HiPhase block

The next qualifying target is HiPhase's **26th-largest block**,
**chr20:54,453,020–55,309,794 (856,775 bp)**. Pgphase stopped at 55,309,789,
five bases short. All larger remaining open seams fail HiPhase's 80% rule
on original primary truth-scorable overlapping reads. Ranking and screening
are preserved in `selection.json`, `ranked-blocks.json` and `screened-seams.json`.

## Cause and repair

The endpoint snarl contains an insertion and deletion, **55,309,794 T>TT /
TTT>T**, excluded as repeat indels. The verified BAM source already phases
**55,309,789 C>CT / CTT>C** on opposite haplotypes. These are the same two
complete T-repeat alleles at shifted positions. The source edits share physical
position 55,309,790; padded catalog keys place the insertion at 55,309,798
and deletion at 55,309,796. Sequence equivalence is checked against the common
reference background, rather than requiring identical CIGAR coordinates.

Terminal recovery accepts a single binary insertion and vetoes intervening
alternatives. That rejects both the verified source pair and the other graph
ALT even though they describe the same two haplotypes. The repair recognizes
this joint endpoint only when both catalog edits are sequence-equivalent to
an MSA- and alignment-verified complementary source pair in the same phase set.
Only those equivalent alternatives are exempt from the veto. Unrelated
intervening candidates still reject attachment. Graph/FASTA contig ID 0 and
original BAM contig ID 11 must be normalized before comparison; block extents
now use contig names so both channels contribute.

Original primary, nonduplicate MAPQ30 molecules must carry exactly one catalog
ALT, its matching complete physical repeat class, and agreeing Q30 SNP calls
at **55,309,475 G>A** and **55,308,854 C>T**. Both catalog alternative rows
must lie inside the molecule's graph profile. One read with only the insertion
row represented is excluded from joint calibration. Summed physical base and
mapping error is at most 1%. Deterministic name hashing separates the two
independent cohorts. Both cohorts must pass the existing diploid association
and 20% Wilson/error rule, agree on the source gauge, and jointly meet the
stricter 10% bound.

Independent original-CIGAR reconstruction (`audit_physical.py`) recovers the
production matrices **[[7,0],[0,10]]** and **[[9,0],[0,8]]**, 34 molecules.
Their separate error bounds are **14.7299%** and the combined bound is
**8.3709%**. Parental truth never enters production.

The two endpoint rows become phased nonanchors, retaining original counts:
insertion DP41, 16 REF / 25 ALT; deletion DP42, 16 REF / 26 ALT.
Two existing read-only rescues enter the core only after their own original
MAPQ30 physical deletion calls meet the 1% error bound, agree with the source
haplotype, and have no contradictory phased clean SNP. Their HP values are
unchanged. Unsupported rescues keep their original assignments.

## Results

| Endpoint reads, including abstentions | Before | Fixed | HiPhase |
|---|---:|---:|---:|
| Original primary truth-scorable | 67 | 67 | 67 |
| Correct | 66 | **66 (98.5%)** | 66 (98.5%) |
| Wrong | 1 | 1 | 0 |
| Unphased | 0 | 0 | 1 |
| Correct in one connected core | 64 | **66** | 66 |

The competitor's original sequence, CIGAR, qualities, flags and mapping
qualities are checked identically. The fix matches HiPhase's total and
connected-core correct counts. The one pre-existing low-MAPQ wrong rescue
remains; HiPhase abstains on that read.

The full **856,775 bp** span now shares one phase set with 747 phased rows,
up from 745. One-thread and four-thread owning candidates, VCF records and
read tags match exactly; the preceding 54–56 Mb continuation retains the
same endpoint and source gauge. Disjoint original molecule cohorts preserve
parental orientation, including the two held-out terminal-only reads.

Every previously correct and phased assignment, existing variant's evidence,
and block extent is preserved. Only two reads change PS, keeping HP.
Whole chr20 retains **254 blocks**, **N50 955,496 bp**, and largest block
**3,012,193 bp**. Among 256,610 primary reads, 230,765 remain correct,
6,658 wrong and 19,187 unphased. Exactly the two validated catalog endpoint
variants are added; existing alleles and counts do not change.

The endpoint satisfies the acceptance contract. Across the entire target,
pgphase still trails HiPhase by **13 correct reads / 35 connected-core correct
reads**; completing geometry does not claim that separate whole-block deficit
is fixed.

## Verification

The committed panel adds exact measured `spans=1`, required source/endpoint
alleles and both SNP witnesses, the certified 80%/HiPhase contract, an owning
regression and preceding-chunk continuation. All prior fixture rows and
HiPhase measurements are preserved. The old binary fails five assertions of the new regression;
the fix passes **655 assertions** in **0.60 seconds** from saved states.

Build and all unit tests, 1,551 predicate assertions, HiFi/ONT golden gates,
and **17,240 gap assertions in four cases**, plus four cache-helper tests,
pass with no new warnings. HiPhase measurements independently cover all
**128** committed windows on identical original alignments. All **234** previous native labels / **122** independent requests are
audited under the final binary: one label changes, with no lost correct/phased
assignment, variant evidence or block extent. The complete suite and cache
helpers finish in **39.15 seconds** from saved states. Native-label
preservation and timings are recorded in `panel-audit.json` and
`test-runtime.json`. The full audits are in `results.json` and
`full-preservation.json`; thread and continuation results are in
`owner-results.json` and `continuation-results.json`.

Final binary SHA256:
`ccd9f822ba30b7202a1cc140eb1dd80d333afd85b253e092c224c8b29b3ecd70`.
Run `bash replay.sh [output-directory]` for chr20, or
`make gap-owner-check GAP=55.309` for the focused saved-state check.
Use `LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu` and a Python with pysam for audits.
