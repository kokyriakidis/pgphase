# Close the seam in HiPhase’s fourth-largest chr20 block

The target is 52,696,940–52,711,825, a 14,885 bp boundary distance inside
HiPhase’s 52,696,410–54,419,033 block (1,722,624 bp). The baseline production
binary is `7e50015b328beeed83607d502a1e6d6a928d2fea8dbbfb77328c3e5610e76be7`;
the fixed optimized binary is
`9c55c8ac6b67240284980560f308b104fee1dc776b10fd0b89988a10aef251c9`.
The target ranking and unchanged original-read competitor audit are in
`../2026-10-05-fourth-largest-block-target/`.

## Cause and fix

The two left anchors are recovered BAM calls: GGA>G at 52,696,410 and G>A at
52,696,940. Their complete, uncut BAM source passes its own path check, but
there is no left graph SNP path. The right source mixes graph and BAM rows.
It starts with a short A insertion at 52,711,825, whose physical bridge calls
have REF support but no reliable ALT witness. Retry nomination required a
BAM-only source or a long insertion, and the deletion fallback stopped at that
insertion’s coordinate. Requiring a left graph path also rejected the valid
BAM-only left path.

The retry now admits a supported mixed source, then uses its unique verified
pure deletion before the first clean right SNP when the insertion bridge
abstains. TTGTG>T at 52,715,881 is such a marker, with its BAM consensus and
alignment verification retained. Its exact key deletes four bases starting at
52,715,882. The complete right graph path and original uncut deletion source
must still pass; another phased block cannot intervene.

Original primary molecule `m84031_231217_062403_s3/163778931/ccs` calls the left
ALT SNP at Q27 and the exact deletion with surviving matching flanks at Q10
and Q17. Both mapping qualities are conservatively counted. The union bound
for that pair is **0.1219498855**. A second overlapping deletion molecule has a
Q3 surviving flank and abstains. Disjoint molecules calibrate the deletion
against physical right clean SNPs: **7 ALT→hap1, 4 REF→hap2, no contradictions**.
Both classes have independent support and the association p-value is
**0.0009765625**. The summed actual base/mapping error across the accepted
bridge and all 11 calibration molecules is **0.1351929490 <=0.20**.
`audit_bridge.py` reproduces this certificate without parental truth.

This is a fixed-molecule quality bound, not a population discordance interval;
the significant association check is separate from the quality bound. The
fallback requires complete independently supported source paths. Existing
sparse-indel bridge predicates retain their previous limits. Contradictory
bridges/calibration, missing classes, weak association, unknown qualities and
an excessive joint quality bound abstain.

The union is deferred until existing rescue passes finish. All captured
anchors must still be present in one consistent final orientation. A second
bug prevented complete source rescues from entering the joined core: those
later rescue calls are absent from the initial BAM solver’s membership list,
and requiring each read to reach the deletion misses reads ending earlier.
Attachment now checks the whole uncut source, its unchanged pre-join gauge,
an agreeing recorded BAM observation and absence of a contrary clean phased
SNP. It retains each saved haplotype and creates no new HP assignment. Truth
is used only to evaluate the result.

## Native owner and permanent regression

The original 136 primary truth-scorable seam overlaps are the denominator,
including abstentions. Rescue PS >=1e9 is excluded from the connected core.

| Seam metric | Before pgphase | Fixed pgphase | HiPhase |
|---|---:|---:|---:|
| Correct | 119 | 119 | 117 |
| Discordant | 8 | 8 | 0 |
| Unphased | 9 | 9 | 19 |
| Largest correct connected core | 47 | 119 | 117 |
| Correct / all 136 original reads | 87.50% | 87.50% | 86.03% |

The two original local cores contain 47+12 correct reads. All 60 correct
rescues now enter their union. The closure passes the >=80% all-read contract
and both total/core HiPhase floors. It preserves the existing eight erroneous
haplotype calls; the connection does not claim to correct them.

The native 52,000,001–53,000,000 owner preserves all **4002 output reads**,
**3652 correct /89 discordant /261 unphased**, and all **719 variant records**.
Alleles, counts, filters and every previous correct/tagged assignment remain.
265 variant phase labels and 1278 read tags change through union/attachment.
One/four-thread candidate TSVs, VCFs and read assignments are identical.
The development owner took **9.9 seconds**; final one/four-thread replays
took **10.7/11.0 seconds** under concurrent verification. The complete repeated
owner regression takes **0.57 seconds** from its saved state, rerunning assertions.

The closure is added to the panel, native replay manifest, strict
certification/HiPhase manifests, required markers and read floors. Its owning
regression checks marker gauges, a retained original bridge witness, read/core
parity and disjoint parental cohorts. Only two original reads qualify for the
50 kb disjoint left flank, so that flank is not used to claim orientation.
Instead, left-SNP-bearing reads ending before the deletion independently
agree **54:1**, and downstream-SNP reads starting beyond the left SNP agree
**1206:0**. The bridge molecule belongs to neither cohort. The old binary
fails **eight assertions**; the fixed focused case passes **1295 assertions**.

## Full chromosome result

The frozen 67-chunk replay closes the complete **52,696,410–54,419,033** span,
**1,722,624 bp**, with 2127 phased heterozygous pgphase rows. The seam again
scores **119 correct /8 discordant /9 unphased**, with **119 correct core
reads**, versus HiPhase’s 117. All 7030 competitor alignments overlapping the
full target match original geometry, CIGAR and sequence. Every old pgphase
block extent remains covered; the number of blocks decreases **261→260**.
N50 stays **904,351 bp**, and the largest block stays **3,012,193 bp**.

All 256,610 output reads and all 64,484 variant records remain. Chromosome
scores are unchanged: **230,612 correct /6,732 discordant /19,266 unphased**.
No previously correct or phased assignment is lost; allele/count/filter and
genotype evidence is preserved. The connection changes 2125 variant phase
labels and 6872 read tags.

**Full-block read parity is not achieved by this seam fix.** Across the 7030
original target-block reads, pgphase still has **6935 correct /18 discordant
/77 unphased**, while HiPhase has **6963 /10 /57**. The connected correct core
improves **6799→6911**, versus HiPhase’s 6963. Thus the target still has a
separate **28-correct-read /52-core-read deficit** even though its endpoints,
internal connection and local seam contract now pass.

## Validation

The optimized build has zero warnings. Unit tests, **1551 predicate assertions
in 47 cases**, gap development checks, HiFi/ONT TSV/VCF goldens and HiFi
one/four-thread determinism pass. The complete window suite passes **15,966
assertions in four cases**, and the four replay-cache helper tests pass.
An initial run found a formatting error in the newly added fraction floor;
it now stores the measured decimal 0.875. No older expectation changed.

All **221 task-start native output labels**, covering **116 independent replay
requests**, have verified current-binary states and unchanged results. Every
old correct/phased read and variant allele/count/filter is preserved, with no
increase in discordant reads. Cache warmers share the existing locks and input
fingerprints; all assertions execute on every test invocation. No production
truth input, new option or coordinate-specific exception was introduced.

## Reproduction

```bash
make -j8
make gap-dev-check
make gap-owner-check GAP='fourth largest'
make unit-tests predicate-tests check
make window-tests

/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-calibrated-source-deletion/audit_bridge.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-calibrated-source-deletion/audit_owner.py

evaluations/2026-10-05-calibrated-source-deletion/replay.sh
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-calibrated-source-deletion/audit_block.py
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-05-calibrated-source-deletion/audit_panel.py
```
