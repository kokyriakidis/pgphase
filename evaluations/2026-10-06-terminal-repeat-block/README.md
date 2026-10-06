# Terminal insertion completes the eleventh-largest HiPhase block

The next qualifying open interval is **12,954,878–12,955,838 (960 bp)**,
at the right end of HiPhase's eleventh-largest chr20 block:
**11,796,969–12,955,838 (1,158,870 bp)**. The larger remaining seams fail
HiPhase's 80% all-original-read correctness contract; see `selection.json`.

The preceding fix already assigns all 84 original primary truth-scorable
reads overlapping this terminal interval correctly to one connected core.
The geometry remains open because the terminal catalog insertion
**12,955,838 T>TC** is classified `REP_HET`, excluded from the clean solve,
and not nominated by recovery that requires a bounded phase-set seam.
HiPhase phases the insertion with the same orientation as the preceding
SNPs, despite its low GQ=6. Its geometry does not require new information:
pgphase already has the catalog allele, GAF calls, original CIGAR and qualities.

## Repair and independent evidence

A final attachment pass checks an unphased binary one-base graph repeat
insertion beyond an existing core's final anchor. It requires two surviving
clean SNPs in that core within 10 kb, separated by at least 100 bp, with no
intervening anchor or repeat alternative. The insertion remains a nonanchor.
The pass changes no read tags or existing phase sets.

Original primary, nonduplicate Q30 BAM molecules must have an exact or
sequence-equivalent insertion call, agree with their GAF allele, physically
call both SNPs consistently, and already belong to the core. Nearby competing
indels, unknown qualities, missing calls and disagreements abstain. The summed
base and mapping error must be at most 1%.

Read-name hashing forms two disjoint molecule cohorts. Both must observe both
alleles and SNP haplotypes, pass the existing p<=0.01 association check and
agree on orientation. One-sided 95% Wilson discordance bounds plus a 1% call
error allowance are capped at 20% per cohort and 10% combined. The target
has **38 accepted molecules, zero conflicts**, with cohorts
`[[14,0],[0,7]]` and `[[10,0],[0,7]]` in the native owner. The measured
bounds plus allowance are **12.41%, 14.73% and 7.65%**, respectively.
`audit_physical.py` independently reconstructs these calls from original
alignments and recorded GAF observations. No parental truth enters production.

The insertion retains its catalog T/TC alleles and original GAF depth/counts
**73 / 44 REF / 29 ALT**. It now shares the core's SNP gauge in both the native
12–13 Mb owner and the 11–13 Mb stitched continuation. Single-thread and
four-thread outputs are identical. All 4,276 owning-chunk primary assignments
are preserved: 4,199 correct, 15 discordant and 62 unphased. Only the terminal
insertion is added to the native VCF; common variant evidence is identical.

## Verification

Build and all unit tests pass with no new warnings. The complete gap suite
passes **16,962 assertions in four test cases**, plus four cache-helper tests.
All 1,551 in-memory phasing predicate assertions and HiFi/ONT golden-output
and thread-determinism gates pass. HiPhase measurements cover all 125 panel
windows on identical original alignments.

The gap has **84/84 original reads correct and 84 connected-core correct**,
matching HiPhase on identical original alignments, with zero discordant reads
and no abstentions. The committed panel has an independently measured HiPhase
row, exact `spans=1`, connected-core parity and disjoint parental orientation
checks. The owner regression also checks the preceding chunk's stitched gauge.
The old binary fails the new regression; the new binary passes **621
assertions in 0.60 s** from saved replay states. All **228 previous native
labels (119 independent requests)** are audited under the final binary. Three
labels gain the target terminal insertion. One short 1.5 Mb replay also
gains an independently supported terminal 1,573,397 A>AT call outside its
panel gap; this does not change that gap's span or full-chromosome output,
where the core continues beyond the insertion. HiPhase also calls this
collateral insertion `1|0` (GQ=12, DP=70, AD=35/34), matching the local
pgphase orientation and catalog DP=70, AD=35/35. Every read tag and prior
variant is preserved. Prior required-site rows and read floors are unchanged.

Full-chromosome verification adds exactly this one variant and preserves
every read tag, every prior variant and all common variant evidence. Pgphase
now covers the complete HiPhase geometry with a **11,792,550–12,955,838
(1,163,289 bp)** block. There are still 256 blocks; **N50 is 944,265 bp** and
the largest block is 3,012,193 bp.

The gap contract applies to this interval. The previous full-block comparison
has a separate 19-read connected-core deficit across the entire 1.16 Mb block;
extending its variant endpoint does not itself fix that read-level deficit.

Evidence: `owner-results.json`, `physical-evidence.json`, `results.json`,
`full-preservation.json`, `panel-audit.json` and `test-runtime.json`.

Reproduce the full chromosome with `bash replay.sh [output-directory]`.
Use `LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu` for pipeline/tests and a Python
with pysam for the audit scripts. `make gap-owner-check GAP=terminal-insertion`
runs the cached owner and continuation checks.
