# Certified insertion reads reach HiPhase core coverage

The preceding short-insertion source-path fix closed the 33,747,591–33,749,688
VCF gap and correctly phased 64/67 overlapping reads, but four correct reads
remained in an output-only rescue phase set. HiPhase put all 64 correct reads
into its connected core. The user's stronger acceptance rule therefore left
that change provisional.

## Final behavior

After a short insertion's independent complete BAM-only source and physical
SNP-to-insertion bridge certify the block gauge, unassigned reads can join
the core using an independent physical insertion call. Primary BAM MAPQ must
be at least 30; base-quality floor is 20; the combined mapping and measured
allele/base error bound must not exceed 0.01. Imported and primary profile
calls must agree with the physical allele. Existing core assignments remain
untouched. Other informative loci must agree in HP and phase set; opposing
complementary descriptions at one locus provide no independent vote.

REF confidence uses actual anchor qualities rather than the insertion caller's
returned configured floor. A disjoint nearby deletion can prevent that caller
from returning any insertion allele. The new read-admission check can still
establish REF when the entire insertion-length footprint matches the reference,
every footprint base is callable and included in the error sum, and no inserted
sequence appears within 64 bases. The ordinary physical insertion caller and
bridge confidence thresholds are unchanged.

An excluded graph description remains contrary evidence. Normalize its selected
REF/ALT edit and test sequence equivalence within 64 bases, using the existing
repeat-equivalence helper. If its observed graph allele disagrees with the
physical call, reject core admission even when that description has no phase
label. Mere coordinate equality misses shifted repeat representations.

The four diagnosed reads are independently callable: one exact TTCT insertion
with Q40 inserted/flanking bases and three REF calls with local minima Q27,
Q27 and Q40. Read 58000534 has an unrelated one-base deletion 11 bases before
the insertion; the exact REF footprint check recovers its otherwise rejected
call. The four source rows have HP0 and no heterozygous observations counted
by the source solve. Both old core refresh paths accepted SNP evidence, while
rescue always wrote the phase-set offset; neither could make this assignment.
Original per-read evidence is retained in the preceding evaluation's
`core-deficit-investigation.json`.

## Controls and measurements

The owning 33–34 Mb replay changes exactly the four diagnosed PS tags from
1033320808 to 33320808, retaining each HP. All output molecules, parental
correctness statuses and VCF bytes remain identical to the provisional output.
Owner counts remain 3,655 scored, 3,523 correct and 132 discordant. The new
checks do not phase either of the two conflicting/unassigned molecules pinned
by the regression.

A first trial using the configured REF floor only admitted the ALT read.
Measuring actual anchors admitted three reads, but the nearby-deletion veto
still excluded the fourth. The footprint fallback admitted all four and also
incorrectly admitted 181539002. That read's graph ALT observation conflicts
with physical/MSA REF, at a graph repeat insertion described 12 bases to the
right of the BAM edit. Exact-position conflict matching did not catch it;
sequence-equivalent edit matching does. The final owning audit confirms that
it stays unassigned, with no new discordance. No parental label or read name
is used in production conditions.

Matched HiPhase comparison verifies identical original alignment coordinates,
CIGARs and sequences for all 67 molecules. pgphase now has 64 correct / 1
incorrect / 2 unphased and all 64 correct reads in one core block. HiPhase has
64 correct / 2 incorrect / 1 unphased, also with 64 correct in one block.
Both correctness and core coverage meet parity; 64/67 = 95.5% exceeds the 80%
floor. `compare_hiphase.py` asserts both counts and no greater discordance.

The owning regression requires core separation >=64/67, pins the four HP/PS
assignments and preserves abstention for the two conflicting reads. The panel's
measured separation floor rises from 0.89 to 0.95. All previous read-count and
accuracy floors remain unchanged. `audit_core.py` rejects changes to any other
read tag, any read's parental correctness, or any VCF byte.

## Full chromosome

The final full-output audit confirms exactly the same four PS-only changes
from the provisional binary. All 256,610 output molecules, 237,334 scored
reads, 230,566 correct / 6,768 incorrect reads, 645 read phase sets, 64,188
VCF keys, 325 VCF blocks and span N50 856,770 remain unchanged. Every read's
parental correctness is preserved; VCF bytes are identical. No other core
assignment changes. The final matched comparison verifies 64/67 correct in
one block, 1 wrong and 2 unphased, versus HiPhase's 64/67, 2 wrong and 1
unphased. Disjoint flank orientation votes are 25:0 and 24:7 in the same
parental orientation, both binomial-supported at p<=0.05.

Relative to the accepted output before the original source-path closure,
exactly that one gap closes, no old span reopens, and all old keys, genotype
alleles, nonphase fields and candidate-block gauges remain. Overall correct
reads increase by three and incorrect reads by one; that earlier outside-gap
rescued-read change remains recorded in the preceding evaluation. The core
admission fix adds no further error. Final full panel spans remain 97/110.

Reproduce with the bench-phasers Python and
`LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu`:

```sh
python evaluations/2026-10-05-certified-insertion-core/audit_core.py \
  --before test_data/tmp_gap_fix61/full-bam-only-final/0 \
  --after test_data/tmp_gap_fix62/full-certified \
  --output /tmp/insertion-core-audit.json
python evaluations/2026-10-05-certified-insertion-core/evaluate_gap.py \
  --before test_data/tmp_gap_fix60/full-final \
  --after test_data/tmp_gap_fix62/full-certified \
  --output /tmp/insertion-gap-audit.json
python evaluations/2026-10-05-certified-insertion-core/compare_hiphase.py
```

## Validation

The frozen final production binary passes all 84 native window cases / 10,779
assertions, including all 110 panel rows and existing required-site witnesses.
The updated test source checks exact core assignments and conflicting-read
abstentions in the same sweep. Build has no warnings; `make unit-tests`,
`make predicate-tests` (47 cases / 1,536 assertions) and `make check` pass.
Full core-change, gap-geometry and matched HiPhase audits all pass.
`manifest.json` records input, source and binary hashes; `validation.json`
records shard results and log hashes. The final binary SHA256 is
`db17374ae2f37a68e803ce397b09f90a75705b0bff23b40aea1498aeccd662a0`.
