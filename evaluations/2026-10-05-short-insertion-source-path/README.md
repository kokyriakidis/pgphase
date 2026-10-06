# Short insertion source path at 33,747,591–33,749,688

Two rules excluded a supported BAM-only block: newly exposed insertion seams
were nominated only for long insertions, and a short insertion's right block
needed a graph SNP path even when it contained only independently recovered
BAM sites. The final change nominates source-supported short insertions and
allows a complete original BAM source to certify a wholly covered BAM block.
Every phased heterozygous anchor must be BAM-injected, have exactly one
adoptable claim from that same source, and agree with its original allele
gauge. The source has no weak or quality cut and has a distinct anchor.
Mixed blocks retain the graph-path requirement. Long-insertion behavior is
unchanged. No coordinates or parental labels enter production conditions.

Physical short-insertion support still requires primary MAPQ/base-quality-30
calls, at least two REF and two ALT molecules, a quality-weighted wrong-parity
bound of 0.001, and the complete left graph path. A diagnostic replay found
24 REF and 10 ALT insertion calls, weighted log odds -164.860097239,
consistent with the existing relative gauge. The endpoints are a clean
33,747,591 T>C SNP and a verified four-base 33,749,688 C>CTTCT insertion
(unanchored candidate coordinate 33,749,689).

## Owner replay and parental truth

The owning 33–34 Mb replay reduces 11 VCF blocks to 10. Scored reads change
3,651→3,655; correct 3,520→3,523; discordant 131→132. Four formerly unphased
reads become correct and one previously correct rescued read becomes wrong.
All old keys, genotype alleles, nonphase fields and uniform old candidate-block
gauges survive. The gap itself has 67 truth-scorable primary molecules:
60 correct / 1 wrong / 6 unphased becomes 64 correct / 1 wrong / 2 unphased.
Thus 95.5% of all eligible reads are correct, including abstentions in the
denominator, exceeding the user's 80% rule. Dominant correctly separated
reads increase 59→60; four new correct reads remain an output-only rescue
cohort. Disjoint left and right flank parental votes are 25:0 and 20:7 in
the same orientation, each with a one-sided binomial tail below 0.05.

The new committed panel row has `spans=1`, concordance >=0.96 and separated
>=0.89. Its owner replay checks absolute correct/discordant counts, parental
orientation and both endpoint allele gauges. All existing required sites
remain; the two new endpoint witnesses are appended.

## Rejected broader source certificate

Allowing any complete source to certify a short insertion's block also joined
36,620,864–36,623,545. That mixed graph/BAM block had only 52/75 correct reads
(69.3%) and weak left-flank orientation. The full trial improved global counts
but failed the local acceptance rule; it is not the accepted implementation.
`rejected-broad-full-audit.json`, `rejected-broad-full-parity.json` and
`rejected-mixed-gap-audit.json` record that trial. Restricting the certificate
to all anchors covered by the same independent source preserves the target
owner output and leaves the entire 36–37 Mb VCF identical to baseline.
The existing owning 36 Mb regression now explicitly checks that these two
mixed-block endpoints keep distinct phase sets.

## Matched HiPhase comparison

All 67 HiPhase input alignment coordinates, CIGARs and sequences match the
original BAM. HiPhase also has 64 correct reads, with 2 wrong and 1 unphased,
and places all 64 correct reads in one core block. pgphase has 1 wrong and
2 unphased, with four correct insertion-only reads in a separate output-only
rescue phase set. Those four begin after the clean left SNP and have no
original source HP assignment: the difference is admission to the core,
not missing input evidence.

Of the two reads correct only in HiPhase, one has a pre-existing graph REF
call conflicting with a Q40 BAM ALT call and original source HP2. The physical
bridge preserves existing core assignments rather than rephasing that read.
The other observes only the recovered insertion, with graph ALT versus MSA
REF and no qualifying base quality, so pgphase abstains. Two reads correct
only in pgphase end before that insertion and have Q40 clean-SNP calls
consistent with the original BAM source. These observations explain our
confidence decisions; they do not establish HiPhase's internal algorithmic
choice without an internal trace. See `hiphase-comparison.json` and
`read-disagreements.json`, reproduced by the accompanying scripts.

`evaluate.py` rejects lost calls, altered genotype alleles/nonphase fields,
mixed old candidate-block gauges, extra new gaps and reopened old spans. It
checks the 80% rule using whole-output parental orientation and disjoint flank
orientation independently. `manifest.json` records the frozen binaries,
input hashes, outputs and reproduction commands.

## Final full-chromosome audit

The frozen final binary closes exactly the target gap and reopens no tracked
gap. All 64,188 variant keys, genotype alleles, nonphase fields and uniform old
candidate-block gauges are preserved. Tagged output reads remain 256,610;
truth-scored reads change 237,330→237,334, correct 230,563→230,566,
discordant 6,767→6,768. VCF blocks change 326→325 and read phase sets 647→645;
span N50 remains 856,770. Only 35 read tags and three VCF rows change. The
newly discordant rescued molecule is
`m84031_231217_062403_s3/113902771/ccs`, outside the accepted gap.
The panel retains all 96 baseline full-output spans and gains the target,
yielding 97/110. The matched HiPhase comparison uses this final whole-chromosome
output, not a short replay.

Reproduce the independent full audit with the bench-phasers Python and
`LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu`:

```sh
python evaluations/2026-10-05-short-insertion-source-path/evaluate.py \
  --before test_data/tmp_gap_fix60/full-final \
  --after test_data/tmp_gap_fix61/full-bam-only-final/0 \
  --output /tmp/short-insertion-full-audit.json
python evaluations/2026-10-05-short-insertion-source-path/compare_hiphase.py
```

The before binary is `3c2cb7ef7a5070416ff98a04e2c7bf9ef6b9c5501a66d9814bcd75cd8a091924`;
final is `3374fb69e745c152bbc862fbab65fc62888a58260fc6a7ac61040bd423a0c00b`.

## Validation

The final production binary passes all 84 native window cases / 10,752
assertions, including all 110 panel rows and existing required sites. The
final test source's added mixed-block guard is checked in a separate owning
rerun alongside the new target (652 assertions / two cases). No existing
read-accuracy floor or unrelated negative span is relaxed. Build is clean;
`make unit-tests`, 47 phase-predicate cases / 1,536 assertions, and `make check`
pass. `validation.json` records shard counts, return codes and log hashes.

## Acceptance superseded by the HiPhase parity requirement

The user subsequently requires every closure to match or beat HiPhase,
including connected core coverage. This change is therefore provisional:
64/67 total correct ties HiPhase, but 60 correct in the dominant core block
falls short of its 64. The four correctly rescued insertion-only reads need
supported core assignments before this gap meets the new acceptance rule.
Passing the earlier 80% audit is insufficient under the stronger requirement.

## Detailed core admission investigation

`core-deficit-investigation.json` traces all four reads. Each has one effective
heterozygous insertion marker at 33,749,689. Three also call REF on the
complementary AAA/AAAA insertion descriptions at 33,762,265; their opposite
haplotype gauges cancel that locus in the rescue scorer. Other SNP calls are
homozygous or noninformative. Both core refresh paths accept SNP evidence;
the rescue path unconditionally writes separate `gap_haps`/offset phase sets.
There is no physically certified insertion-to-core promotion branch.

Original primary BAMs provide MAPQ60 evidence: one exact TTCT insertion with
Q40 inserted/flanking bases and three reference alignments with local minimum
qualities Q27, Q27 and Q40. An insertion-aware core assignment can therefore
be investigated using existing physical information; a strict Q30 floor alone
would miss two. Merely renaming rescue phase sets would bypass the missing
confidence check. This establishes pgphase's admission mechanism, without
claiming an unobserved internal HiPhase decision path.

The subsequent insertion-core fix resolves this deficit; see
[`2026-10-05-certified-insertion-core`](../2026-10-05-certified-insertion-core/README.md)
for the final HiPhase parity measurements and admission checks.
