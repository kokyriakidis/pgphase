# Recovery representation audit

Date: 2026-09-22

## Scope

This audit starts from the retained SNP-first MEC chr20 output and examines the
16 current HiPhase-correct target gaps that remained open. Focused recovery runs
exported the graph/BAM allele matrix for every detected seam. Sites were compared
against the frozen HiPhase 1.6 VCF, then the retained change was rerun over full
chr20 with the same 500 kb chunks, 16 workers, MAPQ 5 graph input, recovery BAM,
and parental truth input as the baseline.

## Reproduced representation defects

At `chr20:14,264,549-14,272,741`, the BAM sub-solve keeps two deletion rows from
a multi-allelic event separate. Their row-wise AFs are 0.657 and 0.314, so the
centered-site gate rejected both. The direct deletion-to-SNP tables are pure:
`23 cross / 0 conflict` and `12 cross / 0 conflict`. A selected alignment-verified
boundary indel may now enter when no centered boundary site exists. The two rows
remain separate. The source-specific graph/BAM gauge is still required. A
focused run joins the gap at 35/35 truth-correct crossing reads. A later
full-chunk replay showed that the source gauge agrees with the local relation:
275 shared reads and 59 sequence-identical candidates all request the same
cross orientation. The full-chunk failure instead came from expanding the exact
MEC problem across the complete neighboring block, which pulled in unrelated
unphased variables and exceeded the 20-variable search bound. The guarded local
retry and its whole-chromosome result are recorded in
`evaluations/2026-09-23-local-edge-retry/`.

At `chr20:14,679,241-14,679,247`, the BAM and graph boundaries are two different
SNPs six bases apart. Their 2x2 allele table is pure (`76 support / 0 conflict`),
but the old validator required one sequence-identical candidate and rejected the
edge. The nearest centered SNP pair may now validate the source read gauge when
its full read set and both deterministic halves independently choose the same
parity at one-sided binomial `p <= 0.05`. The full chr20 run closes the
gap with 76/76 local truth purity.

## Safety experiment

Allowing the same recovered block to use the new relation for a second attachment
was rejected. It merged phase sets across 54.49 Mb and changed 342 previously
correct reads to discordant. Locally stable boundary SNPs cannot rule out an
older polarity change elsewhere in an atomic BAM block. Production therefore
retains one trusted fallback attachment per original block and never uses the
new validation to join BAM/BAM blocks.

With that invariant restored, no previously evaluated read changes between
concordant and discordant. Seven previously skipped reads become evaluable: six
concordant and one discordant. The phased BAM has the same 225,005 tagged reads.
The VCF has 15 fewer phase blocks, while the tracked panel gains the 14.679 Mb
closure and loses none. Full metrics are in `results.tsv`.

## Remaining failure classes

Of the other current target gaps inspected:

- 0.865 Mb has tied full-read MEC and a disconnected read half.
- 21.669 Mb has no two-site read row in the trusted component.
- 21.823 Mb and 46.727 Mb choose opposite parities in the two read halves.
- 57.854 Mb has physical bridge reads, but their shifted repeat deletions are
  unknown at the exact injected rows; no read observes both an eligible injected
  indel and the right SNP. Net-length genotyping shows several competing deletion
  lengths, so treating every nonmatching deletion as reference would invent a
  binary allele and was not adopted.
- Several low-support targets require joining two independent BAM blocks. The
  retained fallback deliberately abstains because local evidence cannot exclude
  an older switch inside either block.
