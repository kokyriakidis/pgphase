# One mode: the hybrid

The catalog's sites phase, the alignment recovers the gaps. There is no second
configuration to choose between.

## What was removed

| removed | why |
|---|---|
| `collect-graph-variation --bam` and its recovery | it was a second mode for the same job; measured equal to the hybrid on one window and 118 correctly-phaseable reads short on the other |
| `--graph-first` / `--no-graph-first` | graph-first IS the mode now, so a flag for it is a flag for "the mode", and its negation is the second mode |
| the `nographfirst` window-test arm | nothing left to select |

`recover_unphased_windows_from_bam` stays -- it is the hybrid's own recovery
entry, and it keeps the `contig_name` and `allow_import` parameters because the
tid translation and the index-safety they exist for are real constraints, now
documented at one caller rather than two.

Behaviour is unchanged: on both panel windows the emitted VCF is **identical**
to the pre-consolidation default, byte for byte over the records. Window tests
66 assertions, unit 4/4.

## The efficiency question, measured rather than assumed

The obvious target was the alignment discovery the first solve appears to throw
away: under graph-first every alignment-derived candidate is withheld before the
catalog's sites are injected, so the discovery pass looks like wasted work.

It is not. Skipping `collect_candidate_sites_from_records` and
`collect_allele_counts_from_records` outright on
`chr20:25,979,591-26,138,679`, single-threaded:

| arm | wall time | candidates | phased hets |
|---|---:|---:|---:|
| as shipped | 4.05 s | 448 | **447** |
| discovery skipped | **2.78 s** | 137 | **137** |

31% faster for a collapse from 447 phased heterozygotes to 137. The candidates
are discarded but their by-product is not: `pre_process_noisy_regs_pgphase` and
the post-classification passes build `chunk.noisy_regions` from them, that model
scopes the noisy-region MSA, and the MSA supplies most of the chunk's phased
sites. Discovery feeds the phasing even though its candidate rows do not reach
the solve.

So the redesign is a consolidation, not a restructuring: discovery runs, its
candidates are withheld from the first solve, the catalog's sites phase, and the
alignment returns inside the windows the solve could not phase. The comment at
the withholding site records the measurement so the next reader does not retry
the same optimisation.

## What a real efficiency gain would need

Separating the noisy-region model from candidate discovery -- deriving
`noisy_regions` from the read digars alone, without building and classifying a
candidate table that is then dropped. `pre_process_noisy_regs_pgphase` already
works from reads; the post-classification refinement is what needs candidates.
Whether the refinement is load-bearing is the measurement to take, and it is not
taken here.

Not measured: chromosome-wide.

## Removing what one mode made dead

Three leftovers were flagged after the consolidation. Two were real.

**`--retry-unphased-with-bam` was redundant** -- the recovery is the default, so
the positive flag selected the default. Removed. `--no-retry-unphased-with-bam`
stays: it is the `noretry` window-test arm, which is how a regression gets
attributed to the recovery rather than the first pass.

**The whitelist mode was dead or destructive, not merely unused.** Under one mode
the alignment channel's candidates are cleared before injection, and the
retention filter ran *after* that clear -- on an empty table, so it could retain
nothing. `--bam-authoritative-bed` excluded catalog sites inside its intervals so
the alignment would own them; with the alignment's candidates already withheld,
that left the interval with no sites at all, so the option only deleted
evidence. And `--private-sites` reached a phasing branch that set
`max_noisy_reg_len = 0`, silently disabling the noisy-region model that supplies
most of the chunk's phased sites.

Removed: `--private-sites`, `--private-msa`,
`--private-msa-admit-all-in-region`, `--private-msa-snp-first`,
`--bam-authoritative-bed`, and with them `load_private_variant_keys`,
`retain_private_bam_candidates`, `is_bam_authoritative_position`,
`load_bam_authority_intervals`, the `BamAuthorityIntervals` type, the
`private_keys` and `bam_authority` parameters threaded through the chunk path,
three CLI validation blocks whose only job was rejecting combinations of these
options, and the unit-test block covering the whitelist. **284 lines removed, 29
added.** The two MSA booleans the removal left behind were default-false, so
their use sites now read `false` with the reason recorded.

`--private-msa-margin` is **renamed `--msa-ambiguity-margin`**, not removed: it
governs the MSA consensus rescoring in `align.cpp` for every noisy region and was
never specific to the whitelist.

**The third leftover was not a defect.** Five flags appeared twice in the help --
`--ref`, `--bam`, `--gaf`, `--graph-sites`, `--phased-vcf-out`. The second
occurrence of each is the usage example at the end of the help text, which is
correct; the duplication was an artifact of the grep that found it. Withdrawn.

Behaviour is unchanged: both panel windows emit a VCF identical to the
pre-removal default. Window tests 66 assertions, unit 4/4.
