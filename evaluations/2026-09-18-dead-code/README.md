# What the removals left behind

`-Wall -Wextra` is already on. After the mode consolidation it reported twelve
functions defined but never called and one unused parameter -- all of them left
by this session's own removals, including a commit whose message claimed "no
unused function in any touched unit".

**485 lines removed, 12 added, behaviour identical on both panel windows.**

## Orphaned by `1ef17ab` (the recovery-path removal)

| function | what used it |
|---|---|
| `gap_fill_unphased_reads` | `--gap-fill` and `--msa-verified-refine` |
| `clean_snp_has_bam_observation` | the three-SNP gap bridge |
| `msa_indel_has_confident_bam_observation` | the in-gap indel evidence backfill |
| `build_cached_gap_proposal` | the gap evidence cache |

Each was checked against history before deletion rather than assumed dead:
`git log -S` dates every one of them to that commit, and none had a caller
after it. `gap_fill_unphased_reads` was the one worth checking carefully -- it
was the shipped stage-2 read recovery -- and its only callers were the two
removed modes, so its loss is by design and not a silent regression.

## Orphaned by `4abf46f` (the whitelist removal)

`resolve_private_vcf_tid`, `authority_bed_tid`, `authority_bed_contig`.

## Orphaned by removing the graph pipeline's recovery entry

`recover_unphased_windows_from_bam` had no caller once the graph arm went -- it
has external linkage, so the compiler stayed silent about it.

**Correction: 2a9a33d did not remove it.** Three attempts asserted out before
their write, and the fourth deleted only the header declaration -- leaving an
externally-linked definition with no declaration, and the six-parameter
signature of `recover_windows_with_targeted_solve` untouched. The commit message
and the first version of this document claimed both were done. What made the
error look plausible was a second wave of warnings appearing in the same cell:
`clone_cached_read`, `merge_cached_allele` and `remap_cached_allele` had in fact
been orphaned by removing `build_cached_gap_proposal` in the *previous* cell, not
by removing this entry point. Actually removed in the follow-up commit, with the
token check now run after the write rather than before it.

`bam_base_quality_at` and `equivalent_shifted_insertion` were genuine
second-wave orphans. Two rounds of build-and-remove were needed to reach zero
warnings, which is the argument for looping rather than sweeping once.

With one caller left, `recover_windows_with_targeted_solve` no longer needs its
`solve_tid` or `allow_import` parameters -- both took a single constant value.
(Also landed only in the follow-up commit, for the same reason.)
The reasons they existed are kept as comments where the region is built: a caller
whose header is not the BAM's must translate by contig **name**, and appending
candidates to a table that is index-parallel to per-site metadata mispairs it.

## The seam recovery was checked, not assumed

Removing the entry point could have taken block-seam recovery with it, since
that entry was where the two window sources were combined. It did not: the
hybrid's own retry path calls `collect_unphased_windows` **and**
`collect_block_seams` directly. Verified by reference, not by inspection of the
diff.

## One unused parameter, and what it exposed

`wfa_collect_noisy_aln_str_no_ps_hap` took an `unassigned` out-parameter and
ignored it, while its with-phase-set sibling fills one. The sibling's list is
built by rescoring reads the clustering left unplaced against the two haplotype
consensuses; this variant runs when no phase set carries reads from both
haplotypes, so there is no haplotype labelling to rescore against. Its
declaration defaulted the parameter to `nullptr` and no caller ever passed one,
so there was no silent data loss through it. Removed.

What the check did surface is a real limitation, measured rather than inferred:
that branch builds its consensuses from **full-cover reads only**, so partial
reads are unused there. On the two panel windows it runs 12 and 6 times against
the sibling's 262 and 255, leaving 8 and 2 partial reads unused. Recovering them
would mean rescoring against the branch's own two consensuses when it produces
two -- an addition with a measurable target, not a bug fix, and not attempted
here.

Zero warnings in our translation units now; the remaining one is in vendored
abPOA. Window tests 66 assertions, unit 4/4.
