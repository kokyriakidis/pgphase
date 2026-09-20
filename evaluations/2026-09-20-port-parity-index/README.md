# Parity index: which functions carry longcallD's name, and where

Generated for side-by-side checking. A shared name means the two sources can
be read against each other directly; it is not by itself evidence that the
bodies agree, which is what the per-function tests are for.

## Functions that now share a name with upstream

| function | ours | upstream |
|---|---|---|
| `abpoa_partial_aln_msa_cons` | align.cpp | align.c |
| `add_phase_set` | align.cpp | align.c |
| `assign_hap_based_on_germline_het_vars_kmeans` | collect_phase.cpp | assign_hap.c |
| `calc_read_error_rate` | align.cpp | seq.c |
| `cal_wfa_partial_aln_beg_end` | align.cpp | align.c |
| `check_agree_haps` | collect_phase.cpp | assign_hap.c |
| `collect_aln_beg_end` | align.cpp | align.c |
| `collect_full_msa_digars` | align.cpp | align.c |
| `collect_left_digars` | align.cpp | align.c |
| `collect_left_msa_digars` | align.cpp | align.c |
| `collect_noisy_read_info` | align.cpp | align.c |
| `collect_noisy_reg_aln_strs` | align.cpp | align.c |
| `collect_noisy_vars1` | collect_phase_noisy.cpp | collect_var.c |
| `collect_partial_aln_beg_end` | align.cpp | align.c |
| `collect_phase_set_with_both_haps` | align.cpp | align.c |
| `collect_read_var_profile` | collect_var.cpp | collect_var.c |
| `collect_reg_ref_bseq` | collect_phase_noisy.cpp | seq.c |
| `collect_right_digars` | align.cpp | align.c |
| `collect_right_msa_digars` | align.cpp | align.c |
| `collect_var_main` | collect_var.cpp | collect_var.c |
| `cr_extend_noisy_regs_with_low_comp` | collect_var.cpp | collect_var.c |
| `edlib_xgaps` | align.cpp | align.c |
| `exact_comp_cand_var` | collect_var.cpp | collect_var.c |
| `exact_comp_var_site` | collect_var.cpp | collect_var.c |
| `exact_comp_var_site_ins` | collect_var.cpp | collect_var.c |
| `full_cover_cmp` | align.cpp | align.c |
| `get_cs_from_digar` | collect_bam_output.cpp | bam_utils.c |
| `get_digar_ave_qual` | collect_var.cpp | bam_utils.c |
| `get_full_cover_from_cons_aln_str` | collect_phase_noisy.cpp | collect_var.c |
| `get_full_cover_from_ref_cons_aln_str` | collect_phase_noisy.cpp | collect_var.c |
| `get_md_from_digar` | collect_bam_output.cpp | bam_utils.c |
| `get_nm_from_digar` | collect_bam_output.cpp | bam_utils.c |
| `get_var_allele_i_from_cons_aln_str` | collect_phase_noisy.cpp | collect_var.c |
| `get_var_init_max_cov_allele` | collect_phase.cpp | assign_hap.c |
| `get_var_site_start` | collect_var.cpp | bam_utils.c |
| `init_assign_read_hap_based_on_cons_alle` | collect_phase.cpp | assign_hap.c |
| `is_collectible_var_digar` | collect_var.cpp | collect_var.c |
| `is_cover_aln_str` | collect_phase_noisy.cpp | collect_var.c |
| `is_homopolymer` | align.cpp | align.c |
| `is_match_aln_str` | collect_phase_noisy.cpp | collect_var.c |
| `is_match_aln_str_del` | collect_phase_noisy.cpp | collect_var.c |
| `iter_update_var_hap_cons_phase_set` | collect_phase.cpp | assign_hap.c |
| `iter_update_var_hap_to_cons_alle` | collect_phase.cpp | assign_hap.c |
| `low_comp_cr_start_end` | collect_var.cpp | collect_var.c |
| `main` | main.cpp | main.c |
| `make_cand_vars_from_baln0` | collect_phase_noisy.cpp | collect_var.c |
| `make_cand_vars_from_msa` | collect_phase_noisy.cpp | collect_var.c |
| `make_cons_read_aln_str` | align.cpp | align.c |
| `make_ref_read_aln_str` | align.cpp | align.c |
| `make_vars_from_msa_cons_aln` | collect_phase_noisy.cpp | collect_var.c |
| `merge_read_var_profile_entries` | collect_phase_noisy.cpp | collect_var.c |
| `merge_var_profile` | collect_phase_noisy.cpp | collect_var.c |
| `push_digar0` | align.cpp | bam_utils.c |
| `push_digar_alt_seq` | align.cpp | bam_utils.c |
| `read_init_hap_phase_set` | collect_phase.cpp | assign_hap.c |
| `read_to_cons_allele_score` | collect_phase.cpp | assign_hap.c |
| `refine_bam1` | collect_bam_output.cpp | bam_utils.c |
| `select_init_var` | collect_phase.cpp | assign_hap.c |
| `sort_noisy_region_reads` | align.cpp | align.c |
| `sort_noisy_regs` | collect_phase_noisy.cpp | collect_var.c |
| `update_bam1_tags` | collect_bam_output.cpp | bam_utils.c |
| `update_cand_var_profile_from_cons_aln_str` | collect_phase_noisy.cpp | collect_var.c |
| `update_cand_var_profile_from_cons_aln_str1` | collect_phase_noisy.cpp | collect_var.c |
| `update_cand_var_profile_from_cons_aln_str2` | collect_phase_noisy.cpp | collect_var.c |
| `update_cand_var_profile_from_cons_aln_str21` | collect_phase_noisy.cpp | collect_var.c |
| `update_digars_from_aln_str` | align.cpp | align.c |
| `update_digars_from_msa1` | align.cpp | align.c |
| `update_read_phase_set` | collect_phase.cpp | assign_hap.c |
| `update_read_var_profile_with_allele` | collect_phase_noisy.cpp | bam_utils.c |
| `update_var_hap_profile_based_on_read_hap` | collect_phase.cpp | assign_hap.c |
| `update_var_hap_profile_cons_alle_based_on_read_hap` | collect_phase.cpp | assign_hap.c |
| `update_var_hap_to_cons_alle` | collect_phase.cpp | assign_hap.c |
| `var_init_hap_profile_cons_allele` | collect_phase.cpp | assign_hap.c |
| `var_init_hap_to_alle_profile` | collect_phase.cpp | assign_hap.c |
| `var_is_homopolymer_indel` | collect_phase_noisy.cpp | collect_var.c |
| `var_noisy_reads_ratio` | collect_var.cpp | collect_var.c |
| `wfa_collect_aln_str` | align.cpp | align.c |
| `wfa_collect_noisy_aln_str_no_ps_hap` | align.cpp | align.c |
| `wfa_collect_noisy_aln_str_with_ps_hap` | align.cpp | align.c |
| `wfa_collect_pretty_alignment` | align.cpp | align.c |
| `wfa_end2end_aln` | align.cpp | align.c |
| `wfa_trim_aln_str` | align.cpp | align.c |

Shared names: **82**.

## Renamed in this pass

Four of ours were the same function under a different name; all four were
checked against the upstream signature before renaming, and the rename was
gated on byte-identical output over `chr20:5,000,000-9,000,000`.

| was | now | upstream |
|---|---|---|
| `init_assign_read_hap` | `init_assign_read_hap_based_on_cons_alle` | `assign_hap.c:151` |
| `update_var_hap_profile` | `update_var_hap_profile_based_on_read_hap` | `assign_hap.c:292` |
| `update_var_hap_profile_cons_alle` | `update_var_hap_profile_cons_alle_based_on_read_hap` | `assign_hap.c:270` |
| `check_agree_alleles` | **reverted** | see below |

## One rename rejected, and why it matters

`check_agree_alleles` looked like `check_agree_haps` (`assign_hap.c:307`) and was
renamed to it -- but our source **already had** a faithful `check_agree_haps`
with upstream's own signature including the `hap` argument, and the two are
deliberately different: `check_agree_haps` tests a read against a given
haplotype, while `check_agree_alleles` consults `chunk.haps[read_i]` and so
ignores unassigned reads. The rename silently turned them into C++ overloads,
which compiles and reads as one function. Reverted, and the header now states
the distinction so the next reader does not repeat it.

## Rejected candidate pairings

A sweep that mined each function's comment for an upstream name produced ten
candidates; six were the name of a function merely MENTIONED in the comment,
not the counterpart, and were rejected on signature:
`build_digars_and_events`/`has_equal_X_in_bam_cigar`,
`classify_cand_vars_pgphase`/`cr_overlap`,
`cr_extend_noisy_regs_with_low_comp`/`cr_merge`, `intervals_to_cr`/`cr_index`,
`passes_lcd_write_var_alt_ref_base_gate`/`write_var_to_vcf`, and
`phase_matrix_var_weight`/`read_to_cons_allele_score` (ours takes a variant and
returns a weight; upstream scores one read's allele against a consensus).
