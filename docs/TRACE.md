# Stage trace

`PGPHASE_TRACE=DIR` makes `pgphase` write `DIR/trace.tsv`: one line for each decision that drops, merges, re-keys or
re-genotypes a candidate, a site or an output record, with a reason code. It is a diagnostic, off by default, and
writing it leaves every output unchanged. With it, a variant can be followed through the pipeline: pgphase-dev's
harness (`scripts/dev/harness`, `pg harness`) joins it to the truth and to the MSA dump.

Columns: `stage`, `event`, `beg`, `end` (the reference span, 1-based inclusive, from the record's VCF POS to its last
REF base), `key` (`TYPE:REF_LEN:ALT` of the candidate's key, or `REF>ALT` from a catalog row's VCF metadata: a catalog
row's key holds its site id, not a sequence), `reason`, `detail`.

The union EM runs twice per chunk: a later `em` event for a candidate supersedes an earlier one.

| stage | event | reason | where | meaning |
|---|---|---|---|---|
| graph | drop | `pair_low_depth`, `pair_high_af`, `pair_low_af` | graph_bam_adapter.cpp `decompose_sites` | a catalog site's REF/ALT pair below the depth or outside the allele-fraction bounds (a homozygous pair is `pair_high_af`) |
| graph | merge | `duplicate_snarl_routed` | same | the pair duplicates another snarl's (minimal VCF identity): its reads go to that one |
| graph | regenotype | `repeat_het_indel_not_phased` | `apply_graph_noise_filter` | a catalog indel in a homopolymer or repeat demoted to RepeatHetIndel: never phased, never written |
| bam | drop | `pruned_category` | collect_var.cpp `prune_not_candidate_variants` | an alignment candidate deleted as LowCoverage, StrandBias or NonVariant (LowCoverage traced with 3+ ALT reads only) |
| msa | region | `resolved`, `no_consensus`, `too_long_for_msa`, `too_many_reads`, `no_reference` | collect_phase_noisy.cpp `collect_noisy_vars1` | a noisy region's MSA outcome (`no_consensus` repeats while it is retried) |
| msa | skip | `mismatch_beside_gap` | `make_cand_vars_from_baln0` | a consensus mismatch beside a gap, left to the gap as its anchor letter |
| msa | call | `admitted`, `replaced_existing` | `merge_var_profile` | an MSA call entering the chunk, alone or in place of a candidate with its key |
| msa | drop | `existing_candidate_kept` | same | an MSA call dropped for an existing candidate with its key |
| msa | merge | `two_allele_insertion`, `two_allele_deletion` | `merge_msa_insertion_alleles`, `merge_msa_colocated_deletions` | two co-located MSA alleles merged into one two-allele record |
| msa | drop | `superseded_by_two_allele_record` | collect_var.cpp `drop_superseded_colocated_records` | a single-allele record describing an allele of a merged two-allele record |
| inject | decision | `injected`, `complement`, `replaces`, `promoted`, `hom`, `island:*`, `depth_gate`, `few_confident_reads`, `category`, `not_germline`, `outside_chunk` | union_inject.cpp `select_sites` | each alignment candidate's injection decision (as the `inj` dump) |
| inject | drop | `msa_reads_as_reference`, `msa_lacks_change` | `suppress_rows_msa_reads_as_reference`, `genotype_substitution_rows` | a catalog row written no record: the MSA region reads it as reference, or no consensus carries its change |
| inject | regenotype | `msa_homozygous` | `genotype_substitution_rows` | a catalog row written 1/1: both consensus sequences carry it |
| inject | drop | `row_replaced_by_msa_call` | `select_sites` | a catalog row an MSA het replaces (its REF cleared) |
| inject | merge / drop | `alias_folded` / `island_duplicate_dropped` | `fold_aliases` | a selection describing an allele another already does |
| inject | drop | `unusable_vcf_meta` | `synthesize_alignment_site_meta` | an alignment call whose VCF form does not round-trip: no record |
| em | not_site | `one_allele_called`, `minor_allele_below_fraction` | union_em.cpp `select_sites` | a germline candidate the EM does not phase |
| em | site | `phased`, `phased_written_unphased`, `unreliable`, `msa_het_unphasable` | `write_site_phases` | an EM site's outcome, with its block and learned error |
| em | site | `skewed_msa_het_unphased` | `write_skewed_msa_hets` | an MSA het its reads split too unevenly to phase, written unphased |
| output | rephase | `lone_het_rephased` | graph_collect.cpp `graph_chunks_to_candidate_table` | a lone het moved into the phase set holding most of its labelled reads |
| output | skip | `no_record`, `bam_alt_deletion_no_ref`, `injected_no_het_consensus`, `no_meta`, `allele_index_out_of_range`, `allele_empty_or_star`, `no_ref` | same | a row or allele the rebuild writes nothing for |
| output | record | `catalog`, `injected`, `mnp_split` | same | an output record built, with its genotype and phase set |
| output | skip | `injected_not_het_category` | same | an injected allele whose re-derived category is not a het's |
| output | drop | `duplicate_collapsed` | same | the losing copy of two records with one key (the kept one in the detail) |
| output | regenotype | `colocated_moved_to_other_hap`, `colocated_split_both_none` | collect_var.cpp `make_colocated_alleles_complementary` | a co-located record moved to the other haplotype |
| output | drop | `conflict_contained`, `conflict_weaker` | `drop_conflicting_haplotype_alleles` | a record dropped for another on its haplotype at its POS |
| output | regenotype | `lone_locus_unphased` | `run_collect_graph_variation` | a phase set of one locus: its hets written unphased |
| writer | drop | `category_not_written`, `no_alt_genotype` | collect_output.cpp `project_vcf_candidates` | a record the VCF writer leaves out |
