# Stage trace

`PGPHASE_TRACE=DIR` makes `pgphase` write `DIR/trace.tsv`: one line for each decision that drops, merges, re-keys or
re-genotypes a candidate, a site or an output record, with a reason code. It is a diagnostic, off by default, and
writing it leaves every output unchanged. With it, a variant can be followed through the pipeline: pgphase-dev's
harness (`scripts/dev/harness`, `pg harness`) joins it to the truth and to the MSA dump.

Columns:

| column | meaning |
|---|---|
| `contig` | the contig being phased (`-` before one is) |
| `chunk` | the chunk (`RegionChunk::chunk_id`) whose worker made the decision; -1 for steps over a whole batch (the VCF records after chunks are merged). A site over a chunk edge is seen by both chunks: one line each |
| `pass` | the union EM's solve: 1 (sites only), 2 (with locus windows); 0 outside the EM. The second solve supersedes the first |
| `beg`, `end` | the reference span, 1-based inclusive: a record's from its VCF POS to its last REF base; a window's, a region's or a read's as named; a chunk boundary's, the two chunks' overlap (or the gap between them); a GAF site's, its position |
| `key` | `TYPE:REF_LEN:ALT` of the candidate's key, or `REF>ALT` from a catalog row's VCF metadata (a catalog row's key holds its site id, not a sequence), a catalog site's boundary walk (`gaf`), or `-` |
| `reason` | the reason code, one of the table's |
| `detail` | space-separated tokens, each `key=value` (no spaces in the value) or a bare lower-case flag (`msa_hom`, `flipped`); a candidate's genotype reads `gt=H1/H2 ps= cat=` (the alleles haplotypes 1 and 2 carry, -1 none), another candidate's the same with a prefix (`to_gt=`, `kept_gt=`) |

`scripts/check_trace.py` (part of `make check`) holds a run's trace to this: every line has the columns, every detail
parses, every (stage, event, reason) is in a table below.

An EM block is named by its first site's position, per chunk. Stitching renames blocks (`stitch rename`:
`merged_block` from the downstream block to the upstream one's id, `colliding_block_id` to a fresh id), and those
renames map every (chunk, block) to the phase set the VCF and BAM write.

`DIR/reads.tsv` holds the reads: one row per read as the phased BAM writes it -- `contig`, `beg` and `end` (the span
of the sites the read was called at, 1-based inclusive), `read`, `hap` (0 unphased), `ps` (-1 unphased), `copies`
(chunks that saw it), and `chunks`: per chunk, separated by ` | `, the label it had there
(`primary`, a `fill` / `fill_only` label applied after stitching, or `none`) and how it came about:

| part | meaning |
|---|---|
| `em hap= block= p= blocks= calls=` | the union EM's label (final solve): the block, the posterior, how many blocks its calls fall in |
| `unlabelled best_p= block= blocks= mixed_blocks= calls=` | no EM label: the best block's posterior (under the label posterior), blocks whose calls point both ways |
| `; deferred` / `; unverified` `hap= block=` | a last-resort label from its calls in one block all pointing one way, or from unverified calls |
| `; island kind= site= err=` | an island label: the site (E: an EM-rejected site, R: a repeat het indel, M: an island-only MSA het, W: a length window) and its error |
| `; from_next_chunk hap= ps=` | taken from the next chunk over a stitched boundary |
| `; fill rescue obs= hap= block= margin= votes=` / `; fill deferred_or_island` | the fill label applied: the rescue's (graph or alignment observations) or the deferred / island one |
| `; alignment_only filled=` | added from the alignment solve, with the calls realignment gave it |
| `; realigned= graded=` | calls realignment set at injected indels, and MSA calls it graded |
| `; snp_replaced=` | graph SNP calls the alignment's base replaced |
| `; aln hap= ps=` | its label in the alignment solve (k-means haplotype, 0 none, and phase set) |

| stage | event | reason | where | meaning |
|---|---|---|---|---|
| graph | drop | `pair_low_depth`, `pair_high_af`, `pair_low_af` | graph_bam_adapter.cpp `decompose_sites` | a catalog site's REF/ALT pair below the depth or outside the allele-fraction bounds (a homozygous pair is `pair_high_af`) |
| graph | merge | `duplicate_snarl_routed` | same | the pair duplicates another snarl's (minimal VCF identity): its reads go to that one |
| graph | regenotype | `repeat_het_indel_not_phased` | `apply_graph_noise_filter` | a catalog indel in a homopolymer or repeat demoted to RepeatHetIndel: never phased, never written |
| bam | drop | `pruned_category` | collect_var.cpp `prune_not_candidate_variants` | an alignment candidate deleted as LowCoverage, StrandBias or NonVariant (LowCoverage traced with 3+ ALT reads only) |
| msa | region | `resolved`, `no_consensus`, `too_long_for_msa`, `too_many_reads`, `no_reference`, `few_confident_reads`, `confident_reads_only` | collect_phase_noisy.cpp `collect_noisy_vars1` | a noisy region's MSA outcome (`confident_reads_only`: under half its reads at MAPQ 30 or more, the MSA takes those alone; `few_confident_reads`: a piece of a split region with fewer than 6 of them, skipped; `no_consensus` repeats while it is retried, with why: `too_few_full_reads N of M (min K)`, `no_full_reads`, `no_phased_reads`, `read_segment_too_long L`, `plegma_clusters N`, `plegma_empty_cluster`, `plegma_error: ...`, `abpoa_clusters N`, `no_reads`) |
| msa | region | `split` | collect_var.cpp `split_long_noisy_regions` | a noisy region longer than the MSA takes, cut into pieces between candidates (plegma only): `pieces`, the `last` piece |
| msa | reads | `hap`, `nohap`, `skipped` | align.cpp `collect_noisy_reg_aln_strs` | the reads one MSA call saw and which its clusters took: the branch (hap: the PS-aware path; nohap: clustered without a phase set; skipped: too few spanning reads), `reads`, `spanning` (both ends), `sampled`, the phase set and its reads per haplotype (`ps_hap1`, `ps_hap2`, `ps_unlabelled`), `cluster1`/`cluster2` sizes, `unclustered` and up to 20 `unclustered_reads` (name:cover, cover 3 spanning both ends); the clusters' own reads are in the MSA dump |
| msa | skip | `mismatch_beside_gap` | `make_cand_vars_from_baln0` | a consensus mismatch beside a gap, left to the gap as its anchor letter |
| msa | call | `admitted`, `replaced_existing` | `merge_var_profile` | an MSA call entering the chunk, alone or in place of a candidate with its key |
| msa | drop | `existing_candidate_kept`, `not_whitelisted` | same | an MSA call dropped for an existing candidate with its key, or absent from the site whitelist |
| msa | merge | `two_allele_insertion`, `two_allele_deletion`, `two_allele_insertion_deletion`, `two_allele_nested_deletion`, `two_allele_snp_deletion`, `two_allele_anchor` | `merge_msa_insertion_alleles`, `merge_msa_colocated_deletions`, `merge_msa_colocated_insertion_and_deletion`, `merge_msa_nested_deletions`, `merge_msa_snp_and_anchored_deletion`, `merge_msa_anchor_alleles` | two co-located MSA alleles merged into one two-allele record (a SNP with the other consensus's deletion anchored on its base: one record from that base) |
| refine | chunk | `msa_alignments` | `inject_alignment_private_sites`, `take_refined_alignments` | with `--refine-aln` and an aligned output, per alignment chunk: the reads whose alignment the MSA rewrote and the output takes (`refined`), the MSA left as it was (`unchanged`), or whose rewrite failed the check against the read's bases, moved its start too far, or met a hard clip or a mate (`rejected`, written as input) |
| msa | drop | `superseded_by_two_allele_record` | collect_var.cpp `drop_superseded_colocated_records` | a single-allele record describing an allele of a merged two-allele record |
| inject | decision | `injected`, `injected:adds_to_row`, `complement`, `replaces`, `promoted`, `hom`, `aln_hom`, `aln_hom:spelled`, `hom:few_confident_reads`, `hom:depth_gate`, `hom:graph_row`, `island:depth_gate`, `island:few_confident_reads`, `depth_gate`, `few_confident_reads`, `category`, `not_germline`, `outside_chunk` | union_inject.cpp `select_sites` | each alignment candidate's injection decision (as the `inj` dump); `aln_hom` a hom call of the alignment's (or the MSA's in a region whose consensus sequences agree) written 1/1 (`:spelled`: admitted because no read spells REF there), `injected:adds_to_row` a call that adds to a catalog row's allele rather than describing it again, written itself, `hom:*` one held back (too few or too small a share of confident reads, the depth gate, a graph row over it) |
| inject | drop | `msa_reads_as_reference`, `msa_lacks_change`, `msa_lacks_indel`, `aln_reads_lack_row`, `msa_edge_reads_lack_row` | `suppress_rows_msa_reads_as_reference`, `genotype_substitution_rows`, `unwrite_indel_rows_msa_lacks`, `unwrite_rows_reads_lack`, `unwrite_edge_rows_reads_lack` | a catalog row written no record: the MSA region reads it as reference, or no consensus carries its change (an indel: anywhere along its repeat tract); or, under a clean alignment het (or across an MSA region's edge), no alignment candidate describes it, fewer than 3 reads spell it and another non-reference sequence 3 reads and 20% of them do |
| inject | regenotype | `msa_homozygous` | `genotype_substitution_rows` | a catalog row written 1/1: both consensus sequences carry it |
| inject | drop | `row_replaced_by_msa_call` | `select_sites` | a catalog row an MSA het replaces (its REF cleared) |
| inject | merge / drop | `alias_folded` / `island_duplicate_dropped` | `fold_aliases` | a selection describing an allele another already does |
| inject | drop | `unusable_vcf_meta` | `synthesize_alignment_site_meta` | an alignment call whose VCF form does not round-trip: no record |
| em | not_site | `one_allele_called`, `minor_allele_below_fraction` | union_em.cpp `select_sites` | a germline candidate the EM does not phase |
| em | site | `phased`, `phased_written_unphased`, `unreliable`, `msa_het_unphasable`, `phased_msa_confirmed`, `phased_by_msa`, `phased_by_alignment_solve` | `write_site_phases`, `confirm_phase_by_msa`, `phase_unphasable_by_msa`, `phase_by_alignment_solve` | an EM site's outcome, with its block and learned error (`join_deep=`: why the site holds no join, `+`-separated: `mapq` poorly mapped, `depth` called by more reads than the chunk's confidently mapped sites, `dup` a duplicated haplotype) |
| em | msa_check | `agrees`, `disagrees` | `confirm_phase_by_msa` | an EM-phased site of a het MSA region against the region's other sites in its block: whether the EM's cis/trans relation is the MSA consensus sequences' (`agree=`/`disagree=` counts, learned error) |
| em | cut | `margin`, `few_spanning_reads`, `weak_bridge`, `hinge` | `cut_blocks` | a block cut between two EM sites and the rule that cut it (`delta`/`strong` margins, `spanning`/`strong_spanning` reads), then what they stand on: `crossing` reads with calls on both sides at any sites (`spanning` counts join-counted calls only), the reliable-reads vote (`reliable_same`/`reliable_other`), and the run of sites around the cut that hold no join (`held_sites`; of them `held_mapq`, `held_depth`, `held_dup` by `join_deep` reason, `held_excluded`) |
| em | join | `weak_bridge_joined`, `reliable_reads_joined` | `cut_blocks` | a join a weak bridge carries, or a weak-bridge or hinge cut undone by two or more crossing reads agreeing at sites of learned error under 0.05, with the same fields as a cut |
| em | locus_window | `phased`, `unreliable` | `write_site_phases` | a locus (`locus`) or catalog allele (`allele`) window of the second solve: an EM site with no VCF record of its own, with its learned error, block and reads |
| gaf | site | `walk_matches_no_allele`, `walk_matches_several_alleles` | graph_query.cpp `scan_indexed_gaf_chunk` | a catalog site some reads crossed (both boundary nodes on the walk) without matching exactly one allele walk: those reads give it no call (`matched=`, `ambiguous=`, `unmatched=` reads and a few names) |
| gaf | chunk | `alignments_below_mapq` | same | the chunk's GAF alignments under the MAPQ floor |
| aln | link | `linked`, `linked_flipped`, `break`, `first_or_no_link` | collect_phase.cpp `assign_hap_based_on_germline_het_vars_kmeans` | the alignment solve's phase-set builder, last pass: each het's link to the earlier het (`to=`), the spanning reads agreeing and conflicting, the phase set it took and its parity (`flags=` tells the solve's rounds apart) |
| reads | chunk | `alignment_only_reads_added` | union_reads.cpp `add_alignment_only_reads` | reads the graph chunk lacked, added from the alignment solve |
| reads | realign | `calls_rewritten` | `realign_indel_observations` | an injected site's reads re-called by realignment: calls set where the MSA left the read uncalled, MSA calls graded agreeing / undecided / contradicted |
| reads | fill | `alignment_only_reads_called` | `fill_missing_observations` | calls given to alignment-only reads at a site |
| reads | snp_call | `graph_call_vs_alignment`, `graph_call_vs_read_alignment`, `recounted` | `confirm_catalog_snp_calls` | a catalog SNP's graph calls the alignment's base replaced, and graph-only calls down-weighted; where the alignment has no SNP, each read's base from its own alignment record (reads absent, not read there, read; replaced); a site whose calls were replaced, its allele counts recounted from them (`before=`, `after=`), as the writer genotypes it |
| window | locus | `built`, `too_long`, `few_reads`, `no_split`, `homozygous`, `one_sided` | union_windows.cpp `build_locus_window_sites` | every noisy locus considered for a haplotype window (spanning reads, seeded reads per side, the two haplotypes' edit distance, reads per side) |
| window | allele | `built`, `collapsed_duplication`, `in_resolved_msa`, `too_long`, `overlaps_previous_window`, `observed_by_em`, `same_or_n_haplotypes`, `one_sided` | `build_allele_window_sites` | every allele window considered (key: `alignment`, `pair`, `catalog`), with its allele and reads per side |
| window | length | `built`, `too_long`, `observed_by_em`, `few_reads`, `one_length`, `one_sided` | `build_length_window_sites` | every repeat tract considered for an island-only length window (period, tract, the two lengths, the het score, reads per side) |
| stitch | boundary | `merged`, `merged_flipped`, `not_merged`, `no_overlap_reads`, `no_candidates` | collect_phase.cpp `flip_chunk_hap` | two adjacent chunks' stitch over their overlap: the votes of the block pair chosen (`same`, `crossed`), the blocks, block pairs linked, overlap reads |
| stitch | rename | `merged_block` | `flip_chunk_hap` | the downstream block of a merged pair takes the upstream block's id (`from=`, `to=`; `flipped` when its haplotypes were swapped) |
| stitch | rename | `colliding_block_id` | `separate_colliding_blocks` | a downstream block renamed off an upstream block's id it was not merged with (`from=`, `to=`) |
| em | window | `repeat_indel_window_phased` | union_phase.cpp `phase_chunk_with_alignment_sites` | a catalog repeat indel written with the genotype and block its allele window was phased with |
| em | site | `skewed_msa_het_unphased` | `write_skewed_msa_hets` | an MSA het its reads split too unevenly to phase, written unphased |
| output | rephase | `lone_het_rephased` | graph_collect.cpp `graph_chunks_to_candidate_table` | a lone het moved into the phase set holding most of its labelled reads, with the vote behind it (`h1=`/`h2=` reads on the winning allele / labelled reads per haplotype there, `ps_n` phase sets voting, `labelled` reads in all) |
| output | skip | `no_record`, `bam_alt_deletion_no_ref`, `injected_no_het_consensus`, `no_meta`, `allele_index_out_of_range`, `allele_empty_or_star`, `no_ref` | same | a row or allele the rebuild writes nothing for |
| output | record | `catalog`, `injected`, `mnp_split` | same | an output record built, with its genotype and phase set |
| output | skip | `injected_not_het_category` | same | an injected allele whose re-derived category is not a het's |
| output | drop | `duplicate_collapsed` | same | the losing copy of two records with one key (the kept one in the detail) |
| output | regenotype | `colocated_moved_to_other_hap`, `colocated_split_both_none` | collect_var.cpp `make_colocated_alleles_complementary` | a co-located record moved to the other haplotype |
| output | merge | `deletion_insertion_one_haplotype` | collect_var.cpp `merge_adjacent_deletion_insertion` | an insertion beside a deletion on one haplotype, written as one record with it |
| output | drop | `walk_contradicted` | collect_var.cpp `drop_contradicted_walks` | a catalog walk (deletion or replacement over 10 bases) a record inside it contradicts on one haplotype |
| output | drop | `conflict_contained`, `conflict_weaker` | `drop_conflicting_haplotype_alleles` | a record dropped for another on its haplotype at its POS |
| output | regenotype | `lone_locus_unphased` | `run_collect_graph_variation` | a phase set of one locus: its hets written unphased |
| writer | drop | `category_not_written`, `no_alt_genotype` | collect_output.cpp `project_vcf_candidates` | a record the VCF writer leaves out |
