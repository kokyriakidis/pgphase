#ifndef PGPHASE_TEST_UPSTREAM_PHASE_BRIDGE_H
#define PGPHASE_TEST_UPSTREAM_PHASE_BRIDGE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Wrappers call original longcallD source compiled into the test binary. */
int upstream_msa_homopolymer(const char *ref, int pos_offset, int is_insertion,
                            int ref_len, const char *alt);
void upstream_var_init(int is_ont, int homopolymer, int category, int n_alleles,
                       const int *coverage, int had_profile, int *consensus,
                       int profile[3][3]);
void upstream_update_cons(int is_ont, int homopolymer, int category, int n_alleles,
                          int hap, int profile[3][3], int *consensus);
int upstream_read_score(int hap, int category, int allele, int *consensus);
void upstream_full_kmeans(int n_vars, int n_reads,
                          const int read_alleles[4][3], int is_ont,
                          int consensus[3][3], int64_t *candidate_phase_sets,
                          int *read_haps, int64_t *read_phase_sets);
int upstream_phase_link(int n_vars, int n_reads, const int *read_haps,
                        const int read_alleles[4][3], int consensus[3][3],
                        int *phase_sets);
int upstream_init_read_hap(int n_vars, const int *categories, const int *is_snp,
                           const int *homopolymer, const int *alleles,
                           int consensus[3][3], int target_category,
                           int *clean_agree, int *clean_conflict);

#ifdef __cplusplus
}
#endif

#endif
