#include "test_upstream_phase_bridge.h"

#include "call_var_main.h"
#include "collect_var.h"
#include "bam_utils.h"
#include "assign_hap.h"

#include <stdlib.h>
#include <string.h>

void var_init_hap_profile_cons_allele(const call_var_opt_t *, cand_var_t *, int *,
                                     int, int *);
int update_var_hap_to_cons_alle(const call_var_opt_t *, cand_var_t *, int, int);
int read_to_cons_allele_score(int, int, cand_var_t *, int, int);
int init_assign_read_hap_based_on_cons_alle(bam_chunk_t *, int, cand_var_t *,
                                            read_var_profile_t *, int *, int);
int iter_update_var_hap_cons_phase_set(bam_chunk_t *, int *, read_var_profile_t *,
                                       cand_var_t *, int, int *);
int var_is_homopolymer_indel(bam_chunk_t *, hts_pos_t, int, int, int, uint8_t *);
int LONGCALLD_VERBOSE = 0;

int upstream_msa_homopolymer(const char *ref, int pos_offset, int is_insertion,
                            int ref_len, const char *alt) {
    bam_chunk_t chunk = {0};
    chunk.ref_beg = 100;
    chunk.ref_seq = (char *)ref;
    uint8_t alt_codes[16] = {0};
    const int alt_len = (int)strlen(alt);
    if (alt_len > 16) return -1;
    static const char bases[] = "ACGT";
    for (int i = 0; i < alt_len; ++i) {
        const char *base = strchr(bases, alt[i]);
        if (base == NULL) return -1;
        alt_codes[i] = (uint8_t)(base - bases);
    }
    return var_is_homopolymer_indel(&chunk, chunk.ref_beg + pos_offset,
                                    is_insertion ? BAM_CINS : BAM_CDEL,
                                    ref_len, alt_len, alt_codes);
}

void upstream_var_init(int is_ont, int homopolymer, int category, int n_alleles,
                       const int *coverage, int had_profile, int *consensus,
                       int profile[3][3]) {
    call_var_opt_t opt = {0};
    opt.is_ont = is_ont;
    cand_var_t var = {0};
    var.is_homopolymer_indel = homopolymer;
    var.n_uniq_alles = n_alleles;
    var.alle_covs = (int *)coverage;
    if (had_profile) {
        var.hap_to_alle_profile = malloc(3 * sizeof(int *));
        for (int h = 0; h < 3; ++h) {
            var.hap_to_alle_profile[h] = calloc((size_t)n_alleles, sizeof(int));
            memcpy(var.hap_to_alle_profile[h], profile[h],
                   (size_t)n_alleles * sizeof(int));
        }
        var.hap_to_cons_alle = malloc(3 * sizeof(int));
        memcpy(var.hap_to_cons_alle, consensus, 3 * sizeof(int));
    }
    int index = 0;
    var_init_hap_profile_cons_allele(&opt, &var, &index, 1, &category);
    memcpy(consensus, var.hap_to_cons_alle, 3 * sizeof(int));
    for (int h = 0; h < 3; ++h) {
        memcpy(profile[h], var.hap_to_alle_profile[h],
               (size_t)n_alleles * sizeof(int));
        free(var.hap_to_alle_profile[h]);
    }
    free(var.hap_to_alle_profile);
    free(var.hap_to_cons_alle);
}

void upstream_update_cons(int is_ont, int homopolymer, int category, int n_alleles,
                          int hap, int profile[3][3], int *consensus) {
    call_var_opt_t opt = {0};
    opt.is_ont = is_ont;
    cand_var_t var = {0};
    int *rows[3] = {profile[0], profile[1], profile[2]};
    var.is_homopolymer_indel = homopolymer;
    var.n_uniq_alles = n_alleles;
    var.hap_to_alle_profile = rows;
    var.hap_to_cons_alle = consensus;
    update_var_hap_to_cons_alle(&opt, &var, category, hap);
}

int upstream_read_score(int hap, int category, int allele, int *consensus) {
    cand_var_t var = {0};
    var.hap_to_cons_alle = consensus;
    return read_to_cons_allele_score(0, hap, &var, category, allele);
}


int upstream_init_read_hap(int n_vars, const int *categories, const int *is_snp,
                           const int *homopolymer, const int *alleles,
                           int consensus[3][3], int target_category,
                           int *clean_agree, int *clean_conflict) {
    cand_var_t vars[3] = {{0}};
    int category_copy[3] = {0};
    for (int i = 0; i < n_vars; ++i) {
        vars[i].var_type = is_snp[i] ? BAM_CDIFF : BAM_CINS;
        vars[i].is_homopolymer_indel = homopolymer[i];
        vars[i].hap_to_cons_alle = consensus[i];
        category_copy[i] = categories[i];
    }
    read_var_profile_t profile = {0};
    profile.start_var_idx = 0;
    profile.end_var_idx = n_vars - 1;
    profile.alleles = (int *)alleles;
    bam_chunk_t chunk = {0};
    chunk.n_clean_agree_snps = clean_agree;
    chunk.n_clean_conflict_snps = clean_conflict;
    return init_assign_read_hap_based_on_cons_alle(
        &chunk, 0, vars, &profile, category_copy, target_category);
}


int upstream_phase_link(int n_vars, int n_reads, const int *read_haps,
                        const int read_alleles[4][3], int consensus[3][3],
                        int *phase_sets) {
    cand_var_t vars[3] = {{0}};
    int categories[3] = {LONGCALLD_CLEAN_HET_SNP,
                         LONGCALLD_CLEAN_HET_SNP,
                         LONGCALLD_CLEAN_HET_SNP};
    int var_index[3] = {0, 1, 2};
    for (int i = 0; i < n_vars; ++i) {
        vars[i].pos = 100 * (i + 1);
        vars[i].var_type = BAM_CDIFF;
        vars[i].hap_to_cons_alle = consensus[i];
    }
    bam_chunk_t chunk = {0};
    chunk.n_reads = n_reads;
    int haps[4] = {0};
    uint8_t skipped[4] = {0};
    int ordered_ids[4] = {0, 1, 2, 3};
    read_var_profile_t profiles[4] = {{0}};
    chunk.haps = haps;
    chunk.is_skipped = skipped;
    chunk.ordered_read_ids = ordered_ids;
    chunk.read_var_cr = cr_init();
    for (int i = 0; i < n_reads; ++i) {
        haps[i] = read_haps[i];
        profiles[i].start_var_idx = 0;
        profiles[i].end_var_idx = n_vars - 1;
        profiles[i].alleles = (int *)read_alleles[i];
        cr_add(chunk.read_var_cr, "cr", 0, n_vars, i);
    }
    cr_index(chunk.read_var_cr);
    const int changed = iter_update_var_hap_cons_phase_set(
        &chunk, var_index, profiles, vars, n_vars, categories);
    for (int i = 0; i < n_vars; ++i)
        phase_sets[i] = (int)vars[i].phase_set;
    cr_destroy(chunk.read_var_cr);
    return changed;
}


void upstream_full_kmeans(int n_vars, int n_reads,
                          const int read_alleles[4][3], int is_ont,
                          int consensus[3][3], int64_t *candidate_phase_sets,
                          int *read_haps, int64_t *read_phase_sets) {
    call_var_opt_t opt = {0};
    opt.is_ont = is_ont;
    cand_var_t vars[3] = {{0}};
    int categories[3] = {LONGCALLD_CLEAN_HET_SNP,
                         LONGCALLD_CLEAN_HET_SNP,
                         LONGCALLD_CLEAN_HET_SNP};
    int coverage[3][2] = {{0}};
    for (int vi = 0; vi < n_vars; ++vi) {
        vars[vi].pos = 100 * (vi + 1);
        vars[vi].var_type = BAM_CDIFF;
        vars[vi].n_uniq_alles = 2;
        vars[vi].alle_covs = coverage[vi];
        vars[vi].total_cov = n_reads;
        for (int ri = 0; ri < n_reads; ++ri)
            ++coverage[vi][read_alleles[ri][vi]];
    }
    bam_chunk_t chunk = {0};
    chunk.n_reads = n_reads;
    chunk.n_cand_vars = n_vars;
    chunk.cand_vars = vars;
    chunk.var_i_to_cate = categories;
    int phase_scores[4] = {0};
    int clean_agree[4] = {0}, clean_conflict[4] = {0};
    uint8_t skipped[4] = {0};
    int ordered_ids[4] = {0, 1, 2, 3};
    read_var_profile_t profiles[4] = {{0}};
    chunk.haps = read_haps;
    chunk.phase_sets = read_phase_sets;
    chunk.phase_scores = phase_scores;
    chunk.n_clean_agree_snps = clean_agree;
    chunk.n_clean_conflict_snps = clean_conflict;
    chunk.is_skipped = skipped;
    chunk.ordered_read_ids = ordered_ids;
    chunk.read_var_profile = profiles;
    chunk.read_var_cr = cr_init();
    for (int ri = 0; ri < n_reads; ++ri) {
        profiles[ri].read_id = ri;
        profiles[ri].start_var_idx = 0;
        profiles[ri].end_var_idx = n_vars - 1;
        profiles[ri].alleles = (int *)read_alleles[ri];
        cr_add(chunk.read_var_cr, "cr", 0, n_vars, ri);
    }
    cr_index(chunk.read_var_cr);
    assign_hap_based_on_germline_het_vars_kmeans(
        &opt, &chunk, LONGCALLD_CAND_GERMLINE_CLEAN_VAR_CATE);
    for (int vi = 0; vi < n_vars; ++vi) {
        memcpy(consensus[vi], vars[vi].hap_to_cons_alle, 3 * sizeof(int));
        candidate_phase_sets[vi] = vars[vi].phase_set;
        for (int h = 0; h < 3; ++h)
            free(vars[vi].hap_to_alle_profile[h]);
        free(vars[vi].hap_to_alle_profile);
        free(vars[vi].hap_to_cons_alle);
    }
    cr_destroy(chunk.read_var_cr);
}
