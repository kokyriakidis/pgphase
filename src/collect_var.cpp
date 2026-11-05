/**
 * @file collect_var.cpp
 * Candidate sites, allele counting, noisy regions, and classification.
 */

#include "trace.hpp"
#include "collect_var.hpp"
#include "collect_phase.hpp"
#include "collect_phase_noisy.hpp"
#include "fisher_exact.hpp"
#include "noise_filter.hpp"
#include "plegma_msa.hpp"

#include "sdust.h"

#include <algorithm>
#include <cctype>
#include <set>
#include <map>
#include <array>
#include <cfloat>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

namespace pgphase_collect {

namespace {

// nt4 encoding table: maps ASCII bases to 0-3 (ACGT), 4 for ambiguous.
// `exact_comp_var_site` / `exact_comp_var_site_ins` (`collect_var.c`).
alignas(64) static constexpr uint8_t kLcdNstNt4Table[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

// Compare ALT sequences using nt4-encoded bytes for deterministic ordering.
int cmp_alt_seq_nt4(const std::string& a, const std::string& b, int len) {
    for (int i = 0; i < len; ++i) {
        const uint8_t ca = kLcdNstNt4Table[static_cast<unsigned char>(a[static_cast<size_t>(i)])];
        const uint8_t cb = kLcdNstNt4Table[static_cast<unsigned char>(b[static_cast<size_t>(i)])];
        if (ca != cb) return (ca < cb) ? -1 : 1;
    }
    return 0;
}

static int cmp_alt_bytes_nt4(const uint8_t* a, const uint8_t* b, int len) {
    for (int i = 0; i < len; ++i) {
        const uint8_t ca = kLcdNstNt4Table[static_cast<unsigned char>(a[i])];
        const uint8_t cb = kLcdNstNt4Table[static_cast<unsigned char>(b[i])];
        if (ca != cb) return (ca < cb) ? -1 : 1;
    }
    return 0;
}

// Compact variant-site representation for merge/dedup comparisons.
struct LcdVarSiteMerge {
    hts_pos_t pos = 0;
    int var_type = 0;
    int ref_len = 0;
    int alt_len = 0;
    const uint8_t* alt_seq = nullptr;
};

static LcdVarSiteMerge merge_site_from_candidate_key(const VariantKey& k) {
    LcdVarSiteMerge s{};
    s.pos = k.pos;
    if (k.type == VariantType::Snp) {
        s.var_type = BAM_CDIFF;
        s.ref_len = 1;
        s.alt_len = static_cast<int>(k.alt.size());
        s.alt_seq = reinterpret_cast<const uint8_t*>(k.alt.data());
    } else if (k.type == VariantType::Insertion) {
        s.var_type = BAM_CINS;
        s.ref_len = 0;
        s.alt_len = static_cast<int>(k.alt.size());
        s.alt_seq = reinterpret_cast<const uint8_t*>(k.alt.data());
    } else {
        s.var_type = BAM_CDEL;
        s.ref_len = k.ref_len;
        s.alt_len = 0;
        s.alt_seq = nullptr;
    }
    return s;
}

static int digar_to_bam_var_type(DigarType t) {
    switch (t) {
        case DigarType::Snp:
            return BAM_CDIFF;
        case DigarType::Insertion:
            return BAM_CINS;
        case DigarType::Deletion:
            return BAM_CDEL;
        case DigarType::SoftClip:
            return BAM_CSOFT_CLIP;
        case DigarType::HardClip:
            return BAM_CHARD_CLIP;
        case DigarType::RefSkip:
            return BAM_CREF_SKIP;
        default:
            return BAM_CEQUAL;
    }
}

// Build a merge-site from a digar op for fuzzy insertion comparison.
static LcdVarSiteMerge merge_site_from_digar_op(const DigarOp& op) {
    LcdVarSiteMerge s{};
    s.pos = op.pos;
    switch (op.type) {
        case DigarType::Snp:
            s.var_type = BAM_CDIFF;
            s.ref_len = 1;
            s.alt_len = static_cast<int>(op.alt.size());
            s.alt_seq = reinterpret_cast<const uint8_t*>(op.alt.data());
            return s;
        case DigarType::Insertion:
            s.var_type = BAM_CINS;
            s.ref_len = 0;
            s.alt_len = static_cast<int>(op.alt.size());
            s.alt_seq = reinterpret_cast<const uint8_t*>(op.alt.data());
            return s;
        case DigarType::Deletion:
            s.var_type = BAM_CDEL;
            s.ref_len = op.len;
            s.alt_len = 0;
            s.alt_seq = nullptr;
            return s;
        default:
            s.var_type = digar_to_bam_var_type(op.type);
            s.ref_len = 1;
            s.alt_len = op.len;
            s.alt_seq = nullptr;
            return s;
    }
}

// Fuzzy insertion comparison on merge-shaped sites (no tid).
static int exact_comp_var_site_ins_merge(const LcdVarSiteMerge* var1, const LcdVarSiteMerge* var2, int min_sv_len) {
    const hts_pos_t var1_pos = var1->var_type == BAM_CDIFF ? var1->pos : var1->pos - 1;
    const hts_pos_t var2_pos = var2->var_type == BAM_CDIFF ? var2->pos : var2->pos - 1;
    if (var1_pos < var2_pos) return -1;
    if (var1_pos > var2_pos) return 1;
    if (var1->var_type < var2->var_type) return -1;
    if (var1->var_type > var2->var_type) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    if (var1->var_type == BAM_CDIFF) {
        if (var1->alt_len < var2->alt_len) return -1;
        if (var1->alt_len > var2->alt_len) return 1;
        return cmp_alt_bytes_nt4(var1->alt_seq, var2->alt_seq, var1->alt_len);
    }
    if (var1->var_type == BAM_CINS) {
        if (var1->alt_len < min_sv_len) {
            if (var1->alt_len < var2->alt_len) return -1;
            if (var1->alt_len > var2->alt_len) return 1;
            return cmp_alt_bytes_nt4(var1->alt_seq, var2->alt_seq, var1->alt_len);
        }
        const int min_len = var1->alt_len < var2->alt_len ? var1->alt_len : var2->alt_len;
        const int max_len = var1->alt_len > var2->alt_len ? var1->alt_len : var2->alt_len;
        if (static_cast<double>(min_len) >= static_cast<double>(max_len) * 0.8) return 0;
        return var1->alt_len - var2->alt_len;
    }
    return 0;
}

} // namespace

// ════════════════════════════════════════════════════════════════════════════
// Variant site comparison
// ════════════════════════════════════════════════════════════════════════════

/**
 * Effective alternate length for sorting and comparison.
 *
 * Deletions use length 0; SNPs and insertions use `alt.size()`.
 * @param v Candidate key.
 * @return Alternate length used in comparators.
 */
static int var_site_alt_len(const VariantKey& v) {
    if (v.type == VariantType::Deletion) return 0;
    return static_cast<int>(v.alt.size());
}

// VCF POS/REF/ALT to the internal candidate key.
VariantKey vcf_to_variant_key(int tid, hts_pos_t vcf_pos,
                              const std::string& vcf_ref,
                              const std::string& vcf_alt) {
    VariantKey key;
    key.tid = tid;
    size_t substitution_prefix = 0;
    size_t substitution_suffix = 0;
    if (vcf_ref.size() == vcf_alt.size()) {
        while (substitution_prefix < vcf_ref.size() &&
               vcf_ref[substitution_prefix] == vcf_alt[substitution_prefix])
            ++substitution_prefix;
        while (substitution_suffix + substitution_prefix < vcf_ref.size() &&
               vcf_ref[vcf_ref.size() - substitution_suffix - 1] ==
                   vcf_alt[vcf_alt.size() - substitution_suffix - 1])
            ++substitution_suffix;
    }
    const size_t substitution_size =
        vcf_ref.size() - substitution_prefix - substitution_suffix;
    if (vcf_ref.size() == vcf_alt.size() && substitution_size == 1) {
        key.type = VariantType::Snp;
        key.pos = vcf_pos + static_cast<hts_pos_t>(substitution_prefix);
        key.ref_len = 1;
        key.alt = vcf_alt.substr(substitution_prefix, 1);
    } else {
        // Match BAM's placement first: consume the full common prefix, even
        // inside a repeat. Remaining common suffix bases are catalog context,
        // not part of the edit. Leaving them in the key makes a padded pure
        // graph indel look like a replacement and blocks its exact BAM match.
        // Trimming the suffix before the prefix would instead shift repeat
        // indels away from the existing BAM/MSA convention.
        size_t shared = 0;
        const size_t max_shared = std::min(vcf_ref.size(), vcf_alt.size());
        while (shared < max_shared && vcf_ref[shared] == vcf_alt[shared])
            ++shared;
        size_t ref_end = vcf_ref.size();
        size_t alt_end = vcf_alt.size();
        if (vcf_ref.size() != vcf_alt.size()) {
            while (ref_end > shared && alt_end > shared &&
                   vcf_ref[ref_end - 1] == vcf_alt[alt_end - 1]) {
                --ref_end;
                --alt_end;
            }
        }
        key.type = vcf_alt.size() > vcf_ref.size()
            ? VariantType::Insertion : VariantType::Deletion;
        key.pos = vcf_pos + static_cast<hts_pos_t>(shared);
        key.ref_len = static_cast<int>(ref_end - shared);
        key.alt = vcf_alt.substr(shared, alt_end - shared);
    }
    return key;
}

/**
 * Total order for candidate keys.
 *
 * Compares tid, sort position, type, `ref_len`, and alternate sequence length / bytes.
 * Return \<0 if \a var1 \< \a var2, 0 if equal, \>0 if \a var1 \> \a var2.
 *
 * @param var1 First variant key.
 * @param var2 Second variant key.
 * @return \<0, 0, or \>0 according to total order.
 */
int exact_comp_var_site(const VariantKey* var1, const VariantKey* var2) {
    if (var1->tid != var2->tid) return var1->tid < var2->tid ? -1 : 1;
    const hts_pos_t p1 = var1->type == VariantType::Snp ? var1->pos : var1->pos - 1;
    const hts_pos_t p2 = var2->type == VariantType::Snp ? var2->pos : var2->pos - 1;
    if (p1 < p2) return -1;
    if (p1 > p2) return 1;
    const int t1 = static_cast<int>(var1->type);
    const int t2 = static_cast<int>(var2->type);
    if (t1 < t2) return -1;
    if (t1 > t2) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    const int a1 = var_site_alt_len(*var1);
    const int a2 = var_site_alt_len(*var2);
    if (a1 < a2) return -1;
    if (a1 > a2) return 1;
    if (var1->type == VariantType::Snp || var1->type == VariantType::Insertion) {
        return cmp_alt_seq_nt4(var1->alt, var2->alt, a1);
    }
    return 0;
}

int exact_comp_cand_var(const CandidateVariant* var1, const CandidateVariant* var2) {
    return exact_comp_var_site(&var1->key, &var2->key);
}

/**
 * Variant comparison with large-insertion fuzzy merge.
 *
 * Same ordering keys as `exact_comp_var_site` for SNPs and deletions. For insertions with
 * `alt` length ≥ \a min_sv_len, treats two alleles as equal when
 * `min(len)/max(len) ≥ 0.8`; shorter insertions still require exact sequence match.
 *
 * @param var1 First variant key.
 * @param var2 Second variant key.
 * @param min_sv_len Minimum insertion length to apply the length-ratio rule.
 * @return \<0, 0, or \>0; zero means merge-equivalent under fuzzy INS rule.
 */
int exact_comp_var_site_ins(const VariantKey* var1, const VariantKey* var2, int min_sv_len) {
    if (var1->tid != var2->tid) return var1->tid < var2->tid ? -1 : 1;
    const hts_pos_t p1 = var1->type == VariantType::Snp ? var1->pos : var1->pos - 1;
    const hts_pos_t p2 = var2->type == VariantType::Snp ? var2->pos : var2->pos - 1;
    if (p1 < p2) return -1;
    if (p1 > p2) return 1;
    const int t1 = static_cast<int>(var1->type);
    const int t2 = static_cast<int>(var2->type);
    if (t1 < t2) return -1;
    if (t1 > t2) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    if (var1->type == VariantType::Snp) {
        const int a1 = var_site_alt_len(*var1);
        const int a2 = var_site_alt_len(*var2);
        if (a1 < a2) return -1;
        if (a1 > a2) return 1;
        return cmp_alt_seq_nt4(var1->alt, var2->alt, a1);
    }
    if (var1->type == VariantType::Insertion) {
        const int a1 = var_site_alt_len(*var1);
        const int a2 = var_site_alt_len(*var2);
        if (a1 < min_sv_len) {
            if (a1 < a2) return -1;
            if (a1 > a2) return 1;
            return cmp_alt_seq_nt4(var1->alt, var2->alt, a1);
        }
        const int min_len = a1 < a2 ? a1 : a2;
        const int max_len = a1 > a2 ? a1 : a2;
        if (static_cast<double>(min_len) >= static_cast<double>(max_len) * 0.8) return 0;
        return a1 - a2;
    }
    return 0;
}

/**
 * @brief Strict-weak ordering for `VariantKey` using `exact_comp_var_site`.
 * @return True if `lhs` is ordered before `rhs`.
 */
bool VariantKeyLess::operator()(const VariantKey& lhs, const VariantKey& rhs) const {
    return exact_comp_var_site(&lhs, &rhs) < 0;
}

// ════════════════════════════════════════════════════════════════════════════
// Candidate site collection
// ════════════════════════════════════════════════════════════════════════════

/**
 * Build a VariantKey from one digar op.
 *
 * SNP: `ref_len` 1; insertion: `ref_len` 0; deletion: `ref_len` = `op.len`, empty `alt`.
 *
 * @param tid Contig index for the read.
 * @param op Parsed SNP/indel operation.
 * @return Normalized `VariantKey` for the candidate table.
 */
VariantKey variant_key_from_digar(int tid, const DigarOp& op) {
    if (op.type == DigarType::Snp) {
        return VariantKey{tid, op.pos, VariantType::Snp, 1, op.alt};
    }
    if (op.type == DigarType::Insertion) {
        return VariantKey{tid, op.pos, VariantType::Insertion, 0, op.alt};
    }
    return VariantKey{tid, op.pos, VariantType::Deletion, op.len, ""};
}

/**
 * @brief Whether a digar op is a non–low-quality SNP/indel inside the region bounds.
 *
 * Use reg_beg or reg_end == -1 to disable that bound.
 *
 * @param digar Parsed operation.
 * @param reg_beg Minimum 1-based position, or -1 for no lower bound.
 * @param reg_end Maximum 1-based position, or -1 for no upper bound.
 * @return True if SNP/indel, not low-quality, and inside bounds.
 */
static bool is_collectible_var_digar(const DigarOp& digar, hts_pos_t reg_beg, hts_pos_t reg_end) {
    const hts_pos_t digar_pos = digar.pos;
    if (reg_beg != -1 && digar_pos < reg_beg) return false;
    if (reg_end != -1 && digar_pos > reg_end) return false;
    if (digar.low_quality) return false;
    return digar.type == DigarType::Snp || digar.type == DigarType::Insertion ||
           digar.type == DigarType::Deletion;
}

/**
 * Lower-bound index into sorted candidates for a read span.
 * @param v Sorted candidates by `sort_pos` / `pos`.
 * @param start Read start (1-based).
 * @param n_total `static_cast<int>(v.size())`.
 * @return First index to consider in the merge walk.
 */
int get_var_site_start(const CandidateTable& v, hts_pos_t start, int n_total) {
    hts_pos_t target = start > 0 ? start - 1 : start;
    int left = 0;
    int right = n_total;
    while (left < right) {
        const int mid = left + (right - left) / 2;
        const hts_pos_t mid_pos = v[static_cast<size_t>(mid)].key.sort_pos();
        if (mid_pos < target) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    while (left < n_total && v[static_cast<size_t>(left)].key.pos < start) {
        ++left;
    }
    return left;
}

/**
 * Deduplicates sorted candidates using fuzzy insertion equivalence.
 *
 * Sorts by `exact_comp_var_site`, then keeps one row per equivalence class where fuzzy large
 * insertions may collapse.
 *
 * @param variants In/out table; replaced with unique keys in sort order.
 */
void collapse_fuzzy_large_insertions(CandidateTable& variants, int min_sv_len) {
    std::sort(variants.begin(), variants.end(), [](const CandidateVariant& a, const CandidateVariant& b) {
        return exact_comp_var_site(&a.key, &b.key) < 0;
    });
    if (variants.empty()) return;
    size_t write_i = 1;
    for (size_t read_i = 1; read_i < variants.size(); ++read_i) {
        if (exact_comp_var_site_ins(
                &variants[write_i - 1].key, &variants[read_i].key, min_sv_len) == 0) {
            continue;
        }
        if (write_i != read_i) variants[write_i] = std::move(variants[read_i]);
        ++write_i;
    }
    variants.resize(write_i);
}

/**
 * @brief Builds the raw candidate table from per-read digars for one chunk.
 *
 * Emits one `CandidateVariant` per collectible SNP/indel overlapping `chunk`, then runs
 * collapse_fuzzy_large_insertions.
 * `collect_all_cand_var_sites`.
 *
 * @param chunk Genomic interval (1-based inclusive bounds).
 * @param reads Parsed alignments with digar ops.
 * @param variants Cleared and filled with deduplicated candidate keys.
 */
void collect_candidate_sites_from_records(const RegionChunk& chunk,
                                          const std::vector<ReadRecord>& reads,
                                          CandidateTable& variants,
                                          int min_sv_len) {
    variants.clear();
    size_t max_var_sites = 0;
    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        for (const DigarOp& digar : read.digars) {
            if (!is_collectible_var_digar(digar, chunk.beg, chunk.end)) continue;
            ++max_var_sites;
        }
    }
    if (max_var_sites == 0) return;

    variants.reserve(max_var_sites);
    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        for (const DigarOp& digar : read.digars) {
            if (!is_collectible_var_digar(digar, chunk.beg, chunk.end)) continue;
            variants.push_back(CandidateVariant{variant_key_from_digar(read.tid, digar), VariantCounts{}});
        }
    }
    collapse_fuzzy_large_insertions(variants, min_sv_len);
}

// ════════════════════════════════════════════════════════════════════════════
// Allele counting
// ════════════════════════════════════════════════════════════════════════════

/**
 * Mean base quality over bases implicated by a digar op.
 * @param d Digar (uses `qi`, `len`, deletion semantics).
 * @param qual Per-read QUAL array.
 * @return Rounded mean Phred, or 0 if not applicable.
 */
int get_digar_ave_qual(const DigarOp& d, const std::vector<uint8_t>& qual) {
    if (d.low_quality) return 0;
    if (d.qi < 0) return 0;
    if (qual.empty()) return 0;
    int q_start = 0;
    int q_end = 0;
    if (d.type == DigarType::Deletion) {
        if (d.qi == 0) {
            q_start = 0;
            q_end = 0;
        } else {
            q_start = d.qi - 1;
            q_end = d.qi;
        }
    } else {
        q_start = d.qi;
        q_end = d.qi + d.len - 1;
    }
    if (q_start < 0) return 0;
    if (q_end >= static_cast<int>(qual.size())) return 0;
    int sum = 0;
    for (int i = q_start; i <= q_end; ++i) sum += static_cast<int>(qual[i]);
    return sum / (q_end - q_start + 1);
}

/**
 * @brief Updates strand-aware depth counters for one read observation at a site.
 *
 * Low-quality observations increment only `low_qual_cov`.
 * Strand is always counted regardless of quality.
 *
 * @param counts Counters to update.
 * @param reverse Read reverse flag.
 * @param alt True if observation supports alternate allele.
 * @param low_quality True if counted toward `low_qual_cov` only.
 */
static void add_coverage(VariantCounts& counts, bool reverse, bool alt, bool low_quality) {
    if (low_quality) {
        counts.low_qual_cov++;
        return;
    }
    counts.total_cov++;
    if (alt) {
        counts.alt_cov++;
        reverse ? counts.reverse_alt++ : counts.forward_alt++;
    } else {
        counts.ref_cov++;
        reverse ? counts.reverse_ref++ : counts.forward_ref++;
    }
}

/**
 * @brief Sweeps reads against sorted candidates to fill allele and strand counts.
 *
 * Initialize candidate allele counts from digars, then merge-update
 * (`bam_utils.c`): merge walk skipping only `BAM_CEQUAL` digars, comparing candidate vs digar
 * with `exact_comp_var_site_ins` on merge-shaped `var_site_t` (clips participate in ordering).
 * Low-quality alt when `digar.low_quality` or mean BQ \< \a min_bq; implicit ref for trailing
 * sites before `read.end`.
 *
 * @param reads Parsed chunk reads.
 * @param variants Sorted candidate table (same order as site collection).
 * @param min_bq Minimum base quality threshold for counting an alt as high-quality.
 */
void collect_allele_counts_from_records(const std::vector<ReadRecord>& reads,
                                        CandidateTable& variants,
                                        int min_bq,
                                        int min_sv_len) {
    const int n = static_cast<int>(variants.size());
    if (n == 0) return;

    for (const ReadRecord& read : reads) {
        if (read.is_skipped) continue;
        if (read.tid < 0) continue;
        const hts_pos_t pos_start = read.beg;
        const hts_pos_t pos_end = read.end;
        int site_i = get_var_site_start(variants, pos_start, n);
        size_t digar_i = 0;
        const std::vector<DigarOp>& dig = read.digars;

        while (site_i < n && digar_i < dig.size()) {
            if (const int stid = variants[static_cast<size_t>(site_i)].key.tid; stid != read.tid) {
                if (stid < read.tid) {
                    ++site_i;
                } else {
                    break;
                }
                continue;
            }
            if (dig[digar_i].type == DigarType::Equal) {
                ++digar_i;
                continue;
            }
            const DigarOp& d = dig[digar_i];
            const int ave = get_digar_ave_qual(d, read.qual);
            const LcdVarSiteMerge cand_site = merge_site_from_candidate_key(variants[static_cast<size_t>(site_i)].key);
            const LcdVarSiteMerge dig_site = merge_site_from_digar_op(d);
            const int ret = exact_comp_var_site_ins_merge(&cand_site, &dig_site, min_sv_len);
            if (ret < 0) {
                VariantCounts& c = variants[static_cast<size_t>(site_i)].counts;
                add_coverage(c, read.reverse, false, false);
                ++site_i;
            } else if (ret == 0) {
                const bool lq = d.low_quality || ave < min_bq;
                add_coverage(
                    variants[static_cast<size_t>(site_i)].counts, read.reverse, true, lq);
                ++site_i;
            } else {
                ++digar_i;
            }
        }

        for (; site_i < n; ++site_i) {
            if (const int stid = variants[static_cast<size_t>(site_i)].key.tid; stid != read.tid) {
                if (stid < read.tid) continue;
                break;
            }
            if (variants[static_cast<size_t>(site_i)].key.pos > pos_end) break;
            add_coverage(variants[static_cast<size_t>(site_i)].counts, read.reverse, false, false);
        }
    }
}


static bool is_not_candidate_category(VariantCategory c) {
    return c == VariantCategory::NonVariant ||
           c == VariantCategory::LowCoverage ||
           c == VariantCategory::StrandBias;
}

void prune_not_candidate_variants(PhasingChunk& chunk) {
    CandidateTable kept;
    kept.reserve(chunk.candidates.size());
    for (CandidateVariant& cv : chunk.candidates) {
        if (is_not_candidate_category(cv.counts.category)) {
            // Traced with real support only: most pruned candidates are a read's sequencing error.
            if (cv.counts.alt_cov >= 3 || cv.counts.category != VariantCategory::LowCoverage)
                pgphase_trace::candidate("bam", "drop", cv, "pruned_category",
                                     "cat=" + std::to_string(static_cast<int>(cv.counts.category)) + " init=" +
                                     std::to_string(static_cast<int>(cv.counts.candvarcate_initial)) + " alt=" +
                                         std::to_string(cv.counts.alt_cov) + "/" + std::to_string(cv.counts.total_cov));
            continue;
        }
        kept.push_back(std::move(cv));
    }
    chunk.candidates = std::move(kept);
}


// ════════════════════════════════════════════════════════════════════════════
// Variant classification — math utilities
// ════════════════════════════════════════════════════════════════════════════


/**
 * @brief Maps ASCII base to 0–3 (ACGT) or 4 for other.
 */
static int nt4_from_ref_char(char ch) {
    static const std::array<int8_t, 256> table = [] {
        std::array<int8_t, 256> t{};
        t.fill(4);
        t['A'] = t['a'] = 0;
        t['C'] = t['c'] = 1;
        t['G'] = t['g'] = 2;
        t['T'] = t['t'] = 3;
        return t;
    }();
    return table[static_cast<unsigned char>(ch)];
}

/**
 * @brief Encoded base at absolute `abs_pos` within a loaded reference substring starting at `ref_beg`.
 */
static int ref_nt4_at(const std::string& ref_seq, hts_pos_t ref_beg, hts_pos_t abs_pos) {
    if (ref_seq.empty() || abs_pos < ref_beg) return 4;
    const hts_pos_t rel = abs_pos - ref_beg;
    if (rel < 0 || static_cast<size_t>(rel) >= ref_seq.size()) return 4;
    return nt4_from_ref_char(ref_seq[static_cast<size_t>(rel)]);
}

/**
 * @brief True if a short indel sits inside or immediately beside a homopolymer or tandem repeat run.
 *
 * Only applied to indels with span ≤ `xid` (default 5 bp from `noisy_reg_max_xgaps`); larger
 * indels return false immediately. SNPs are also gated: `ref_len` must be ≤ `xid`.
 *
 * **Algorithm:** Checks up to 6 bp **downstream** of the variant for a repeat run, then up to
 * 6 bp **upstream**. For each direction, tests repeat unit lengths 1–6 bp; a unit qualifies if
 * at least **3 consecutive copies** are present in the reference. The first matching unit in
 * either direction makes the variant a homopolymer candidate and returns `true`.
 *
 * Example — homopolymer (unit = 1 bp):
 * ```
 * reference: ...C A A A A A A A G...
 *                  ^ deletion of one A
 * downstream: AAAAAA — 6 copies of "A" → homopolymer
 * ```
 *
 * Example — short tandem repeat (unit = 2 bp):
 * ```
 * reference: ...G C A C A C A C A T...
 *                    ^ deletion of "CA"
 * downstream: CACACA — 3 copies of "CA" → tandem repeat
 * ```
 *
 * @param var Variant to test (type must be Insertion or Deletion; SNPs are gated by xid).
 * @param ref_seq Reference substring for the chunk.
 * @param ref_beg Absolute 1-based coordinate of `ref_seq[0]`.
 * @param ref_end Unused (reserved for future use).
 * @param xid Maximum indel span to check; variants larger than this return false.
 */
bool var_is_homopolymer_pg(const VariantKey& var,
                                  const std::string& ref_seq,
                                  hts_pos_t ref_beg,
                                  hts_pos_t ref_end,
                                  int xid) {
    if (ref_seq.empty()) return false;
    hts_pos_t start_pos = 0;
    hts_pos_t end_pos = 0;
    if (var.type == VariantType::Snp) {
        start_pos = var.pos - 1;
        end_pos = var.pos + 1;
    } else if (var.type == VariantType::Insertion) {
        const int alt_len = static_cast<int>(var.alt.size());
        if (alt_len > xid) return false;
        start_pos = var.pos - 1;
        end_pos = var.pos;
    } else {
        if (var.ref_len > xid) return false;
        start_pos = var.pos + var.ref_len - 1;
        end_pos = var.pos;
    }
    (void)ref_end;
    int is_homopolymer = 1;
    constexpr int max_unit_len = 6;
    constexpr int n_check_copy_num = 3;
    uint8_t ref_bases[6];
    for (int i = 0; i < 6; ++i) {
        ref_bases[i] = static_cast<uint8_t>(ref_nt4_at(ref_seq, ref_beg, end_pos + i));
    }
    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_homopolymer = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (ref_nt4_at(ref_seq, ref_beg, end_pos + i * rep_unit_len + j) != ref_bases[j]) {
                    is_homopolymer = 0;
                    break;
                }
            }
            if (is_homopolymer == 0) break;
        }
        if (is_homopolymer) break;
    }
    if (is_homopolymer) return true;
    for (int i = 0; i < 6; ++i) {
        ref_bases[i] = static_cast<uint8_t>(ref_nt4_at(ref_seq, ref_beg, start_pos - i));
    }
    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_homopolymer = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (ref_nt4_at(ref_seq, ref_beg, start_pos - i * rep_unit_len - j) != ref_bases[j]) {
                    is_homopolymer = 0;
                    break;
                }
            }
            if (is_homopolymer == 0) break;
        }
        if (is_homopolymer) break;
    }
    return is_homopolymer != 0;
}

/**
 * @brief True if the reference flanking an indel is a tandem repeat of the deleted/inserted motif.
 *
 * Only applied to indels with span ≤ `xid` (default 5 bp); larger indels return false immediately.
 * Complementary to `var_is_homopolymer_pg`: catches micro-homology / STR contexts where the
 * alignment can be shifted left or right with equal edit cost, making the exact breakpoint ambiguous.
 *
 * **Algorithm:**
 * - **Deletion** of length L: compares the deleted sequence (fetched from `ref_seq` at `pos`)
 *   against the `3 × L` bases immediately following the deletion endpoint. If they match the
 *   deleted motif repeated 3 times, returns `true`.
 * - **Insertion** of length L: constructs a 3-copy repeat of the inserted sequence and compares
 *   it against `3 × L` bases of the reference starting at `pos`. Returns `true` on match.
 *
 * Example — deletion in a tandem repeat:
 * ```
 * delete "AT" at pos 100 (reference: ...ATATAT...)
 * check ref[102..107] = "ATATAT" == "AT"×3 → true
 * ```
 *
 * Example — insertion consistent with STR:
 * ```
 * insert "GC" at pos 200 (reference: ...GCGCGC...)
 * check ref[200..205] = "GCGCGC" == "GC"×3 → true
 * ```
 *
 * @param var Variant to test (must be Insertion or Deletion).
 * @param ref_seq Reference substring for the chunk.
 * @param ref_beg Absolute 1-based coordinate of `ref_seq[0]`.
 * @param ref_end Last valid absolute coordinate in `ref_seq`.
 * @param xid Maximum indel span to check; variants larger than this return false.
 * @return True if the flanking reference matches 3 tandem copies of the indel motif.
 */
bool var_is_repeat_region_pg(const VariantKey& var,
                                    const std::string& ref_seq,
                                    hts_pos_t ref_beg,
                                    hts_pos_t ref_end,
                                    int xid) {
    if (ref_seq.empty()) return false;
    const hts_pos_t pos = var.pos;
    if (var.type == VariantType::Deletion) {
        const int del_len = var.ref_len;
        if (del_len > xid) return false;
        const int len = del_len * 3;
        if (pos < ref_beg || pos + del_len + len > ref_end) return false;
        const size_t off = static_cast<size_t>(pos - ref_beg);
        const size_t off2 = static_cast<size_t>(pos + del_len - ref_beg);
        if (off + static_cast<size_t>(len) > ref_seq.size() ||
            off2 + static_cast<size_t>(len) > ref_seq.size()) return false;
        // Compare through nt4 rather than with memcmp: the reference is
        // soft-masked in exactly the tandem repeats this asks about, and the two
        // windows can straddle a mask boundary. An ambiguous base is not a
        // match, so a run of N cannot pass as a tandem repeat.
        for (int i = 0; i < len; ++i) {
            const int a = ref_nt4_at(ref_seq, ref_beg, pos + i);
            const int b = ref_nt4_at(ref_seq, ref_beg, pos + del_len + i);
            if (a > 3 || b > 3 || a != b) return false;
        }
        return true;
    }
    if (var.type == VariantType::Insertion) {
        const int ins_len = static_cast<int>(var.alt.size());
        if (ins_len > xid) return false;
        const int len = ins_len * 3;
        if (pos < ref_beg || pos + len > ref_end) return false;
        const size_t off = static_cast<size_t>(pos - ref_beg);
        if (off + static_cast<size_t>(len) > ref_seq.size()) return false;
        // Both sides are normalised to nt4 before comparison. A byte-exact
        // `ref_b == alt_b` compared an UPPERCASE alt (candidate bases are
        // uppercased at bam_digar.cpp:323) against raw reference bytes, so an
        // insertion inside a soft-masked tandem repeat -- which is most of them
        // in CHM13 -- never matched and was never classified RepeatHetIndel.
        // The sibling predicate var_is_homopolymer_pg reads the reference through
        // nt4 and is case-insensitive, so the OR at the call site gave different
        // answers for the same locus depending on which branch reached it.
        std::vector<int> ref_b(static_cast<size_t>(len));
        for (int j = 0; j < len; ++j)
            ref_b[static_cast<size_t>(j)] = ref_nt4_at(ref_seq, ref_beg, pos + j);
        std::vector<int> alt_b = ref_b;
        for (int j = ins_len; j < len; ++j)
            alt_b[static_cast<size_t>(j)] = alt_b[static_cast<size_t>(j - ins_len)];
        for (int j = 0; j < ins_len; ++j)
            alt_b[static_cast<size_t>(j)] = nt4_from_ref_char(var.alt[static_cast<size_t>(j)]);
        for (int j = 0; j < len; ++j) {
            if (ref_b[static_cast<size_t>(j)] > 3 || alt_b[static_cast<size_t>(j)] > 3) return false;
            if (ref_b[static_cast<size_t>(j)] != alt_b[static_cast<size_t>(j)]) return false;
        }
        return true;
    }
    return false;
}

// ════════════════════════════════════════════════════════════════════════════
// Noisy reads ratio cache: interval trees for coverage and error-rate queries.
// ════════════════════════════════════════════════════════════════════════════

/**
 * @brief Lazily builds read coverage and merged-XID interval trees for `var_noisy_reads_ratio`.
 * @param chunk Chunk whose `var_noisy_read_*` fields are populated once per classification pass.
 */
static void build_var_noisy_read_cache(PhasingChunk& chunk) {
    if (chunk.var_noisy_read_cov_cr) return;

    chunk.var_noisy_read_cov_cr.reset(cr_init());
    chunk.var_noisy_read_err_cr.reset(cr_init());
    chunk.var_noisy_read_marks.assign(chunk.reads.size(), 0);
    chunk.var_noisy_read_mark_id = 0;

    for (size_t read_i = 0; read_i < chunk.reads.size(); ++read_i) {
        const ReadRecord& r = chunk.reads[read_i];
        if (r.is_skipped || r.beg > r.end) continue;

        cr_add(chunk.var_noisy_read_cov_cr.get(), "cr",
               static_cast<int32_t>(r.beg - 1), static_cast<int32_t>(r.end),
               static_cast<int32_t>(read_i));

        // Merge adjacent XID events per read → noisy intervals.
        hts_pos_t noisy_start = -1;
        hts_pos_t noisy_end = -1;
        bool has_noisy = false;
        for (const DigarOp& d : r.digars) {
            if (d.type != DigarType::Snp && d.type != DigarType::Insertion &&
                d.type != DigarType::Deletion) continue;
            const hts_pos_t cur_start = d.pos - 1;
            hts_pos_t cur_end = d.pos;
            if (d.type == DigarType::Snp || d.type == DigarType::Deletion) cur_end += d.len - 1;
            if (!has_noisy) {
                noisy_start = cur_start;
                noisy_end = cur_end;
                has_noisy = true;
                continue;
            }
            if (cur_start < noisy_end) {
                if (cur_end > noisy_end) noisy_end = cur_end;
                continue;
            }
            cr_add(chunk.var_noisy_read_err_cr.get(), "cr",
                   static_cast<int32_t>(noisy_start), static_cast<int32_t>(noisy_end),
                   static_cast<int32_t>(read_i));
            noisy_start = cur_start;
            noisy_end = cur_end;
        }
        if (has_noisy) {
            cr_add(chunk.var_noisy_read_err_cr.get(), "cr",
                   static_cast<int32_t>(noisy_start), static_cast<int32_t>(noisy_end),
                   static_cast<int32_t>(read_i));
        }
    }
    cr_index(chunk.var_noisy_read_cov_cr.get());
    cr_index(chunk.var_noisy_read_err_cr.get());
}

/**
 * @brief Fraction of overlapping reads that carry a noisy digar span inside `[var_start,var_end]`.
 * @param chunk Chunk with reads (cache built on first use).
 * @param var_start Variant span start (1-based).
 * @param var_end Variant span end (1-based).
 * @return Ratio in [0,1], or 0 if no overlapping reads.
 */
static double var_noisy_reads_ratio(PhasingChunk& chunk, hts_pos_t var_start, hts_pos_t var_end) {
    build_var_noisy_read_cache(chunk);

    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    const int64_t total_n = cr_overlap(chunk.var_noisy_read_cov_cr.get(), "cr",
                                       static_cast<int32_t>(var_start - 1),
                                       static_cast<int32_t>(var_end), &ovlp_b, &max_b);
    const int total_reads = static_cast<int>(total_n);
    int noisy_reads = 0;

    if (total_reads > 0) {
        ++chunk.var_noisy_read_mark_id;
        if (chunk.var_noisy_read_mark_id == INT_MAX) {
            std::fill(chunk.var_noisy_read_marks.begin(), chunk.var_noisy_read_marks.end(), 0);
            chunk.var_noisy_read_mark_id = 1;
        }
        int64_t* noisy_b = nullptr;
        int64_t noisy_max_b = 0;
        const int64_t noisy_n = cr_overlap(chunk.var_noisy_read_err_cr.get(), "cr",
                                           static_cast<int32_t>(var_start - 1),
                                           static_cast<int32_t>(var_end), &noisy_b, &noisy_max_b);
        for (int64_t i = 0; i < noisy_n; ++i) {
            const int ri =
                static_cast<int>(cr_label(chunk.var_noisy_read_err_cr.get(), noisy_b[i]));
            if (ri < 0 || static_cast<size_t>(ri) >= chunk.var_noisy_read_marks.size()) continue;
            if (chunk.var_noisy_read_marks[static_cast<size_t>(ri)] == chunk.var_noisy_read_mark_id)
                continue;
            chunk.var_noisy_read_marks[static_cast<size_t>(ri)] = chunk.var_noisy_read_mark_id;
            ++noisy_reads;
        }
        std::free(noisy_b);
    }
    std::free(ovlp_b);
    return total_reads > 0 ? static_cast<double>(noisy_reads) / static_cast<double>(total_reads) : 0.0;
}

/**
 * @brief Adds a variant’s reference span (expanded by low-complexity overlap) to `var_cr` when filters pass.
 * @param var_cr Accumulator for dense-variant / repeat-het intervals.
 * @param low_comp_cr Low-complexity tree or nullptr.
 * @param var Locus to add.
 * @param chunk Chunk for noisy-read ratio when \a check_noisy_reads_ratio is true.
 * @param check_noisy_reads_ratio If true, requires `var_noisy_reads_ratio` ≥ `opts.min_af`.
 * @param opts Thresholds (`min_af`).
 */
static void cr_add_var_to_noisy_cr(cgranges_t* var_cr,
                                   cgranges_t* low_comp_cr,
                                   const VariantKey& var,
                                   PhasingChunk& chunk,
                                   bool check_noisy_reads_ratio,
                                   const Options& opts) {
    hts_pos_t var_start = 0;
    hts_pos_t var_end = 0;
    variant_genomic_span(var, var_start, var_end);
    if (low_comp_cr != nullptr && low_comp_cr->n_r > 0) {
        int64_t* low_comp_b = nullptr;
        int64_t max_low_comp_n = 0;
        const int64_t low_comp_n = cr_overlap(low_comp_cr, "cr",
                                              static_cast<int32_t>(var_start - 1),
                                              static_cast<int32_t>(var_end),
                                              &low_comp_b, &max_low_comp_n);
        for (int64_t j = 0; j < low_comp_n; ++j) {
            const int32_t start = cr_start(low_comp_cr, low_comp_b[j]) + 1;
            const int32_t end = cr_end(low_comp_cr, low_comp_b[j]);
            if (start < var_start) var_start = start;
            if (end > var_end) var_end = end;
        }
        std::free(low_comp_b);
    }
    if (!check_noisy_reads_ratio || var_noisy_reads_ratio(chunk, var_start, var_end) >= opts.min_af) {
        cr_add(var_cr, "cr",
               static_cast<int32_t>(var_start - 1), static_cast<int32_t>(var_end), 1);
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Variant classification
// ════════════════════════════════════════════════════════════════════════════

/**
 * First-pass category from depth, ONT strand bias, AF band, and repeat context.
 *
 * Sets `counts.allele_fraction`. ONT runs Fisher exact on forward vs reverse alt counts.
 *
 * @param key Variant identity.
 * @param counts Allele and strand tallies (in/out).
 * @param ref_slice Loaded reference sequence for the chunk.
 * @param ref_beg Absolute coordinate of `ref_slice[0]`.
 * @param ref_end Last valid absolute coordinate for repeat/homopolymer checks.
 * @param opts Depth, AF, strand p-value, and `noisy_reg_max_xgaps`.
 * @return Initial `VariantCategory` before noisy-overlap second pass.
 */
VariantCategory classify_variant_initial(const VariantKey& key,
                                         VariantCounts& counts,
                                         const std::string& ref_slice,
                                         hts_pos_t ref_beg,
                                         hts_pos_t ref_end,
                                         const Options& opts) {
    const int depth_with_low_quality = counts.total_cov + counts.low_qual_cov;
    counts.allele_fraction =
        counts.total_cov == 0
            ? 0.0
            : static_cast<double>(counts.alt_cov) / static_cast<double>(counts.total_cov);

    if (depth_with_low_quality < opts.min_depth) return VariantCategory::LowCoverage;
    if (counts.alt_cov < opts.min_alt_depth) return VariantCategory::LowCoverage;

    if (opts.is_ont()) {
        const int fa = counts.forward_alt;
        const int ra = counts.reverse_alt;
        const int expected = (fa + ra) / 2;
        if (expected > 0) {
            const double p = fisher_exact_two_tail(fa, ra, expected, expected);
            if (p < opts.strand_bias_pval) return VariantCategory::StrandBias;
        }
    }

    if (counts.allele_fraction < opts.min_af) return VariantCategory::LowAlleleFraction;
    if (counts.allele_fraction > opts.max_af) return VariantCategory::CleanHom;

    if ((key.type == VariantType::Insertion || key.type == VariantType::Deletion) &&
        (var_is_homopolymer_pg(key, ref_slice, ref_beg, ref_end, opts.noisy_reg_max_xgaps) ||
         var_is_repeat_region_pg(key, ref_slice, ref_beg, ref_end, opts.noisy_reg_max_xgaps))) {
        return VariantCategory::RepeatHetIndel;
    }
    return key.type == VariantType::Snp ? VariantCategory::CleanHetSnp : VariantCategory::CleanHetIndel;
}

/** Two-pass candidate classification: two passes by necessity, not two separate classifiers.
 *
 * **Pass 1** — Only *local* rules: `classify_variant_initial` fills `cats[]` and
 * `candvarcate_initial` (coverage, allele fraction, strand bias, repeat het, clean categories).
 * Every candidate that is not `LOW_COV` (and not skipped for ONT strand bias when building the
 * index) inserts its reference interval into `var_pos_cr`, which is then indexed. Density checks
 * need the full set of surviving loci in that structure before any per-site `n_ov` query.
 *
 * **Pass 2** — *Context* rules: overlap with existing `chunk_noisy` → `NON_VAR` (HiFi);
 * `RepeatHetIndel` handling adds to `noisy_var_cr`; `cr_overlap` on `var_pos_cr` with `n_ov > 1`
 * flags dense clusters for `noisy_var_cr`; then `LOW_AF` is folded to `LOW_COV` (
 * `var_i_to_cate` after this loop). Merged ranges are written back to `chunk.noisy_regions`.
 *
 * Note:  `NOISY_CAND_HET` / `NOISY_CAND_HOM` (`e` / `h`) are **not** produced here; they
 * come from MSA-based recall inside noisy regions. This pass matches  `classify_cand_vars`
 * (candidates remain `CLEAN_*` / `REP_HET_INDEL` until `NON_VAR` from containment, etc.).
 *
 * @param chunk Chunk with `candidates`, `ref_seq`, `low_complexity_regions`, and `noisy_regions` inputs.
 * @param opts Platform and merge options.
 */
static void classify_cand_vars_pgphase(PhasingChunk& chunk, const Options& opts) {
    CandidateTable& variants = chunk.candidates;
    if (variants.empty()) return;

    CrangesOwner low_own;
    low_own.adopt(intervals_to_cr(noisy_extension_intervals(chunk)));
    cgranges_t* low_comp = low_own.cr;
    if (low_comp != nullptr) cr_index(low_comp);

    // Classification does not modify cand.alt_ref_base (only the init pass
    // sets 4 and noisy `make_cand_vars_from_baln0` sets consensus column bytes). SNP ref_base is filled
    // here for classifier parity with existing pgphase paths.
    for (auto& cv : variants) {
        if (cv.key.type == VariantType::Snp) {
            cv.ref_base = static_cast<uint8_t>(ref_nt4_at(chunk.ref_seq, chunk.ref_beg, cv.key.pos));
        }
    }

    // Pass 1: local labels + interval index of all non–LOW_COV candidates (see function comment).
    std::vector<VariantCategory> cats;
    cats.reserve(variants.size());
    for (auto& cv : variants) {
        const VariantCategory init = classify_variant_initial(
            cv.key, cv.counts, chunk.ref_seq, chunk.ref_beg, chunk.ref_end, opts);
        cv.counts.candvarcate_initial = init;
        cats.push_back(init);
    }

    cgranges_t* var_pos_cr = cr_init();
    for (size_t i = 0; i < variants.size(); ++i) {
        const VariantCategory c = cats[i];
        if (c == VariantCategory::LowCoverage) continue;
        if (opts.is_ont() && c == VariantCategory::StrandBias) continue;
        const VariantKey& k = variants[i].key;
        if (k.type == VariantType::Insertion) {
            cr_add(var_pos_cr, "cr", static_cast<int32_t>(k.pos - 1), static_cast<int32_t>(k.pos), 1);
        } else {
            cr_add(var_pos_cr, "cr",
                   static_cast<int32_t>(k.pos - 1),
                   static_cast<int32_t>(k.pos + k.ref_len - 1), 1);
        }
    }
    cr_index(var_pos_cr);

    cgranges_t* chunk_noisy = intervals_to_cr(chunk.noisy_regions);
    if (chunk_noisy == nullptr) {
        chunk_noisy = cr_init();
    } else {
        cr_index(chunk_noisy);
    }

    cgranges_t* noisy_var_cr = cr_init();
    int64_t* ovlp_b = nullptr;
    int64_t max_b = 0;
    const hts_pos_t reg_beg = chunk.region.beg;
    const hts_pos_t reg_end = chunk.region.end;

    // Pass 2: noisy mask, dense overlaps via var_pos_cr, noisy_var_cr merge, LOW_AF→LOW_COV.
    for (size_t i = 0; i < variants.size(); ++i) {
        VariantCategory c = cats[i];
        // Pass 2: skip NON_VAR and STRAND_BIAS before overlap checks.
        if (c == VariantCategory::NonVariant || c == VariantCategory::StrandBias) continue;

        if (chunk_noisy->n_r > 0) {
            int64_t n = 0;
            if (variants[i].key.type == VariantType::Insertion) {
                n = cr_overlap(chunk_noisy, "cr",
                               static_cast<int32_t>(variants[i].key.pos - 1),
                               static_cast<int32_t>(variants[i].key.pos),
                               &ovlp_b, &max_b);
            } else {
                n = cr_overlap(chunk_noisy, "cr",
                               static_cast<int32_t>(variants[i].key.pos - 1),
                               static_cast<int32_t>(variants[i].key.pos + variants[i].key.ref_len - 1),
                               &ovlp_b, &max_b);
            }
            if (n > 0) {
                cats[i] = VariantCategory::NonVariant;
                continue;
            }
        }

        if (c == VariantCategory::LowCoverage) continue;

        if (c == VariantCategory::RepeatHetIndel) {
            if (variants[i].key.pos >= reg_beg && variants[i].key.pos <= reg_end) {
                cr_add_var_to_noisy_cr(noisy_var_cr, low_comp, variants[i].key, chunk, false, opts);
            }
            continue;
        }

        int64_t n_ov = 0;
        if (variants[i].key.type == VariantType::Insertion) {
            n_ov = cr_overlap(var_pos_cr, "cr",
                              static_cast<int32_t>(variants[i].key.pos - 1),
                              static_cast<int32_t>(variants[i].key.pos),
                              &ovlp_b, &max_b);
        } else {
            n_ov = cr_overlap(var_pos_cr, "cr",
                              static_cast<int32_t>(variants[i].key.pos - 1),
                              static_cast<int32_t>(variants[i].key.pos + variants[i].key.ref_len - 1),
                              &ovlp_b, &max_b);
        }
        if (n_ov > 1) {
            if (variants[i].key.pos >= reg_beg && variants[i].key.pos <= reg_end) {
                cr_add_var_to_noisy_cr(noisy_var_cr, low_comp, variants[i].key, chunk, true, opts);
            }
        }

        // LOW_AF is rewritten to LOW_COV in var_i_to_cate after this loop.
        if (c == VariantCategory::LowAlleleFraction) cats[i] = VariantCategory::LowCoverage;
    }

    std::free(ovlp_b);
    cr_destroy(var_pos_cr);

    // A catalog allele with a pileup call within kCatalogPileupGap bp of it that no catalog allele is: two descriptions
    // of one locus, resolved by the region's MSA (union gap phasing's alignment solve).
    if (opts.catalog_loci != nullptr && !opts.catalog_loci->empty()) {
        constexpr hts_pos_t kCatalogPileupGap = 10;
        constexpr hts_pos_t kMaxPileupSpan = 1000;  // a longer pileup call starting further back is not looked for
        using KeyTuple = std::tuple<VariantType, hts_pos_t, int, std::string>;
        const auto tuple_of = [](const VariantKey& k) { return KeyTuple{k.type, k.pos, k.ref_len, k.alt}; };
        std::vector<KeyTuple> catalog_keys;
        for (const CatalogLocus& locus : *opts.catalog_loci)
            if (locus.end >= reg_beg && locus.beg <= reg_end) catalog_keys.push_back(tuple_of(locus.key));
        std::sort(catalog_keys.begin(), catalog_keys.end());
        std::vector<std::pair<hts_pos_t, size_t>> by_pos;  // pileup candidates by start
        for (size_t i = 0; i < variants.size(); ++i) {
            if (cats[i] == VariantCategory::LowCoverage || cats[i] == VariantCategory::NonVariant ||
                cats[i] == VariantCategory::StrandBias) continue;
            if (std::binary_search(catalog_keys.begin(), catalog_keys.end(), tuple_of(variants[i].key))) continue;
            by_pos.emplace_back(variants[i].key.pos, i);
        }
        std::sort(by_pos.begin(), by_pos.end());
        for (const CatalogLocus& locus : *opts.catalog_loci) {
            if (locus.end < reg_beg || locus.beg > reg_end) continue;
            hts_pos_t lo = locus.beg, hi = locus.end;
            bool differs = false;
            for (auto it = std::lower_bound(by_pos.begin(), by_pos.end(),
                                            std::make_pair(locus.beg - kCatalogPileupGap - kMaxPileupSpan, size_t{0}));
                 it != by_pos.end() && it->first <= locus.end + kCatalogPileupGap; ++it) {
                const VariantKey& k = variants[it->second].key;
                hts_pos_t vs = 0, ve = 0;
                variant_genomic_span(k, vs, ve);
                if (ve < locus.beg - kCatalogPileupGap) continue;
                differs = true;
                lo = std::min(lo, vs);
                hi = std::max(hi, ve);
            }
            if (!differs) continue;
            VariantKey span;
            span.type = VariantType::Deletion;
            span.pos = lo;
            span.ref_len = static_cast<int>(hi - lo + 1);
            cr_add_var_to_noisy_cr(noisy_var_cr, low_comp, span, chunk, false, opts);
        }
    }

    if (noisy_var_cr->n_r > 0) {
        cr_index(noisy_var_cr);
        cgranges_t* merged =
            cr_merge2(chunk_noisy, noisy_var_cr, -1, opts.noisy_reg_merge_dis, opts.min_sv_len);
        cr_destroy(chunk_noisy);
        cr_destroy(noisy_var_cr);
        chunk_noisy = merged;
    } else {
        cr_destroy(noisy_var_cr);
    }

    intervals_from_cr(chunk_noisy, chunk.noisy_regions);
    cr_destroy(chunk_noisy);

    for (size_t i = 0; i < variants.size(); ++i) {
        variants[i].counts.category = cats[i];
        variants[i].lcd_var_i_to_cate = category_to_flag(cats[i]);
    }
}

/**
 * Classify all chunk candidates (two-pass logic).
 * @param chunk Chunk with allele counts populated.
 * @param opts Classification options.
 * @param header Reserved (unused).
 */
void classify_chunk_candidates(PhasingChunk& chunk, const Options& opts, const bam_hdr_t* header) {
    (void)header;
    classify_cand_vars_pgphase(chunk, opts);
}



/**
 * Per-chunk candidate collection and phasing scaffold.
 *
 * BAM/FASTA I/O is intentionally outside this function: callers prepare `chunk.reads`,
 * reference slices, and initial noisy/low-complexity state before entering the numbered
 * steps below.
 *
 * **Phase 1 — site collection (steps 1.1–1.2)**
 * Candidate sites (SNP/INS/DEL) are extracted from read digars, then each site is visited
 * again to accumulate ref/alt/low-qual counts and optional per-read support rows.
 *
 * **Phase 2 — classification (steps 2.1–2.6)**
 * Read-level noisy intervals are merged and filtered (`pre_process_noisy_regs`), candidates
 * are labelled (CLEAN_HET_SNP, LOW_COV, STRAND_BIAS, REP_HET_INDEL, …), noisy regions are
 * widened using classified context (`post_process_noisy_regs`), and — for non-ONT reads —
 * candidates fully contained in noisy spans are demoted to NON_VAR.
 *
 * **Guard** — if both `chunk.candidates` and `chunk.noisy_regions` are empty after
 * classification the function returns early: there is nothing to phase and nothing to
 * propagate to downstream callers.
 *
 * **Phase 3 — phasing scaffold (steps 3.1–3.2)**
 * `collect_read_var_profile` builds a sparse per-read allele table: for every candidate site
 * each read covers, it records 0 (ref), 1 (alt), or -2 (low-quality alt).
 * `assign_hap_based_on_germline_het_vars_kmeans` runs **once** here with **`kCandGermlineClean`**
 * (kCandCleanHetSnp | kCandCleanHetIndel | kCandCleanHom).
 * The second germline k-means with `kCandGermlineVarCate` happens after noisy-region MSA
 * (step 4), which `collect_noisy_vars_step4` runs when MSA-recalled noisy candidates are
 * merged.
 * The somatic path is not yet implemented.
 *
 * **Example — two het sites phased by shared reads, one hom site**
 * @code
 *   Candidate sites:   A (pos 100, CleanHetSnp)   B (pos 200, CleanHetSnp)   C (pos 300, CleanHom)
 *
 *   Read 1: ----[A=alt]--------[B=ref]----     Agrees with (A=1, B=0) consensus
 *   Read 2: ----[A=alt]--------[B=ref]----     Reinforces (A=1, B=0)
 *   Read 3:            --------[B=alt]--[C=alt] Conflicting read if assigned to same hap; implies opposite
 *
 *   K-Means Output Phase set anchored at A:
 *     A → phase_set=100, hap_alt=1, hap_ref=2   (anchor: alt on hap1)
 *     B → phase_set=100, hap_alt=2, hap_ref=1   (polarized relative to A: alt on hap2)
 *     C → phase_set=300, hap_alt=3, hap_ref=0   (CleanHom: alt on both, no K-Means iteration involvement)
 * @endcode
 *
 * @param chunk Chunk with reads and region pre-filled by the caller.
 * @param opts Thresholds (min BQ, AF, strand p-value, SV length, verbosity, …).
 * @param header BAM header for contig names used in verbose logging.
 */
// True when a variant's VCF position falls inside the chunk's active region.
static bool active_region_contains(hts_pos_t active_reg_beg,
                                                     hts_pos_t active_reg_end,
                                                     const VariantKey& key) {
    const hts_pos_t vcf_pos = (key.type == VariantType::Deletion || key.type == VariantType::Insertion)
                                  ? key.pos - 1
                                  : key.pos;
    return !(vcf_pos < active_reg_beg || vcf_pos > active_reg_end);
}

// A co-located deletion that duplicates an allele of a merged multiallelic
// record is a second description of the same locus, and the wrong one. The merge
// runs per noisy region, so a record produced by a different region or an
// earlier pass survives it: at chr20:55,883,019 the merged record reads
// AATATAT -> AAT,A at 1|2, matching the competitor, while a leftover 4 bp
// deletion reported 1 reference against 36 alt at allele fraction 0.973 -- it
// scores the other haplotype's reads against its own allele -- and was emitted
// 1|1, a homozygous call at a locus whose two haplotypes are the alleles beside
// it. Demote the leftover so the prune pass drops it; its allele is not lost, it
// is allele 1 or 2 of the merged record.
//
// The same holds for insertions, where the allele is compared directly rather
// than through a retained length. At chr20:55,919,945 the MSA merges the locus
// into one record carrying alleles A and AA, emitted C>CA,CAA at 2|1 with 11
// reference, 26 and 19 reads. Injection has already added the catalog's C>CAA
// as a separate single-allele candidate with alt "AA", because injection runs
// BEFORE the MSA -- at that moment the alignment record is still single-allele
// and there is no allele set to match against, which is why fixing this in
// find_matching_candidate was measured inert. The leftover can only be
// recognised here, once the merged record exists.
static void drop_superseded_colocated_records(PhasingChunk& chunk) {
    for (const CandidateVariant& merged : chunk.candidates) {
        if (merged.msa_insertion_alts.empty()) continue;
        for (CandidateVariant& other : chunk.candidates) {
            if (&other == &merged) continue;
            if (other.key.type != merged.key.type ||
                other.key.pos != merged.key.pos ||
                !other.msa_insertion_alts.empty()) continue;
            bool duplicate = false;
            if (merged.key.type == VariantType::Deletion) {
                const int retained = merged.key.ref_len - other.key.ref_len;
                if (retained < 0) continue;
                duplicate = retained == 0;
                for (const std::string& allele : merged.msa_insertion_alts)
                    if (static_cast<size_t>(retained) == allele.size()) duplicate = true;
            } else {
                // Only an exact allele match is a duplicate description. A
                // record naming an allele the merge did not carry is a
                // different event and must survive.
                for (const std::string& allele : merged.msa_insertion_alts)
                    if (allele == other.key.alt) duplicate = true;
            }
            if (!duplicate) continue;
            // A record imported from a recovery sub-solve that carries a DECIDED
            // heterozygous consensus is not a duplicate description of the
            // merged record -- it describes the other haplotype. At
            // chr20:55,336,460 the alignment emits the locus as a complementary
            // pair, C>CGT at 0|1 and CGTGT>C at 1|0, and demoting the deletion
            // here left the graph arm with half the locus.
            if (other.bam_injected && other.hap_to_cons_alle.size() > 2 &&
                other.hap_to_cons_alle[1] >= 0 && other.hap_to_cons_alle[2] >= 0 &&
                other.hap_to_cons_alle[1] != other.hap_to_cons_alle[2])
                continue;
            pgphase_trace::candidate("msa", "drop", other, "superseded_by_two_allele_record",
                                     "by=" + pgphase_trace::key_text(merged.key));
            other.counts.category = VariantCategory::LowCoverage;
            other.lcd_var_i_to_cate = 0;
        }
    }
}

/// Discover candidates from read alignments, count their alleles, and classify
/// them. Read-level noisy regions are built first; classification adds candidate
/// spans, then containment and pruning remove calls inside noisy intervals.
// Diagnostic (Options::diag_dir): every candidate as classified, before
// pruning, and the noisy regions. Columns: VCF-convention position (indels at
// their anchor), key position, type S/I/D, REF length, ALT (first 40 bases),
// initial and final category codes (VariantCategory order), REF / ALT / total
// coverage, and 1 when kept from a noisy region too long for MSA.
static void write_classify_diag(const PhasingChunk& chunk, const std::string& dir) {
    const std::string tag = std::to_string(chunk.region.tid) + "." + std::to_string(chunk.region.beg) + ".tsv";
    if (FILE* fh = std::fopen((dir + "/aln." + tag).c_str(), "w")) {
        for (const CandidateVariant& c : chunk.candidates) {
            const char type = c.key.type == VariantType::Snp ? 'S' : c.key.type == VariantType::Insertion ? 'I' : 'D';
            std::fprintf(fh, "%lld\t%lld\t%c\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
                         static_cast<long long>(c.key.sort_pos()), static_cast<long long>(c.key.pos), type,
                         c.key.ref_len, c.key.alt.substr(0, 40).c_str(),
                         static_cast<int>(c.counts.candvarcate_initial), static_cast<int>(c.counts.category),
                         c.counts.ref_cov, c.counts.alt_cov, c.counts.total_cov, c.unresolved_noisy ? 1 : 0);
        }
        std::fclose(fh);
    }
    if (FILE* fh = std::fopen((dir + "/noisy." + tag).c_str(), "w")) {
        for (const Interval& r : chunk.noisy_regions)
            std::fprintf(fh, "%lld\t%lld\n", static_cast<long long>(r.beg), static_cast<long long>(r.end));
        std::fclose(fh);
    }
}

// The clean hom SNPs inside the chunk's noisy regions that are no candidate any longer -- suppressed for their
// region's MSA to recall (apply_noisy_containment_filter), or masked by classification's per-read noisy windows
// before those were merged into the regions -- set aside for restore_failed_noisy_homs, before the prune drops them.
static void hold_noisy_homs(PhasingChunk& chunk) {
    chunk.held_noisy_homs.clear();
    if (chunk.noisy_regions.empty()) return;
    std::vector<std::pair<hts_pos_t, hts_pos_t>> spans;
    for (const Interval& r : chunk.noisy_regions) spans.emplace_back(r.beg, r.end);
    std::sort(spans.begin(), spans.end());
    std::vector<hts_pos_t> reach(spans.size());  // the furthest end of the spans up to each
    for (size_t i = 0; i < spans.size(); ++i) reach[i] = std::max(spans[i].second, i ? reach[i - 1] : spans[i].second);
    for (const CandidateVariant& cv : chunk.candidates) {
        if (cv.counts.category != VariantCategory::NonVariant ||
            cv.counts.candvarcate_initial != VariantCategory::CleanHom || cv.key.type != VariantType::Snp)
            continue;
        const auto it = std::upper_bound(spans.begin(), spans.end(), std::make_pair(cv.key.pos, cv.key.pos),
                                         [](const auto& a, const auto& b) { return a.first < b.first; });
        if (it == spans.begin() || reach[static_cast<size_t>(it - spans.begin()) - 1] < cv.key.pos) continue;
        chunk.held_noisy_homs.push_back(cv);
    }
}

void collect_var_classify(PhasingChunk& chunk,
                          const Options& opts,
                          const bam_hdr_t* header) {
    // Parsed X/I/D events propose sites; read profiles later score the same sites.
    collect_candidate_sites_from_records(
        chunk.region, chunk.reads, chunk.candidates, opts.min_sv_len);

    collect_allele_counts_from_records(chunk.reads,
                                       chunk.candidates,
                                       opts.min_bq,
                                       opts.min_sv_len);

    // Read-level noise gives classification its initial exclusion intervals.
    pre_process_noisy_regs_pgphase(chunk, opts);

    // Classification uses allele counts and augments the noisy-region model.
    if (opts.verbose >= 2) {
        std::fprintf(stderr,
                     "PRE_CLASSIFY noisy_regions (%zu intervals):\n",
                     chunk.noisy_regions.size());
        for (const Interval& iv : chunk.noisy_regions) {
            std::fprintf(stderr,
                         "  %" PRId64 "-%" PRId64 " (label=%d)\n",
                         static_cast<int64_t>(iv.beg),
                         static_cast<int64_t>(iv.end),
                         iv.label);
        }
    }
    classify_chunk_candidates(chunk, opts, header);

    // Merge noisy spans after candidate context has been added.
    post_process_noisy_regs_pgphase(chunk, chunk.candidates);

    // A contained candidate is suppressed before profiles and phasing are built.
    split_long_noisy_regions(chunk, opts);
    apply_noisy_containment_filter(chunk, opts.keep_calls_in_unresolved_noisy ? opts.max_noisy_reg_len : 0);
    if (!opts.diag_dir.empty()) write_classify_diag(chunk, opts.diag_dir);
    if (opts.restore_failed_noisy_homs) hold_noisy_homs(chunk);
    prune_not_candidate_variants(chunk);

    dump_all_noisy_regions(chunk, opts, header);
}


/// Phase the clean candidates, recall sites in noisy regions with MSA, and
/// re-phase when recall adds sites. Then fill MSA strand counts, remove records
/// superseded by merged alleles, and mark candidates inside the active region.
void collect_var_run_phasing(PhasingChunk& chunk, const Options& opts,
                             const VariantKeySet* noisy_site_whitelist) {
    if (chunk.candidates.empty() && chunk.noisy_regions.empty()) return;

    if (!chunk.candidates.empty()) {
        // co-phasing using clean-region SNPs + indels + clean hom.
        assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineClean);
    }

    // MSA can add candidates and rerun k-means on the expanded site set.
    collect_noisy_vars_step4(chunk, opts, noisy_site_whitelist);
    // The MSA path derives ref_cov/alt_cov from alle_covs and never writes the
    // strand tallies, so its candidates reach the strand-bias screen with 0+0
    // and are exempted from it. One sweep fills them from the read profiles.
    derive_msa_candidate_strand_counts(chunk);
    // Step 4 is what creates the merged multiallelic records, so the superseded
    // duplicates can only be identified after it has run.
    if (opts.collapse_colocated_alleles)
        drop_superseded_colocated_records(chunk);

    const hts_pos_t active_reg_beg = chunk.region.beg;
    const hts_pos_t active_reg_end = chunk.region.end;
    for (CandidateVariant& cand : chunk.candidates) {
        cand.lcd_make_variants_region_pass =
            active_region_contains(active_reg_beg, active_reg_end, cand.key);
    }
}

/// Build read profiles and solve the BAM chunk. Graph recovery calls the two
/// stages separately so it can add catalog sites between them.
void collect_var_phase(PhasingChunk& chunk,
                       const Options& opts) {
    collect_var_build_profiles(chunk, opts);
    collect_var_run_phasing(chunk, opts);
}

} // namespace pgphase_collect
