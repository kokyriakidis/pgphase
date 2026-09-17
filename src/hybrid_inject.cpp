/// @file hybrid_inject.cpp
/// @brief Augment BAM-derived PhasingChunk with graph snarl observations.

#include "hybrid_inject.hpp"

#include "collect_phase.hpp"
#include "collect_var.hpp"
#include "fisher_exact.hpp"
#include "noise_filter.hpp"

#include <htslib/vcf.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <unordered_set>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

static int resolve_private_vcf_tid(const std::string& contig,
                                   const bam_hdr_t* bam_header) {
    int matched_tid = -1;
    for (int tid = 0; tid < bam_header->n_targets; ++tid) {
        const std::string name = bam_header->target_name[tid];
        if (name == contig) return tid;
        const size_t hash = name.rfind('#');
        const std::string suffix = hash == std::string::npos
                                       ? name
                                       : name.substr(hash + 1);
        if (suffix != contig) continue;
        if (matched_tid >= 0) return -1;
        matched_tid = tid;
    }
    return matched_tid;
}

VariantKeySet load_private_variant_keys(const std::string& path,
                                        const bam_hdr_t* bam_header) {
    using BcfFilePtr = std::unique_ptr<htsFile, decltype(&hts_close)>;
    using BcfHeaderPtr = std::unique_ptr<bcf_hdr_t, decltype(&bcf_hdr_destroy)>;
    using BcfRecordPtr = std::unique_ptr<bcf1_t, decltype(&bcf_destroy)>;

    BcfFilePtr fp(bcf_open(path.c_str(), "r"), hts_close);
    if (!fp) throw std::runtime_error("failed to open private-sites VCF: " + path);
    BcfHeaderPtr header(bcf_hdr_read(fp.get()), bcf_hdr_destroy);
    if (!header)
        throw std::runtime_error("failed to read private-sites VCF header: " + path);
    BcfRecordPtr record(bcf_init(), bcf_destroy);
    if (!record)
        throw std::runtime_error("failed to allocate private-sites VCF record");

    VariantKeySet keys;
    int read_status = 0;
    while ((read_status = bcf_read(fp.get(), header.get(), record.get())) == 0) {
        bcf_unpack(record.get(), BCF_UN_STR);
        if (record->n_allele != 2) continue;
        const char* contig = bcf_hdr_id2name(header.get(), record->rid);
        if (contig == nullptr) continue;
        const int tid = resolve_private_vcf_tid(contig, bam_header);
        if (tid < 0) continue;
        keys.insert(vcf_to_variant_key(
            tid, record->pos + 1, record->d.allele[0], record->d.allele[1]));
    }
    if (read_status < -1)
        throw std::runtime_error("failed while reading private-sites VCF: " + path);
    return keys;
}

size_t retain_private_bam_candidates(PhasingChunk& chunk,
                                     const VariantKeySet& private_keys) {
    CandidateTable retained;
    retained.reserve(std::min(chunk.candidates.size(), private_keys.size()));
    for (CandidateVariant& candidate : chunk.candidates) {
        if (private_keys.find(candidate.key) != private_keys.end())
            retained.push_back(std::move(candidate));
    }
    chunk.candidates = std::move(retained);
    return chunk.candidates.size();
}

void clear_bam_evidence_at_graph_candidates(
        PhasingChunk& chunk,
        const std::unordered_set<int>& graph_only_candidates) {
    for (const int candidate_i : graph_only_candidates) {
        if (candidate_i < 0 || candidate_i >= static_cast<int>(chunk.candidates.size()))
            continue;
        CandidateVariant& candidate = chunk.candidates[static_cast<size_t>(candidate_i)];
        candidate.counts.ref_cov = 0;
        candidate.counts.alt_cov = 0;
        candidate.counts.total_cov = 0;
        candidate.counts.low_qual_cov = 0;
        candidate.counts.forward_ref = 0;
        candidate.counts.reverse_ref = 0;
        candidate.counts.forward_alt = 0;
        candidate.counts.reverse_alt = 0;
        candidate.counts.alle_covs.clear();
    }

    for (ReadVariantProfile& profile : chunk.read_var_profile) {
        if (profile.start_var_idx < 0) continue;
        for (const int candidate_i : graph_only_candidates) {
            if (candidate_i < profile.start_var_idx || candidate_i > profile.end_var_idx)
                continue;
            const size_t offset = static_cast<size_t>(candidate_i - profile.start_var_idx);
            if (offset < profile.alleles.size()) profile.alleles[offset] = -1;
        }
    }
}

// ────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ────────────────────────────────────────────────────────────────────────────

// Declared in hybrid_inject.hpp.
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
    } else if (vcf_alt.size() > vcf_ref.size()) {
        // Strip the full shared prefix, mirroring the deletion branch and the
        // BAM convention (variant_key_from_digar: ref_len=0, alt=inserted
        // bases).  A clean left-anchored insertion consumes the entire REF
        // (shared == ref.size()) and yields ref_len=0; a multi-base anchor like
        // TA->TAAA reduces to pos after the shared run, alt="A", ref_len=0.  The
        // previous single-base strip mis-encoded such sites (pos and alt off by
        // the extra anchor bases).  For a single-base anchor this reduces to the
        // old result.
        size_t shared = 0;
        const size_t max_shared = std::min(vcf_ref.size(), vcf_alt.size());
        while (shared < max_shared && vcf_ref[shared] == vcf_alt[shared]) ++shared;
        key.type = VariantType::Insertion;
        key.pos = vcf_pos + static_cast<hts_pos_t>(shared);
        key.ref_len = static_cast<int>(vcf_ref.size() - shared);
        key.alt = vcf_alt.substr(shared);
    } else {
        // Strip the full shared prefix so deletions use the same normalized form
        // as the BAM path (variant_key_from_digar) and the standalone graph path
        // (graph_collect.cpp): pos = first deleted base, ref_len = deleted span,
        // alt = "" for a pure deletion.  The previous single-base anchor strip
        // left a residual base in alt and an inflated ref_len for homopolymer
        // deletions (e.g. TAA->TA yielded ref_len=2, alt="A"), which (a) blocked
        // bridging to BAM deletions in find_matching_candidate (alt never matched
        // BAM's "") and (b) collided with the BAM deletion in
        // merge_chunk_candidates, silently overwriting the verified BAM call.
        // For a single-base anchor (shared == 1) this reduces to the old result.
        size_t shared = 0;
        const size_t max_shared = std::min(vcf_ref.size(), vcf_alt.size());
        while (shared < max_shared && vcf_ref[shared] == vcf_alt[shared]) ++shared;
        key.type = VariantType::Deletion;
        key.pos = vcf_pos + static_cast<hts_pos_t>(shared);
        key.ref_len = static_cast<int>(vcf_ref.size() - shared);
        key.alt = vcf_alt.substr(shared);
    }
    return key;
}

/// Find a BAM candidate matching a VariantKey by exact position + type + allele.
///
/// The binary search requires candidates[0, n) to be sorted by sort_pos().
/// Only the original BAM candidates satisfy this (collect_var_classify sorts
/// them); graph-only candidates are appended unsorted, so callers must pass
/// n = original BAM candidate count to keep the search valid.
static int find_matching_candidate(const CandidateTable& candidates,
                                   const VariantKey& target,
                                   int n) {
    if (n > static_cast<int>(candidates.size()))
        n = static_cast<int>(candidates.size());
    const hts_pos_t target_sort = target.sort_pos();
    int lo = 0;
    int hi = n - 1;
    while (lo <= hi) {
        const int mid = lo + (hi - lo) / 2;
        const hts_pos_t mid_sort = candidates[static_cast<size_t>(mid)].key.sort_pos();
        if (mid_sort < target_sort) lo = mid + 1;
        else if (mid_sort > target_sort) hi = mid - 1;
        else {
            // Scan for exact key match at this sort_pos.
            int start = mid;
            while (start > 0 &&
                   candidates[static_cast<size_t>(start - 1)].key.sort_pos() == target_sort)
                --start;
            for (int i = start;
                 i < n &&
                 candidates[static_cast<size_t>(i)].key.sort_pos() == target_sort;
                 ++i) {
                const VariantKey& k = candidates[static_cast<size_t>(i)].key;
                if (k.type == target.type &&
                    k.ref_len == target.ref_len &&
                    k.alt == target.alt) {
                    return i;
                }
            }
            return -1;
        }
    }
    return -1;
}

/// Add a new CandidateVariant from a graph site.
///
/// The candidate is added UNCLASSIFIED (category LowCoverage, flag 0) so it is
/// excluded from k-means until its allele counts have been accumulated from
/// BAM and graph reads.  classify_graph_only_candidates() then applies the same
/// depth/AF/het gates the BAM pipeline uses before any graph site can become a
/// CleanHet phasing anchor.  Stamping CleanHet here (before counts exist) let
/// homozygous and low-support graph sites flood k-means and degrade phasing.
static int add_graph_only_candidate(PhasingChunk& chunk,
                                    const GraphSite& site,
                                    const std::string& vcf_alt,
                                    int tid) {
    const int idx = static_cast<int>(chunk.candidates.size());
    CandidateVariant cand;
    cand.key = vcf_to_variant_key(tid, site.pos, site.ref, vcf_alt);
    cand.graph_site = true;
    if (cand.key.type == VariantType::Snp && cand.key.pos >= site.pos &&
        static_cast<size_t>(cand.key.pos - site.pos) < site.ref.size()) {
        switch (site.ref[static_cast<size_t>(cand.key.pos - site.pos)]) {
            case 'A': case 'a': cand.ref_base = 0; break;
            case 'C': case 'c': cand.ref_base = 1; break;
            case 'G': case 'g': cand.ref_base = 2; break;
            case 'T': case 't': cand.ref_base = 3; break;
            default: break;
        }
    }
    cand.counts.n_uniq_alles = 2;
    cand.counts.category = VariantCategory::LowCoverage;
    cand.counts.candvarcate_initial = VariantCategory::LowCoverage;
    cand.lcd_var_i_to_cate = 0;  // excluded from k-means until gated
    cand.lcd_make_variants_region_pass = true;
    chunk.candidates.push_back(std::move(cand));
    return idx;
}

// ────────────────────────────────────────────────────────────────────────────
// Phase A: site augmentation
// ────────────────────────────────────────────────────────────────────────────

SiteToCandidateMap inject_graph_sites(
        PhasingChunk& chunk,
        const GraphSiteCatalogView& graph_sites,
        const std::unordered_map<std::string, std::string>& chrom_remap,
        const Options& opts,
        int* sites_bridged_out,
        int* sites_added_out,
        std::unordered_set<int>* graph_only_candidates_out,
        GraphOnlyVcfAlleles* graph_only_vcf_alleles_out,
        std::unordered_set<int>* all_graph_candidates_out) {
    (void)opts;
    (void)chrom_remap;
    SiteToCandidateMap site_to_candidate;
    int bridged = 0, added = 0;

    // Track pre-sort indices of graph-only candidates and their original VCF
    // (ref, alt) strings so the noise filter can screen on the catalog
    // representation, matching the standalone graph pipeline.
    std::vector<int> graph_only_pre_sort;
    std::unordered_set<int> all_graph_pre_sort;
    std::unordered_map<int, GraphOnlyVcfAllele> pre_sort_vcf_alleles;

    if (graph_sites.empty()) {
        if (sites_bridged_out) *sites_bridged_out = 0;
        if (sites_added_out) *sites_added_out = 0;
        return site_to_candidate;
    }

    const int chunk_tid = chunk.region.tid;
    const int orig_count = static_cast<int>(chunk.candidates.size());

    // find_matching_candidate binary-searches candidates[0, orig_count) and
    // requires them sorted by sort_pos().  collect_var_classify guarantees
    // this; assert it so a future change to candidate ordering fails loudly
    // here rather than silently missing bridges.
    assert(std::is_sorted(
        chunk.candidates.begin(),
        chunk.candidates.begin() + orig_count,
        [](const CandidateVariant& a, const CandidateVariant& b) {
            return a.key.sort_pos() < b.key.sort_pos();
        }));

    for (size_t si = 0; si < graph_sites.size(); ++si) {
        const GraphSite& site = graph_sites[si];
        if (!site.eligible) continue;
        if (site.alts.empty()) continue;

        // Use the same key format as GraphReadAllele.site_id.
        const std::string site_key = graph_site_key_str(site);

        // Try every non-spanning ALT for a match before falling back to adding
        // the first as graph-only. Stopping at the first ALT loses the catalog's
        // confirmation whenever the allele actually present is a later one: at
        // chr20:48,149,567 the catalog carries 11 ALTs of a 5-mer repeat
        // (+5, +20, +10, -5, +15, +35, +30, +25, -15, -10, -20) and 74 of the
        // 1178 catalog sites across this window are multi-allelic, so the first
        // ALT is frequently not the one a read supports.
        int fallback_ai = -1;
        for (size_t ai = 0; ai < site.alts.size(); ++ai) {
            const std::string& vcf_alt = site.alts[ai];
            if (vcf_alt == "*") continue;
            if (fallback_ai < 0) fallback_ai = static_cast<int>(ai);

            const VariantKey target = vcf_to_variant_key(
                chunk_tid, site.pos, site.ref, vcf_alt);

            const int match_idx = find_matching_candidate(
                chunk.candidates, target, orig_count);

            if (match_idx < 0) continue;  // try the next ALT before adding

            if (match_idx >= 0 && match_idx < orig_count) {
                chunk.candidates[match_idx].graph_site = true;
                site_to_candidate[site_key] = match_idx;
                all_graph_pre_sort.insert(match_idx);
                pre_sort_vcf_alleles[match_idx] =
                    GraphOnlyVcfAllele{site.pos, site.ref, vcf_alt};
                ++bridged;
            }
            fallback_ai = -1;  // matched, nothing to add
            break;
        }
        // No ALT matched an existing candidate exactly. Before adding one, look
        // for an indel the alignment channel already called at this position:
        // the catalog's claim is that the locus varies, and which length is
        // present there is a read measurement the alignment has already made.
        // Adding the catalog's length as a second record puts two descriptions
        // of one event in the table, and the added one is then counted as though
        // the reads carried it -- at chr20:55,903,460 the alignment calls a 1 bp
        // insertion with 58 of 63 reads behind it and read truth agrees (58 of
        // 61 covering reads at +1, exactly one at +2), yet the catalog's 2 bp
        // claim was added beside it and accumulated 59 of 63 as alt support, and
        // both were emitted.
        //
        // So the claim is honoured on the record that holds the measurement:
        // graph_site is set, which is what lets classify_graph_only_candidates
        // promote the locus, while the allele and counts stay the alignment's.
        // pre_sort_vcf_alleles is deliberately NOT set here -- it carries the
        // catalog's own ref/alt for consumers that screen on it, and this
        // candidate's allele is not the catalog's.
        if (fallback_ai >= 0) {
            const VariantKey target = vcf_to_variant_key(
                chunk_tid, site.pos, site.ref,
                site.alts[static_cast<size_t>(fallback_ai)]);
            if (target.type != VariantType::Snp) {
                for (int ci = 0; ci < orig_count; ++ci) {
                    CandidateVariant& cand = chunk.candidates[static_cast<size_t>(ci)];
                    if (cand.graph_site || cand.key.type != target.type ||
                        cand.key.pos != target.pos) continue;
                    cand.graph_site = true;
                    site_to_candidate[site_key] = ci;
                    all_graph_pre_sort.insert(ci);
                    ++bridged;
                    fallback_ai = -1;
                    break;
                }
            }
        }
        if (fallback_ai >= 0) {
            const std::string& vcf_alt = site.alts[static_cast<size_t>(fallback_ai)];
            const int new_idx = add_graph_only_candidate(
                chunk, site, vcf_alt, chunk_tid);
            site_to_candidate[site_key] = new_idx;
            graph_only_pre_sort.push_back(new_idx);
            all_graph_pre_sort.insert(new_idx);
            pre_sort_vcf_alleles[new_idx] =
                GraphOnlyVcfAllele{site.pos, site.ref, vcf_alt};
            ++added;
        }
    }

    // Re-sort candidates if new ones were added.
    if (added > 0) {
        std::vector<size_t> order(chunk.candidates.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            const hts_pos_t pa = chunk.candidates[a].key.sort_pos();
            const hts_pos_t pb = chunk.candidates[b].key.sort_pos();
            if (pa != pb) return pa < pb;
            return static_cast<int>(chunk.candidates[a].key.type) <
                   static_cast<int>(chunk.candidates[b].key.type);
        });

        std::vector<int> old_to_new(chunk.candidates.size());
        for (size_t i = 0; i < order.size(); ++i)
            old_to_new[order[i]] = static_cast<int>(i);

        CandidateTable sorted;
        sorted.reserve(chunk.candidates.size());
        for (size_t i = 0; i < order.size(); ++i)
            sorted.push_back(std::move(chunk.candidates[order[i]]));
        chunk.candidates = std::move(sorted);

        for (auto& [key, idx] : site_to_candidate)
            idx = old_to_new[static_cast<size_t>(idx)];

        // Remap graph-only indices through the permutation.
        if (graph_only_candidates_out) {
            for (int pre_idx : graph_only_pre_sort)
                graph_only_candidates_out->insert(
                    old_to_new[static_cast<size_t>(pre_idx)]);
        }
        if (all_graph_candidates_out) {
            for (int pre_idx : all_graph_pre_sort)
                all_graph_candidates_out->insert(
                    old_to_new[static_cast<size_t>(pre_idx)]);
        }
        if (graph_only_vcf_alleles_out) {
            for (const auto& [pre_idx, alleles] : pre_sort_vcf_alleles)
                (*graph_only_vcf_alleles_out)[
                    old_to_new[static_cast<size_t>(pre_idx)]] = alleles;
        }
    } else {
        if (graph_only_candidates_out) {
            for (int idx : graph_only_pre_sort)
                graph_only_candidates_out->insert(idx);
        }
        if (all_graph_candidates_out)
            *all_graph_candidates_out = std::move(all_graph_pre_sort);
        if (graph_only_vcf_alleles_out) {
            for (const auto& [idx, alleles] : pre_sort_vcf_alleles)
                (*graph_only_vcf_alleles_out)[idx] = alleles;
        }
    }

    if (sites_bridged_out) *sites_bridged_out = bridged;
    if (sites_added_out) *sites_added_out = added;
    return site_to_candidate;
}

// ────────────────────────────────────────────────────────────────────────────
// Phase B: graph read injection (after BAM profiles are built)
// ────────────────────────────────────────────────────────────────────────────

/// Extend an existing BAM read's profile at graph-only sites and record when
/// its GAF walk confirms the same allele at an existing BAM candidate.
static bool extend_bam_profile_with_graph_obs(
        PhasingChunk& chunk,
        int read_i,
        const std::vector<std::pair<int,int>>& graph_obs,  // (candidate_idx, allele)
        const std::unordered_set<int>& graph_only_candidates) {
    ReadVariantProfile& prof = chunk.read_var_profile[static_cast<size_t>(read_i)];
    const int bam_start = prof.start_var_idx;
    auto bam_alleles = std::move(prof.bam_alleles);
    auto bam_qi = std::move(prof.bam_qi);

    // Track which observations are actually applied so we only update
    // allele counts for slots that were filled (not already occupied).
    std::vector<std::pair<int,int>> applied;  // graph-only (candidate_idx, allele)
    bool confirmed = false;

    for (const auto& [cand_idx, allele] : graph_obs) {
        const bool graph_only = graph_only_candidates.count(cand_idx) != 0;
        const auto& key = chunk.candidates[static_cast<size_t>(cand_idx)].key;
        if (key.type == VariantType::Snp) {
            bool deleted = false;
            for (const auto& op : chunk.reads[static_cast<size_t>(read_i)].digars) {
                if (op.pos > key.pos) break;
                if ((op.type == DigarType::Deletion || op.type == DigarType::RefSkip) &&
                    key.pos < op.pos + op.len) {
                    deleted = true;
                    break;
                }
            }
            // A graph traversal cannot resolve the nucleotide of a SNP that
            // the BAM alignment explicitly deletes or skips.
            if (deleted) continue;
        }

        if (prof.start_var_idx < 0) {
            if (!graph_only) continue;
            // Empty profile — initialize with this single observation.
            prof.start_var_idx = cand_idx;
            prof.end_var_idx = cand_idx;
            prof.alleles = {allele};
            prof.alt_qi = {kGraphConfirmedAltQi};
            prof.graph_alleles = {allele};
            applied.emplace_back(cand_idx, allele);
            continue;
        }

        if (cand_idx >= prof.start_var_idx && cand_idx <= prof.end_var_idx) {
            const size_t offset = static_cast<size_t>(cand_idx - prof.start_var_idx);
            if (offset >= prof.alleles.size()) continue;
            const int previous_allele = prof.alleles[offset];
            if (prof.graph_alleles.size() < prof.alleles.size())
                prof.graph_alleles.resize(prof.alleles.size(), -1);
            prof.graph_alleles[offset] = allele;
            if (previous_allele == allele || previous_allele < 0) {
                if (prof.alt_qi.size() < prof.alleles.size())
                    prof.alt_qi.resize(prof.alleles.size(), -1);
                prof.alleles[offset] = allele;
                prof.alt_qi[offset] = kGraphConfirmedAltQi;
                confirmed = true;
                if (!graph_only) continue;
                if (previous_allele < 0) applied.emplace_back(cand_idx, allele);
            }
        } else if (cand_idx < prof.start_var_idx) {
            if (!graph_only) continue;
            // Extend left: prepend slots.
            const int gap = prof.start_var_idx - cand_idx;
            std::vector<int> new_alleles(static_cast<size_t>(gap), -1);
            std::vector<int> new_qi(static_cast<size_t>(gap), -1);
            std::vector<int> new_graph(static_cast<size_t>(gap), -1);
            new_alleles[0] = allele;
            new_qi[0] = kGraphConfirmedAltQi;
            new_graph[0] = allele;
            new_alleles.insert(new_alleles.end(),
                               prof.alleles.begin(), prof.alleles.end());
            new_qi.insert(new_qi.end(),
                          prof.alt_qi.begin(), prof.alt_qi.end());
            if (prof.graph_alleles.size() < prof.alleles.size())
                prof.graph_alleles.resize(prof.alleles.size(), -1);
            new_graph.insert(new_graph.end(), prof.graph_alleles.begin(),
                             prof.graph_alleles.end());
            prof.alleles = std::move(new_alleles);
            prof.alt_qi = std::move(new_qi);
            prof.graph_alleles = std::move(new_graph);
            prof.start_var_idx = cand_idx;
            applied.emplace_back(cand_idx, allele);
        } else {
            if (!graph_only) continue;
            // Extend right: append slots.
            const size_t new_span = static_cast<size_t>(cand_idx - prof.start_var_idx + 1);
            prof.alleles.resize(new_span, -1);
            prof.alt_qi.resize(new_span, -1);
            prof.graph_alleles.resize(new_span, -1);
            prof.alleles[new_span - 1] = allele;
            prof.alt_qi[new_span - 1] = kGraphConfirmedAltQi;
            prof.graph_alleles[new_span - 1] = allele;
            prof.end_var_idx = cand_idx;
            applied.emplace_back(cand_idx, allele);
        }
    }

    // Update allele counts only for observations that were actually applied.
    for (const auto& [cand_idx, allele] : applied) {
        CandidateVariant& cand =
            chunk.candidates[static_cast<size_t>(cand_idx)];
        if (allele == 0) {
            ++cand.counts.ref_cov;
            ++cand.counts.total_cov;
        } else {
            ++cand.counts.alt_cov;
            ++cand.counts.total_cov;
        }
    }

    prof.bam_alleles.assign(prof.alleles.size(), -1);
    prof.bam_qi.assign(prof.alleles.size(), -1);
    for (size_t i = 0; i < bam_alleles.size(); ++i) {
        const int offset = bam_start + static_cast<int>(i) - prof.start_var_idx;
        if (offset < 0 || static_cast<size_t>(offset) >= prof.alleles.size()) continue;
        prof.bam_alleles[offset] = bam_alleles[i];
        if (i < bam_qi.size()) prof.bam_qi[offset] = bam_qi[i];
    }
    return !applied.empty() || confirmed;
}

int inject_graph_reads(
        PhasingChunk& chunk,
        const std::vector<GraphReadAllele>& graph_rows,
        const SiteToCandidateMap& site_to_candidate,
        const std::unordered_set<int>& graph_only_candidates,
        const Options& opts,
        int* reads_extended_out) {
    for (auto& profile : chunk.read_var_profile) {
        if (!profile.bam_alleles.empty()) continue;
        profile.bam_alleles = profile.alleles;
        profile.bam_qi = profile.alt_qi;
    }
    if (graph_rows.empty() || site_to_candidate.empty()) {
        if (reads_extended_out) *reads_extended_out = 0;
        return 0;
    }

    const int chunk_tid = chunk.region.tid;

    // Build name → read index map for existing BAM reads.
    std::unordered_map<std::string, int> existing_read_idx;
    existing_read_idx.reserve(chunk.reads.size());
    for (size_t i = 0; i < chunk.reads.size(); ++i)
        existing_read_idx[chunk.reads[i].qname] = static_cast<int>(i);

    // Group graph rows by read name.
    struct ReadObs {
        int candidate_idx;
        int allele;
        int mapq;
    };
    std::unordered_map<std::string, std::vector<ReadObs>> read_observations;

    for (const GraphReadAllele& row : graph_rows) {
        // site_to_candidate represents REF and the first ALT.  A different ALT
        // is not evidence for that binary candidate.
        if (row.allele < 0 || row.allele > 1) continue;
        if (row.mapq < opts.min_mapq) continue;

        auto it = site_to_candidate.find(row.site_id);
        if (it == site_to_candidate.end()) continue;

        const int cand_idx = it->second;
        const int allele = row.allele;
        read_observations[row.read_name].push_back(
            ReadObs{cand_idx, allele, row.mapq});
    }

    int injected = 0;
    int extended = 0;

    for (auto& [read_name, obs_vec] : read_observations) {
        if (obs_vec.empty()) continue;

        auto existing_it = existing_read_idx.find(read_name);
        if (existing_it != existing_read_idx.end()) {
            // Doubly-mapped read: extend graph-only sites and retain exact
            // graph confirmation at shared BAM/graph sites.
            std::vector<std::pair<int,int>> graph_obs;
            graph_obs.reserve(obs_vec.size());
            for (const ReadObs& obs : obs_vec)
                graph_obs.emplace_back(obs.candidate_idx, obs.allele);
            if (extend_bam_profile_with_graph_obs(
                    chunk, existing_it->second, graph_obs,
                    graph_only_candidates)) {
                ++extended;
            }
            continue;
        }

        // Graph-only read: create synthetic read + profile.
        std::sort(obs_vec.begin(), obs_vec.end(),
                  [](const ReadObs& a, const ReadObs& b) {
                      return a.candidate_idx < b.candidate_idx;
                  });

        const int first_idx = obs_vec.front().candidate_idx;
        const int last_idx = obs_vec.back().candidate_idx;
        const int max_mapq = std::max_element(obs_vec.begin(), obs_vec.end(),
            [](const ReadObs& a, const ReadObs& b) {
                return a.mapq < b.mapq;
            })->mapq;

        const int read_id = static_cast<int>(chunk.reads.size());

        ReadRecord read;
        read.tid = chunk_tid;
        read.input_index = 0;
        read.qname = read_name;
        read.mapq = max_mapq;
        read.is_skipped = false;
        read.beg = chunk.candidates[static_cast<size_t>(first_idx)].key.sort_pos();
        read.end = chunk.candidates[static_cast<size_t>(last_idx)].key.sort_pos();
        if (read.end < read.beg) read.end = read.beg;

        ReadVariantProfile profile;
        profile.read_id = read_id;
        profile.start_var_idx = first_idx;
        profile.end_var_idx = last_idx;
        const size_t span = static_cast<size_t>(last_idx - first_idx + 1);
        profile.alleles.assign(span, -1);
        profile.alt_qi.assign(span, kGraphConfirmedAltQi);
        profile.graph_alleles.assign(span, -1);
        profile.bam_alleles.assign(span, -1);
        profile.bam_qi.assign(span, -1);

        for (const ReadObs& obs : obs_vec) {
            const int offset = obs.candidate_idx - first_idx;
            if (offset >= 0 && static_cast<size_t>(offset) < span) {
                profile.alleles[static_cast<size_t>(offset)] = obs.allele;
                profile.graph_alleles[static_cast<size_t>(offset)] = obs.allele;
            }
        }

        chunk.reads.push_back(std::move(read));
        chunk.read_var_profile.push_back(std::move(profile));

        // Update allele counts on candidates.
        for (const ReadObs& obs : obs_vec) {
            CandidateVariant& cand =
                chunk.candidates[static_cast<size_t>(obs.candidate_idx)];
            if (obs.allele == 0) {
                ++cand.counts.ref_cov;
                ++cand.counts.total_cov;
            } else {
                ++cand.counts.alt_cov;
                ++cand.counts.total_cov;
            }
        }

        ++injected;
    }

    // Rebuild the read_var_cr interval tree with all reads (extended profiles
    // may have changed span boundaries).
    if (injected > 0 || extended > 0) {
        cgranges_t* cr = cr_init();
        if (cr == nullptr)
            throw std::runtime_error("failed to allocate hybrid read_var_cr");
        for (const ReadVariantProfile& p : chunk.read_var_profile) {
            if (p.start_var_idx < 0 || p.end_var_idx < p.start_var_idx)
                continue;
            cr_add(cr, "cr", p.start_var_idx, p.end_var_idx + 1, p.read_id);
        }
        cr_index(cr);
        chunk.read_var_cr.reset(cr);
    }

    if (reads_extended_out) *reads_extended_out = extended;
    return injected;
}

// ────────────────────────────────────────────────────────────────────────────
// Noise filter for graph-only candidates
// ────────────────────────────────────────────────────────────────────────────

void apply_hybrid_noise_filter(
        PhasingChunk& chunk,
        const std::string& ref_seq,
        hts_pos_t ref_beg,
        hts_pos_t ref_end,
        const std::unordered_set<int>& graph_only_candidates,
        int max_xgaps,
        const GraphOnlyVcfAlleles* graph_only_vcf_alleles,
        bool trim_minimal) {
    if (ref_seq.empty() || graph_only_candidates.empty()) return;

    const std::vector<Interval> lc = find_low_complexity_intervals(ref_seq, ref_beg);

    for (int ci : graph_only_candidates) {
        if (ci < 0 || static_cast<size_t>(ci) >= chunk.candidates.size()) continue;
        CandidateVariant& cand = chunk.candidates[static_cast<size_t>(ci)];
        if (cand.counts.category != VariantCategory::CleanHetIndel) continue;

        const VariantKey& k = cand.key;
        if (k.type == VariantType::Snp) {
            // SNPs in low-complexity regions are NOT demoted: both the BAM
            // pipeline (recalls them as NOISY_CAND_HET via noisy-region MSA)
            // and the standalone graph pipeline (CLEAN_HET_SNP) keep these as
            // real het calls, so the hybrid pipeline must keep them too.
            continue;
        }

        // Prefer the original VCF (pos, ref, alt) from the snarl catalog so the
        // homopolymer/repeat verdict matches the standalone graph pipeline
        // (apply_graph_noise_filter screens on the catalog representation).
        // Reconstructing from the normalized VariantKey can yield a different
        // string/position for the same site and over-demote graph het indels
        // that the graph pipeline keeps as phasing anchors.  Fall back to
        // reconstruction only when the original is unavailable.
        hts_pos_t noisy_pos = k.pos;
        std::string vcf_ref, vcf_alt;
        bool have_alleles = false;
        if (graph_only_vcf_alleles) {
            auto it = graph_only_vcf_alleles->find(ci);
            if (it != graph_only_vcf_alleles->end()) {
                noisy_pos = it->second.pos;
                vcf_ref = it->second.ref;
                vcf_alt = it->second.alt;
                have_alleles = true;
            }
        }
        if (!have_alleles) {
            if (k.type == VariantType::Insertion) {
                const hts_pos_t anchor = k.pos - 1;
                if (anchor < ref_beg || anchor > ref_end) continue;
                const char anchor_base =
                    ref_seq[static_cast<size_t>(anchor - ref_beg)];
                vcf_ref = std::string(1, anchor_base);
                vcf_alt = anchor_base + k.alt;
            } else {
                const hts_pos_t anchor = k.pos - 1;
                if (anchor < ref_beg || k.pos + k.ref_len - 1 > ref_end) continue;
                vcf_ref = ref_seq.substr(
                    static_cast<size_t>(anchor - ref_beg),
                    static_cast<size_t>(k.ref_len + 1));
                const char anchor_base =
                    ref_seq[static_cast<size_t>(anchor - ref_beg)];
                vcf_alt = std::string(1, anchor_base);
                if (!k.alt.empty()) vcf_alt += k.alt;
            }
        }

        // Experimental: trim to minimal VCF form before the noise check, the
        // same normalization apply_graph_noise_filter performs. Catalog alleles
        // carry the full repeat run on both flanks, so the derived indel can be
        // longer than max_xgaps and escape repeat detection; trimming the shared
        // suffix then prefix exposes the true indel span.
        if (trim_minimal) {
            trim_to_minimal_vcf(noisy_pos, vcf_ref, vcf_alt);
        }

        if (is_noisy_site(noisy_pos, vcf_ref, vcf_alt, ref_seq, ref_beg, ref_end,
                          lc, max_xgaps)) {
            cand.counts.category = VariantCategory::RepeatHetIndel;
            cand.counts.candvarcate_initial = VariantCategory::RepeatHetIndel;
            cand.lcd_var_i_to_cate = kLongcalldRepHetVar;
        }
    }
}

// ────────────────────────────────────────────────────────────────────────────
// BAM count backfill for graph-only candidates
// ────────────────────────────────────────────────────────────────────────────

void backfill_graph_candidate_counts(
        PhasingChunk& chunk,
        const std::unordered_set<int>& graph_only_candidates) {
    if (graph_only_candidates.empty()) return;

    for (const ReadVariantProfile& prof : chunk.read_var_profile) {
        if (prof.start_var_idx < 0) continue;

        for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
            if (!graph_only_candidates.count(vi)) continue;

            const size_t offset = static_cast<size_t>(vi - prof.start_var_idx);
            if (offset >= prof.alleles.size()) continue;
            const int allele = prof.alleles[offset];

            if (allele < 0) continue;  // -1 (uninformative) or -2 (low qual)

            CandidateVariant& cand = chunk.candidates[static_cast<size_t>(vi)];
            if (allele == 0) {
                ++cand.counts.ref_cov;
            } else {
                ++cand.counts.alt_cov;
            }
            ++cand.counts.total_cov;
        }
    }

    // Publish the allele fraction these counts imply. Accumulating ref_cov,
    // alt_cov and total_cov without it left allele_fraction at its initial 0,
    // and every downstream gate reads that field rather than the counts -- so a
    // graph-only candidate with textbook heterozygous coverage was judged as
    // having no alternate support and dropped. chr20:48,173,317
    // (TGGGGATG>T, a 7 bp deletion hiphase phases) backfilled to ref_cov 23,
    // alt_cov 42, total_cov 65 and allele_fraction 0.000, so it never reached
    // the candidate table at all. Same expression the BAM and graph-only paths
    // use (gap_evidence.cpp:256, graph_bam_adapter.cpp:693).
    for (const int vi : graph_only_candidates) {
        if (vi < 0 || static_cast<size_t>(vi) >= chunk.candidates.size()) continue;
        VariantCounts& c = chunk.candidates[static_cast<size_t>(vi)].counts;
        c.allele_fraction = c.total_cov > 0
                ? static_cast<double>(c.alt_cov) / static_cast<double>(c.total_cov)
                : 0.0;
    }
}

// ────────────────────────────────────────────────────────────────────────────
// Quality gate for graph-only candidates
// ────────────────────────────────────────────────────────────────────────────

int classify_graph_only_candidates(
        PhasingChunk& chunk,
        const std::unordered_set<int>& graph_only_candidates,
        const Options& opts) {
    if (graph_only_candidates.empty()) return 0;

    int promoted = 0;
    for (int ci : graph_only_candidates) {
        if (ci < 0 || static_cast<size_t>(ci) >= chunk.candidates.size())
            continue;
        CandidateVariant& cand = chunk.candidates[static_cast<size_t>(ci)];
        VariantCounts& c = cand.counts;

        // Classify graph-only candidates with the SAME depth/AF/het/hom logic
        // the standalone graph pipeline uses (classify_graph_candidates), NOT
        // the BAM classify_variant_initial.  The two differ in one place: the
        // BAM classifier demotes homopolymer/repeat het indels to
        // RepeatHetIndel inline, keyed on the normalized VariantKey.  Doing
        // that here over-demotes graph het indels (the graph's main phasing
        // advantage) and pre-empts apply_hybrid_noise_filter, which screens on
        // the original catalog ref/alt.  So homopolymer screening is left
        // entirely to apply_hybrid_noise_filter (called next), mirroring the
        // graph pipeline's classify_graph_candidates + apply_graph_noise_filter
        // ordering.  Only CleanHet{Snp,Indel} enter het k-means.
        c.allele_fraction =
            c.total_cov == 0
                ? 0.0
                : static_cast<double>(c.alt_cov) /
                      static_cast<double>(c.total_cov);

        VariantCategory cat;
        if (c.total_cov < opts.min_depth || c.alt_cov < opts.min_alt_depth) {
            cat = VariantCategory::LowCoverage;
        } else if (opts.is_ont() &&
                   [&] {
                       const int fa = c.forward_alt;
                       const int ra = c.reverse_alt;
                       const int expected = (fa + ra) / 2;
                       if (expected <= 0) return false;
                       return fisher_exact_two_tail(fa, ra, expected, expected) <
                              opts.strand_bias_pval;
                   }()) {
            cat = VariantCategory::StrandBias;
        } else if (c.allele_fraction < opts.min_af) {
            // Fold LowAlleleFraction into LowCoverage so it is pruned from
            // output, matching the BAM pipeline's classify_chunk_candidates
            // pass 2 (which rewrites LOW_AF to LOW_COV before prune).
            cat = VariantCategory::LowCoverage;
        } else if (c.allele_fraction > opts.max_af) {
            cat = VariantCategory::CleanHom;
        } else if (cand.key.type == VariantType::Snp) {
            cat = VariantCategory::CleanHetSnp;
        } else if (std::abs(c.allele_fraction - 0.5) <= opts.graph_indel_af_margin &&
                   c.alt_cov >= opts.graph_indel_min_alt) {
            cat = VariantCategory::CleanHetIndel;
        } else {
            // An indel whose AF is off the centre of the graph window but still
            // inside the BAM classifier's own [min_af, max_af] is a het the BAM
            // pipeline would call -- as NoisyCandHet, the class that exists for
            // a het needing MSA verification. Mapping it to LowCoverage instead
            // deleted it, because that is what prune_not_candidate_variants
            // removes: chr20:48,173,317 (TGGGGATG>T) is discovered AND phased
            // 0|1 by collect-bam-variation standalone at ref_cov 23 /
            // alt_cov 41, AF 0.6406, yet vanished from the hybrid run purely
            // because the catalog claimed the site and routed it through this
            // classifier's tighter +/-0.11 window. Same thresholds as the BAM
            // path, so a catalog-claimed site is not judged more harshly than
            // the identical site would be without the catalog.
            //
            // Admitting it here as NoisyCandHet is NOT yet correct, and the
            // reason is a verification asymmetry: `msa_verified` is set in
            // exactly one place, inside the factory that CONSTRUCTS a candidate
            // from the MSA consensus (collect_phase_noisy.cpp:213). It is a
            // property of MSA-created candidates, not a stamp applied to
            // existing ones, so a catalog-claimed candidate -- which exists
            // before that pass -- can never acquire it. Probed over
            // chr20:48,145,000-48,240,000: all 14 BAM-discovered NoisyCandHet
            // sites carry msa_verified = 1, all 5 catalog-claimed ones carry 0.
            // Admitting the latter puts unverified sites, three of them at
            // AF 0.725-0.791, into noisy k-means: read concordance on the arm
            // fell from 100.00% (393/393) to 88.39% (449/508), with both blocks
            // internally inconsistent at 87.06% and 90.45%.
            //
            // The fix is to let the noisy MSA pass construct the verified
            // version of such a site, as it already does for the identical site
            // when the catalog does not claim it -- not to admit it raw, and
            // not to hide the admission behind a flag. Until then the site is
            // still dropped, but now for a named reason.
            cat = VariantCategory::LowCoverage;
        }
        c.category = cat;
        c.candvarcate_initial = cat;
        const bool clean_het = cat == VariantCategory::CleanHetSnp ||
                               cat == VariantCategory::CleanHetIndel;
        const bool centred_anchor =
            std::abs(c.allele_fraction - 0.5) <= opts.anchor_af_margin;
        // Preserve off-centre graph hets as calls while preventing them from
        // voting. This matches the standalone graph pipeline and protects the
        // joint model when several paralogous loci collapse onto one graph
        // interval, producing confident but non-diploid AF near 0.25 or 0.75.
        cand.lcd_var_i_to_cate =
            clean_het && !centred_anchor
                ? kCandNonAnchorHet
                : category_to_flag(cat);
        if (clean_het && centred_anchor) {
            ++promoted;
        }
    }
    return promoted;
}

}  // namespace pgphase_collect
