#include "collect_var.hpp"
#include <map>
#include "graph_bam_adapter.hpp"

#include "collect_phase.hpp"
#include "fisher_exact.hpp"
#include "noise_filter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

extern "C" {
#include "cgranges.h"
#include "htslib/sam.h"
}

namespace pgphase_collect {

namespace {

// A single read's observation at one candidate site: which site and which allele.
struct GraphProfileObservation {
    int site_index = -1;
    int allele = -1;
};

// Minimal-VCF identity of an allele pair, expressed as views into the original
// ref/alt strings (no allocation). Two allele pairs that denote the same physical
// variant — even when emitted by different overlapping snarls — share the same
// (pos, ref-range, alt-range), so comparing these views is an exact dedup test.
struct MinimalVcfId {
    hts_pos_t pos = 0;
    const char* ref = nullptr;
    size_t ref_len = 0;
    const char* alt = nullptr;
    size_t alt_len = 0;

    bool operator==(const MinimalVcfId& o) const {
        return pos == o.pos && ref_len == o.ref_len && alt_len == o.alt_len &&
               std::memcmp(ref, o.ref, ref_len) == 0 &&
               std::memcmp(alt, o.alt, alt_len) == 0;
    }
};

// Reduce (pos, ref, alt) to minimal VCF form as offset ranges into the inputs:
// trim the shared suffix, then the shared prefix (advancing pos), keeping ≥1 base
// in each. Allocation-free — mirrors trim_to_minimal_vcf without copying.
MinimalVcfId minimal_vcf_id(hts_pos_t pos, const std::string& ref,
                            const std::string& alt) {
    size_t rb = 0, re = ref.size();
    size_t ab = 0, ae = alt.size();
    while (re - rb > 1 && ae - ab > 1 && ref[re - 1] == alt[ae - 1]) {
        --re;
        --ae;
    }
    while (re - rb > 1 && ae - ab > 1 && ref[rb] == alt[ab]) {
        ++rb;
        ++ab;
        ++pos;
    }
    return MinimalVcfId{pos, ref.data() + rb, re - rb, alt.data() + ab, ae - ab};
}

// 64-bit FNV-1a fingerprint of a minimal-VCF identity. Used as the hash-map key;
// exact MinimalVcfId comparison guards against collisions, so a clash can never
// cause an incorrect merge.
uint64_t minimal_vcf_fingerprint(const MinimalVcfId& id) {
    uint64_t h = 14695981039346656037ULL;
    auto mix = [&h](uint64_t v) {
        for (int b = 0; b < 8; ++b) {
            h ^= static_cast<uint8_t>(v & 0xFF);
            h *= 1099511628211ULL;
            v >>= 8;
        }
    };
    mix(static_cast<uint64_t>(id.pos));
    for (size_t k = 0; k < id.ref_len; ++k) {
        h ^= static_cast<uint8_t>(id.ref[k]);
        h *= 1099511628211ULL;
    }
    h ^= 0xFFULL;  // delimiter so "AB|C" and "A|BC" differ
    h *= 1099511628211ULL;
    for (size_t k = 0; k < id.alt_len; ++k) {
        h ^= static_cast<uint8_t>(id.alt[k]);
        h *= 1099511628211ULL;
    }
    return h;
}

// True when a site was eligible and had ≥2 allele walks but the walk vectors
// have been cleared (CompactGraphSiteIndex took ownership of the walk data).
bool graph_site_has_released_walk_storage(const GraphSite& site) {
    if (!site.eligible || site.allele_walks.size() < 2) return false;
    return std::all_of(site.allele_walks.begin(), site.allele_walks.end(),
                       [](const GraphWalk& walk) { return walk.empty(); });
}

// Create a CandidateVariant from a GraphSite and append it to the chunk.
// Allele counts are initially zero; they are filled after parent-gated
// observation counting.
void add_graph_candidate(GraphChunkBuildResult& out,
                         const GraphSite& site,
                         const std::string& key,
                         const std::vector<int>& allele_counts,
                         int tid) {
    CandidateVariant candidate;
    candidate.key.tid = tid;
    candidate.key.pos = site.order_pos() > 0 ? site.order_pos() : 1;
    candidate.key.type = VariantType::Snp;
    candidate.key.ref_len = 1;
    candidate.key.alt = key;
    candidate.counts.n_uniq_alles = static_cast<int>(allele_counts.size());
    candidate.counts.alle_covs = allele_counts;
    candidate.counts.ref_cov = allele_counts.empty() ? 0 : allele_counts[0];
    candidate.counts.alt_cov = 0;
    for (size_t allele = 1; allele < allele_counts.size(); ++allele) {
        candidate.counts.alt_cov += allele_counts[allele];
    }
    candidate.counts.total_cov = candidate.counts.ref_cov + candidate.counts.alt_cov;
    candidate.counts.forward_ref = candidate.counts.ref_cov;
    candidate.counts.forward_alt = candidate.counts.alt_cov;
    candidate.counts.allele_fraction =
        candidate.counts.total_cov > 0
            ? static_cast<double>(candidate.counts.alt_cov) / candidate.counts.total_cov
            : 0.0;
    candidate.counts.category = VariantCategory::CleanHetSnp;
    candidate.counts.candvarcate_initial = VariantCategory::CleanHetSnp;
    candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
    candidate.hap_to_cons_alle[0] = -1;
    candidate.hap_to_cons_alle[1] = -1;
    candidate.hap_to_cons_alle[2] = -1;
    // Match longcallD candidate initialization. Reads use a separate -1
    // sentinel, but an as-yet unphased candidate starts at 0.
    candidate.phase_set = kUnsetCandidatePhaseSet;
    out.chunk.candidates.push_back(std::move(candidate));
    out.site_ids.push_back(key);
    out.site_meta.push_back({
        site.ref_contig.empty() ? site.chrom : site.ref_contig,
        site.pos,
        site.ref,
        site.alts
    });
}

// Build a ReadRecord + ReadVariantProfile from a read's site observations
// and append them to the chunk.  The profile's allele vector spans from the
// first to last observed site index, with -1 for unobserved positions.
void add_read_profile(GraphChunkBuildResult& out,
                      const std::string& read_name,
                      const std::vector<GraphProfileObservation>& observations,
                      int mapq) {
    if (observations.empty()) return;

    ReadRecord read;
    read.tid = 0;
    read.input_index = 0;
    read.qname = read_name;
    read.mapq = mapq;
    read.is_skipped = false;
    read.beg = out.chunk.candidates[static_cast<size_t>(observations.front().site_index)].key.sort_pos();
    read.end = out.chunk.candidates[static_cast<size_t>(observations.back().site_index)].key.sort_pos();
    if (read.end < read.beg) read.end = read.beg;

    ReadVariantProfile profile;
    profile.read_id = static_cast<int>(out.chunk.reads.size());
    profile.start_var_idx = observations.front().site_index;
    profile.end_var_idx = observations.back().site_index;
    profile.alleles.assign(static_cast<size_t>(profile.end_var_idx - profile.start_var_idx + 1), -1);
    profile.alt_qi.assign(profile.alleles.size(), 0);
    for (const GraphProfileObservation& obs : observations) {
        const int offset = obs.site_index - profile.start_var_idx;
        if (offset < 0 || static_cast<size_t>(offset) >= profile.alleles.size()) continue;
        profile.alleles[static_cast<size_t>(offset)] = obs.allele;
    }

    out.chunk.reads.push_back(std::move(read));
    out.chunk.read_var_profile.push_back(std::move(profile));
}

// Build the cgranges interval index over read profiles so that
// k-means phasing can quickly find reads overlapping each candidate.

// Record a read's allele at a site.  If the same site was already seen with
// a different allele (from an overlapping chunk), mark it conflicted (-1).
void merge_phase_read_observation(PhaseReadOutputRow& row,
                                  const std::string& site_id,
                                  int allele) {
    auto [it, inserted] = row.allele_by_site.emplace(site_id, allele);
    if (!inserted && it->second != allele) it->second = -1;
}

// Fold a chunk's hap/PS assignment into a read's running output row. A later
// phased assignment owns the read, matching downstream chunk ownership. An
// unphased overlap cannot erase an earlier valid assignment: it contributes no
// contrary haplotype evidence and previously discarded truth-consistent tags.
void merge_phase_read_assignment(PhaseReadOutputRow& row,
                                 int chunk_id,
                                 int hap,
                                 hts_pos_t phase_set,
                                 bool is_primary,
                                 bool fill_only) {
    if (row.copies == 0) {
        row.chunk_id = chunk_id;
    } else if (row.chunk_id != chunk_id) {
        row.chunk_id = -1;
    }
    ++row.copies;

    const bool phased = (hap == 1 || hap == 2) && phase_set > 0;
    const bool may_replace = !fill_only || !row.has_phased_assignment;
    if (may_replace &&
        ((phased && (is_primary || !row.has_primary_assignment)) ||
         !row.has_phased_assignment)) {
        row.hap = hap;
        row.phase_set = phase_set;
        row.has_phased_assignment = phased;
    }
    if (phased && is_primary) row.has_primary_assignment = true;
}

// Classify graph candidates using count-based rules matching the BAM pipeline's
// classify_variant_initial (collect_var.cpp), minus reference-sequence-dependent
// checks (homopolymer/repeat detection).  Runs in a single pass over candidates.
void classify_graph_candidates(PhasingChunk& chunk, const Options& opts) {
    for (auto& cand : chunk.candidates) {
        VariantCounts& c = cand.counts;

        // Low coverage / low alt depth.
        if (c.total_cov < opts.min_depth || c.alt_cov < opts.min_alt_depth) {
            c.category = VariantCategory::LowCoverage;
            c.candvarcate_initial = VariantCategory::LowCoverage;
            cand.lcd_var_i_to_cate = category_to_flag(VariantCategory::LowCoverage);
            continue;
        }

        // ONT strand-bias filter (Fisher exact test on alt forward vs reverse).
        if (opts.is_ont()) {
            const int fa = c.forward_alt;
            const int ra = c.reverse_alt;
            const int expected = (fa + ra) / 2;
            if (expected > 0) {
                const double p = fisher_exact_two_tail(fa, ra, expected, expected);
                if (p < opts.strand_bias_pval) {
                    c.category = VariantCategory::StrandBias;
                    c.candvarcate_initial = VariantCategory::StrandBias;
                    cand.lcd_var_i_to_cate = category_to_flag(VariantCategory::StrandBias);
                    continue;
                }
            }
        }

        if (c.n_uniq_alles > 2) {
            int max_cov = 0;
            for (int cov : c.alle_covs) max_cov = std::max(max_cov, cov);
            const double top_fraction =
                c.total_cov > 0 ? static_cast<double>(max_cov) / c.total_cov : 0.0;
            if (top_fraction > opts.max_af) {
                c.category = VariantCategory::CleanHom;
                c.candvarcate_initial = VariantCategory::CleanHom;
                cand.lcd_var_i_to_cate = category_to_flag(VariantCategory::CleanHom);
                continue;
            }
            c.category = VariantCategory::CleanHetIndel;
            c.candvarcate_initial = VariantCategory::CleanHetIndel;
            cand.lcd_var_i_to_cate = kCandCleanHetIndel;
            continue;
        }

        // Low allele fraction → folded to LowCoverage (matches BAM pipeline pass 2).
        if (c.allele_fraction < opts.min_af) {
            c.category = VariantCategory::LowCoverage;
            c.candvarcate_initial = VariantCategory::LowAlleleFraction;
            cand.lcd_var_i_to_cate = category_to_flag(VariantCategory::LowCoverage);
            continue;
        }

        // Homozygous (high AF).
        if (c.allele_fraction > opts.max_af) {
            c.category = VariantCategory::CleanHom;
            c.candvarcate_initial = VariantCategory::CleanHom;
            cand.lcd_var_i_to_cate = category_to_flag(VariantCategory::CleanHom);
            continue;
        }

        // Surviving het — type-aware classification.
        //
        // Allele fraction gates k-means participation for every site, not just
        // indels.  Where the graph collapses two paralogous loci into one, reads
        // from both pile up together: a site that is het on one copy and hom on
        // the other lands near AF 0.25 or 0.75, which clears the 0.20/0.80
        // depth filters and then votes as if it were a haplotype marker.  The
        // pipeline consequently separates paralogs rather than haplotypes, and
        // does so *confidently* -- on HG002 chr20:65.99-66.21 Mb the discordant
        // reads carry a higher evidence margin (24) than the concordant ones
        // (16).  In that window 32.9% of anchors sit outside AF 0.35-0.65
        // against 11.8% chromosome-wide.
        //
        // Only lcd_var_i_to_cate changes: the site is still emitted as a call,
        // it just stops voting.
        const bool af_centred =
            std::abs(c.allele_fraction - 0.5) <= opts.anchor_af_margin;
        // The comment above says the site "is still emitted as a call, it just
        // stops voting" -- but kLongcalldLowAfVar is outside kCandGermlineClean,
        // which is what gates VCF output, so the site was in fact dropped
        // entirely.  kCandNonAnchorHet expresses the stated intent: outside the
        // anchor mask, inside the germline mask.
        const uint32_t non_anchor_cate =
            opts.emit_nonanchor_hets ? kCandNonAnchorHet : kLongcalldLowAfVar;
        if (cand.key.type == VariantType::Snp) {
            c.category = VariantCategory::CleanHetSnp;
            c.candvarcate_initial = VariantCategory::CleanHetSnp;
            cand.lcd_var_i_to_cate =
                af_centred ? kCandCleanHetSnp : non_anchor_cate;
        } else if (!af_centred) {
            c.category = VariantCategory::CleanHetIndel;
            c.candvarcate_initial = VariantCategory::CleanHetIndel;
            cand.lcd_var_i_to_cate = non_anchor_cate;
        } else {
            c.category = VariantCategory::CleanHetIndel;
            c.candvarcate_initial = VariantCategory::CleanHetIndel;
            cand.lcd_var_i_to_cate = kCandCleanHetIndel;
        }
    }
}

} // namespace

// Defined at namespace scope so callers that reorder candidates (the in-chunk
// recovery merge) can rebuild the index; it was file-local until then.
void rebuild_read_var_cr(PhasingChunk& chunk) {
    cgranges_t* cr = cr_init();
    if (cr == nullptr) throw std::runtime_error("failed to allocate graph read profile cgranges");
    for (const ReadVariantProfile& profile : chunk.read_var_profile) {
        if (profile.start_var_idx < 0 || profile.end_var_idx < profile.start_var_idx) continue;
        cr_add(cr, "cr", profile.start_var_idx, profile.end_var_idx + 1, profile.read_id);
    }
    cr_index(cr);
    chunk.read_var_cr.reset(cr);
}

size_t promote_link_supported_repeat_indels(GraphChunkBuildResult& result,
                                            const Options& opts) {
    if (!opts.link_earned_repeat_indels) return 0;
    PhasingChunk& chunk = result.chunk;
    const size_t n = chunk.candidates.size();
    if (n == 0) return 0;

    // Per-candidate read alleles, binarised: 0 reference, 1 any alt. The
    // profile carries an allele index, and which alt a read carries does not
    // matter for agreement with a biallelic neighbour.
    std::vector<std::vector<std::pair<int, int>>> obs(n);  // (read_id, 0|1)
    for (const ReadVariantProfile& prof : chunk.read_var_profile) {
        if (prof.start_var_idx < 0) continue;
        for (size_t k = 0; k < prof.alleles.size(); ++k) {
            const size_t ci = static_cast<size_t>(prof.start_var_idx) + k;
            if (ci >= n) break;
            const int a = prof.alleles[k];
            if (a < 0) continue;
            obs[ci].emplace_back(prof.read_id, a == 0 ? 0 : 1);
        }
    }

    std::vector<size_t> trusted;
    for (size_t i = 0; i < n; ++i)
        if (chunk.candidates[i].counts.category == VariantCategory::CleanHetSnp)
            trusted.push_back(i);
    if (trusted.empty()) return 0;

    size_t promoted = 0;
    for (size_t i = 0; i < n; ++i) {
        CandidateVariant& cand = chunk.candidates[i];
        if (cand.counts.category != VariantCategory::RepeatHetIndel) continue;
        if (obs[i].size() < static_cast<size_t>(opts.link_earned_min_reads)) continue;
        std::map<int, int> mine;
        for (const auto& o : obs[i]) mine[o.first] = o.second;

        // Nearest trusted neighbours on either side, by candidate order.
        const auto it = std::lower_bound(trusted.begin(), trusted.end(), i);
        std::vector<size_t> near;
        for (auto j = it; j != trusted.end() && near.size() < 3; ++j) near.push_back(*j);
        for (auto j = it; j != trusted.begin() && near.size() < 6;) near.push_back(*--j);

        bool earned = false;
        for (size_t j : near) {
            int same = 0, cross = 0;
            for (const auto& o : obs[j]) {
                const auto f = mine.find(o.first);
                if (f == mine.end()) continue;
                if (f->second == o.second) ++same;
                else ++cross;
            }
            const int total = same + cross;
            if (total < opts.link_earned_min_reads) continue;
            const double purity = static_cast<double>(std::max(same, cross)) / total;
            if (purity >= opts.link_earned_min_purity) { earned = true; break; }
        }
        if (!earned) continue;
        cand.counts.category = VariantCategory::CleanHetIndel;
        cand.counts.candvarcate_initial = VariantCategory::CleanHetIndel;
        cand.lcd_var_i_to_cate = kCandCleanHetIndel;
        ++promoted;
    }
    return promoted;
}

void apply_graph_noise_filter(GraphChunkBuildResult& result,
                              const std::string& ref_seq,
                              hts_pos_t ref_beg,
                              hts_pos_t ref_end,
                              int max_xgaps) {
    if (ref_seq.empty()) return;

    const std::vector<Interval> lc = find_low_complexity_intervals(ref_seq, ref_beg);

    for (size_t ci = 0; ci < result.chunk.candidates.size(); ++ci) {
        CandidateVariant& cand = result.chunk.candidates[ci];
        // Only reclassify surviving het indels.
        if (cand.counts.category != VariantCategory::CleanHetIndel) continue;

        // Get VCF-style ref/alt from the parallel site_meta vector.
        if (ci >= result.site_meta.size()) continue;
        const GraphSiteMeta& meta = result.site_meta[ci];
        if (meta.alts.empty()) continue;

        // Check every alt allele of the (possibly multiallelic) site. A het indel
        // sitting in a homopolymer / STR context is noise regardless of which alt
        // is eventually phased, so demote the site if ANY alt is noisy.
        //
        // Graph alleles are not minimal: a small indel inside a repeat is emitted
        // with the full repeat run on both sides (e.g. CAAAAAAAAAAA > CAAAAAAAA for a
        // 3 bp deletion). is_noisy_site derives the indel length from the allele
        // sizes, so without trimming the length exceeds max_xgaps and the
        // homopolymer/repeat scan is skipped. Trim each allele pair to its minimal
        // VCF representation first.
        bool noisy = false;
        for (const std::string& raw_alt : meta.alts) {
            std::string vcf_ref = meta.ref;
            std::string vcf_alt = raw_alt;
            hts_pos_t vcf_pos = meta.pos;
            trim_to_minimal_vcf(vcf_pos, vcf_ref, vcf_alt);
            if (is_noisy_site(vcf_pos, vcf_ref, vcf_alt, ref_seq, ref_beg, ref_end,
                              lc, max_xgaps)) {
                noisy = true;
                break;
            }
        }
        if (noisy) {
            cand.counts.category = VariantCategory::RepeatHetIndel;
            cand.counts.candvarcate_initial = VariantCategory::RepeatHetIndel;
            cand.lcd_var_i_to_cate = kLongcalldRepHetVar;
        }
    }
}

// Main entry point: convert graph-space allele observations into a PhasingChunk.
//
// Pipeline within this function:
//   1. Enumerate eligible catalog sites → initial multi-allelic candidates
//   2. Wire parent→child snarl relationships for gating
//   3. Single-pass over rows: intern read names, collect per-read allele obs
//   4. Parent-gated observation counting (child obs only counted when the
//      read also observes a qualifying parent allele)
//   5. Phase 1: drop individual alt alleles below min_alt_depth
//   6. Phase 2: decompose multi-allelic sites into biallelic ref/alt pairs,
//      apply depth + AF filters per pair
//   7. Classify surviving candidates (het/hom/low-cov/strand-bias)
//   8. Phase 3: remap read observations to biallelic pair space, drop pruned
//   9. Build ReadVariantProfile + cgranges index for k-means
GraphChunkBuildResult build_graph_chunk(const GraphSiteCatalogView& catalog,
                                               const std::vector<GraphReadAllele>& rows,
                                               const std::string& /*contig*/,
                                               hts_pos_t beg,
                                               hts_pos_t end,
                                               int chunk_id,
                                               const Options& opts) {
    const int chunk_tid = 0;

    GraphChunkBuildResult out;
    out.chunk.region.chunk_id = chunk_id;
    out.chunk.region.reg_chunk_i = chunk_id;
    out.chunk.region.tid = chunk_tid;
    out.chunk.region.beg = beg + 1;
    out.chunk.region.end = end;
    out.chunk.ref_beg = beg + 1;
    out.chunk.ref_end = end;
    out.chunk.chunk_min_qual = 60;
    out.chunk.chunk_first_quar_qual = 60;
    out.chunk.chunk_median_qual = 60;
    out.chunk.chunk_third_quar_qual = 60;
    out.chunk.chunk_max_qual = 60;
    out.chunk.up_ovlp_read_i.assign(1, {});
    out.chunk.down_ovlp_read_i.assign(1, {});
    out.chunk.n_up_ovlp_skip_reads.assign(1, 0);
    out.chunk.n_down_ovlp_skip_reads.assign(1, 0);

    std::unordered_map<std::string, int> site_to_candidate;
    std::vector<std::vector<int>> allele_counts;
    // Per-site, per-allele strand counts — accumulated from raw rows before dedup.
    std::vector<std::vector<int>> fwd_strand_counts;
    std::vector<std::vector<int>> rev_strand_counts;
    std::vector<int> parent_candidate;
    std::vector<std::vector<int>> conditional_parent_alleles;
    for (size_t site_i = 0; site_i < catalog.size(); ++site_i) {
        const GraphSite& site = catalog[site_i];
        const std::string sid = graph_site_key_str(site);
        if (!site.eligible) {
            out.filtered_sites.push_back({sid, site.ref_contig.empty() ? site.chrom : site.ref_contig,
                                          site.pos, 0, 0, 0, 0.0, "precandidate_ineligible"});
            continue;
        }
        const bool released_walk_storage = graph_site_has_released_walk_storage(site);
        if (!released_walk_storage && !graph_site_is_queryable(site)) {
            out.filtered_sites.push_back({sid, site.ref_contig.empty() ? site.chrom : site.ref_contig,
                                          site.pos, 0, 0, 0, 0.0, "precandidate_not_queryable"});
            continue;
        }
        const int n_alleles = static_cast<int>(site.allele_walks.size());
        if (n_alleles < 2) {
            out.filtered_sites.push_back({sid, site.ref_contig.empty() ? site.chrom : site.ref_contig,
                                          site.pos, 0, 0, 0, 0.0, "precandidate_monoallelic"});
            continue;
        }
        site_to_candidate.emplace(sid, static_cast<int>(out.chunk.candidates.size()));
        allele_counts.emplace_back(static_cast<size_t>(n_alleles), 0);
        fwd_strand_counts.emplace_back(static_cast<size_t>(n_alleles), 0);
        rev_strand_counts.emplace_back(static_cast<size_t>(n_alleles), 0);
        parent_candidate.push_back(-1);
        conditional_parent_alleles.push_back(site.conditional_parent_alleles);
        add_graph_candidate(out, site, sid, allele_counts.back(), chunk_tid);
    }

    // Wire parent→child snarl relationships.  A nested snarl's observations
    // are only counted when the read also traverses a qualifying allele of the
    // parent snarl (conditional_parent_alleles, from the VCF PA field).
    for (size_t site_i = 0; site_i < catalog.size(); ++site_i) {
        const GraphSite& site = catalog[site_i];
        if (site.parent.empty()) continue;
        const std::string key = graph_site_key_str(site);
        auto child_it = site_to_candidate.find(key);
        auto parent_it = site_to_candidate.find(site.parent);
        if (child_it == site_to_candidate.end() || parent_it == site_to_candidate.end()) continue;
        parent_candidate[static_cast<size_t>(child_it->second)] = parent_it->second;
    }

    // Intern read names into contiguous integer IDs to replace string-keyed maps.
    std::unordered_map<std::string, uint32_t> read_name_to_id;
    std::vector<std::string> read_id_to_name;
    read_name_to_id.reserve(rows.size());
    auto intern_read = [&](const std::string& name) -> uint32_t {
        auto [it, inserted] = read_name_to_id.emplace(name, static_cast<uint32_t>(read_id_to_name.size()));
        if (inserted) read_id_to_name.push_back(name);
        return it->second;
    };

    // Pack allele + reverse into a single int: low 16 bits = allele, bit 16 = reverse.
    // Conflict sentinel is -1 (allele < 0).
    constexpr int kRevBit = 0x10000;
    auto pack_allele_rev = [](int allele, bool rev) -> int {
        return allele | (rev ? kRevBit : 0);
    };
    auto unpack_allele = [](int packed) -> int { return packed & 0xFFFF; };
    auto unpack_reverse = [&](int packed) -> bool { return (packed & kRevBit) != 0; };

    // Single pass over rows: intern read names, track max mapq, and collect
    // per-read site allele observations. Replaces three separate passes.
    struct SiteAlleleEntry { int site_i; int packed; };
    // Temporary per-row (rid, site_i, packed) triples; sorted into per-read
    // vectors after the pass when n_reads is known.
    struct RowTriple { uint32_t rid; int site_i; int packed; int mapq; };
    std::vector<RowTriple> row_triples;
    row_triples.reserve(rows.size());
    // Counted because a row dropped here leaves its site with no evidence at all,
    // which surfaces much later as a "ref_only" filter that names neither cause.
    int64_t row_site_missing = 0, row_allele_oob = 0, row_kept = 0;
    for (const GraphReadAllele& row : rows) {
        const uint32_t rid = intern_read(row.read_name);
        auto it = site_to_candidate.find(row.site_id);
        if (it == site_to_candidate.end()) { ++row_site_missing; continue; }
        const int site_i = it->second;
        CandidateVariant& candidate = out.chunk.candidates[static_cast<size_t>(site_i)];
        if (row.allele < 0 || row.allele >= candidate.counts.n_uniq_alles) { ++row_allele_oob; continue; }
        ++row_kept;
        row_triples.push_back({rid, site_i, pack_allele_rev(row.allele, row.reverse), row.mapq});
    }
    if (opts.verbose >= 1) {
        std::fprintf(stderr, "[rows] total %zu kept %" PRId64 " site-unknown %" PRId64
                             " allele-out-of-range %" PRId64 "\n",
                     rows.size(), row_kept, row_site_missing, row_allele_oob);
    }
    const uint32_t n_reads = static_cast<uint32_t>(read_id_to_name.size());

    // Build max mapq and per-read allele vectors from the collected triples.
    std::vector<int> read_max_mapq(n_reads, 0);
    std::vector<std::vector<SiteAlleleEntry>> allele_by_read_site(n_reads);
    for (const RowTriple& t : row_triples) {
        if (t.mapq > read_max_mapq[t.rid]) read_max_mapq[t.rid] = t.mapq;
        allele_by_read_site[t.rid].push_back({t.site_i, t.packed});
    }
    // Sort each read's entries by site_i, then deduplicate (mark conflicts as -1).
    for (auto& entries : allele_by_read_site) {
        if (entries.empty()) continue;
        std::sort(entries.begin(), entries.end(),
                  [](const SiteAlleleEntry& a, const SiteAlleleEntry& b) {
                      return a.site_i < b.site_i;
                  });
        size_t out_i = 0;
        for (size_t j = 1; j < entries.size(); ++j) {
            if (entries[j].site_i == entries[out_i].site_i) {
                // Same site: mark conflict if alleles differ.
                if (entries[out_i].packed >= 0 &&
                    unpack_allele(entries[out_i].packed) != unpack_allele(entries[j].packed)) {
                    entries[out_i].packed = -1;
                }
            } else {
                entries[++out_i] = entries[j];
            }
        }
        entries.resize(out_i + 1);
    }
    // A read whose walk re-enters the same snarl (tandem repeats do this) yields
    // several observations for one site; if they disagree the site is marked
    // conflicted and contributes nothing.  Counted because the result is a site
    // with zero evidence, reported downstream only as "ref_only".
    if (opts.verbose >= 1) {
        int64_t conflicted = 0, kept_obs = 0;
        for (const auto& entries : allele_by_read_site)
            for (const auto& e : entries) (e.packed < 0 ? conflicted : kept_obs)++;
        std::fprintf(stderr, "[dedup] observations kept %" PRId64 " conflicted %" PRId64 "\n",
                     kept_obs, conflicted);
    }

    // Binary search helper for sorted SiteAlleleEntry vectors.
    auto find_site_entry = [](const std::vector<SiteAlleleEntry>& v, int site_i)
        -> const SiteAlleleEntry* {
        auto it = std::lower_bound(v.begin(), v.end(), site_i,
                                   [](const SiteAlleleEntry& e, int s) { return e.site_i < s; });
        if (it != v.end() && it->site_i == site_i) return &*it;
        return nullptr;
    };

    // Build per-read observation lists, applying parent gating.
    // Counted so the cost of the gate is visible: dropping an observation here
    // leaves the site with no allele evidence at all, which surfaces downstream
    // as a "ref_only" filter rather than as anything mentioning nesting.
    int64_t gate_no_parent_cand = 0, gate_no_parent_obs = 0, gate_allele_mismatch = 0, gate_kept = 0;
    std::vector<std::vector<GraphProfileObservation>> read_obs(n_reads);
    for (uint32_t rid = 0; rid < n_reads; ++rid) {
        const auto& by_site = allele_by_read_site[rid];
        for (const auto& entry : by_site) {
            const int site_i = entry.site_i;
            const int packed = entry.packed;
            if (packed < 0) continue;
            const int allele = unpack_allele(packed);
            const bool rev = unpack_reverse(packed);
            const std::vector<int>& conditional = conditional_parent_alleles[static_cast<size_t>(site_i)];
            const int parent_i = parent_candidate[static_cast<size_t>(site_i)];
            if (!conditional.empty()) {
                if (parent_i < 0) { ++gate_no_parent_cand; continue; }
                const SiteAlleleEntry* parent_obs = find_site_entry(by_site, parent_i);
                if (!parent_obs || parent_obs->packed < 0) { ++gate_no_parent_obs; continue; }
                const int parent_allele = unpack_allele(parent_obs->packed);
                if (std::find(conditional.begin(), conditional.end(), parent_allele) ==
                    conditional.end()) {
                    ++gate_allele_mismatch;
                    continue;
                }
                ++gate_kept;
            }
            ++allele_counts[static_cast<size_t>(site_i)][static_cast<size_t>(allele)];
            auto& sc = rev ? rev_strand_counts[static_cast<size_t>(site_i)]
                           : fwd_strand_counts[static_cast<size_t>(site_i)];
            ++sc[static_cast<size_t>(allele)];
            read_obs[rid].push_back(GraphProfileObservation{site_i, allele});
        }
    }
    if (opts.verbose >= 1) {
        const int64_t dropped = gate_no_parent_cand + gate_no_parent_obs + gate_allele_mismatch;
        std::fprintf(stderr,
                     "[parent-gate] kept %" PRId64 " dropped %" PRId64
                     " (no parent candidate %" PRId64 ", parent unobserved %" PRId64
                     ", parent allele mismatch %" PRId64 ")\n",
                     gate_kept, dropped, gate_no_parent_cand, gate_no_parent_obs,
                     gate_allele_mismatch);
    }

    for (size_t site_i = 0; site_i < out.chunk.candidates.size(); ++site_i) {
        CandidateVariant& candidate = out.chunk.candidates[site_i];
        const std::vector<int>& counts = allele_counts[site_i];
        const std::vector<int>& fwd = fwd_strand_counts[site_i];
        const std::vector<int>& rev = rev_strand_counts[site_i];
        candidate.counts.alle_covs = counts;
        candidate.counts.ref_cov = counts.empty() ? 0 : counts[0];
        candidate.counts.alt_cov = 0;
        for (size_t allele = 1; allele < counts.size(); ++allele) candidate.counts.alt_cov += counts[allele];
        candidate.counts.total_cov = candidate.counts.ref_cov + candidate.counts.alt_cov;
        candidate.counts.forward_ref = fwd.empty() ? 0 : fwd[0];
        candidate.counts.reverse_ref = rev.empty() ? 0 : rev[0];
        candidate.counts.forward_alt = 0;
        candidate.counts.reverse_alt = 0;
        for (size_t allele = 1; allele < fwd.size(); ++allele) candidate.counts.forward_alt += fwd[allele];
        for (size_t allele = 1; allele < rev.size(); ++allele) candidate.counts.reverse_alt += rev[allele];
        candidate.counts.allele_fraction =
            candidate.counts.total_cov > 0
                ? static_cast<double>(candidate.counts.alt_cov) / candidate.counts.total_cov
                : 0.0;
    }

    // Candidate filtering — three phases mirroring BAM §13 het classification.
    //
    // Phase 1: drop individual alt alleles below min_alt_depth (noise walks).
    // Build allele_remap[site_i][old_allele] = new_allele (or -1 if dropped).
    // Recompute per-site counts from surviving alleles only.
    const size_t n_cands = out.chunk.candidates.size();
    std::vector<std::vector<int>> allele_remap(n_cands);
    for (size_t i = 0; i < n_cands; ++i) {
        CandidateVariant& cand = out.chunk.candidates[i];
        std::vector<int>& ac = allele_counts[i];
        std::vector<int>& fc = fwd_strand_counts[i];
        std::vector<int>& rc = rev_strand_counts[i];
        allele_remap[i].assign(ac.size(), -1);
        allele_remap[i][0] = 0;  // ref walk always kept at index 0
        std::vector<int> new_ac = {ac[0]};
        std::vector<int> new_fc = {fc[0]};
        std::vector<int> new_rc = {rc[0]};
        int next_allele = 1;
        for (size_t a = 1; a < ac.size(); ++a) {
            if (ac[a] >= opts.min_alt_depth) {
                allele_remap[i][a] = next_allele++;
                new_ac.push_back(ac[a]);
                new_fc.push_back(a < fc.size() ? fc[a] : 0);
                new_rc.push_back(a < rc.size() ? rc[a] : 0);
            }
        }
        {
            static const char* dbg = getenv("PGPHASE_DEBUG_SITE");
            static const long dbg_pos = dbg ? atol(dbg) : -1;
            if (dbg_pos >= 0 && i < out.site_meta.size() &&
                std::llabs(static_cast<long long>(out.site_meta[i].pos) - dbg_pos) <= 2) {
                std::string before, after;
                for (int x : ac) before += std::to_string(x) + ",";
                for (int x : new_ac) after += std::to_string(x) + ",";
                std::fprintf(stderr, "[compact %ld] cate=0x%x min_alt_depth=%d counts before=[%s] after=[%s]\n",
                             static_cast<long>(out.site_meta[i].pos), cand.lcd_var_i_to_cate,
                             opts.min_alt_depth, before.c_str(), after.c_str());
            }
        }
        allele_counts[i] = new_ac;
        fwd_strand_counts[i] = new_fc;
        rev_strand_counts[i] = new_rc;
        cand.counts.alle_covs = new_ac;
        cand.counts.n_uniq_alles = static_cast<int>(new_ac.size());
        cand.counts.ref_cov = new_ac[0];
        cand.counts.alt_cov = 0;
        for (size_t a = 1; a < new_ac.size(); ++a) cand.counts.alt_cov += new_ac[a];
        cand.counts.total_cov = cand.counts.ref_cov + cand.counts.alt_cov;
        cand.counts.forward_ref = new_fc[0];
        cand.counts.reverse_ref = new_rc[0];
        cand.counts.forward_alt = 0;
        cand.counts.reverse_alt = 0;
        for (size_t a = 1; a < new_fc.size(); ++a) cand.counts.forward_alt += new_fc[a];
        for (size_t a = 1; a < new_rc.size(); ++a) cand.counts.reverse_alt += new_rc[a];
        cand.counts.allele_fraction = cand.counts.total_cov > 0
            ? static_cast<double>(cand.counts.alt_cov) / cand.counts.total_cov
            : 0.0;
    }

    // Build inverse allele mapping: for each surviving new allele index, record the original
    // allele-walk index so callers can map back to GraphSite.alts[orig-1].
    std::vector<std::vector<int>> allele_orig_idx(n_cands);
    for (size_t i = 0; i < n_cands; ++i) {
        allele_orig_idx[i].push_back(0);  // new index 0 = ref = orig index 0
        for (size_t old_a = 1; old_a < allele_remap[i].size(); ++old_a) {
            if (allele_remap[i][old_a] >= 0) {
                allele_orig_idx[i].push_back(static_cast<int>(old_a));
            }
        }
    }

    // Phase 2: decompose each surviving site into biallelic (ref vs alt_i) pairs,
    // mirroring the BAM path where every candidate is a single ref/alt pair.
    // AF and depth filters are applied per pair; filtered pairs go to out.filtered_sites.
    // A read observing ref contributes allele 0 to every pair from that site;
    // with --snarl-allele-phasing, a read observing alt_i contributes allele 1
    // to its own pair and allele 0 to the other pairs from that site.
    struct NewPairEntry { int new_idx; int old_alt_phase1; };
    std::vector<std::vector<NewPairEntry>> old_site_to_pairs(n_cands);
    std::vector<bool> source_is_multi(n_cands, false);

    // Overlapping/nested snarls can emit the same physical variant from different
    // sites. Keyed by snarl order_pos() these look distinct, but in minimal VCF
    // form they are identical and would enter k-means as separate anchors, giving
    // one variant 2-4x the votes it should cast. Map each minimal-VCF identity to
    // the first candidate that produced it; later duplicates pool their counts
    // into that canonical candidate instead of creating a new one. Per-read
    // observation dedup (Phase 3) then collapses a read seen at both snarls into a
    // single vote.
    //
    // Keyed by a 64-bit fingerprint (cache-friendly int hashing, no per-pair
    // string allocation); the bucket stores the exact MinimalVcfId + canonical
    // index so fingerprint collisions are resolved without ever merging distinct
    // variants. Reserved to the pair upper bound to avoid rehashing.
    // old_i records the source site so the two alts of a single multiallelic
    // snarl are never collapsed into each other: they are distinct variants by
    // construction even if their minimal-VCF strings coincide degenerately
    // (e.g. when node sequences are unavailable). Only cross-site matches dedup.
    struct CanonEntry { MinimalVcfId id; int new_idx; size_t old_i; };
    std::unordered_map<uint64_t, std::vector<CanonEntry>> canonical_by_fp;
    canonical_by_fp.reserve(n_cands * 2);

    std::vector<CandidateVariant> new_cands;
    std::vector<std::string>      new_ids;
    std::vector<GraphSiteMeta>    new_meta;
    std::vector<std::vector<int>> new_allele_counts;
    std::vector<std::vector<int>> new_fwd_strand;
    std::vector<std::vector<int>> new_rev_strand;
    std::vector<std::vector<int>> new_orig_idx;
    // Most sites are biallelic (1 surviving alt); reserve n_cands as a lower bound.
    new_cands.reserve(n_cands);
    new_ids.reserve(n_cands);
    new_meta.reserve(n_cands);
    new_allele_counts.reserve(n_cands);
    new_fwd_strand.reserve(n_cands);
    new_rev_strand.reserve(n_cands);
    new_orig_idx.reserve(n_cands);

    for (size_t i = 0; i < n_cands; ++i) {
        const std::vector<int>& ac = allele_counts[i];
        if (ac.size() < 2) {
            // Chunks overlap, so a site near a boundary is built in two chunks and
            // sees reads in only one of them.  The empty copy is not a filtered
            // site -- it is an artifact of chunking -- and labelling it "ref_only"
            // makes the dump read as though the site had no alt evidence anywhere.
            // Distinguish the two so the dump can be deduplicated honestly.
            const char* reason = (ac[0] == 0) ? "no_reads_in_chunk" : "ref_only";
            out.filtered_sites.push_back({out.site_ids[i], out.site_meta[i].chrom,
                                          out.site_meta[i].pos, ac[0], 0, ac[0], 0.0,
                                          reason});
            continue;
        }
        const std::vector<int>& fc = fwd_strand_counts[i];
        const std::vector<int>& rc_s = rev_strand_counts[i];
        const int ref_c = ac[0];
        // Copied out for the filtered-site dump: the per-pair drop below runs
        // before `meta` is bound, and site_meta[i] is stable for the whole loop.
        const std::string& meta_chrom = out.site_meta[i].chrom;
        const hts_pos_t meta_pos = out.site_meta[i].pos;
        // Total depth across every observed allele at this site.  The default
        // denominator, ref_c + alt_c, treats each alt as if the site were
        // biallelic against the graph's reference allele.  Where no read carries
        // that reference allele -- 34,698 sites on HG002 chr18, median alt depth
        // 66 -- every alt then scores AF = 1.0 and the site is discarded as
        // homozygous, however the reads actually split between the alts.  A het
        // between two non-reference alleles is invisible to that model.
        // Measuring against site depth makes such a site read as, say, 40/58 and
        // 12/58 rather than 1.0 and 1.0.  For a biallelic site the two
        // denominators are equal, so only multi-allele sites change.
        int site_total = 0;
        for (int c : ac) site_total += c;
        // A multi-allelic snarl is not a set of independent ref/alt questions.
        // For pair_i the meaningful contrast is "carries alt_i" vs "carries
        // something else" -- a read on alt_j is direct evidence against alt_i,
        // not a missing observation.  Scoring it that way makes an alt_1/alt_2
        // heterozygote read as 0.5/0.5 instead of 1.0/1.0, which is what made it
        // invisible.  (--af-vs-site-depth changed only the denominator and left
        // the read profile monomorphic, which is why it degraded phasing.)
        const bool multi = opts.snarl_allele_phasing && ac.size() > 2;
        source_is_multi[i] = multi;
        // Keep a multi-allelic snarl whole rather than splitting it into
        // independent ref/alt questions.  k-means already works on integer
        // allele indices -- hap_to_cons_alle[hap] is "which allele this
        // haplotype carries", and any non-zero index projects as ALT -- so a
        // snarl with n alleles can be a single anchor where two reads agree iff
        // they carry the same allele.  That is well defined no matter how the
        // reads distribute across the alleles, which binary contrasts are not:
        // in a repeat snarl with 8+ alleles no single alt reaches an
        // informative allele fraction, and the site is lost either way.
        if (multi && opts.snarl_keep_whole) {
            int alt_total = 0;
            for (size_t a = 1; a < ac.size(); ++a) alt_total += ac[a];
            // Need two alleles with real support for the locus to say anything.
            int n_supported = 0;
            for (int c : ac) if (c >= opts.min_alt_depth) ++n_supported;
            if (site_total < opts.min_depth || n_supported < 2) {
                out.filtered_sites.push_back({out.site_ids[i], meta_chrom, meta_pos,
                                              ac[0], alt_total, site_total, 0.0,
                                              site_total < opts.min_depth ? "low_depth"
                                                                          : "multiallelic_unsupported"});
                continue;
            }
            // The sample is diploid: however many alleles the graph offers here,
            // its reads should concentrate on at most two.  When they spread
            // wider, the snarl is not resolving two haplotypes -- measured on
            // HG002 chr20, within-haplotype allele purity falls from 99% at two
            // alleles to 38% at sixteen -- and using it as an anchor assigns
            // reads confidently to the wrong haplotype.  Nesting level does not
            // predict this (LV0/LV1/LV2 purity is flat); allele spread does.
            std::vector<int> desc(ac);
            std::sort(desc.begin(), desc.end(), std::greater<int>());
            const int top2 = desc[0] + (desc.size() > 1 ? desc[1] : 0);
            const double top2_frac =
                site_total > 0 ? static_cast<double>(top2) / site_total : 0.0;
            if (top2_frac < opts.snarl_top2_frac) {
                // Too spread to trust as a single n-allelic anchor -- but dropping
                // it outright would discard signal the biallelic decomposition
                // does extract, so fall through to that instead of skipping the
                // site.  (Dropping made the gate non-monotonic: tightening it
                // removed anchors faster than it removed noise.)
                goto decompose_site;
            }
            const int new_idx = static_cast<int>(new_cands.size());
            CandidateVariant whole = out.chunk.candidates[i];
            whole.counts.alle_covs       = ac;
            whole.counts.n_uniq_alles    = static_cast<int>(ac.size());
            whole.counts.ref_cov         = ac[0];
            whole.counts.alt_cov         = alt_total;
            whole.counts.total_cov       = site_total;
            whole.counts.allele_fraction =
                site_total > 0 ? static_cast<double>(alt_total) / site_total : 0.0;
            new_cands.push_back(whole);
            new_ids.push_back(out.site_ids[i]);
            new_meta.push_back(out.site_meta[i]);
            new_allele_counts.push_back(ac);
            new_fwd_strand.push_back(fc);
            new_rev_strand.push_back(rc_s);
            new_orig_idx.push_back(allele_orig_idx[i]);
            // old_alt_phase1 = -1 marks "identity mapping": observations keep
            // their allele index instead of collapsing to 0/1.
            old_site_to_pairs[i].push_back({new_idx, -1});
            continue;
        }
    decompose_site:
        for (size_t a = 1; a < ac.size(); ++a) {
            const int alt_c   = ac[a];
            const int pair_ref_c = multi ? (site_total - alt_c) : ref_c;
            const int total_c = multi ? site_total
                                      : (opts.af_vs_site_depth ? site_total : ref_c + alt_c);
            const double af   = total_c > 0 ? static_cast<double>(alt_c) / total_c : 0.0;
            const int orig_alt = allele_orig_idx[i][a];

            const std::string pair_id = (ac.size() > 2)
                ? out.site_ids[i] + ":" + std::to_string(orig_alt)
                : out.site_ids[i];

            std::string drop_reason;
            if (total_c < opts.min_depth) {
                drop_reason = "low_depth";
            } else if (af > opts.max_af) {
                drop_reason = "high_af";
            } else if (af < opts.min_af) {
                drop_reason = "low_af";
            }
            if (!drop_reason.empty()) {
                out.filtered_sites.push_back({pair_id, meta_chrom, meta_pos,
                                              pair_ref_c, alt_c, total_c, af,
                                              std::move(drop_reason)});
                continue;
            }

            int fwd_ref = fc.empty() ? 0 : fc[0];
            int rev_ref = rc_s.empty() ? 0 : rc_s[0];
            if (multi) {
                fwd_ref = 0; rev_ref = 0;
                for (size_t b = 0; b < fc.size(); ++b) if (b != a) fwd_ref += fc[b];
                for (size_t b = 0; b < rc_s.size(); ++b) if (b != a) rev_ref += rc_s[b];
            }
            const int fwd_alt = a < fc.size() ? fc[a] : 0;
            const int rev_alt = a < rc_s.size() ? rc_s[a] : 0;

            // Compute the minimal-VCF identity (pos, ref, alt) of this pair so
            // duplicates from overlapping snarls share a key. Identity views point
            // into meta.ref / meta.alts[alt_idx], which outlive this loop.
            const GraphSiteMeta& meta = out.site_meta[i];
            const size_t alt_idx = static_cast<size_t>(orig_alt - 1);
            const std::string& id_alt_src =
                alt_idx < meta.alts.size() ? meta.alts[alt_idx] : meta.ref;
            const MinimalVcfId identity = minimal_vcf_id(meta.pos, meta.ref, id_alt_src);
            const uint64_t fp = minimal_vcf_fingerprint(identity);

            int canon = -1;
            std::vector<CanonEntry>& bucket = canonical_by_fp[fp];
            for (const CanonEntry& e : bucket) {
                if (e.old_i != i && e.id == identity) { canon = e.new_idx; break; }
            }
            if (canon >= 0) {
                // Duplicate variant from an overlapping snarl: pool counts into the
                // canonical candidate and route this site's observations there. Do
                // not create a new k-means anchor.
                // Route this snarl's observations to the canonical anchor so the
                // variant is not double-counted in k-means, but keep the canonical
                // candidate's original counts: pooling coverage across overlapping
                // snarls over-weights these anchors and net-regresses Hamming. The
                // per-read observation dedup in Phase 3 still collapses a read seen
                // at both snarls into one vote.
                const size_t canon_idx = static_cast<size_t>(canon);
                old_site_to_pairs[i].push_back(
                    {static_cast<int>(canon_idx), static_cast<int>(a)});
                continue;
            }

            const int new_idx = static_cast<int>(new_cands.size());
            bucket.push_back({identity, new_idx, i});
            old_site_to_pairs[i].push_back({new_idx, static_cast<int>(a)});

            CandidateVariant pair_cand = out.chunk.candidates[i];

            // Derive variant type from VCF REF/ALT sequence lengths.
            {
                const size_t ref_len = meta.ref.size();
                const size_t alt_len = alt_idx < meta.alts.size() ? meta.alts[alt_idx].size() : ref_len;
                if (ref_len < alt_len) {
                    pair_cand.key.type = VariantType::Insertion;
                    pair_cand.key.ref_len = 1;
                } else if (ref_len > alt_len) {
                    pair_cand.key.type = VariantType::Deletion;
                    pair_cand.key.ref_len = static_cast<int>(ref_len);
                } else {
                    pair_cand.key.type = VariantType::Snp;
                    pair_cand.key.ref_len = static_cast<int>(ref_len);
                }
            }


            pair_cand.counts.alle_covs       = {pair_ref_c, alt_c};
            pair_cand.counts.n_uniq_alles    = 2;
            pair_cand.counts.ref_cov         = pair_ref_c;
            pair_cand.counts.alt_cov         = alt_c;
            pair_cand.counts.total_cov       = total_c;
            pair_cand.counts.forward_ref     = fwd_ref;
            pair_cand.counts.reverse_ref     = rev_ref;
            pair_cand.counts.forward_alt     = fwd_alt;
            pair_cand.counts.reverse_alt     = rev_alt;
            pair_cand.counts.allele_fraction = af;
            new_cands.push_back(pair_cand);
            new_ids.push_back(pair_id);
            new_meta.push_back(out.site_meta[i]);
            new_allele_counts.push_back({ref_c, alt_c});
            new_fwd_strand.push_back({fwd_ref, fwd_alt});
            new_rev_strand.push_back({rev_ref, rev_alt});
            new_orig_idx.push_back({0, orig_alt});
        }
    }

    out.chunk.candidates       = std::move(new_cands);
    out.site_ids               = std::move(new_ids);
    out.site_meta              = std::move(new_meta);
    allele_counts              = std::move(new_allele_counts);
    fwd_strand_counts          = std::move(new_fwd_strand);
    rev_strand_counts          = std::move(new_rev_strand);
    allele_orig_idx            = std::move(new_orig_idx);

    // Classify candidates before building read profiles so that pruned sites
    // (LowCoverage, StrandBias) never enter the profile / cr_overlap index.
    classify_graph_candidates(out.chunk, opts);

    // Build a fast lookup for pruned candidates.
    const size_t n_final_cands = out.chunk.candidates.size();
    std::vector<bool> cand_pruned(n_final_cands, false);
    for (size_t ci = 0; ci < n_final_cands; ++ci) {
        const VariantCategory cat = out.chunk.candidates[ci].counts.category;
        if (cat == VariantCategory::LowCoverage || cat == VariantCategory::StrandBias ||
            cat == VariantCategory::NonVariant) {
            cand_pruned[ci] = true;
        }
    }

    // Phase 3: remap read observations to the new biallelic pair space.
    // Allele 0 (ref) fans out to all new pairs from its original site (allele 0 in each).
    // Allele j (alt) maps to allele 1 in the one pair that holds alt j.
    // Observations at pruned candidates are dropped.
    for (uint32_t rid = 0; rid < n_reads; ++rid) {
        auto& obs = read_obs[rid];
        std::vector<GraphProfileObservation> new_obs;
        new_obs.reserve(obs.size());
        for (const auto& o : obs) {
            if (o.site_index < 0 || static_cast<size_t>(o.site_index) >= n_cands) continue;
            const size_t old_si = static_cast<size_t>(o.site_index);
            const auto& pairs = old_site_to_pairs[old_si];
            if (pairs.empty()) continue;
            // Was the *original* snarl multi-allelic?  If so, a read carrying one
            // alt is evidence against every other alt of that snarl, so it must
            // appear in those pairs as allele 0 rather than be omitted.  Omitting
            // it leaves each pair seeing only its own carriers -- monomorphic, and
            // useless to k-means -- which is how alt-vs-alt heterozygotes became
            // invisible even after their allele fractions were corrected.
            const bool multi_src =
                old_si < source_is_multi.size() && source_is_multi[old_si];
            // Fast path: biallelic site (one surviving pair) — no loop needed.
            if (pairs.size() == 1) {
                if (cand_pruned[static_cast<size_t>(pairs[0].new_idx)]) continue;
                if (pairs[0].old_alt_phase1 < 0) {
                    // Whole multi-allelic snarl: carry the allele index through.
                    if (o.allele >= 0 &&
                        static_cast<size_t>(o.allele) < allele_remap[old_si].size()) {
                        const int mapped = allele_remap[old_si][static_cast<size_t>(o.allele)];
                        if (mapped >= 0) new_obs.push_back({pairs[0].new_idx, mapped});
                    }
                    continue;
                }
                if (o.allele == 0) {
                    new_obs.push_back({pairs[0].new_idx, 0});
                } else if (o.allele > 0 &&
                           static_cast<size_t>(o.allele) < allele_remap[old_si].size()) {
                    if (allele_remap[old_si][static_cast<size_t>(o.allele)] == pairs[0].old_alt_phase1)
                        new_obs.push_back({pairs[0].new_idx, 1});
                    else if (multi_src)
                        new_obs.push_back({pairs[0].new_idx, 0});
                }
                continue;
            }
            // General path: multiallelic site with multiple surviving pairs.
            if (o.allele == 0) {
                for (const NewPairEntry& pe : pairs)
                    if (!cand_pruned[static_cast<size_t>(pe.new_idx)])
                        new_obs.push_back({pe.new_idx, 0});
            } else if (o.allele > 0 &&
                       static_cast<size_t>(o.allele) < allele_remap[old_si].size()) {
                const int phase1_alt = allele_remap[old_si][static_cast<size_t>(o.allele)];
                if (phase1_alt >= 0) {
                    for (const NewPairEntry& pe : pairs) {
                        if (cand_pruned[static_cast<size_t>(pe.new_idx)]) continue;
                        if (pe.old_alt_phase1 == phase1_alt) new_obs.push_back({pe.new_idx, 1});
                        else if (multi_src) new_obs.push_back({pe.new_idx, 0});
                    }
                }
            }
        }
        obs = std::move(new_obs);
    }

    // Build sorted read order by name for deterministic output.
    std::vector<uint32_t> read_order(n_reads);
    std::iota(read_order.begin(), read_order.end(), 0);
    std::sort(read_order.begin(), read_order.end(),
              [&](uint32_t a, uint32_t b) { return read_id_to_name[a] < read_id_to_name[b]; });

    for (uint32_t rid : read_order) {
        auto& observations = read_obs[rid];
        if (observations.empty()) continue;
        std::sort(observations.begin(), observations.end(),
                  [](const GraphProfileObservation& lhs, const GraphProfileObservation& rhs) {
                      if (lhs.site_index != rhs.site_index) return lhs.site_index < rhs.site_index;
                      return lhs.allele < rhs.allele;
                  });
        std::vector<GraphProfileObservation> dedup;
        dedup.reserve(observations.size());
        for (const GraphProfileObservation& obs : observations) {
            if (!dedup.empty() && dedup.back().site_index == obs.site_index) {
                if (dedup.back().allele != obs.allele) dedup.back().allele = -1;
                continue;
            }
            dedup.push_back(obs);
        }
        dedup.erase(std::remove_if(dedup.begin(), dedup.end(),
                                   [](const GraphProfileObservation& obs) { return obs.allele < 0; }),
                    dedup.end());
        const int mapq = read_max_mapq[rid];
        add_read_profile(out, read_id_to_name[rid], dedup, mapq);
    }

    out.chunk.haps.assign(out.chunk.reads.size(), 0);
    out.chunk.phase_sets.assign(out.chunk.reads.size(), kUnphasedReadPhaseSet);
    rebuild_read_var_cr(out.chunk);
    out.site_allele_orig_idx = std::move(allele_orig_idx);

    return out;
}

// Merge-intersect overlap detection between adjacent chunks.
// Reads are inserted in sorted name order by build_graph_chunk,
// so we can merge-intersect directly in O(n+m) without sorting.
static void populate_graph_chunk_pair_overlap_impl(PhasingChunk& pre, PhasingChunk& cur) {
    pre.down_ovlp_read_i.assign(1, {});
    cur.up_ovlp_read_i.assign(1, {});
    pre.n_down_ovlp_skip_reads.assign(1, 0);
    cur.n_up_ovlp_skip_reads.assign(1, 0);

    size_t pi = 0, ci = 0;
    while (pi < pre.reads.size() && ci < cur.reads.size()) {
        const int cmp = pre.reads[pi].qname.compare(cur.reads[ci].qname);
        if (cmp < 0) { ++pi; }
        else if (cmp > 0) { ++ci; }
        else {
            pre.down_ovlp_read_i[0].push_back(static_cast<int>(pi));
            cur.up_ovlp_read_i[0].push_back(static_cast<int>(ci));
            ++pi; ++ci;
        }
    }
}

void populate_graph_chunk_overlaps(std::vector<GraphChunkBuildResult>& graph_chunks) {
    for (size_t i = 1; i < graph_chunks.size(); ++i)
        populate_graph_chunk_pair_overlap_impl(graph_chunks[i - 1].chunk, graph_chunks[i].chunk);
}

void phase_graph_chunks(std::vector<GraphChunkBuildResult>& graph_chunks,
                            const Options& opts) {
    for (GraphChunkBuildResult& graph_chunk : graph_chunks) {
        assign_hap_based_on_germline_het_vars_kmeans(graph_chunk.chunk, opts, kCandGermlineClean);
    }
    populate_graph_chunk_overlaps(graph_chunks);
    std::vector<PhasingChunk> chunks;
    chunks.reserve(graph_chunks.size());
    for (GraphChunkBuildResult& graph_chunk : graph_chunks) {
        chunks.push_back(std::move(graph_chunk.chunk));
    }
    stitch_chunk_haps(chunks, &opts, nullptr);
    for (size_t i = 0; i < chunks.size(); ++i) {
        graph_chunks[i].chunk = std::move(chunks[i]);
        rescue_unphased_graph_reads(graph_chunks[i].chunk);
    }
}

namespace {

constexpr double kRescueSiteOrientationPValue = 0.01;
constexpr double kIndependentBlockAssociationPValue = 0.01;
constexpr double kIndependentBlockMaxDiscordance = 0.10;
constexpr double kRescueSingletonMaxDiscordance = 0.15;
constexpr double kOneSided95PercentZ = 1.6448536269514722;

struct RescueSiteVote {
    // [read haplotype - 1][observed allele]
    std::array<std::array<int, 2>, 2> counts{};
    // Read-only rescue assignments may extend a two-locus chain, but they
    // cannot establish the stronger evidence needed for a singleton rescue.
    std::array<std::array<int, 2>, 2> primary_counts{};
};

struct RescueMarker {
    hts_pos_t phase_set = kUnphasedReadPhaseSet;
    std::array<int, 2> allele_to_hap{};
    bool is_snp = false;
    bool is_direct = false;
    bool singleton_safe = false;
};

// Exact one-sided P(X >= successes), X ~ Binomial(trials, 0.5). Sites enter
// rescue only when their allele/haplotype association is unlikely under an
// unlinked null, avoiding a fixed read-count threshold that changes meaning
// with local depth.
double rescue_binomial_tail(int successes, int trials) {
    if (successes <= 0) return 1.0;
    if (successes > trials) return 0.0;

    const long double log_term =
        std::lgamma(static_cast<long double>(trials + 1)) -
        std::lgamma(static_cast<long double>(successes + 1)) -
        std::lgamma(static_cast<long double>(trials - successes + 1)) -
        static_cast<long double>(trials) * std::log(2.0L);
    long double term = std::exp(log_term);
    long double tail = term;
    for (int k = successes; k < trials; ++k) {
        term *= static_cast<long double>(trials - k) /
                static_cast<long double>(k + 1);
        tail += term;
    }
    return static_cast<double>(std::min(1.0L, tail));
}

bool is_oriented_biallelic_candidate(const CandidateVariant& candidate) {
    const int hap1 = candidate.hap_to_cons_alle[1];
    const int hap2 = candidate.hap_to_cons_alle[2];
    return candidate.phase_set > 0 &&
           (hap1 == 0 || hap1 == 1) &&
           (hap2 == 0 || hap2 == 1) && hap1 != hap2;
}

double one_sided_wilson_upper_bound(int discordant, int total) {
    if (total <= 0 || discordant < 0 || discordant > total) return 1.0;
    const double n = static_cast<double>(total);
    const double rate = static_cast<double>(discordant) / n;
    const double z2 = kOneSided95PercentZ * kOneSided95PercentZ;
    const double center = rate + z2 / (2.0 * n);
    const double spread = kOneSided95PercentZ * std::sqrt(
        rate * (1.0 - rate) / n + z2 / (4.0 * n * n));
    return (center + spread) / (1.0 + z2 / n);
}

}  // namespace

bool independent_bam_block_is_supported(
        const std::vector<IndependentBamBlockLink>& links) {
    int total = 0;
    int discordant = 0;
    size_t informative_links = 0;

    for (const IndependentBamBlockLink& link : links) {
        const auto& counts = link.counts;
        const int bam_hap1 = counts[0][0] + counts[0][1];
        const int bam_hap2 = counts[1][0] + counts[1][1];
        const int graph_hap1 = counts[0][0] + counts[1][0];
        const int graph_hap2 = counts[0][1] + counts[1][1];
        if (bam_hap1 == 0 || bam_hap2 == 0 ||
            graph_hap1 == 0 || graph_hap2 == 0) {
            continue;
        }

        const int same = counts[0][0] + counts[1][1];
        const int cross = counts[0][1] + counts[1][0];
        const int link_total = same + cross;
        const int link_concordant = std::max(same, cross);
        const double association_p = std::min(
            1.0, 2.0 * rescue_binomial_tail(link_concordant, link_total));
        if (association_p > kIndependentBlockAssociationPValue) return false;

        total += link_total;
        discordant += std::min(same, cross);
        ++informative_links;
    }
    if (informative_links == 0 || total == 0) return false;

    // A p-value against random association alone becomes permissive at high
    // depth. The one-sided Wilson bound also limits the plausible block-wide
    // discordance while accounting for the amount of supporting evidence.
    return one_sided_wilson_upper_bound(discordant, total) <=
           kIndependentBlockMaxDiscordance;
}

size_t apply_independent_bam_read_blocks(PhasingChunk& chunk) {
    chunk.haps.resize(chunk.reads.size(), 0);
    chunk.phase_sets.resize(chunk.reads.size(), kUnphasedReadPhaseSet);
    chunk.gap_haps.resize(chunk.reads.size(), 0);
    chunk.gap_phase_sets.resize(chunk.reads.size(), kUnphasedReadPhaseSet);

    const size_t n = std::min(chunk.reads.size(),
                              chunk.bam_fallback_haps.size());
    size_t applied = 0;
    for (size_t read_i = 0; read_i < n; ++read_i) {
        if ((chunk.haps[read_i] == 1 || chunk.haps[read_i] == 2) &&
            chunk.phase_sets[read_i] > 0) {
            continue;
        }
        if ((chunk.gap_haps[read_i] == 1 || chunk.gap_haps[read_i] == 2) &&
            chunk.gap_phase_sets[read_i] > 0) {
            continue;
        }
        if ((chunk.bam_fallback_haps[read_i] != 1 &&
             chunk.bam_fallback_haps[read_i] != 2) ||
            read_i >= chunk.bam_fallback_phase_sets.size() ||
            chunk.bam_fallback_phase_sets[read_i] <= 0) {
            continue;
        }
        chunk.gap_haps[read_i] = chunk.bam_fallback_haps[read_i];
        chunk.gap_phase_sets[read_i] =
            chunk.bam_fallback_phase_sets[read_i];
        ++applied;
    }
    return applied;
}

static int rescue_observed_allele(const ReadVariantProfile& profile,
                                  size_t offset,
                                  bool use_bam_observations) {
    if (offset < profile.alleles.size() &&
        (profile.alleles[offset] == 0 || profile.alleles[offset] == 1)) {
        return profile.alleles[offset];
    }
    if (use_bam_observations && offset < profile.bam_alleles.size() &&
        (profile.bam_alleles[offset] == 0 ||
         profile.bam_alleles[offset] == 1)) {
        return profile.bam_alleles[offset];
    }
    return -1;
}

static size_t rescue_unphased_graph_read_layer(
        PhasingChunk& chunk, bool use_bam_observations) {
    if (chunk.candidates.empty() || chunk.read_var_profile.empty()) return 0;

    chunk.haps.resize(chunk.reads.size(), 0);
    chunk.phase_sets.resize(chunk.reads.size(), kUnphasedReadPhaseSet);
    chunk.gap_haps.resize(chunk.reads.size(), 0);
    chunk.gap_phase_sets.resize(chunk.reads.size(), kUnphasedReadPhaseSet);
    chunk.gap_from_bam_observation.resize(chunk.reads.size(), false);

    // Each unphased candidate normally overlaps only a few phase sets, so a
    // sparse map avoids allocating candidates * phase_sets dense counters.
    using VotesByPhaseSet = std::unordered_map<hts_pos_t, RescueSiteVote>;
    std::vector<VotesByPhaseSet> site_votes(chunk.candidates.size());
    for (const ReadVariantProfile& profile : chunk.read_var_profile) {
        if (profile.read_id < 0 ||
            static_cast<size_t>(profile.read_id) >= chunk.reads.size()) {
            continue;
        }
        const size_t read_i = static_cast<size_t>(profile.read_id);
        int hap = chunk.haps[read_i];
        hts_pos_t phase_set = chunk.phase_sets[read_i];
        const bool primary_assignment =
            (hap == 1 || hap == 2) && phase_set > 0;
        if ((hap != 1 && hap != 2) &&
            (chunk.gap_haps[read_i] == 1 || chunk.gap_haps[read_i] == 2)) {
            hap = chunk.gap_haps[read_i];
            phase_set = chunk.gap_phase_sets[read_i] - kGapFillPsOffset;
        }
        if ((hap != 1 && hap != 2) || phase_set <= 0) continue;

        for (int site_i = profile.start_var_idx;
             site_i <= profile.end_var_idx; ++site_i) {
            if (site_i < 0 ||
                static_cast<size_t>(site_i) >= chunk.candidates.size()) {
                continue;
            }
            const size_t offset = static_cast<size_t>(
                site_i - profile.start_var_idx);
            const int allele = rescue_observed_allele(
                profile, offset, use_bam_observations);
            if (allele < 0) continue;
            RescueSiteVote& vote =
                site_votes[static_cast<size_t>(site_i)][phase_set];
            ++vote.counts[static_cast<size_t>(hap - 1)]
                         [static_cast<size_t>(allele)];
            if (primary_assignment) {
                ++vote.primary_counts[static_cast<size_t>(hap - 1)]
                                     [static_cast<size_t>(allele)];
            }
        }
    }

    std::vector<RescueMarker> markers(chunk.candidates.size());
    for (size_t site_i = 0; site_i < chunk.candidates.size(); ++site_i) {
        const CandidateVariant& candidate = chunk.candidates[site_i];
        RescueMarker& marker = markers[site_i];
        marker.is_snp = candidate.key.type == VariantType::Snp;

        if (is_oriented_biallelic_candidate(candidate)) {
            marker.phase_set = candidate.phase_set;
            marker.is_direct = true;
            marker.allele_to_hap[static_cast<size_t>(
                candidate.hap_to_cons_alle[1])] = 1;
            marker.allele_to_hap[static_cast<size_t>(
                candidate.hap_to_cons_alle[2])] = 2;
            continue;
        }

        // An excluded site is usable only when exactly one established phase
        // set gives it a significant orientation. Evidence from different
        // blocks is never pooled, because their numeric HP labels need not use
        // the same gauge.
        RescueMarker supported;
        int supported_phase_sets = 0;
        for (const auto& entry : site_votes[site_i]) {
            const auto& counts = entry.second.counts;
            const int hap1 = counts[0][0] + counts[0][1];
            const int hap2 = counts[1][0] + counts[1][1];
            const int allele0 = counts[0][0] + counts[1][0];
            const int allele1 = counts[0][1] + counts[1][1];
            if (hap1 == 0 || hap2 == 0 || allele0 == 0 || allele1 == 0) {
                continue;
            }

            const int same = counts[0][0] + counts[1][1];
            const int cross = counts[0][1] + counts[1][0];
            if (same == cross) continue;
            const int concordant = std::max(same, cross);
            const int total = same + cross;
            if (rescue_binomial_tail(concordant, total) >
                kRescueSiteOrientationPValue) {
                continue;
            }

            ++supported_phase_sets;
            supported.phase_set = entry.first;
            supported.is_snp = marker.is_snp;
            supported.allele_to_hap =
                same > cross ? std::array<int, 2>{1, 2}
                             : std::array<int, 2>{2, 1};

            // A single inferred locus may tag a read only when primary graph
            // assignments independently establish both its association and a
            // low discordance rate. This is deliberately stronger than the
            // test used when two independent inferred loci agree.
            const auto& primary = entry.second.primary_counts;
            const int primary_hap1 = primary[0][0] + primary[0][1];
            const int primary_hap2 = primary[1][0] + primary[1][1];
            const int primary_allele0 = primary[0][0] + primary[1][0];
            const int primary_allele1 = primary[0][1] + primary[1][1];
            const int primary_same = primary[0][0] + primary[1][1];
            const int primary_cross = primary[0][1] + primary[1][0];
            const int primary_total = primary_same + primary_cross;
            const int primary_concordant =
                std::max(primary_same, primary_cross);
            const bool same_orientation =
                primary_same != primary_cross &&
                ((same > cross) == (primary_same > primary_cross));
            supported.singleton_safe =
                primary_hap1 > 0 && primary_hap2 > 0 &&
                primary_allele0 > 0 && primary_allele1 > 0 &&
                same_orientation &&
                rescue_binomial_tail(primary_concordant, primary_total) <=
                    kRescueSiteOrientationPValue &&
                one_sided_wilson_upper_bound(
                    std::min(primary_same, primary_cross), primary_total) <=
                    kRescueSingletonMaxDiscordance;
        }
        if (supported_phase_sets == 1) marker = supported;
    }

    struct PhaseSetReadScore {
        std::array<int, 2> snp{};
        std::array<int, 2> indel{};
        std::array<int, 2> direct_snp{};
        std::array<int, 2> direct_indel{};
        std::array<int, 2> singleton_safe_snp{};
        std::array<int, 2> singleton_safe_indel{};
    };

    struct LocusVote {
        int hap = 0;
        bool is_snp = false;
        bool is_direct = false;
        bool singleton_safe = false;
    };

    size_t rescued = 0;
    for (const ReadVariantProfile& profile : chunk.read_var_profile) {
        if (profile.read_id < 0 ||
            static_cast<size_t>(profile.read_id) >= chunk.reads.size()) {
            continue;
        }
        const size_t read_i = static_cast<size_t>(profile.read_id);
        if (chunk.haps[read_i] == 1 || chunk.haps[read_i] == 2 ||
            chunk.gap_haps[read_i] == 1 || chunk.gap_haps[read_i] == 2) {
            continue;
        }
        // A decisive whole-chunk BAM assignment already carries more than one
        // clean-site equivalent of evidence. The augmented singleton pass is a
        // fill path and must not replace that staged assignment.
        if (use_bam_observations &&
            read_i < chunk.bam_fallback_haps.size() &&
            read_i < chunk.bam_fallback_phase_sets.size() &&
            (chunk.bam_fallback_haps[read_i] == 1 ||
             chunk.bam_fallback_haps[read_i] == 2) &&
            chunk.bam_fallback_phase_sets[read_i] > 0) {
            continue;
        }

        // Collapse co-located candidate rows before scoring. Split or nested
        // representations at one coordinate provide one observation; if they
        // disagree on the proposed haplotype, that locus contributes nothing.
        using LocusKey = std::pair<hts_pos_t, hts_pos_t>;  // PS, position
        std::map<LocusKey, LocusVote> locus_votes;
        for (int site_i = profile.start_var_idx;
             site_i <= profile.end_var_idx; ++site_i) {
            if (site_i < 0 ||
                static_cast<size_t>(site_i) >= markers.size()) {
                continue;
            }
            const size_t offset = static_cast<size_t>(
                site_i - profile.start_var_idx);
            const int allele = rescue_observed_allele(
                profile, offset, use_bam_observations);
            if (allele < 0) continue;

            const RescueMarker& marker = markers[static_cast<size_t>(site_i)];
            if (marker.phase_set <= 0) continue;
            const int proposed_hap =
                marker.allele_to_hap[static_cast<size_t>(allele)];
            if (proposed_hap != 1 && proposed_hap != 2) continue;

            const hts_pos_t pos =
                chunk.candidates[static_cast<size_t>(site_i)].key.sort_pos();
            const LocusKey key{marker.phase_set, pos};
            auto [it, inserted] = locus_votes.emplace(
                key, LocusVote{proposed_hap, marker.is_snp,
                               marker.is_direct, marker.singleton_safe});
            if (!inserted && it->second.hap != proposed_hap) {
                it->second.hap = 0;
            } else if (!inserted) {
                // A direct SNP is the strongest form when equivalent rows agree.
                it->second.is_snp = it->second.is_snp || marker.is_snp;
                it->second.is_direct =
                    it->second.is_direct || marker.is_direct;
                it->second.singleton_safe =
                    it->second.singleton_safe || marker.singleton_safe;
            }
        }

        std::unordered_map<hts_pos_t, PhaseSetReadScore> scores;
        for (const auto& entry : locus_votes) {
            const LocusVote& vote = entry.second;
            if (vote.hap != 1 && vote.hap != 2) continue;
            PhaseSetReadScore& score = scores[entry.first.first];
            std::array<int, 2>& tier = vote.is_snp ? score.snp : score.indel;
            std::array<int, 2>& direct =
                vote.is_snp ? score.direct_snp : score.direct_indel;
            std::array<int, 2>& singleton_safe =
                vote.is_snp ? score.singleton_safe_snp
                            : score.singleton_safe_indel;
            const size_t hap_i = static_cast<size_t>(vote.hap - 1);
            ++tier[hap_i];
            if (vote.is_direct) ++direct[hap_i];
            if (vote.singleton_safe) ++singleton_safe[hap_i];
        }

        hts_pos_t best_phase_set = kUnphasedReadPhaseSet;
        int best_hap = 0;
        int best_margin = 0;
        int best_total = 0;
        int best_direct = 0;
        int best_singleton_safe = 0;
        bool tied = false;
        for (const auto& entry : scores) {
            const PhaseSetReadScore& score = entry.second;
            const bool use_snps = score.snp[0] + score.snp[1] > 0;
            const std::array<int, 2>& tier = use_snps ? score.snp
                                                      : score.indel;
            const std::array<int, 2>& direct =
                use_snps ? score.direct_snp : score.direct_indel;
            const std::array<int, 2>& singleton_safe =
                use_snps ? score.singleton_safe_snp
                         : score.singleton_safe_indel;
            if (tier[0] == tier[1]) continue;
            const int margin = std::abs(tier[0] - tier[1]);
            const int total = tier[0] + tier[1];
            const int hap = tier[0] > tier[1] ? 1 : 2;
            if (margin > best_margin ||
                (margin == best_margin && total > best_total)) {
                best_phase_set = entry.first;
                best_hap = hap;
                best_margin = margin;
                best_total = total;
                best_direct = direct[0] + direct[1];
                best_singleton_safe =
                    singleton_safe[static_cast<size_t>(hap - 1)];
                tied = false;
            } else if (margin == best_margin && total == best_total) {
                tied = true;
            }
        }

        // One allele is enough when the ordinary solve directly phased the
        // site, or when primary graph assignments gave an inferred site the
        // stronger statistical guarantees above.
        if (best_hap == 0 || tied ||
            (best_total == 1 && best_direct == 0 &&
             best_singleton_safe == 0)) {
            continue;
        }
        chunk.gap_haps[read_i] = best_hap;
        chunk.gap_phase_sets[read_i] = best_phase_set + kGapFillPsOffset;
        chunk.gap_from_bam_observation[read_i] =
            use_bam_observations;
        ++rescued;
    }
    return rescued;
}

size_t rescue_unphased_graph_reads(PhasingChunk& chunk) {
    size_t total_rescued = 0;
    const auto run_to_fixed_point = [&](bool use_bam_observations) {
        size_t rescued = 0;
        while (true) {
            // Grow from established blocks toward the middle of a gap. Every
            // next layer independently passes the site-orientation test.
            const size_t layer_rescued =
                rescue_unphased_graph_read_layer(
                    chunk, use_bam_observations);
            if (layer_rescued == 0) break;
            rescued += layer_rescued;
        }
        return rescued;
    };

    // Preserve every assignment supported by the original graph profiles.
    // Only after that fixed point is stable may exact BAM observations fill
    // missing alleles for reads that are still unphased.
    total_rescued += run_to_fixed_point(false);
    total_rescued += run_to_fixed_point(true);
    return total_rescued;
}

void merge_graph_chunk_into_read_rows(
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const GraphChunkBuildResult& gc,
    int min_read_hap_margin) {
    const PhasingChunk& chunk = gc.chunk;
    for (const ReadVariantProfile& profile : chunk.read_var_profile) {
        const size_t read_i = static_cast<size_t>(profile.read_id);
        const ReadRecord& read = chunk.reads[read_i];
        PhaseReadOutputRow& row = rows_by_read[read.qname];
        if (row.read_name.empty()) row.read_name = read.qname;

        int hap = read_i < chunk.haps.size() ? chunk.haps[read_i] : 0;
        hts_pos_t phase_set =
            read_i < chunk.phase_sets.size()
                ? chunk.phase_sets[read_i]
                : kUnphasedReadPhaseSet;
        bool is_primary = (hap == 1 || hap == 2) && phase_set > 0;
        if (!is_primary && read_i < chunk.gap_haps.size() &&
            read_i < chunk.gap_phase_sets.size() &&
            (chunk.gap_haps[read_i] == 1 || chunk.gap_haps[read_i] == 2) &&
            chunk.gap_phase_sets[read_i] > 0) {
            hap = chunk.gap_haps[read_i];
            phase_set = chunk.gap_phase_sets[read_i];
        }
        // Read-confidence gate.  init_assign_read_hap_based_on_cons_alle commits a read to a
        // haplotype on any non-zero score, with no minimum evidence: a read
        // agreeing with a single clean het SNP and contradicting none is
        // assigned as confidently as one agreeing with thirty.  On HG002 chr20
        // those one-SNP-margin reads are 11% of reads and 69% of all phasing
        // errors, at an 18x elevated error rate.  Leaving them unphased is
        // better than guessing.
        if (min_read_hap_margin > 0 &&
            read.n_clean_agree_snps - read.n_clean_conflict_snps < min_read_hap_margin) {
            hap = 0;
        }

        for (int site_i = profile.start_var_idx; site_i <= profile.end_var_idx; ++site_i) {
            const int offset = site_i - profile.start_var_idx;
            if (offset < 0 || static_cast<size_t>(offset) >= profile.alleles.size()) continue;
            const int allele = profile.alleles[static_cast<size_t>(offset)];
            if (allele < 0) continue;
            if (site_i < 0 || static_cast<size_t>(site_i) >= gc.site_ids.size()) continue;
            merge_phase_read_observation(row,
                                         gc.site_ids[static_cast<size_t>(site_i)],
                                         allele);
        }
        const bool fill_only =
            !is_primary &&
            read_i < chunk.gap_from_bam_observation.size() &&
            chunk.gap_from_bam_observation[read_i];
        merge_phase_read_assignment(row, chunk.region.chunk_id, hap, phase_set,
                                    is_primary, fill_only);
    }

    // Reads without GAF catalog observations have no ReadVariantProfile and
    // therefore never enter the loop above. Their independent BAM assignment
    // is output-only; a primary graph assignment from an overlapping chunk
    // still takes precedence in merge_phase_read_assignment.
    for (const ReadPhaseAssignment& assignment :
         chunk.bam_output_fallback_reads) {
        PhaseReadOutputRow& row = rows_by_read[assignment.qname];
        if (row.read_name.empty()) row.read_name = assignment.qname;
        merge_phase_read_assignment(row, chunk.region.chunk_id, assignment.hap,
                                    assignment.phase_set, false, false);
    }
}

void write_graph_phase_sites_tsv_header(std::ostream& out) {
    out << "CHUNK_ID\tSITE_INDEX\tSITE_ID\tPOS\tN_ALLELES\tDEPTH\tALLELE_COUNTS\tPHASE_SET\tHAP1_ALLELE\tHAP2_ALLELE\n";
}

void write_graph_phase_sites_tsv_rows(std::ostream& out,
                                          const GraphChunkBuildResult& gc) {
    for (size_t i = 0; i < gc.chunk.candidates.size(); ++i) {
        const CandidateVariant& candidate = gc.chunk.candidates[i];
        out << gc.chunk.region.chunk_id << '\t'
            << i << '\t'
            << gc.site_ids[i] << '\t'
            << candidate.key.pos << '\t'
            << candidate.counts.n_uniq_alles << '\t'
            << candidate.counts.total_cov << '\t';
        for (size_t allele = 0; allele < candidate.counts.alle_covs.size(); ++allele) {
            if (allele > 0) out << ',';
            out << candidate.counts.alle_covs[allele];
        }
        out << '\t'
            << candidate.phase_set << '\t'
            << candidate.hap_to_cons_alle[1] << '\t'
            << candidate.hap_to_cons_alle[2] << '\n';
    }
}

void write_graph_phase_sites_tsv(std::ostream& out,
                                     const std::vector<GraphChunkBuildResult>& graph_chunks) {
    write_graph_phase_sites_tsv_header(out);
    for (const auto& gc : graph_chunks)
        write_graph_phase_sites_tsv_rows(out, gc);
}

// Write a minimal unmapped BAM record carrying HP and PS aux tags.
// Used by flush_graph_phase_bam_after_merge to emit phased read assignments.
static void write_bam_record(samFile* out_sam, sam_hdr_t* hdr,
                             const std::string& name, int hap, hts_pos_t ps) {
    bam1_t* rec = bam_init1();
    if (!rec) throw std::runtime_error("bam_init1 failed");
    struct RecGuard { bam1_t* r; ~RecGuard() { bam_destroy1(r); } } rg{rec};

    if (bam_set_qname(rec, name.c_str()) < 0)
        throw std::runtime_error("bam_set_qname failed for: " + name);
    rec->core.flag = BAM_FUNMAP;
    rec->core.tid  = -1;
    rec->core.pos  = -1;
    rec->core.mtid = -1;
    rec->core.mpos = -1;
    rec->core.qual = 255;

    if (hap > 0) {
        const int32_t h = static_cast<int32_t>(hap);
        bam_aux_append(rec, "HP", 'i', sizeof(int32_t),
                       reinterpret_cast<const uint8_t*>(&h));
    }
    if (ps >= 0) {
        const int32_t p = static_cast<int32_t>(ps);
        bam_aux_append(rec, "PS", 'i', sizeof(int32_t),
                       reinterpret_cast<const uint8_t*>(&p));
    }
    if (sam_write1(out_sam, hdr, rec) < 0)
        throw std::runtime_error("failed to write BAM record for: " + name);
}

void flush_graph_phase_bam_after_merge(
    samFile* phase_bam_out,
    sam_hdr_t* phase_bam_hdr,
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const std::unordered_set<std::string>* next_chunk_qnames,
    std::unordered_set<std::string>& emitted_read_names) {
    std::vector<std::pair<std::string, PhaseReadOutputRow>> drained;
    drained.reserve(rows_by_read.size());
    for (auto it = rows_by_read.begin(); it != rows_by_read.end(); ) {
        if (next_chunk_qnames && next_chunk_qnames->count(it->first)) {
            ++it;
            continue;
        }
        std::string key        = std::move(it->first);
        PhaseReadOutputRow row = std::move(it->second);
        it                     = rows_by_read.erase(it);
        drained.emplace_back(std::move(key), std::move(row));
    }
    std::sort(drained.begin(), drained.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    for (auto& kv : drained) {
        emitted_read_names.insert(kv.first);
        if (phase_bam_out && phase_bam_hdr) {
            const int hap =
                kv.second.has_phased_assignment ? kv.second.hap : 0;
            const hts_pos_t ps =
                kv.second.has_phased_assignment
                    ? kv.second.phase_set
                    : static_cast<hts_pos_t>(-1);
            write_bam_record(phase_bam_out, phase_bam_hdr, kv.first, hap, ps);
        }
    }
}

} // namespace pgphase_collect
