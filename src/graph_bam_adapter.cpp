#include "collect_var.hpp"
#include <map>
#include "graph_bam_adapter.hpp"

#include "collect_phase.hpp"
#include "collect_phase_noisy.hpp"
#include "fisher_exact.hpp"
#include "noise_filter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
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
    candidate.graph_site = true;
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

bool core_dominates_rescue_coverage(const PhasingChunk& chunk, const hts_pos_t phase_set) {
    std::map<hts_pos_t, int> coverage;
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
        if (chunk.reads[ri].is_skipped || ri >= chunk.haps.size()) continue;
        int delta = 0;
        if ((chunk.haps[ri] == 1 || chunk.haps[ri] == 2) &&
            ri < chunk.phase_sets.size() && chunk.phase_sets[ri] == phase_set)
            delta = 1;
        else if (chunk.haps[ri] == 0 && ri < chunk.gap_haps.size() &&
                 (chunk.gap_haps[ri] == 1 || chunk.gap_haps[ri] == 2) &&
                 ri < chunk.gap_phase_sets.size() &&
                 chunk.gap_phase_sets[ri] == phase_set + kGapFillPsOffset)
            delta = -1;
        coverage[chunk.reads[ri].beg] += delta;
        coverage[chunk.reads[ri].end] -= delta;
    }
    int depth = 0;
    for (const auto& [pos, delta] : coverage) {
        (void)pos;
        depth += delta;
        if (depth < 0) return false;
    }
    return true;
}

int masked_snp_deletion_haplotype(const PhasingChunk& chunk,
                                const ReadVariantProfile& profile, size_t deletion_i) {
    if (deletion_i >= chunk.candidates.size() || profile.start_var_idx < 0 ||
        deletion_i < static_cast<size_t>(profile.start_var_idx)) return 0;
    const CandidateVariant& deletion = chunk.candidates[deletion_i];
    if (!deletion.bam_injected || !deletion.msa_verified ||
        deletion.key.type != VariantType::Deletion || deletion.key.ref_len <= 0 ||
        !deletion.key.alt.empty() || !is_phase_set_anchor(deletion) ||
        deletion.hap_to_cons_alle[1] < 0 || deletion.hap_to_cons_alle[1] > 1 ||
        deletion.hap_to_cons_alle[2] != 1 - deletion.hap_to_cons_alle[1]) return 0;
    const size_t offset = deletion_i - profile.start_var_idx;
    if (offset >= profile.alleles.size() || offset >= profile.bam_alleles.size() ||
        profile.alleles[offset] != 1 || profile.bam_alleles[offset] != 1) return 0;
    const int hap = deletion.hap_to_cons_alle[1] == 1 ? 1 : 2;
    bool masked_snp = false;
    for (size_t off = 0; off < profile.alleles.size(); ++off) {
        const size_t ci = static_cast<size_t>(profile.start_var_idx) + off;
        if (ci >= chunk.candidates.size()) break;
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (!is_phase_set_anchor(candidate)) continue;
        if (profile.alleles[off] >= 0 &&
            (candidate.phase_set != deletion.phase_set ||
             profile.alleles[off] != candidate.hap_to_cons_alle[hap])) return 0;
        if (off < profile.bam_alleles.size() && profile.bam_alleles[off] >= 0 &&
            (candidate.phase_set != deletion.phase_set ||
             profile.bam_alleles[off] != candidate.hap_to_cons_alle[hap])) return 0;
        if (candidate.msa_verified && candidate.key.type == VariantType::Snp &&
            candidate.phase_set == deletion.phase_set &&
            candidate.key.pos >= deletion.key.pos &&
            candidate.key.pos < deletion.key.pos + deletion.key.ref_len &&
            candidate.hap_to_cons_alle[hap] == 0 &&
            candidate.hap_to_cons_alle[3 - hap] == 1 &&
            profile.alleles[off] < 0 && off < profile.bam_alleles.size() &&
            profile.bam_alleles[off] < 0) masked_snp = true;
    }
    return masked_snp ? hap : 0;
}

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

bool bam_site_has_only_weak_links(const std::vector<CandidateVariant>& candidates,
                                  size_t site_index,
                                  const std::vector<hts_pos_t>& weak_cuts) {
    const CandidateVariant& site = candidates[site_index];
    if (!is_phase_set_anchor(site) || weak_cuts.empty()) return false;
    const hts_pos_t pos = site.key.sort_pos();
    hts_pos_t left = 0;
    hts_pos_t right = std::numeric_limits<hts_pos_t>::max();
    for (size_t ci = 0; ci < candidates.size(); ++ci) {
        if (ci == site_index) continue;
        const CandidateVariant& other = candidates[ci];
        if (other.phase_set != site.phase_set || !is_phase_set_anchor(other)) continue;
        const hts_pos_t other_pos = other.key.sort_pos();
        if (other_pos == pos) return false;
        if (other_pos < pos) left = std::max(left, other_pos);
        else right = std::min(right, other_pos);
    }
    const bool has_left = left > 0;
    const bool has_right = right != std::numeric_limits<hts_pos_t>::max();
    return (has_left || has_right) &&
        (!has_left || std::binary_search(weak_cuts.begin(), weak_cuts.end(), left)) &&
        (!has_right || std::binary_search(weak_cuts.begin(), weak_cuts.end(), pos));
}

bool adopt_unphased_graph_allele_from_bam(CandidateVariant& graph,
                                         const CandidateVariant& source,
                                         hts_pos_t phase_set) {
    if (!graph.graph_site || graph.bam_injected || graph.phase_set > 0 ||
        graph.counts.category != VariantCategory::RepeatHetIndel ||
        graph.counts.n_uniq_alles != 2 || phase_set <= 0 ||
        source.counts.category != VariantCategory::NoisyCandHet ||
        !source.msa_verified || source.counts.n_uniq_alles != 2 ||
        !is_phase_set_anchor(source) || source.hap_to_cons_alle[1] > 1 ||
        source.hap_to_cons_alle[2] > 1)
        return false;
    // Exact matching has already proved the allele meaning. Keep the single
    // catalog row, but carry the BAM key, genotype and proof together so its
    // original repeat demotion cannot hide the recovered VCF call. Catalog
    // walk metadata and selected ALT indices remain in the parallel arrays.
    graph = source;
    graph.graph_site = true;
    graph.bam_injected = true;
    graph.alignment_verified = true;
    graph.phase_set = phase_set;
    if (graph.counts.alle_covs.empty())
        graph.counts.alle_covs = {graph.counts.ref_cov, graph.counts.alt_cov};
    return true;
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

void reclassify_physically_validated_graph_snps(
        const GraphSiteCatalogView& catalog, GraphChunkBuildResult& out) {
    // Keep the initial solve's gauge on surviving anchors. Only a phase set
    // supported solely by an invalid SNP loses its read labels.
    std::unordered_set<std::string> homozygous_alt_sites;
    for (size_t i = 0; i < catalog.size(); ++i)
        if (catalog[i].bam_homozygous_alt)
            homozygous_alt_sites.insert(graph_site_key_str(catalog[i]));
    std::unordered_set<std::string> ref_absent_sites;
    for (size_t i = 0; i < catalog.size(); ++i)
        if (catalog[i].bam_alt_deletion_no_ref)
            ref_absent_sites.insert(graph_site_key_str(catalog[i]));
    std::unordered_set<std::string> low_fraction_sites;
    for (size_t i = 0; i < catalog.size(); ++i)
        if (catalog[i].bam_low_fraction_snp)
            low_fraction_sites.insert(graph_site_key_str(catalog[i]));
    std::unordered_set<hts_pos_t> retired_phase_sets;
    for (size_t i = 0; i < out.chunk.candidates.size(); ++i) {
        const bool ref_absent = ref_absent_sites.count(out.site_ids[i]) != 0;
        const bool low_fraction = low_fraction_sites.count(out.chunk.candidates[i].key.alt) != 0 &&
            out.chunk.candidates[i].counts.category == VariantCategory::CleanHetSnp;
        if (!ref_absent && !low_fraction && homozygous_alt_sites.count(out.site_ids[i]) == 0) continue;
        CandidateVariant& candidate = out.chunk.candidates[i];
        if (candidate.phase_set > 0) retired_phase_sets.insert(candidate.phase_set);
        if (low_fraction) {
            out.site_meta[i].bam_low_fraction_snp = true;
            candidate.counts.category = VariantCategory::LowAlleleFraction;
            candidate.lcd_var_i_to_cate = kCandNonAnchorHet;
            candidate.hap_to_cons_alle[1] = candidate.hap_to_cons_alle[2] = -1;
        } else if (ref_absent) {
            candidate.lcd_var_i_to_cate = kCandNonAnchorHet;
            candidate.hap_to_cons_alle[1] = candidate.hap_to_cons_alle[2] = -1;
            out.site_meta[i].bam_alt_deletion_no_ref = true;
        } else {
            candidate.counts.category = VariantCategory::CleanHom;
            candidate.counts.candvarcate_initial = VariantCategory::CleanHom;
            candidate.lcd_var_i_to_cate = kCandCleanHom;
            candidate.hap_to_cons_alle[1] = candidate.hap_to_cons_alle[2] = 1;
        }
        candidate.phase_set = kUnsetCandidatePhaseSet;
    }
    for (const CandidateVariant& candidate : out.chunk.candidates)
        if (is_phase_set_anchor(candidate)) retired_phase_sets.erase(candidate.phase_set);
    for (size_t ri = 0; ri < out.chunk.reads.size(); ++ri) {
        if (retired_phase_sets.count(out.chunk.phase_sets[ri]) == 0) continue;
        out.chunk.haps[ri] = 0;
        out.chunk.phase_sets[ri] = kUnphasedReadPhaseSet;
    }
}

size_t supplement_phased_snp_branches(const GraphSiteCatalogView& catalog,
        const std::vector<GraphReadAllele>& rows, GraphChunkBuildResult& graph_chunk,
        const Options& opts) {
    // Full-walk genotyping has already fixed the candidates and their gauges.
    // Other catalog alleles can still carry the exact same local SNP branch;
    // recovering those observations must not repeat genotype classification.
    PhasingChunk& chunk = graph_chunk.chunk;
    std::unordered_map<std::string, const GraphSite*> sites;
    for (size_t si = 0; si < catalog.size(); ++si)
        sites.emplace(graph_site_key_str(catalog[si]), &catalog[si]);
    struct Projection {
        size_t candidate;
        std::vector<int> alleles;
    };
    std::unordered_map<std::string, std::vector<Projection>> projections;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.counts.category != VariantCategory::CleanHetSnp ||
            !is_phase_set_anchor(candidate)) continue;
        const std::string* sequence = selected_graph_candidate_alt(graph_chunk, ci);
        if (sequence == nullptr) continue;
        const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
        const MinimalVcfId id = minimal_vcf_id(meta.pos, meta.ref, *sequence);
        if (id.ref_len != 1 || id.alt_len != 1 || *id.ref == *id.alt) continue;
        std::string sid = graph_chunk.site_ids[ci];
        auto found = sites.find(sid);
        if (found == sites.end()) {
            const size_t separator = sid.rfind(':');
            if (separator == std::string::npos) continue;
            sid.resize(separator);
            found = sites.find(sid);
        }
        if (found == sites.end()) continue;
        const GraphSite& site = *found->second;
        const size_t alt_i = static_cast<size_t>(graph_chunk.site_allele_orig_idx[ci][1]);
        if (alt_i >= site.allele_walks.size() || site.allele_walks.size() <= 2 ||
            !site.conditional_parent_alleles.empty()) continue;
        const GraphWalk& ref = site.allele_walks[0];
        const GraphWalk& alt = site.allele_walks[alt_i];
        if (ref.size() < 3 || ref.size() != alt.size()) continue;
        // An isolated node substitution between identical oriented neighbors
        // identifies the SNP without aligning or interpreting another indel.
        size_t branch = 0;
        while (branch < ref.size() && ref[branch] == alt[branch]) ++branch;
        if (branch == 0 || branch + 1 >= ref.size() ||
            !std::equal(ref.begin() + branch + 1, ref.end(), alt.begin() + branch + 1)) continue;
        Projection projection{ci, std::vector<int>(site.allele_walks.size(), -1)};
        for (size_t ai = 0; ai < site.allele_walks.size(); ++ai) {
            const GraphWalk& walk = site.allele_walks[ai];
            int occurrences = 0;
            int allele = -1;
            for (size_t wi = 1; wi + 1 < walk.size(); ++wi) {
                if (!(walk[wi - 1] == ref[branch - 1]) ||
                    !(walk[wi + 1] == ref[branch + 1])) continue;
                if (walk[wi] == ref[branch]) allele = 0;
                else if (walk[wi] == alt[branch]) allele = 1;
                else continue;
                ++occurrences;
            }
            // Repeated or opposing branches do not locate a unique SNP call.
            if (occurrences == 1) projection.alleles[ai] = allele;
        }
        if (projection.alleles[0] == 0 && projection.alleles[alt_i] == 1)
            projections[sid].push_back(std::move(projection));
    }
    if (projections.empty()) return 0;
    std::unordered_map<std::string, size_t> read_index;
    read_index.reserve(chunk.reads.size());
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri)
        read_index.emplace(chunk.reads[ri].qname, ri);
    std::vector<size_t> profile_index(chunk.reads.size(), chunk.read_var_profile.size());
    for (size_t pi = 0; pi < chunk.read_var_profile.size(); ++pi) {
        const int ri = chunk.read_var_profile[pi].read_id;
        if (ri >= 0 && static_cast<size_t>(ri) < profile_index.size()) profile_index[ri] = pi;
    }
    struct ProposedCall {
        size_t candidate;
        int original_allele;
        int allele;
    };
    std::vector<std::vector<ProposedCall>> proposed(chunk.reads.size());
    for (const GraphReadAllele& row : rows) {
        const auto source = projections.find(row.site_id);
        const auto read = read_index.find(row.read_name);
        if (source == projections.end() || read == read_index.end() ||
            row.allele < 0 || row.mapq < opts.min_mapq || chunk.reads[read->second].is_skipped ||
            profile_index[read->second] >= chunk.read_var_profile.size()) continue;
        for (const Projection& projection : source->second) {
            if (static_cast<size_t>(row.allele) >= projection.alleles.size()) continue;
            proposed[read->second].push_back({projection.candidate, row.allele,
                                            projection.alleles[row.allele]});
        }
    }
    std::vector<std::array<int, 2>> extra_counts(chunk.candidates.size(), {0, 0});
    for (size_t ri = 0; ri < proposed.size(); ++ri) {
        auto& calls = proposed[ri];
        if (calls.empty()) continue;
        std::sort(calls.begin(), calls.end(),
                  [](const ProposedCall& a, const ProposedCall& b) {
                      return a.candidate < b.candidate;
                  });
        ReadVariantProfile& profile = chunk.read_var_profile[profile_index[ri]];
        for (size_t beg = 0; beg < calls.size();) {
            size_t end = beg + 1;
            bool conflicted = false;
            while (end < calls.size() && calls[end].candidate == calls[beg].candidate) {
                conflicted |= calls[end].original_allele != calls[beg].original_allele;
                ++end;
            }
            const int ci = static_cast<int>(calls[beg].candidate);
            const int allele = calls[beg].allele;
            // Preserve the builder's duplicate-observation veto. Re-entering
            // a snarl with different walks cannot create a first-wins call.
            const bool known = profile.start_var_idx >= 0 &&
                ci >= profile.start_var_idx && ci <= profile.end_var_idx &&
                profile.alleles[static_cast<size_t>(ci - profile.start_var_idx)] >= 0;
            if (!conflicted && allele >= 0 && !known) {
                ++extra_counts[ci][allele];
            } else {
                calls[beg].allele = -1;
            }
            for (size_t duplicate = beg + 1; duplicate < end; ++duplicate)
                calls[duplicate].allele = -1;
            beg = end;
        }
    }
    std::vector<bool> eligible(chunk.candidates.size(), false);
    for (size_t ci = 0; ci < extra_counts.size(); ++ci) {
        const auto& counts = chunk.candidates[ci].counts;
        const int ref = counts.ref_cov + extra_counts[ci][0];
        const int alt = counts.alt_cov + extra_counts[ci][1];
        const double af = ref + alt > 0 ? static_cast<double>(alt) / (ref + alt) : 0.0;
        // Do not strengthen a frozen heterozygote when the recovered branch
        // evidence contradicts its genotype under the original AF filters.
        eligible[ci] = af >= opts.min_af && af <= opts.max_af;
    }
    size_t added = 0;
    for (size_t ri = 0; ri < proposed.size(); ++ri) {
        if (proposed[ri].empty()) continue;
        ReadVariantProfile& profile = chunk.read_var_profile[profile_index[ri]];
        for (const ProposedCall& call : proposed[ri]) {
            if (call.allele < 0 || !eligible[call.candidate]) continue;
            update_read_var_profile_with_allele(
                static_cast<int>(call.candidate), call.allele, -1, profile);
            ++added;
        }
    }
    if (added != 0) rebuild_read_var_cr(chunk);
    return added;
}

const std::string* selected_graph_candidate_alt(
        const GraphChunkBuildResult& graph_chunk, size_t candidate_index) {
    if (candidate_index >= graph_chunk.chunk.candidates.size() ||
        candidate_index >= graph_chunk.site_meta.size() ||
        candidate_index >= graph_chunk.site_allele_orig_idx.size() ||
        graph_chunk.chunk.candidates[candidate_index].counts.n_uniq_alles != 2) {
        return nullptr;
    }
    const GraphSiteMeta& meta = graph_chunk.site_meta[candidate_index];
    const std::vector<int>& original =
        graph_chunk.site_allele_orig_idx[candidate_index];
    if (meta.ref.empty() || original.size() != 2) return nullptr;
    const int alt_index = original[1] - 1;
    if (alt_index < 0 || static_cast<size_t>(alt_index) >= meta.alts.size())
        return nullptr;
    const std::string& alt = meta.alts[static_cast<size_t>(alt_index)];
    return alt.empty() || alt == "*" ? nullptr : &alt;
}

bool retained_source_snp_path_anchor_supported(
        const GraphChunkBuildResult& gc, size_t candidate_index) {
    constexpr int kMinQuality = 30;
    constexpr int kUnknownQuality = 255;
    constexpr size_t kMinIndependentLoci = 2;
    constexpr int kMinAlleleSupport = 2;
    constexpr double kMaxWrongParity = 0.001;
    if (candidate_index >= gc.chunk.candidates.size() ||
        candidate_index >= gc.site_meta.size()) return false;
    const CandidateVariant& candidate = gc.chunk.candidates[candidate_index];
    if (candidate.bam_injected || !is_phase_set_anchor(candidate) ||
        candidate.counts.category != VariantCategory::CleanHetSnp) return false;
    const RecoverySourceSite* source = nullptr;
    for (const RecoverySourceSite& site : gc.recovery_source_sites) {
        if (site.candidate_index != candidate_index) continue;
        if (source != nullptr) return false;
        source = &site;
    }
    if (source == nullptr || !source->can_adopt || !source->msa_key ||
        source->phase_set <= 0 || source->msa_key->type != VariantType::Snp ||
        source->msa_key->ref_len != 1 || source->msa_key->alt.size() != 1 ||
        source->hap1_allele < 0 || source->hap1_allele > 1 ||
        source->hap2_allele != 1 - source->hap1_allele) return false;
    const auto path = gc.recovery_source_path_supported.find(source->phase_set);
    const auto weak = gc.recovery_source_weak_cuts.find(source->phase_set);
    const auto quality = gc.recovery_source_quality_cuts.find(source->phase_set);
    if (path == gc.recovery_source_path_supported.end() || !path->second ||
        weak == gc.recovery_source_weak_cuts.end() || !weak->second.empty() ||
        (quality != gc.recovery_source_quality_cuts.end() && !quality->second.empty()))
        return false;
    const std::string* alt = selected_graph_candidate_alt(gc, candidate_index);
    if (alt == nullptr) return false;
    const GraphSiteMeta& meta = gc.site_meta[candidate_index];
    const VariantKey physical = vcf_to_variant_key(source->msa_key->tid, meta.pos, meta.ref, *alt);
    if (exact_comp_var_site(&physical, &*source->msa_key) != 0) return false;
    for (const RecoveryPhaseGauge& gauge : gc.recovery_phase_gauges) {
        std::optional<size_t> marker;
        for (size_t si = 0; si < gauge.bam_sites.size(); ++si) {
            const RecoveryBamSite& saved = gauge.bam_sites[si];
            if (saved.phase_set != source->phase_set ||
                exact_comp_var_site(&saved.key, &physical) != 0) continue;
            if (marker || saved.hap1_allele != source->hap1_allele ||
                saved.hap2_allele != source->hap2_allele) return false;
            marker = si;
        }
        if (!marker) continue;
        std::array<int, 2> support{};
        std::unordered_set<std::string> seen;
        for (const RecoveryBamRead& read : gauge.bam_reads) {
            if (read.mapq < kMinQuality || read.mapq == kUnknownQuality ||
                !seen.insert(read.qname).second) continue;
            int marker_hap = 0;
            std::array<std::set<hts_pos_t>, 3> loci;
            for (size_t oi = 0; oi < read.observations.size(); ++oi) {
                const auto [si, allele] = read.observations[oi];
                if (si >= gauge.bam_sites.size() || oi >= read.base_qualities.size() ||
                    read.base_qualities[oi] < kMinQuality ||
                    read.base_qualities[oi] == kUnknownQuality) continue;
                const RecoveryBamSite& saved = gauge.bam_sites[si];
                if (saved.phase_set != source->phase_set ||
                    (allele != 0 && allele != 1)) continue;
                const int hap = allele == saved.hap1_allele ? 1 :
                    allele == saved.hap2_allele ? 2 : 0;
                if (hap == 0) continue;
                if (si == *marker) marker_hap = hap;
                else if (saved.clean_snp && saved.key.type == VariantType::Snp &&
                         saved.key.ref_len == 1 && saved.key.alt.size() == 1 &&
                         saved.key.pos != physical.pos)
                    loci[hap].insert(saved.key.pos);
            }
            if (marker_hap == 0) continue;
            if (!loci[3 - marker_hap].empty()) return false;
            if (loci[marker_hap].size() >= kMinIndependentLoci) ++support[marker_hap - 1];
        }
        if (support[0] >= kMinAlleleSupport && support[1] >= kMinAlleleSupport &&
            std::ldexp(1.0, -(support[0] + support[1])) <= kMaxWrongParity)
            return true;
    }
    return false;
}

std::optional<VariantKey> retained_shared_deletion_key(
        const GraphChunkBuildResult& gc, size_t candidate_index) {
    if (candidate_index >= gc.chunk.candidates.size() ||
        candidate_index >= gc.site_meta.size()) return std::nullopt;
    const CandidateVariant& candidate = gc.chunk.candidates[candidate_index];
    if (candidate.bam_injected || !candidate.alignment_verified ||
        !is_phase_set_anchor(candidate) || candidate.hap_to_cons_alle[1] > 1 ||
        candidate.hap_to_cons_alle[2] != 1 - candidate.hap_to_cons_alle[1])
        return std::nullopt;
    const RecoverySourceSite* origin = nullptr;
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        if (source.candidate_index != candidate_index) continue;
        if (origin != nullptr) return std::nullopt;
        origin = &source;
    }
    if (origin == nullptr || !origin->can_adopt || !origin->msa_key ||
        origin->msa_key->type != VariantType::Deletion ||
        !origin->msa_key->alt.empty()) return std::nullopt;
    const std::string* alt = selected_graph_candidate_alt(gc, candidate_index);
    if (alt == nullptr) return std::nullopt;
    const GraphSiteMeta& meta = gc.site_meta[candidate_index];
    const VariantKey key = vcf_to_variant_key(
        origin->msa_key->tid, meta.pos, meta.ref, *alt);
    if (key.type != VariantType::Deletion || !key.alt.empty() || key.ref_len <= 0 ||
        key.ref_len != origin->msa_key->ref_len) return std::nullopt;
    return key;
}

int masked_bam_snp_haplotype(const PhasingChunk& chunk,
                            const ReadVariantProfile& profile, size_t witness_i,
                            const std::vector<size_t>& physically_masked) {
    constexpr int kMinBaseQuality = 30;
    constexpr int kUnknownQuality = 255;
    if (profile.start_var_idx < 0 || witness_i >= chunk.candidates.size() ||
        witness_i < static_cast<size_t>(profile.start_var_idx)) return 0;
    const CandidateVariant& witness = chunk.candidates[witness_i];
    const size_t offset = witness_i - profile.start_var_idx;
    if (!is_phase_set_anchor(witness) || witness.bam_injected ||
        witness.counts.category != VariantCategory::CleanHetSnp ||
        offset >= profile.alleles.size() || offset >= profile.graph_alleles.size() ||
        offset >= profile.bam_alleles.size() || offset >= profile.bam_base_qualities.size() ||
        profile.bam_base_qualities[offset] < kMinBaseQuality ||
        profile.bam_base_qualities[offset] == kUnknownQuality ||
        profile.alleles[offset] < 0 || profile.alleles[offset] > 1 ||
        profile.graph_alleles[offset] != profile.alleles[offset] ||
        profile.bam_alleles[offset] != profile.alleles[offset]) return 0;
    const int hap = profile.alleles[offset] == witness.hap_to_cons_alle[1] ? 1 : 2;
    bool masked = false;
    for (size_t off = 0; off < profile.alleles.size(); ++off) {
        const size_t ci = static_cast<size_t>(profile.start_var_idx) + off;
        if (ci >= chunk.candidates.size()) break;
        const CandidateVariant& site = chunk.candidates[ci];
        if (!is_phase_set_anchor(site)) continue;
        if (std::find(physically_masked.begin(), physically_masked.end(), ci) !=
                physically_masked.end()) {
            if (site.phase_set != witness.phase_set || site.bam_injected ||
                site.counts.category != VariantCategory::CleanHetSnp ||
                site.key.type != VariantType::Snp || profile.alleles[off] != 0 ||
                off >= profile.bam_alleles.size() || profile.bam_alleles[off] != 0 ||
                off >= profile.bam_base_qualities.size() ||
                profile.bam_base_qualities[off] != 0) return 0;
            masked = true;
            continue;
        }
        for (const auto* channel : {&profile.alleles, &profile.bam_alleles})
            if (off < channel->size() && (*channel)[off] >= 0 &&
                (site.phase_set != witness.phase_set ||
                 (*channel)[off] != site.hap_to_cons_alle[hap])) return 0;
    }
    return masked ? hap : 0;
}

int physical_snp_call(const bam1_t* aln, hts_pos_t pos,
                      char ref_base, char alt_base,
                      int* base_quality) {
    if (aln == nullptr) return -1;
    hts_pos_t ref_pos = aln->core.pos + 1;
    int query_pos = 0;
    const uint32_t* cigar = bam_get_cigar(aln);
    for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int len = bam_cigar_oplen(cigar[i]);
        const int consumed = bam_cigar_type(op);
        if ((consumed & 2) != 0 && ref_pos <= pos && pos < ref_pos + len) {
            if (op == BAM_CDEL) return 1;
            if (op != BAM_CMATCH && op != BAM_CEQUAL && op != BAM_CDIFF)
                return -1;
            const int qi = query_pos + static_cast<int>(pos - ref_pos);
            if (qi < 0 || qi >= aln->core.l_qseq) return -1;
            const char base = seq_nt16_str[bam_seqi(bam_get_seq(aln), qi)];
            if (base_quality != nullptr)
                *base_quality = bam_get_qual(aln)[qi];
            // Reference FASTA can be soft-masked while BAM bases are encoded
            // uppercase. Compare nucleotides, not their original case.
            if (base == std::toupper(static_cast<unsigned char>(ref_base)))
                return 0;
            if (base == std::toupper(static_cast<unsigned char>(alt_base)))
                return 2;
            return -1;
        }
        if ((consumed & 2) != 0) ref_pos += len;
        if ((consumed & 1) != 0) query_pos += len;
        if (ref_pos > pos) break;
    }
    return -1;
}

uint8_t bam_snp_observation_quality(const bam1_t* alignment, hts_pos_t pos,
                                  char ref_base, char alt_base, int allele) {
    constexpr int kMissingBaseQuality = 255;
    if (allele < 0 || allele > 1 || base_to_nt4(ref_base) > 3 ||
        base_to_nt4(alt_base) > 3 || base_to_nt4(ref_base) == base_to_nt4(alt_base))
        return 0;
    int quality = 0;
    const int physical = physical_snp_call(
        alignment, pos, ref_base, alt_base, &quality);
    // MSA can move an observation relative to the original CIGAR. Its allele
    // remains usable, but cannot borrow the quality of a different BAM base.
    return physical == (allele == 0 ? 0 : 2) && quality != kMissingBaseQuality
        ? static_cast<uint8_t>(quality) : 0;
}

// Re-solving a seam created inside a completed graph gap can reorient an
// established block. Require a callable, directionally consistent SNP pair
// before exposing that new pair to the BAM sub-solve.
bool has_direct_snp_parity_for_retry(
        const GraphChunkBuildResult& graph_chunk, const RecoverySeam& seam,
        WorkerContext& context, int tid,
        std::optional<bool>* graph_parity) {
    if (graph_parity != nullptr) graph_parity->reset();
    const PhasingChunk& chunk = graph_chunk.chunk;
    constexpr int kMinMapq = 30;
    constexpr int kMinBaseq = 30;
    constexpr int kUnknownQuality = 255;
    // High-quality calls do not make a small or one-haplotype cohort a
    // reliable diploid connection. Test molecule counts as well as quality odds.
    constexpr int kMinPairedReads = 4;
    constexpr int kMinSupportPerHap = 2;
    constexpr double kMinDominantFraction = 0.75;
    constexpr double kMaxWrongParity = 0.001;
    if (context.bams.empty() || context.indexes.empty()) return false;
    std::optional<size_t> left_i, right_i;
    std::optional<VariantKey> left_key, right_key;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if ((candidate.counts.category != VariantCategory::CleanHetSnp &&
             !(candidate.counts.category == VariantCategory::NoisyCandHet &&
               candidate.msa_verified && candidate.alignment_verified)) ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1])
            continue;
        if (candidate.phase_set != seam.left_phase_set &&
            candidate.phase_set != seam.right_phase_set)
            continue;
        // Graph keys identify walks, not nucleotides. Project only the selected
        // allele; normalization can also move a padded SNP's reference base.
        VariantKey key = candidate.key;
        if (!candidate.bam_injected) {
            const std::string* alt = selected_graph_candidate_alt(graph_chunk, ci);
            if (alt == nullptr || ci >= graph_chunk.site_meta.size()) continue;
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
            key = vcf_to_variant_key(candidate.key.tid, meta.pos, meta.ref, *alt);
        }
        if (key.type != VariantType::Snp || key.ref_len != 1 || key.alt.size() != 1)
            continue;
        if (candidate.phase_set == seam.left_phase_set &&
            key.pos <= seam.beg &&
            (!left_i || key.pos > left_key->pos)) {
            left_i = ci;
            left_key = std::move(key);
        } else if (candidate.phase_set == seam.right_phase_set &&
            key.pos >= seam.end &&
            (!right_i || key.pos < right_key->pos)) {
            right_i = ci;
            right_key = std::move(key);
        }
    }
    if (!left_i || !right_i) return false;
    const CandidateVariant& left = chunk.candidates[*left_i];
    const CandidateVariant& right = chunk.candidates[*right_i];
    if (left_key->pos >= right_key->pos) return false;
    const char left_ref = context.ref.base(
        tid, left_key->pos, context.primary_header());
    const char right_ref = context.ref.base(
        tid, right_key->pos, context.primary_header());
    if (left_ref == 'N' || right_ref == 'N') return false;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       left_key->pos - 1, right_key->pos), &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    std::unordered_set<std::string> seen;
    std::array<int, 2> votes{};
    std::array<std::array<int, 2>, 2> hap_votes{};
    double log_odds = 0.0;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                         alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq ||
            read->core.qual == kUnknownQuality ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        int left_quality = 0, right_quality = 0;
        const int left_call = physical_snp_call(
            read, left_key->pos, left_ref, left_key->alt[0], &left_quality);
        const int right_call = physical_snp_call(
            read, right_key->pos, right_ref, right_key->alt[0], &right_quality);
        if ((left_call != 0 && left_call != 2) ||
            (right_call != 0 && right_call != 2) ||
            left_quality < kMinBaseq || right_quality < kMinBaseq ||
            left_quality == kUnknownQuality ||
            right_quality == kUnknownQuality)
            continue;
        seen.insert(bam_get_qname(read));
        const bool left_hap1 = (left_call == 2) ==
            (left.hap_to_cons_alle[1] == 1);
        const bool right_hap1 = (right_call == 2) ==
            (right.hap_to_cons_alle[1] == 1);
        const bool flip = left_hap1 != right_hap1;
        ++votes[flip ? 1 : 0];
        ++hap_votes[left_hap1 ? 0 : 1][flip ? 1 : 0];
        const double p =
            std::pow(10.0, -left_quality / 10.0) +
            std::pow(10.0, -right_quality / 10.0) +
            2.0 * std::pow(10.0, -read->core.qual / 10.0);
        if (p > 0.0 && p < 0.5)
            log_odds += (flip ? 1.0 : -1.0) *
                std::log((1.0 - p) / p);
    }
    const int total = votes[0] + votes[1];
    const int winner = std::max(votes[0], votes[1]);
    if (!left.bam_injected || !right.bam_injected) {
        const size_t parity = votes[1] > votes[0] ? 1 : 0;
        // Projected graph anchors do not carry the BAM source's independent
        // MSA phase evidence. Require both haplotypes to support their parity;
        // a larger cohort cannot hide a reversal on the other haplotype.
        for (const auto& hap : hap_votes)
            if (hap[parity] < kMinSupportPerHap ||
                hap[parity] <= hap[1 - parity])
                return false;
        if (binomial_upper_tail(total, winner, 0.5) > kMaxWrongParity)
            return false;
    }
    const bool admitted = total >= kMinPairedReads &&
           static_cast<double>(winner) / total >= kMinDominantFraction &&
           (log_odds > 0.0) == (votes[1] > votes[0]) &&
           std::abs(log_odds) >=
               std::log((1.0 - kMaxWrongParity) / kMaxWrongParity);
    if (admitted && graph_parity != nullptr &&
        (!left.bam_injected || !right.bam_injected))
        *graph_parity = votes[1] > votes[0];
    return admitted;
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

std::optional<bool> calibrated_indel_bridge_flip(
        const IndependentBamBlockLink& gauge, const std::array<int, 2>& parity,
        double log_odds) {
    constexpr double kMaxBridgeWrongParity = 0.05;
    constexpr double kMaxJointError = 0.20;
    constexpr double kMaxGaugeCallError = 0.01;
    const auto& counts = gauge.counts;
    const int same = counts[0][0] + counts[1][1];
    const int cross = counts[0][1] + counts[1][0];
    if (same <= cross || counts[0][0] + counts[0][1] == 0 ||
        counts[1][0] + counts[1][1] == 0 ||
        counts[0][0] + counts[1][0] == 0 || counts[0][1] + counts[1][1] == 0 ||
        2.0 * rescue_binomial_tail(same, same + cross) > kIndependentBlockAssociationPValue ||
        (parity[0] != 0 && parity[1] != 0) || parity[0] + parity[1] == 0)
        return std::nullopt;
    const bool flip = parity[1] != 0;
    const double signed_odds = flip ? log_odds : -log_odds;
    if (signed_odds < std::log((1.0 - kMaxBridgeWrongParity) / kMaxBridgeWrongParity))
        return std::nullopt;
    const double bridge_error = 1.0 / (1.0 + std::exp(signed_odds));
    if (one_sided_wilson_upper_bound(cross, same + cross) +
        kMaxGaugeCallError + bridge_error > kMaxJointError)
        return std::nullopt;
    return flip;
}

std::optional<bool> calibrated_source_deletion_bridge_flip(
        const IndependentBamBlockLink& gauge, const std::array<int, 2>& parity,
        double quality_error_bound) {
    constexpr int kMinAlleleClassSupport = 2;
    constexpr double kMaxJointError = 0.20;
    const auto& counts = gauge.counts;
    const int same = counts[0][0] + counts[1][1];
    if (counts[0][0] < kMinAlleleClassSupport || counts[1][1] < kMinAlleleClassSupport ||
        counts[0][1] != 0 || counts[1][0] != 0 ||
        2.0 * rescue_binomial_tail(same, same) > kIndependentBlockAssociationPValue ||
        (parity[0] != 0 && parity[1] != 0) || parity[0] + parity[1] == 0 ||
        !std::isfinite(quality_error_bound) || quality_error_bound < 0.0 ||
        quality_error_bound > kMaxJointError) return std::nullopt;
    return parity[1] != 0;
}

int complementary_insertion_length_class(int observed, int first, int second) {
    if (observed <= 0 || first <= 0 || second <= 0 || first == second) return -1;
    const int first_distance = std::abs(observed - first);
    const int second_distance = std::abs(observed - second);
    return first_distance == second_distance ? -1 :
        first_distance < second_distance ? 0 : 1;
}

std::optional<bool> calibrated_repeat_insertion_bridge_flip(
        const IndependentBamBlockLink& gauge,
        const std::array<std::array<int, 2>, 2>& parity,
        const std::array<double, 2>& call_errors) {
    constexpr double kMaxJointError = 0.20;
    constexpr double kMaxGaugeCallError = 0.01;
    const auto& counts = gauge.counts;
    const int same = counts[0][0] + counts[1][1];
    const int cross = counts[0][1] + counts[1][0];
    if (same <= cross || counts[0][0] == 0 || counts[1][1] == 0 ||
        2.0 * rescue_binomial_tail(same, same + cross) > kIndependentBlockAssociationPValue)
        return std::nullopt;
    std::optional<bool> flip;
    double joint_error = 1.0;
    for (size_t hap = 0; hap < parity.size(); ++hap) {
        if ((parity[hap][0] != 0 && parity[hap][1] != 0) ||
            parity[hap][0] + parity[hap][1] == 0 ||
            !std::isfinite(call_errors[hap]) || call_errors[hap] < 0.0 ||
            call_errors[hap] >= 0.5) return std::nullopt;
        const bool current = parity[hap][1] != 0;
        if (flip && *flip != current) return std::nullopt;
        flip = current;
        // One incorrect class would give contrary parity and veto the join.
        // A coherent wrong parity needs independent errors in both classes.
        joint_error *= std::min(1.0,
            one_sided_wilson_upper_bound(counts[hap][1 - hap],
                counts[hap][0] + counts[hap][1]) + kMaxGaugeCallError + call_errors[hap]);
    }
    return joint_error <= kMaxJointError ? flip : std::nullopt;
}

std::optional<bool> calibrated_repeat_snp_bridge_flip(
        const std::array<IndependentBamBlockLink, 2>& gauges,
        const std::array<int, 2>& parity, double wrong_parity_bound) {
    constexpr double kMaxWrongParity = 0.001;
    constexpr double kMaxJointError = 0.20;
    constexpr double kMaxGaugeCallError = 0.01;
    if ((parity[0] != 0 && parity[1] != 0) || parity[0] + parity[1] == 0 ||
        !std::isfinite(wrong_parity_bound) || wrong_parity_bound < 0.0 ||
        wrong_parity_bound > kMaxWrongParity)
        return std::nullopt;
    double joint_error = wrong_parity_bound;
    for (const IndependentBamBlockLink& gauge : gauges) {
        const auto& counts = gauge.counts;
        const int same = counts[0][0] + counts[1][1];
        const int cross = counts[0][1] + counts[1][0];
        if (same <= cross || counts[0][0] == 0 || counts[1][1] == 0 ||
            2.0 * rescue_binomial_tail(same, same + cross) >
                kIndependentBlockAssociationPValue)
            return std::nullopt;
        joint_error += one_sided_wilson_upper_bound(cross, same + cross) +
            kMaxGaugeCallError;
    }
    return joint_error <= kMaxJointError ?
        std::optional<bool>(parity[1] != 0) : std::nullopt;
}

std::optional<bool> calibrated_deletion_chain_flip(
        const std::array<IndependentBamBlockLink, 2>& gauges,
        const std::array<int, 2>& parity) {
    constexpr double kMinAgreement = 0.80;
    constexpr double kMaxWrongParity = 0.20;
    constexpr int kMinIndependentBridgeMolecules = 3;
    double error = 0.0;
    for (const IndependentBamBlockLink& gauge : gauges) {
        const auto& counts = gauge.counts;
        const int same = counts[0][0] + counts[1][1];
        const int cross = counts[0][1] + counts[1][0];
        const int total = same + cross;
        if (counts[0][0] == 0 || counts[1][1] == 0 || total == 0 ||
            static_cast<double>(same) / total < kMinAgreement ||
            2.0 * rescue_binomial_tail(same, total) > kIndependentBlockAssociationPValue)
            return std::nullopt;
        // Repeat-length error is measured on separate molecules; base quality
        // alone cannot account for an erroneous but well-aligned repeat length.
        error += static_cast<double>(cross + 1) / (total + 2);
    }
    const int total = parity[0] + parity[1];
    const int winner = std::max(parity[0], parity[1]);
    if (total < kMinIndependentBridgeMolecules ||
        static_cast<double>(winner) / total < kMinAgreement || error >= 0.5)
        return std::nullopt;
    const double log_odds = (2 * winner - total) * std::log((1.0 - error) / error);
    if (log_odds < std::log((1.0 - kMaxWrongParity) / kMaxWrongParity))
        return std::nullopt;
    return parity[1] > parity[0];
}

std::optional<bool> calibrated_repeat_chain_flip(
        const IndependentBamBlockLink& gauge,
        const std::array<std::array<int, 2>, 2>& bridge,
        const std::array<int, 2>& right_parity, double right_log_odds) {
    constexpr int kMinClassSupport = 2;
    constexpr int kMinRightMolecules = 6;
    constexpr double kMinAgreement = 0.80;
    constexpr double kMaxWrongParity = 0.001;
    constexpr double kMaxCalibrationError = 0.20;
    constexpr double kMaxCallError = 0.01;
    const auto& counts = gauge.counts;
    if (counts[0][0] == 0 || counts[1][1] == 0 ||
        counts[0][0] + counts[1][1] <= counts[0][1] + counts[1][0] ||
        2.0 * rescue_binomial_tail(counts[0][0] + counts[1][1],
            counts[0][0] + counts[0][1] + counts[1][0] + counts[1][1]) > kIndependentBlockAssociationPValue ||
        one_sided_wilson_upper_bound(counts[0][1] + counts[1][0],
            counts[0][0] + counts[0][1] + counts[1][0] + counts[1][1]) + kMaxCallError > kMaxCalibrationError ||
        bridge[0][0] < kMinClassSupport || bridge[1][1] < kMinClassSupport ||
        bridge[0][1] != 0 || bridge[1][0] != 0) return std::nullopt;
    const int total = right_parity[0] + right_parity[1];
    const int winner = std::max(right_parity[0], right_parity[1]);
    const bool flip = right_parity[1] > right_parity[0];
    if (total < kMinRightMolecules || static_cast<double>(winner) / total < kMinAgreement ||
        !std::isfinite(right_log_odds) || (flip ? right_log_odds : -right_log_odds) <
            std::log((1.0 - kMaxWrongParity) / kMaxWrongParity)) return std::nullopt;
    return flip;
}

std::optional<int> calibrated_verified_insertion_hap1(
        const std::array<IndependentBamBlockLink, 2>& cohorts) {
    constexpr double kMaxJointError = 0.20;
    constexpr double kMaxCallError = 0.01;
    std::optional<int> hap1;
    int total = 0, discordant = 0;
    for (const IndependentBamBlockLink& cohort : cohorts) {
        const auto& counts = cohort.counts;
        const int same = counts[0][0] + counts[1][1];
        const int cross = counts[0][1] + counts[1][0];
        if (same == cross || counts[0][0] + counts[0][1] == 0 ||
            counts[1][0] + counts[1][1] == 0 ||
            counts[0][0] + counts[1][0] == 0 || counts[0][1] + counts[1][1] == 0 ||
            2.0 * rescue_binomial_tail(std::max(same, cross), same + cross) >
                kIndependentBlockAssociationPValue) return std::nullopt;
        const int allele = same > cross ? 0 : 1;
        if (hap1 && *hap1 != allele) return std::nullopt;
        hap1 = allele;
        total += same + cross;
        discordant += std::min(same, cross);
    }
    return one_sided_wilson_upper_bound(discordant, total) + kMaxCallError <= kMaxJointError ? hap1 : std::nullopt;
}

std::optional<int> calibrated_terminal_insertion_hap1(
        const std::array<IndependentBamBlockLink, 2>& cohorts) {
    constexpr double kMaxCohortError = 0.20;
    constexpr double kMaxCombinedError = 0.10;
    constexpr double kMaxCallError = 0.01;
    std::optional<int> hap1;
    int total = 0;
    int discordant = 0;
    for (const IndependentBamBlockLink& cohort : cohorts) {
        const auto& counts = cohort.counts;
        const int same = counts[0][0] + counts[1][1];
        const int cross = counts[0][1] + counts[1][0];
        if (same == cross || counts[0][0] + counts[0][1] == 0 ||
            counts[1][0] + counts[1][1] == 0 ||
            counts[0][0] + counts[1][0] == 0 || counts[0][1] + counts[1][1] == 0 ||
            2.0 * rescue_binomial_tail(std::max(same, cross), same + cross) >
                kIndependentBlockAssociationPValue ||
            one_sided_wilson_upper_bound(std::min(same, cross), same + cross) +
                kMaxCallError > kMaxCohortError) return std::nullopt;
        const int allele = same > cross ? 0 : 1;
        if (hap1 && *hap1 != allele) return std::nullopt;
        hap1 = allele;
        total += same + cross;
        discordant += std::min(same, cross);
    }
    if (one_sided_wilson_upper_bound(discordant, total) + kMaxCallError >
        kMaxCombinedError) return std::nullopt;
    return hap1;
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
                                  bool use_bam_observations,
                                  bool use_graph_observations = false) {
    if (!use_bam_observations && use_graph_observations) {
        return offset < profile.graph_alleles.size() &&
            (profile.graph_alleles[offset] == 0 || profile.graph_alleles[offset] == 1)
                ? profile.graph_alleles[offset] : -1;
    }
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

// Validate an excluded graph indel against the closest phased clean SNP
// in the proposed block. This compares alleles on the same molecules, without
// trusting a possibly switched read HP label elsewhere in that block. The
// selected pair must agree with the inferred marker in all reads and in both
// deterministic, disjoint read halves.
static bool rescue_indel_has_direct_snp_support(
        const PhasingChunk& chunk, size_t site_i,
        const RescueMarker& marker,
        const std::vector<RecoverySeam>& recovery_windows) {
    constexpr double kDirectSnpPValue = 0.05;
    constexpr double kMaxDirectSnpDiscordance = 0.25;
    constexpr double kMaxUnlinkedSiteFraction = 0.5;
    const CandidateVariant& site = chunk.candidates[site_i];
    if (site.bam_injected || site.phase_set > 0 ||
        site.lcd_var_i_to_cate != kLongcalldRepHetVar ||
        site.key.type == VariantType::Snp || recovery_windows.empty()) {
        return false;
    }
    const hts_pos_t site_pos = site.key.sort_pos();
    const bool in_gap = std::any_of(
        recovery_windows.begin(), recovery_windows.end(),
        [site_pos](const RecoverySeam& window) {
            return window.beg < site_pos && site_pos < window.end;
        });
    if (!in_gap) return false;
    // Two unphased graph alternatives at one coordinate may encode the same
    // event differently. A direct SNP pair cannot validate both rows as
    // independent singletons.
    for (size_t i = 0; i < chunk.candidates.size(); ++i) {
        if (i != site_i && !chunk.candidates[i].bam_injected &&
            chunk.candidates[i].phase_set <= 0 &&
            chunk.candidates[i].key.sort_pos() == site_pos) {
            return false;
        }
    }

    size_t anchor_i = chunk.candidates.size();
    hts_pos_t nearest = std::numeric_limits<hts_pos_t>::max();
    for (size_t i = 0; i < chunk.candidates.size(); ++i) {
        const CandidateVariant& anchor = chunk.candidates[i];
        if (anchor.phase_set != marker.phase_set ||
            anchor.key.type != VariantType::Snp ||
            anchor.lcd_var_i_to_cate != kCandCleanHetSnp ||
            !is_oriented_biallelic_candidate(anchor)) {
            continue;
        }
        const hts_pos_t distance =
            std::llabs(anchor.key.sort_pos() - site_pos);
        if (distance < nearest) {
            anchor_i = i;
            nearest = distance;
        }
    }
    if (anchor_i == chunk.candidates.size()) return false;

    std::array<int, 3> agrees{};
    std::array<int, 3> conflicts{};
    std::array<bool, 2> observed_site_allele{};
    std::array<bool, 2> observed_anchor_hap{};
    int site_observations = 0;
    const CandidateVariant& anchor = chunk.candidates[anchor_i];
    for (const ReadVariantProfile& profile : chunk.read_var_profile) {
        if (profile.read_id < 0 ||
            static_cast<size_t>(profile.read_id) >= chunk.reads.size() ||
            chunk.reads[static_cast<size_t>(profile.read_id)].is_skipped ||
            profile.start_var_idx < 0 ||
            site_i < static_cast<size_t>(profile.start_var_idx) ||
            site_i > static_cast<size_t>(profile.end_var_idx)) {
            continue;
        }
        const int site_allele = rescue_observed_allele(
            profile, site_i - static_cast<size_t>(profile.start_var_idx),
            false);
        if (site_allele < 0) continue;
        ++site_observations;
        if (anchor_i < static_cast<size_t>(profile.start_var_idx) ||
            anchor_i > static_cast<size_t>(profile.end_var_idx)) {
            continue;
        }
        const int anchor_allele = rescue_observed_allele(
            profile, anchor_i - static_cast<size_t>(profile.start_var_idx),
            false);
        if (anchor_allele < 0) continue;
        const int site_hap = marker.allele_to_hap[
            static_cast<size_t>(site_allele)];
        const int anchor_hap = anchor_allele ==
            anchor.hap_to_cons_alle[1] ? 1 : 2;
        observed_site_allele[static_cast<size_t>(site_allele)] = true;
        observed_anchor_hap[static_cast<size_t>(anchor_hap - 1)] = true;

        uint64_t hash = 14695981039346656037ULL;
        const std::string& qname =
            chunk.reads[static_cast<size_t>(profile.read_id)].qname;
        for (const unsigned char byte : qname) {
            hash ^= byte;
            hash *= 1099511628211ULL;
        }
        const size_t fold = static_cast<size_t>(hash & 1ULL) + 1;
        std::array<int, 3>& votes =
            site_hap == anchor_hap ? agrees : conflicts;
        ++votes[0];
        ++votes[fold];
    }
    const int joint_observations = agrees[0] + conflicts[0];
    if (!observed_site_allele[0] || !observed_site_allele[1] ||
        !observed_anchor_hap[0] || !observed_anchor_hap[1] ||
        one_sided_wilson_upper_bound(conflicts[0],
                                     joint_observations) >
            kMaxDirectSnpDiscordance ||
        one_sided_wilson_upper_bound(
            site_observations - joint_observations,
            site_observations) > kMaxUnlinkedSiteFraction) {
        return false;
    }
    for (size_t fold = 0; fold < agrees.size(); ++fold) {
        if (agrees[fold] <= conflicts[fold] ||
            rescue_binomial_tail(agrees[fold],
                                 agrees[fold] + conflicts[fold]) >
                kDirectSnpPValue) {
            return false;
        }
    }
    return true;
}

// A singleton needs primary read support in its proposed orientation. Rescued
// reads cannot certify another singleton, even when their block is now joined.
static bool rescue_singleton_has_primary_support(
        const RescueSiteVote& vote, const RescueMarker& marker) {
    const auto& primary = vote.primary_counts;
    const int hap1 = primary[0][0] + primary[0][1];
    const int hap2 = primary[1][0] + primary[1][1];
    const int allele0 = primary[0][0] + primary[1][0];
    const int allele1 = primary[0][1] + primary[1][1];
    const int same = primary[0][0] + primary[1][1];
    const int cross = primary[0][1] + primary[1][0];
    const int total = same + cross;
    return hap1 > 0 && hap2 > 0 && allele0 > 0 && allele1 > 0 &&
        same != cross &&
        ((same > cross) == (marker.allele_to_hap[0] == 1)) &&
        rescue_binomial_tail(std::max(same, cross), total) <=
            kRescueSiteOrientationPValue &&
        one_sided_wilson_upper_bound(std::min(same, cross), total) <=
            kRescueSingletonMaxDiscordance;
}

static size_t rescue_unphased_graph_read_layer(
        PhasingChunk& chunk, bool use_bam_observations,
        const std::vector<RecoverySeam>& recovery_windows) {
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
                profile, offset, use_bam_observations,
                chunk.candidates[static_cast<size_t>(site_i)].bam_independent_genotype);
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

        // K-means can retain distinct internal alleles on a homozygous row.
        // That row needs independent read association, just like an excluded
        // site; its internal labels alone cannot certify a direct singleton.
        if (is_oriented_biallelic_candidate(candidate) &&
            candidate.counts.category != VariantCategory::CleanHom &&
            candidate.counts.category != VariantCategory::NoisyCandHom) {
            marker.phase_set = candidate.phase_set;
            marker.is_direct = true;
            marker.allele_to_hap[static_cast<size_t>(
                candidate.hap_to_cons_alle[1])] = 1;
            marker.allele_to_hap[static_cast<size_t>(
                candidate.hap_to_cons_alle[2])] = 2;
            // A retained noisy MSA genotype can provide block connectivity
            // without making its allele a reliable singleton read call.
            if (candidate.read_rescue_requires_validation || marker.is_snp) {
                if (candidate.read_rescue_requires_validation) marker.is_direct = false;
                const auto found = site_votes[site_i].find(marker.phase_set);
                marker.singleton_safe = found != site_votes[site_i].end() &&
                    rescue_singleton_has_primary_support(found->second, marker);
            }
            if (marker.is_direct || marker.singleton_safe ||
                !candidate.bam_independent_genotype) continue;
            // A verified local genotype can remain independent of the reads'
            // established blocks. If its own PS has no primary support, infer
            // a rescue marker from independent block associations below; this
            // neither changes the genotype nor joins those phase sets.
            marker = RescueMarker{};
            marker.is_snp = candidate.key.type == VariantType::Snp;
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
            supported.singleton_safe =
                rescue_singleton_has_primary_support(entry.second, supported);
            if (!supported.singleton_safe) {
                supported.singleton_safe =
                    rescue_indel_has_direct_snp_support(
                        chunk, site_i, supported, recovery_windows);
            }
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
        std::array<int, 2> certified_snp{};
    };

    struct LocusVote {
        int hap = 0;
        bool is_snp = false;
        bool is_direct = false;
        bool singleton_safe = false;
        bool certified_snp = false;
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
                profile, offset, use_bam_observations,
                chunk.candidates[static_cast<size_t>(site_i)].bam_independent_genotype);
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
                               marker.is_direct, marker.singleton_safe,
                               marker.is_snp && marker.is_direct && marker.singleton_safe});
            if (!inserted && it->second.hap != proposed_hap) {
                it->second.hap = 0;
            } else if (!inserted) {
                // A direct SNP is the strongest form when equivalent rows agree.
                it->second.is_snp = it->second.is_snp || marker.is_snp;
                it->second.is_direct =
                    it->second.is_direct || marker.is_direct;
                it->second.singleton_safe =
                    it->second.singleton_safe || marker.singleton_safe;
                it->second.certified_snp = it->second.certified_snp ||
                    (marker.is_snp && marker.is_direct && marker.singleton_safe);
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
            if (vote.certified_snp)
                ++score.certified_snp[hap_i];
        }

        hts_pos_t best_phase_set = kUnphasedReadPhaseSet;
        int best_hap = 0;
        int best_margin = 0;
        int best_total = 0;
        int best_direct = 0;
        int best_singleton_safe = 0;
        bool best_has_certified_snp = false;
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
            const bool has_certified_snp = use_snps &&
                score.certified_snp[static_cast<size_t>(hap - 1)] > 0;
            // A directly phased SNP can break equal block support only when
            // primary reads independently certify its allele gauge. An inferred
            // excluded-site marker cannot turn its own association into this
            // certificate. Equal certified scores still abstain; no blocks join.
            if (margin > best_margin ||
                (margin == best_margin && (total > best_total ||
                 (total == best_total && has_certified_snp && !best_has_certified_snp)))) {
                best_phase_set = entry.first;
                best_hap = hap;
                best_margin = margin;
                best_total = total;
                best_direct = direct[0] + direct[1];
                best_singleton_safe =
                    singleton_safe[static_cast<size_t>(hap - 1)];
                best_has_certified_snp = has_certified_snp;
                tied = false;
            } else if (margin == best_margin && total == best_total &&
                       has_certified_snp == best_has_certified_snp) {
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

static bool msa_snp_has_physical_gauge(const PhasingChunk& chunk, size_t site_i) {
    constexpr int kMinMapq = 30;
    constexpr int kUnknownMapq = 255;
    constexpr int kMinBaseQuality = 30;
    constexpr hts_pos_t kIndependentSnpDistance = 100;
    const CandidateVariant& site = chunk.candidates[site_i];
    RescueSiteVote vote;
    for (const auto& profile : chunk.read_var_profile) {
        if (profile.read_id < 0 ||
            chunk.reads[static_cast<size_t>(profile.read_id)].is_skipped ||
            profile.bam_mapq < kMinMapq || profile.bam_mapq == kUnknownMapq ||
            profile.start_var_idx < 0 ||
            site_i < static_cast<size_t>(profile.start_var_idx)) continue;
        const size_t site_offset = site_i - static_cast<size_t>(profile.start_var_idx);
        const auto physical_allele = [&](size_t offset) {
            if (offset >= profile.alleles.size() || offset >= profile.bam_alleles.size() ||
                offset >= profile.bam_base_qualities.size() ||
                profile.bam_base_qualities[offset] < kMinBaseQuality ||
                profile.alleles[offset] != profile.bam_alleles[offset]) return -1;
            return profile.bam_alleles[offset];
        };
        const int allele = physical_allele(site_offset);
        if (allele < 0 || allele > 1) continue;
        int anchor_hap = -1;
        bool conflict = false;
        // A nearest-site lookup discards molecules whose usable SNP is farther
        // into the same block. Use every covered clean SNP, but give each
        // molecule only one vote and abstain if its clean gauge disagrees.
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            const CandidateVariant& anchor = chunk.candidates[
                static_cast<size_t>(profile.start_var_idx) + offset];
            if (anchor.phase_set != site.phase_set ||
                anchor.key.type != VariantType::Snp || anchor.key.ref_len != 1 ||
                anchor.counts.category != VariantCategory::CleanHetSnp ||
                !is_oriented_biallelic_candidate(anchor) ||
                std::llabs(anchor.key.sort_pos() - site.key.sort_pos()) <
                    kIndependentSnpDistance) continue;
            const int observed = physical_allele(offset);
            if (observed < 0 || observed > 1) continue;
            const int hap = observed == anchor.hap_to_cons_alle[1] ? 0 : 1;
            if (anchor_hap >= 0 && anchor_hap != hap) { conflict = true; break; }
            anchor_hap = hap;
        }
        if (conflict || anchor_hap < 0) continue;
        ++vote.primary_counts[static_cast<size_t>(anchor_hap)][static_cast<size_t>(allele)];
    }
    RescueMarker marker;
    marker.phase_set = site.phase_set;
    marker.allele_to_hap = {site.hap_to_cons_alle[1] == 0 ? 1 : 2,
                          site.hap_to_cons_alle[1] == 1 ? 1 : 2};
    // Certify against stable SNP alleles, never the read HP being repaired.
    // Reuse singleton rescue's association test and discordance confidence
    // bound, including support for both haplotypes and both site alleles.
    return rescue_singleton_has_primary_support(vote, marker);
}

size_t refresh_recovered_read_haps_from_bam_snps(PhasingChunk& chunk) {
    constexpr int kMinMapq = 30;
    constexpr int kUnknownMapq = 255;
    constexpr int kMinBaseQuality = 30;
    constexpr hts_pos_t kDistinctSnpDistance = 100;
    size_t refreshed = 0;
    std::vector<int8_t> msa_certified(chunk.candidates.size(), -1);
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
        if (chunk.reads[ri].is_skipped || chunk.haps[ri] == 0 ||
            chunk.phase_sets[ri] <= 0) continue;
        const auto& profile = chunk.read_var_profile[ri];
        if (profile.bam_mapq < kMinMapq || profile.bam_mapq == kUnknownMapq)
            continue;
        int proposed_hap = 0, observations = 0, clean_observations = 0;
        hts_pos_t first = 0, last = 0, last_counted = 0, last_clean = 0;
        hts_pos_t first_clean = 0;
        unsigned weak_clean_haps = 0;
        bool verified_msa_witness = false;
        bool conflict = false;
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            if (offset >= profile.bam_alleles.size() ||
                offset >= profile.bam_base_qualities.size()) break;
            const uint8_t quality = profile.bam_base_qualities[offset];
            if (quality == 0) continue;
            const size_t candidate_i = static_cast<size_t>(profile.start_var_idx) + offset;
            const auto& candidate = chunk.candidates[candidate_i];
            const bool clean = candidate.counts.category == VariantCategory::CleanHetSnp;
            const bool verified_msa = candidate.bam_injected && candidate.msa_verified &&
                candidate.counts.category == VariantCategory::NoisyCandHet;
            if (candidate.phase_set != chunk.phase_sets[ri] ||
                candidate.key.type != VariantType::Snp || candidate.key.ref_len != 1 ||
                (!clean && !verified_msa) ||
                !is_phase_set_anchor(candidate)) continue;
            const int allele = profile.bam_alleles[offset];
            if (allele < 0 || allele != profile.alleles[offset]) continue;
            const int hap = allele == candidate.hap_to_cons_alle[1] ? 1 :
                            allele == candidate.hap_to_cons_alle[2] ? 2 : 0;
            if (hap == 0) continue;
            if (quality < kMinBaseQuality) {
                if (clean) weak_clean_haps |= 1u << hap;
                continue;
            }
            if (verified_msa) {
                if (msa_certified[candidate_i] < 0)
                    msa_certified[candidate_i] = msa_snp_has_physical_gauge(chunk, candidate_i);
                if (msa_certified[candidate_i] == 0) continue;
            }
            if (proposed_hap != 0 && proposed_hap != hap) { conflict = true; break; }
            const hts_pos_t pos = candidate.key.sort_pos();
            if (proposed_hap == 0) first = last = pos;
            first = std::min(first, pos);
            last = std::max(last, pos);
            proposed_hap = hap;
            verified_msa_witness |= verified_msa;
            if (clean && (clean_observations == 0 || pos != last_clean)) {
                if (clean_observations == 0) first_clean = pos;
                ++clean_observations;
                last_clean = pos;
            }
            // Candidate order is coordinate order; duplicate descriptions of
            // one base must not increase read confidence.
            if (observations == 0 || pos != last_counted) {
                ++observations;
                last_counted = pos;
            }
        }
        // Imported HP may still reflect the preceding indel clustering.
        // A source-phased MSA SNP also verified by Q30 BAM bases and an
        // independent clean-SNP gauge can correct an indel-derived HP alone.
        // Graph-only SNPs retain their spaced two-locus rule.
        // Any contrary eligible SNP above vetoes both certificates.
        if (conflict || proposed_hap == 0 || proposed_hap == chunk.haps[ri] ||
            (!verified_msa_witness && last - first < kDistinctSnpDistance)) continue;
        // A low-quality physical clean call cannot assign HP, but it is still
        // contrary evidence against overriding HP with one noisy locus. The
        // original spaced Q30 clean-SNP certificate remains sufficient.
        const bool independent_clean_witness = clean_observations >= 2 &&
            last_clean - first_clean >= kDistinctSnpDistance;
        if (verified_msa_witness && !independent_clean_witness &&
            (weak_clean_haps & (1u << (3 - proposed_hap))) != 0) continue;
        chunk.haps[ri] = proposed_hap;
        auto& read = chunk.reads[ri];
        read.n_clean_agree_snps = clean_observations;
        read.n_bridge_agree_snps = observations - clean_observations;
        read.n_bridge_conflict_snps = 0;
        read.n_clean_conflict_snps = 0;
        read.hap_score_margin = observations;
        read.n_vars_scored = observations;
        ++refreshed;
    }
    return refreshed;
}

size_t rescue_unphased_graph_reads(
        PhasingChunk& chunk,
        const std::vector<RecoverySeam>& recovery_windows) {
    size_t total_rescued = 0;
    const auto run_to_fixed_point = [&](bool use_bam_observations) {
        size_t rescued = 0;
        while (true) {
            // Grow from established blocks toward the middle of a gap. Every
            // next layer independently passes the site-orientation test.
            const size_t layer_rescued =
                rescue_unphased_graph_read_layer(
                    chunk, use_bam_observations, recovery_windows);
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

// A seam stitch can absorb a whole BAM PS even when its source-path audit found
// a cut inside it. Keep graph-to-graph joins already validated by the stitch,
// but remove private BAM sites on the far side of an unsupported cut from the
// near graph block. Those sites retain their exact BAM allele orientation and
// may form or attach to a local block in the subsequent component-aware pass.
void detach_bam_sites_across_weak_cuts(GraphChunkBuildResult& gc) {
    PhasingChunk& chunk = gc.chunk;
    std::map<hts_pos_t, std::vector<const RecoverySourceSite*>> by_source;
    std::map<hts_pos_t, std::vector<const RecoverySourceRead*>> reads_by_source;
    for (const RecoverySourceRead& read : gc.recovery_source_reads)
        if (read.phase_set > 0 && read.read_index < chunk.reads.size())
            reads_by_source[read.phase_set].push_back(&read);
    std::set<hts_pos_t> occupied;
    std::set<hts_pos_t> source_ids;
    for (const CandidateVariant& candidate : chunk.candidates)
        if (candidate.phase_set > 0) occupied.insert(candidate.phase_set);
    for (const hts_pos_t read_ps : chunk.phase_sets)
        if (read_ps > 0) occupied.insert(read_ps);
    for (const hts_pos_t gap_ps : chunk.gap_phase_sets)
        if (gap_ps > 0) occupied.insert(gap_ps);
    for (const RecoverySourceSite& site : gc.recovery_source_sites) {
        if (site.phase_set > 0) source_ids.insert(site.phase_set);
        if (site.phase_set > 0 && site.can_adopt &&
            site.candidate_index < chunk.candidates.size())
            by_source[site.phase_set].push_back(&site);
    }

    for (const auto& [source_ps, sites] : by_source) {
        const auto found = gc.recovery_source_weak_cuts.find(source_ps);
        if (found == gc.recovery_source_weak_cuts.end()) continue;
        const std::vector<hts_pos_t>& cuts = found->second;
        for (size_t cut_i = 0; cut_i < cuts.size(); ++cut_i) {
            const hts_pos_t cut = cuts[cut_i];
            const hts_pos_t next_cut = cut_i + 1 < cuts.size()
                ? cuts[cut_i + 1]
                : std::numeric_limits<hts_pos_t>::max();
            // Several graph blocks can own rows before the same source cut.
            // The closest row's owner cannot certify another owner's far-side
            // island; inspect each owner independently rather than strand it.
            std::map<hts_pos_t, hts_pos_t> near_positions;
            for (const RecoverySourceSite* site : sites) {
                const CandidateVariant& candidate =
                    chunk.candidates[site->candidate_index];
                const hts_pos_t pos = candidate.key.sort_pos();
                if (pos <= cut && candidate.phase_set > 0)
                    near_positions[candidate.phase_set] = std::max(
                        near_positions[candidate.phase_set], pos);
            }
            // Additional owners must reach the component immediately before
            // this cut. A label last seen before an older cut may have a later
            // independently certified join; retain the legacy closest-owner
            // check without reopening that unrelated earlier component. An
            // additional graph owner needs an exact shared clean SNP here;
            // inherited indel-only labels cannot certify its source gauge.
            const hts_pos_t previous_cut = cut_i == 0 ? 0 : cuts[cut_i - 1];
            hts_pos_t closest_pos = 0;
            for (const auto& entry : near_positions)
                closest_pos = std::max(closest_pos, entry.second);
            for (auto owner = near_positions.begin(); owner != near_positions.end();) {
                const bool shared_anchor = std::any_of(
                    sites.begin(), sites.end(), [&](const RecoverySourceSite* site) {
                        const CandidateVariant& row = chunk.candidates[site->candidate_index];
                        const hts_pos_t pos = row.key.sort_pos();
                        return site->clean_shared_snp && site->can_adopt &&
                            !row.bam_injected && is_phase_set_anchor(row) &&
                            row.phase_set == owner->first && previous_cut < pos && pos <= cut;
                    });
                if (owner->second != closest_pos &&
                    (owner->second <= previous_cut || !shared_anchor))
                    owner = near_positions.erase(owner);
                else
                    ++owner;
            }
            for (const auto& [near_ps, near_pos] : near_positions) {
                if (near_ps == source_ps) continue;
                bool graph_anchored = false;
                std::vector<const RecoverySourceSite*> far_sites;
                for (const RecoverySourceSite* site : sites) {
                    const CandidateVariant& candidate =
                        chunk.candidates[site->candidate_index];
                    const hts_pos_t pos = candidate.key.sort_pos();
                    if (pos <= cut || pos > next_cut ||
                        candidate.phase_set != near_ps)
                        continue;
                    if (!candidate.bam_injected)
                        graph_anchored = true;
                    else if (site->hap1_allele >= 0 &&
                             site->hap2_allele >= 0 &&
                             site->hap1_allele != site->hap2_allele)
                        far_sites.push_back(site);
                }
                if (graph_anchored || far_sites.empty()) continue;

                // A lone imported site before the cut can be anchored by a graph
                // block. Detach its far-side site only when no recovery read calls
                // a site on each side; a physical bridge may already certify the
                // graph join despite this weak BAM-source path.
                const bool earlier_bam_site = std::any_of(
                    sites.begin(), sites.end(), [&](const RecoverySourceSite* site) {
                        const CandidateVariant& candidate =
                            chunk.candidates[site->candidate_index];
                        return candidate.phase_set == near_ps &&
                               candidate.bam_injected &&
                               candidate.key.sort_pos() < near_pos &&
                               site->hap1_allele >= 0 &&
                               site->hap2_allele >= 0 &&
                               site->hap1_allele != site->hap2_allele;
                    });
                if (!earlier_bam_site) {
                    // A clean SNP at the cut can still make a certified physical
                    // SNP-to-deletion bridge after this transfer pass.
                    const bool indel_cut = std::any_of(
                        sites.begin(), sites.end(), [&](const RecoverySourceSite* site) {
                            const CandidateVariant& candidate =
                                chunk.candidates[site->candidate_index];
                            return candidate.phase_set == near_ps &&
                                   candidate.key.sort_pos() == near_pos &&
                                   candidate.key.type != VariantType::Snp;
                        });
                    if (!indel_cut) continue;
                    const bool graph_anchor = std::any_of(
                        chunk.candidates.begin(), chunk.candidates.end(),
                        [near_ps, cut](const CandidateVariant& candidate) {
                            return candidate.phase_set == near_ps &&
                                   !candidate.bam_injected &&
                                   candidate.key.sort_pos() <= cut &&
                                   candidate.hap_to_cons_alle[1] >= 0 &&
                                   candidate.hap_to_cons_alle[2] >= 0 &&
                                   candidate.hap_to_cons_alle[1] !=
                                       candidate.hap_to_cons_alle[2];
                        });
                    if (!graph_anchor) continue;
                }
                // An additional owner is an orphan only when no observed read
                // connects its near anchors to this island. Retain any existing
                // physical bridge, even if the BAM source split at this cut.
                if (!earlier_bam_site || near_pos != closest_pos) {
                    std::set<size_t> far_indices;
                    for (const RecoverySourceSite* site : far_sites)
                        far_indices.insert(site->candidate_index);
                    const bool crossing_read = std::any_of(
                        chunk.read_var_profile.begin(),
                        chunk.read_var_profile.end(),
                        [&](const ReadVariantProfile& profile) {
                            if (profile.start_var_idx < 0) return false;
                            bool near_call = false;
                            bool far_call = false;
                            for (size_t offset = 0;
                                 offset < profile.alleles.size(); ++offset) {
                                if (profile.alleles[offset] < 0) continue;
                                const size_t ci = static_cast<size_t>(
                                    profile.start_var_idx) + offset;
                                if (ci >= chunk.candidates.size()) break;
                                const CandidateVariant& candidate =
                                    chunk.candidates[ci];
                                near_call |= candidate.phase_set == near_ps &&
                                    candidate.key.sort_pos() <= cut;
                                far_call |= far_indices.count(ci) != 0;
                                if (near_call && far_call) return true;
                            }
                            return false;
                        });
                    if (crossing_read) continue;
                }
                const auto same_alleles = [&chunk](const RecoverySourceSite* site) {
                    const CandidateVariant& candidate =
                        chunk.candidates[site->candidate_index];
                    return candidate.hap_to_cons_alle[1] == site->hap1_allele &&
                           candidate.hap_to_cons_alle[2] == site->hap2_allele;
                };
                const auto reversed_alleles = [&chunk](
                        const RecoverySourceSite* site) {
                    const CandidateVariant& candidate =
                        chunk.candidates[site->candidate_index];
                    return candidate.hap_to_cons_alle[1] == site->hap2_allele &&
                           candidate.hap_to_cons_alle[2] == site->hap1_allele;
                };
                if (std::any_of(far_sites.begin(), far_sites.end(),
                                [&](const RecoverySourceSite* site) {
                                    return !same_alleles(site) &&
                                           !reversed_alleles(site);
                                }))
                    continue;
                hts_pos_t local_ps = source_ps;
                while (local_ps <= 0 || occupied.count(local_ps) != 0 ||
                       (local_ps != source_ps && source_ids.count(local_ps) != 0))
                    ++local_ps;
                occupied.insert(local_ps);
                const auto flip_hap = [](int hap) {
                    return hap == 1 ? 2 : (hap == 2 ? 1 : hap);
                };
                for (const RecoverySourceSite* site : far_sites) {
                    CandidateVariant& candidate =
                        chunk.candidates[site->candidate_index];
                    if (reversed_alleles(site)) {
                        std::swap(candidate.hap_to_cons_alle[1],
                                  candidate.hap_to_cons_alle[2]);
                        std::swap(candidate.hap_to_alle_profile[1],
                                  candidate.hap_to_alle_profile[2]);
                        candidate.hap_alt = flip_hap(candidate.hap_alt);
                        candidate.hap_ref = flip_hap(candidate.hap_ref);
                    }
                    candidate.phase_set = local_ps;
                    candidate.gap_link_supported = false;
                }
                // The site split must carry its source reads too. Otherwise a
                // read calling only this detached component can retain near_ps
                // through an unrelated homozygous row, creating a false read
                // phase-set span across the weak cut.
                std::set<size_t> detached_indices;
                for (const RecoverySourceSite* site : far_sites)
                    detached_indices.insert(site->candidate_index);
                for (const RecoverySourceRead* read : reads_by_source[source_ps]) {
                    const size_t ri = read->read_index;
                    if (ri >= chunk.phase_sets.size() ||
                        ri >= chunk.read_var_profile.size() ||
                        chunk.phase_sets[ri] != near_ps ||
                        (read->hap != 1 && read->hap != 2))
                        continue;
                    const ReadVariantProfile& profile = chunk.read_var_profile[ri];
                    if (profile.start_var_idx < 0) continue;
                    bool sees_detached = false;
                    bool sees_other = false;
                    const auto has_call = [&profile](size_t offset) {
                        return (offset < profile.alleles.size() &&
                                profile.alleles[offset] >= 0) ||
                               (offset < profile.graph_alleles.size() &&
                                profile.graph_alleles[offset] >= 0) ||
                               (offset < profile.bam_alleles.size() &&
                                profile.bam_alleles[offset] >= 0);
                    };
                    for (const RecoverySourceSite* site : sites) {
                        const size_t ci = site->candidate_index;
                        if (ci < static_cast<size_t>(profile.start_var_idx) ||
                            ci > static_cast<size_t>(profile.end_var_idx))
                            continue;
                        const size_t offset =
                            ci - static_cast<size_t>(profile.start_var_idx);
                        if (!has_call(offset)) continue;
                        if (detached_indices.count(ci) != 0)
                            sees_detached = true;
                        else if (site->hap1_allele >= 0 &&
                                 site->hap2_allele >= 0 &&
                                 site->hap1_allele != site->hap2_allele)
                            sees_other = true;
                    }
                    // Graph observations can also support the near block even
                    // when none of its sites are shared with this BAM source.
                    // Do not move a read that still calls such an anchor.
                    for (size_t offset = 0; offset < profile.alleles.size() &&
                                            !sees_other; ++offset) {
                        const size_t ci =
                            static_cast<size_t>(profile.start_var_idx) + offset;
                        if (ci >= chunk.candidates.size()) break;
                        const CandidateVariant& candidate = chunk.candidates[ci];
                        if (has_call(offset) &&
                            candidate.phase_set == near_ps &&
                            candidate.hap_to_cons_alle[1] >= 0 &&
                            candidate.hap_to_cons_alle[2] >= 0 &&
                            candidate.hap_to_cons_alle[1] !=
                                candidate.hap_to_cons_alle[2])
                            sees_other = true;
                    }
                    if (!sees_detached || sees_other) continue;
                    chunk.phase_sets[ri] = local_ps;
                    chunk.haps[ri] = read->hap;
                }
            }
        }
    }
}

bool bam_source_site_path_supported(const GraphChunkBuildResult& gc,
                                    size_t candidate_index) {
    if (candidate_index >= gc.chunk.candidates.size()) return false;
    const CandidateVariant& candidate = gc.chunk.candidates[candidate_index];
    if (!candidate.bam_injected || !is_phase_set_anchor(candidate)) return false;
    const RecoverySourceSite* source = nullptr;
    for (const RecoverySourceSite& site : gc.recovery_source_sites) {
        if (site.candidate_index != candidate_index) continue;
        if (source != nullptr) return false;
        source = &site;
    }
    if (source == nullptr || !source->can_adopt || source->phase_set <= 0)
        return false;
    const auto path = gc.recovery_source_path_supported.find(source->phase_set);
    if (path == gc.recovery_source_path_supported.end() || !path->second)
        return false;
    const auto weak = gc.recovery_source_weak_cuts.find(source->phase_set);
    const auto quality = gc.recovery_source_quality_cuts.find(source->phase_set);
    if ((weak != gc.recovery_source_weak_cuts.end() && !weak->second.empty()) ||
        (quality != gc.recovery_source_quality_cuts.end() &&
         !quality->second.empty()))
        return false;
    const auto orientation = [](const CandidateVariant& row,
                                 const RecoverySourceSite& origin) {
        if (row.hap_to_cons_alle[1] < 0 || row.hap_to_cons_alle[1] > 1 ||
            row.hap_to_cons_alle[2] != 1 - row.hap_to_cons_alle[1] ||
            origin.hap1_allele < 0 || origin.hap1_allele > 1 ||
            origin.hap2_allele != 1 - origin.hap1_allele)
            return -1;
        return static_cast<int>(row.hap_to_cons_alle[1] != origin.hap1_allele);
    };
    const int parity = orientation(candidate, *source);
    if (parity < 0) return false;
    bool anchored = false;
    for (const RecoverySourceSite& site : gc.recovery_source_sites) {
        if (site.phase_set != source->phase_set || !site.can_adopt ||
            site.candidate_index >= gc.chunk.candidates.size())
            continue;
        const CandidateVariant& row = gc.chunk.candidates[site.candidate_index];
        if (row.phase_set != candidate.phase_set) continue;
        // Classified homozygotes can retain stale unequal allele labels.
        // Unknown source gauges still fail the orientation check below.
        if (!is_phase_set_anchor(row) && row.hap_to_cons_alle[1] >= 0 &&
            row.hap_to_cons_alle[2] >= 0)
            continue;
        if (orientation(row, site) != parity) return false;
        anchored |= row.key.sort_pos() != candidate.key.sort_pos();
    }
    return anchored;
}

bool bam_source_prefix_to_graph_supported(const GraphChunkBuildResult& gc,
                                          size_t candidate_index) {
    const PhasingChunk& chunk = gc.chunk;
    if (candidate_index >= chunk.candidates.size()) return false;
    const CandidateVariant& marker = chunk.candidates[candidate_index];
    if (!marker.bam_injected || !is_phase_set_anchor(marker)) return false;
    const RecoverySourceSite* origin = nullptr;
    for (const RecoverySourceSite& site : gc.recovery_source_sites) {
        if (site.candidate_index != candidate_index) continue;
        if (origin != nullptr) return false;
        origin = &site;
    }
    if (origin == nullptr || origin->phase_set <= 0 || !origin->can_adopt)
        return false;
    const auto orientation = [](const CandidateVariant& row,
                                 const RecoverySourceSite& site) {
        if (!site.can_adopt || row.hap_to_cons_alle[1] < 0 ||
            row.hap_to_cons_alle[1] > 1 ||
            row.hap_to_cons_alle[2] != 1 - row.hap_to_cons_alle[1] ||
            site.hap1_allele < 0 || site.hap1_allele > 1 ||
            site.hap2_allele != 1 - site.hap1_allele) return -1;
        return static_cast<int>(row.hap_to_cons_alle[1] != site.hap1_allele);
    };
    const int gauge = orientation(marker, *origin);
    if (gauge < 0) return false;
    const hts_pos_t first = marker.key.sort_pos();
    hts_pos_t last = std::numeric_limits<hts_pos_t>::max();
    for (const RecoverySourceSite& site : gc.recovery_source_sites) {
        if (site.phase_set != origin->phase_set || !site.clean_shared_snp ||
            site.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& row = chunk.candidates[site.candidate_index];
        if (!row.bam_injected && row.phase_set == marker.phase_set &&
            row.counts.category == VariantCategory::CleanHetSnp &&
            is_phase_set_anchor(row) && row.key.sort_pos() > first)
            last = std::min(last, row.key.sort_pos());
    }
    if (last == std::numeric_limits<hts_pos_t>::max()) return false;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& row = chunk.candidates[ci];
        if (row.phase_set != marker.phase_set || !is_phase_set_anchor(row) ||
            row.key.sort_pos() < first || row.key.sort_pos() > last) continue;
        size_t claims = 0;
        for (const RecoverySourceSite& site : gc.recovery_source_sites) {
            if (site.candidate_index != ci || site.phase_set != origin->phase_set)
                continue;
            if ((!row.bam_injected && !site.clean_shared_snp) ||
                orientation(row, site) != gauge) return false;
            ++claims;
        }
        if (claims != 1) return false;
    }
    const auto weak = gc.recovery_source_weak_cuts.find(origin->phase_set);
    if (weak == gc.recovery_source_weak_cuts.end()) return false;
    const auto crosses = [first, last](const auto& cuts) {
        return std::any_of(cuts.begin(), cuts.end(), [first, last](hts_pos_t cut) {
            return first <= cut && cut < last;
        });
    };
    const auto quality = gc.recovery_source_quality_cuts.find(origin->phase_set);
    return !crosses(weak->second) &&
        (quality == gc.recovery_source_quality_cuts.end() || !crosses(quality->second));
}

// Require both observed allele classes and the same significant orientation
// in deterministic read halves before a boundary may orient a whole block.
std::optional<bool> local_run_boundary_flip(
        const PhasingChunk& chunk, size_t source_i, size_t target_i,
        int min_mapq, double max_p) {
    std::array<std::pair<int, int>, 2> halves{};
    std::array<int, 2> source_alleles{};
    const CandidateVariant& source = chunk.candidates[source_i];
    const CandidateVariant& target = chunk.candidates[target_i];
    const bool bam_pair = source.bam_injected && target.bam_injected;
    for (size_t ri = 0; ri < chunk.read_var_profile.size(); ++ri) {
        if (ri >= chunk.reads.size() || chunk.reads[ri].is_skipped) continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        // An independently mapped BAM pair is not gated by its GAF mapping.
        // Use that channel's calls too, so the MAPQ certifies the observations
        // being counted rather than a different alignment of the same read.
        const int mapq = bam_pair ? profile.bam_mapq : chunk.reads[ri].mapq;
        if (mapq < min_mapq || (min_mapq > 0 && mapq == 255)) continue;
        const auto& alleles = bam_pair ? profile.bam_alleles : profile.alleles;
        if (profile.start_var_idx < 0 ||
            source_i < static_cast<size_t>(profile.start_var_idx) ||
            target_i < static_cast<size_t>(profile.start_var_idx) ||
            source_i > static_cast<size_t>(profile.end_var_idx) ||
            target_i > static_cast<size_t>(profile.end_var_idx))
            continue;
        const size_t source_offset =
            source_i - static_cast<size_t>(profile.start_var_idx);
        const size_t target_offset =
            target_i - static_cast<size_t>(profile.start_var_idx);
        if (source_offset >= alleles.size() ||
            target_offset >= alleles.size())
            continue;
        const int source_allele = alleles[source_offset];
        const int target_allele = alleles[target_offset];
        if ((source_allele != 0 && source_allele != 1) ||
            (target_allele != 0 && target_allele != 1))
            continue;
        ++source_alleles[static_cast<size_t>(source_allele)];
        uint64_t hash = 14695981039346656037ULL;
        for (const unsigned char byte : chunk.reads[ri].qname) {
            hash ^= byte;
            hash *= 1099511628211ULL;
        }
        auto& votes = halves[static_cast<size_t>(hash & 1ULL)];
        if ((source_allele == source.hap_to_cons_alle[1]) ==
            (target_allele == target.hap_to_cons_alle[1]))
            ++votes.first;
        else
            ++votes.second;
    }
    if (source_alleles[0] == 0 || source_alleles[1] == 0) return std::nullopt;
    const int same = halves[0].first + halves[1].first;
    const int cross = halves[0].second + halves[1].second;
    const int winner = std::max(same, cross);
    const int total = same + cross;
    if (winner == 0 || same == cross ||
        halves[0].first == halves[0].second ||
        halves[1].first == halves[1].second)
        return std::nullopt;
    const bool flip = cross > same;
    for (const auto& half : halves)
        if ((half.second > half.first) != flip) return std::nullopt;
    double term = std::exp(
        std::lgamma(static_cast<double>(total + 1)) -
        std::lgamma(static_cast<double>(winner + 1)) -
        std::lgamma(static_cast<double>(total - winner + 1)) -
        static_cast<double>(total) * std::log(2.0));
    double tail = term;
    for (int k = winner; k < total; ++k) {
        term *= static_cast<double>(total - k) /
                static_cast<double>(k + 1);
        tail += term;
    }
    return tail <= max_p ? std::optional<bool>(flip) : std::nullopt;
}

std::optional<hts_pos_t> bam_prefix_before_deletion_pair(
        const GraphChunkBuildResult& gc, size_t first_i, size_t second_i,
        hts_pos_t graph_end) {
    constexpr int kMinMapq = 30;
    constexpr int kUnknownMapq = 255;
    constexpr int kMinAlleleClassSupport = 2;
    constexpr double kMaxWrongParity = 0.001;
    const PhasingChunk& chunk = gc.chunk;
    if (first_i >= chunk.candidates.size() || second_i >= chunk.candidates.size() ||
        first_i == second_i) return std::nullopt;
    const CandidateVariant& a = chunk.candidates[first_i];
    const CandidateVariant& b = chunk.candidates[second_i];
    if (!a.bam_injected || !b.bam_injected || !a.msa_verified || !b.msa_verified ||
        !a.alignment_verified || !b.alignment_verified ||
        !is_phase_set_anchor(a) || !is_phase_set_anchor(b) ||
        a.key.type != VariantType::Deletion || b.key.type != VariantType::Deletion ||
        a.key.pos != b.key.pos || a.key.ref_len == b.key.ref_len ||
        a.phase_set != b.phase_set || a.hap_to_cons_alle[1] == b.hap_to_cons_alle[1])
        return std::nullopt;
    const hts_pos_t pair_pos = a.key.sort_pos();
    std::vector<const RecoverySourceSite*> origin(chunk.candidates.size(), nullptr);
    for (const RecoverySourceSite& site : gc.recovery_source_sites) {
        if (site.candidate_index >= origin.size()) continue;
        if (origin[site.candidate_index] != nullptr) return std::nullopt;
        origin[site.candidate_index] = &site;
    }
    if (origin[first_i] == nullptr || origin[second_i] == nullptr)
        return std::nullopt;
    const hts_pos_t source_ps = origin[first_i]->phase_set;
    const auto orientation = [&](size_t ci) {
        const RecoverySourceSite* site = origin[ci];
        const CandidateVariant& row = chunk.candidates[ci];
        if (site == nullptr || !site->can_adopt || site->phase_set != source_ps ||
            row.hap_to_cons_alle[1] < 0 || row.hap_to_cons_alle[1] > 1 ||
            row.hap_to_cons_alle[2] != 1 - row.hap_to_cons_alle[1] ||
            site->hap1_allele < 0 || site->hap1_allele > 1 ||
            site->hap2_allele != 1 - site->hap1_allele) return -1;
        return static_cast<int>(row.hap_to_cons_alle[1] != site->hap1_allele);
    };
    const int gauge = orientation(first_i);
    if (source_ps <= 0 || gauge < 0 || orientation(second_i) != gauge)
        return std::nullopt;
    std::optional<size_t> snp_i;
    hts_pos_t snp_pos = graph_end;
    hts_pos_t run_beg = pair_pos;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& row = chunk.candidates[ci];
        const hts_pos_t pos = row.key.sort_pos();
        if (row.phase_set != a.phase_set || !is_phase_set_anchor(row) ||
            pos <= graph_end || pos > pair_pos) continue;
        if (pos == pair_pos) {
            if (ci != first_i && ci != second_i) return std::nullopt;
            continue;
        }
        // Source membership, not a numeric PS inherited from attachment,
        // certifies each site of the transferred prefix.
        if (!row.bam_injected || orientation(ci) != gauge) return std::nullopt;
        run_beg = std::min(run_beg, pos);
        if (row.counts.category == VariantCategory::CleanHetSnp &&
            row.key.type == VariantType::Snp && row.key.ref_len == 1 &&
            row.key.alt.size() == 1 && pos > snp_pos) {
            snp_i = ci;
            snp_pos = pos;
        }
    }
    if (!snp_i || run_beg <= graph_end || run_beg >= pair_pos)
        return std::nullopt;
    const auto cuts = gc.recovery_source_weak_cuts.find(source_ps);
    if (cuts == gc.recovery_source_weak_cuts.end()) return std::nullopt;
    const auto crosses = [run_beg, pair_pos](const auto& positions) {
        return std::any_of(positions.begin(), positions.end(),
                          [run_beg, pair_pos](hts_pos_t pos) {
                              return run_beg <= pos && pos < pair_pos;
                          });
    };
    const auto quality = gc.recovery_source_quality_cuts.find(source_ps);
    if (crosses(cuts->second) ||
        (quality != gc.recovery_source_quality_cuts.end() && crosses(quality->second)))
        return std::nullopt;
    std::array<int, 2> support{};
    std::unordered_set<std::string> counted;
    for (size_t ri = 0; ri < chunk.read_var_profile.size() && ri < chunk.reads.size(); ++ri) {
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        const ReadRecord& read = chunk.reads[ri];
        if (read.is_skipped || profile.bam_mapq < kMinMapq ||
            profile.bam_mapq == kUnknownMapq || profile.start_var_idx < 0 ||
            counted.count(read.qname) != 0) continue;
        const auto call = [&](size_t ci) {
            if (ci < static_cast<size_t>(profile.start_var_idx)) return -1;
            const size_t offset = ci - static_cast<size_t>(profile.start_var_idx);
            return offset < profile.bam_alleles.size() ?
                static_cast<int>(profile.bam_alleles[offset]) : -1;
        };
        const int snp_call = call(*snp_i);
        const int a_call = call(first_i);
        const int b_call = call(second_i);
        // Double-REF and double-ALT do not identify a deletion length class.
        if ((snp_call != 0 && snp_call != 1) ||
            !((a_call == 1 && b_call == 0) || (a_call == 0 && b_call == 1))) continue;
        counted.insert(read.qname);
        const CandidateVariant& deletion = a_call == 1 ? a : b;
        if ((snp_call == chunk.candidates[*snp_i].hap_to_cons_alle[1]) !=
            (deletion.hap_to_cons_alle[1] == 1)) return std::nullopt;
        ++support[static_cast<size_t>(snp_call)];
    }
    if (support[0] < kMinAlleleClassSupport || support[1] < kMinAlleleClassSupport ||
        std::ldexp(1.0, -(support[0] + support[1])) > kMaxWrongParity)
        return std::nullopt;
    return run_beg;
}

// A source path certifies its retained component, even after attachment
// changes its PS. An independent physical bridge may admit exact shared
// catalog rows; every anchor still needs the same original source gauge.
bool bam_source_run_supported(const GraphChunkBuildResult& gc,
                              hts_pos_t phase_set,
                              bool include_shared_graph) {
    const PhasingChunk& chunk = gc.chunk;
    std::vector<const RecoverySourceSite*> source_by_candidate(
        chunk.candidates.size(), nullptr);
    for (const RecoverySourceSite& site : gc.recovery_source_sites) {
        if (site.candidate_index >= chunk.candidates.size() ||
            chunk.candidates[site.candidate_index].phase_set != phase_set)
            continue;
        const CandidateVariant& candidate = chunk.candidates[site.candidate_index];
        if (!is_phase_set_anchor(candidate) && candidate.hap_to_cons_alle[1] >= 0 &&
            candidate.hap_to_cons_alle[2] >= 0)
            continue;
        const RecoverySourceSite*& slot =
            source_by_candidate[site.candidate_index];
        if (slot != nullptr) return false;
        slot = &site;
    }

    hts_pos_t source_ps = 0;
    hts_pos_t first = std::numeric_limits<hts_pos_t>::max();
    hts_pos_t last = 0;
    int orientation = -1;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set != phase_set) continue;
        // BAM homozygotes retain a source PS but provide no orientation.
        // Excluding them also keeps the cut interval at the actual anchors.
        if (!is_phase_set_anchor(candidate) && candidate.hap_to_cons_alle[1] >= 0 &&
            candidate.hap_to_cons_alle[2] >= 0)
            continue;
        const RecoverySourceSite* found = source_by_candidate[ci];
        if ((!candidate.bam_injected && !include_shared_graph) ||
            found == nullptr) return false;
        const RecoverySourceSite& site = *found;
        if (!site.can_adopt || site.phase_set <= 0 ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[2] < 0 ||
            candidate.hap_to_cons_alle[1] == candidate.hap_to_cons_alle[2] ||
            site.hap1_allele < 0 || site.hap2_allele < 0 ||
            site.hap1_allele == site.hap2_allele ||
            (source_ps > 0 && source_ps != site.phase_set))
            return false;
        const bool same =
            candidate.hap_to_cons_alle[1] == site.hap1_allele &&
            candidate.hap_to_cons_alle[2] == site.hap2_allele;
        const bool swapped =
            candidate.hap_to_cons_alle[1] == site.hap2_allele &&
            candidate.hap_to_cons_alle[2] == site.hap1_allele;
        if (!same && !swapped) return false;
        const int row_orientation = swapped;
        if (orientation >= 0 && orientation != row_orientation) return false;
        orientation = row_orientation;
        source_ps = site.phase_set;
        const hts_pos_t pos = candidate.key.sort_pos();
        first = std::min(first, pos);
        last = std::max(last, pos);
    }
    if (source_ps <= 0 || first >= last) return false;
    const auto cuts = gc.recovery_source_weak_cuts.find(source_ps);
    if (cuts == gc.recovery_source_weak_cuts.end()) return false;
    const auto crosses = [first, last](const auto& positions) {
        return std::any_of(positions.begin(), positions.end(),
                           [first, last](hts_pos_t cut) {
                               return first <= cut && cut < last;
                           });
    };
    const auto quality = gc.recovery_source_quality_cuts.find(source_ps);
    return !crosses(cuts->second) &&
        (quality == gc.recovery_source_quality_cuts.end() ||
         !crosses(quality->second));
}

bool graph_snp_ref_absence_supported(int ref_count, int alt_count,
        int deletion_count, int other_count, size_t tested_sites) {
    constexpr int kMinDeletionObservations = 10;
    constexpr double kMinDeletionFraction = 0.2;
    constexpr double kFamilywiseError = 0.01;
    const int callable = ref_count + alt_count;
    const int total = callable + deletion_count;
    return ref_count == 0 && other_count == 0 && total > 0 && tested_sites > 0 &&
        (deletion_count == 0 || (deletion_count >= kMinDeletionObservations &&
         static_cast<double>(deletion_count) / total >= kMinDeletionFraction)) &&
        std::ldexp(1.0, -callable) * tested_sites <= kFamilywiseError;
}

bool graph_snp_padded_deletion_supported(int ref_count, int alt_count,
        int deletion_count, int other_count, size_t tested_sites) {
    constexpr int kMinDeletionObservations = 2;
    constexpr double kMinDeletionFraction = 0.2;
    constexpr double kFamilywiseError = 0.01;
    const int total = alt_count + deletion_count;
    return ref_count == 0 && other_count == 0 && tested_sites > 0 &&
        deletion_count >= kMinDeletionObservations && total > 0 &&
        static_cast<double>(deletion_count) / total >= kMinDeletionFraction &&
        std::ldexp(1.0, -alt_count) * tested_sites <= kFamilywiseError;
}

bool graph_snp_low_alt_fraction_supported(int ref_count, int alt_count,
        int deletion_count, int other_count, double min_af, size_t tested_sites) {
    constexpr double kFamilywiseError = 0.01;
    const int total = ref_count + alt_count;
    return ref_count > 0 && alt_count >= 0 && deletion_count == 0 &&
        other_count == 0 && tested_sites > 0 && min_af > 0.0 && min_af < 0.5 &&
        static_cast<double>(alt_count) / total < min_af &&
        rescue_binomial_tail(ref_count, total) * tested_sites <= kFamilywiseError;
}

std::optional<int> physical_deletion_gauge_haplotype(
        const std::array<int, 2>& hap_counts, double wrong_gauge_bound) {
    constexpr int kMinPairs = 2;
    constexpr double kMaxWrongGauge = 0.001;
    if (!std::isfinite(wrong_gauge_bound) || wrong_gauge_bound < 0.0 ||
        wrong_gauge_bound > kMaxWrongGauge || hap_counts[0] < 0 || hap_counts[1] < 0 ||
        (hap_counts[0] > 0 && hap_counts[1] > 0)) return std::nullopt;
    if (hap_counts[0] >= kMinPairs) return 1;
    if (hap_counts[1] >= kMinPairs) return 2;
    return std::nullopt;
}

bool graph_snp_cohort_is_physically_contradicted(
        const std::optional<GraphSnpReferenceEvidence>& terminal,
        const std::optional<GraphSnpReferenceEvidence>& partner) {
    constexpr int kMinContradictedAltReads = 2;
    constexpr double kMaxWrongAlternate = 0.001;
    return terminal && partner && terminal->reference_class_reads > 0 &&
        terminal->alternate_class_deletions >= kMinContradictedAltReads &&
        partner->reference_class_reads > 0 &&
        partner->alternate_class_reference_reads >= kMinContradictedAltReads &&
        partner->wrong_alternate_bound <= kMaxWrongAlternate;
}

} // namespace pgphase_collect
