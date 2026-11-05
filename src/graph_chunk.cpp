#include "trace.hpp"
#include "collect_var.hpp"
#include <map>
#include "graph_chunk.hpp"

#include "collect_phase.hpp"
#include "collect_phase_noisy.hpp"
#include "fisher_exact.hpp"
#include "noise_filter.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <limits>
#include <numeric>
#include <set>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
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
// Allele counts are initially zero; they are filled after observation
// counting.
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
    candidate.graph_site_key = key;
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
        if (cand.key.type == VariantType::Snp) {
            c.category = VariantCategory::CleanHetSnp;
            c.candvarcate_initial = VariantCategory::CleanHetSnp;
            cand.lcd_var_i_to_cate = kCandCleanHetSnp;
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
            if (pgphase_trace::on() && cand.counts.category != VariantCategory::RepeatHetIndel)
                pgphase_trace::emit("graph", "regenotype", meta.pos,
                                    meta.pos + std::max<hts_pos_t>(1, meta.ref.size()) - 1,
                                    meta.ref + ">" + (meta.alts.empty() ? "?" : meta.alts[0].substr(0, 60)),
                                    "repeat_het_indel_not_phased",
                                    "cat=" + std::to_string(static_cast<int>(cand.counts.category)));
            cand.counts.category = VariantCategory::RepeatHetIndel;
            cand.counts.candvarcate_initial = VariantCategory::RepeatHetIndel;
            cand.lcd_var_i_to_cate = kLongcalldRepHetVar;
        }
    }
}

namespace {

// One graph chunk built from catalog sites and GAF rows, stage by stage (see
// build_graph_chunk for the order). The per-site tables are indexed by the
// chunk's candidates and are rebuilt with them when sites are decomposed into
// biallelic pairs.
class GraphChunkBuilder {
public:
    GraphChunkBuilder(const GraphSiteCatalogView& catalog, const std::vector<GraphReadAllele>& rows,
                      hts_pos_t beg, hts_pos_t end, int chunk_id, const Options& opts)
        : catalog_(catalog), rows_(rows), beg_(beg), end_(end), chunk_id_(chunk_id), opts_(opts) {}

    GraphChunkBuildResult build() {
        init_chunk();
        add_eligible_sites();
        const std::vector<std::vector<SiteAlleleEntry>> allele_by_read_site = collect_read_alleles();
        count_observations(allele_by_read_site);
        set_site_counts();
        drop_weak_alleles();
        decompose_sites();
        // Classify candidates before building read profiles so that pruned sites
        // (LowCoverage, StrandBias) never enter the profile / cr_overlap index.
        classify_graph_candidates(out_.chunk, opts_);
        remap_observations();
        add_read_profiles();
        out_.chunk.haps.assign(out_.chunk.reads.size(), 0);
        out_.chunk.phase_sets.assign(out_.chunk.reads.size(), kUnphasedReadPhaseSet);
        rebuild_read_var_cr(out_.chunk);
        out_.site_allele_orig_idx = std::move(allele_orig_idx_);
        return std::move(out_);
    }

private:
    static constexpr int kChunkTid = 0;
    // Allele and strand packed in one int: low 16 bits the allele, bit 16
    // reverse. A conflict is -1 (allele < 0).
    static constexpr int kRevBit = 0x10000;
    static int pack_allele_rev(int allele, bool rev) { return allele | (rev ? kRevBit : 0); }
    static int unpack_allele(int packed) { return packed & 0xFFFF; }
    static bool unpack_reverse(int packed) { return (packed & kRevBit) != 0; }

    struct SiteAlleleEntry { int site_i; int packed; };
    // Each original site's surviving pairs after decomposition: the new
    // candidate and the Phase 1 alt it holds.
    struct NewPairEntry { int new_idx; int old_alt_phase1; };

    const GraphSiteCatalogView& catalog_;
    const std::vector<GraphReadAllele>& rows_;
    const hts_pos_t beg_, end_;
    const int chunk_id_;
    const Options& opts_;

    GraphChunkBuildResult out_;
    std::unordered_map<std::string, int> site_to_candidate_;
    std::vector<std::vector<int>> allele_counts_;
    // Per-site, per-allele strand counts -- accumulated from raw rows before dedup.
    std::vector<std::vector<int>> fwd_strand_counts_;
    std::vector<std::vector<int>> rev_strand_counts_;
    // Reads interned to contiguous ids.
    std::vector<std::string> read_id_to_name_;
    std::vector<int> read_max_mapq_;
    std::vector<std::vector<GraphProfileObservation>> read_obs_;
    // Phase 1-2 bookkeeping over the original sites.
    size_t n_cands_ = 0;
    std::vector<std::vector<int>> allele_remap_;
    std::vector<std::vector<int>> allele_orig_idx_;
    std::vector<std::vector<NewPairEntry>> old_site_to_pairs_;

    void init_chunk() {
        out_.chunk.region.chunk_id = chunk_id_;
        out_.chunk.region.reg_chunk_i = chunk_id_;
        out_.chunk.region.tid = kChunkTid;
        out_.chunk.region.beg = beg_ + 1;
        out_.chunk.region.end = end_;
        out_.chunk.ref_beg = beg_ + 1;
        out_.chunk.ref_end = end_;
        out_.chunk.chunk_min_qual = 60;
        out_.chunk.chunk_first_quar_qual = 60;
        out_.chunk.chunk_median_qual = 60;
        out_.chunk.chunk_third_quar_qual = 60;
        out_.chunk.chunk_max_qual = 60;
        out_.chunk.up_ovlp_read_i.assign(1, {});
        out_.chunk.down_ovlp_read_i.assign(1, {});
        out_.chunk.n_up_ovlp_skip_reads.assign(1, 0);
        out_.chunk.n_down_ovlp_skip_reads.assign(1, 0);
    }

    // 1. Each eligible, queryable catalog site with two or more allele walks
    //    becomes one multi-allelic candidate; the others are filtered sites.
    void add_eligible_sites() {
        for (size_t site_i = 0; site_i < catalog_.size(); ++site_i) {
            const GraphSite& site = catalog_[site_i];
            const std::string sid = graph_site_key_str(site);
            if (!site.eligible) {
                out_.filtered_sites.push_back({sid, site.ref_contig.empty() ? site.chrom : site.ref_contig,
                                               site.pos, 0, 0, 0, 0.0, "precandidate_ineligible"});
                continue;
            }
            const bool released_walk_storage = graph_site_has_released_walk_storage(site);
            if (!released_walk_storage && !graph_site_is_queryable(site)) {
                out_.filtered_sites.push_back({sid, site.ref_contig.empty() ? site.chrom : site.ref_contig,
                                               site.pos, 0, 0, 0, 0.0, "precandidate_not_queryable"});
                continue;
            }
            const int n_alleles = static_cast<int>(site.allele_walks.size());
            if (n_alleles < 2) {
                out_.filtered_sites.push_back({sid, site.ref_contig.empty() ? site.chrom : site.ref_contig,
                                               site.pos, 0, 0, 0, 0.0, "precandidate_monoallelic"});
                continue;
            }
            site_to_candidate_.emplace(sid, static_cast<int>(out_.chunk.candidates.size()));
            allele_counts_.emplace_back(static_cast<size_t>(n_alleles), 0);
            fwd_strand_counts_.emplace_back(static_cast<size_t>(n_alleles), 0);
            rev_strand_counts_.emplace_back(static_cast<size_t>(n_alleles), 0);
            add_graph_candidate(out_, site, sid, allele_counts_.back(), kChunkTid);
        }
    }

    // 2. One pass over the rows: intern read names, track each read's best
    //    MAPQ and collect its allele per site, sorted by site; a read seen at
    //    one site with two different alleles is a conflict there (-1).
    std::vector<std::vector<SiteAlleleEntry>> collect_read_alleles() {
        std::unordered_map<std::string, uint32_t> read_name_to_id;
        read_name_to_id.reserve(rows_.size());
        const auto intern_read = [&](const std::string& name) -> uint32_t {
            auto [it, inserted] = read_name_to_id.emplace(name, static_cast<uint32_t>(read_id_to_name_.size()));
            if (inserted) read_id_to_name_.push_back(name);
            return it->second;
        };
        // Temporary per-row (rid, site_i, packed) triples; sorted into per-read
        // vectors after the pass when n_reads is known.
        struct RowTriple { uint32_t rid; int site_i; int packed; int mapq; };
        std::vector<RowTriple> row_triples;
        row_triples.reserve(rows_.size());
        // Counted because a row dropped here leaves its site with no evidence at all,
        // which surfaces much later as a "ref_only" filter that names neither cause.
        int64_t row_site_missing = 0, row_allele_oob = 0, row_kept = 0;
        for (const GraphReadAllele& row : rows_) {
            const uint32_t rid = intern_read(row.read_name);
            auto it = site_to_candidate_.find(row.site_id);
            if (it == site_to_candidate_.end()) { ++row_site_missing; continue; }
            const int site_i = it->second;
            const CandidateVariant& candidate = out_.chunk.candidates[static_cast<size_t>(site_i)];
            if (row.allele < 0 || row.allele >= candidate.counts.n_uniq_alles) { ++row_allele_oob; continue; }
            ++row_kept;
            row_triples.push_back({rid, site_i, pack_allele_rev(row.allele, row.reverse), row.mapq});
        }
        if (opts_.verbose >= 1) {
            std::fprintf(stderr, "[rows] total %zu kept %" PRId64 " site-unknown %" PRId64
                                 " allele-out-of-range %" PRId64 "\n",
                         rows_.size(), row_kept, row_site_missing, row_allele_oob);
        }
        const uint32_t n_reads = static_cast<uint32_t>(read_id_to_name_.size());
        read_max_mapq_.assign(n_reads, 0);
        std::vector<std::vector<SiteAlleleEntry>> allele_by_read_site(n_reads);
        for (const RowTriple& t : row_triples) {
            if (t.mapq > read_max_mapq_[t.rid]) read_max_mapq_[t.rid] = t.mapq;
            allele_by_read_site[t.rid].push_back({t.site_i, t.packed});
        }
        // Sort each read's entries by site_i, then deduplicate (mark conflicts as -1).
        for (auto& entries : allele_by_read_site) {
            if (entries.empty()) continue;
            std::sort(entries.begin(), entries.end(),
                      [](const SiteAlleleEntry& a, const SiteAlleleEntry& b) { return a.site_i < b.site_i; });
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
        if (opts_.verbose >= 1) {
            int64_t conflicted = 0, kept_obs = 0;
            for (const auto& entries : allele_by_read_site)
                for (const auto& e : entries) (e.packed < 0 ? conflicted : kept_obs)++;
            std::fprintf(stderr, "[dedup] observations kept %" PRId64 " conflicted %" PRId64 "\n",
                         kept_obs, conflicted);
        }
        return allele_by_read_site;
    }

    // 3. Per-read observation lists and per-site allele and strand counts.
    void count_observations(const std::vector<std::vector<SiteAlleleEntry>>& allele_by_read_site) {
        const uint32_t n_reads = static_cast<uint32_t>(read_id_to_name_.size());
        read_obs_.assign(n_reads, {});
        for (uint32_t rid = 0; rid < n_reads; ++rid) {
            for (const auto& entry : allele_by_read_site[rid]) {
                const int site_i = entry.site_i;
                const int packed = entry.packed;
                if (packed < 0) continue;
                const int allele = unpack_allele(packed);
                const bool rev = unpack_reverse(packed);
                ++allele_counts_[static_cast<size_t>(site_i)][static_cast<size_t>(allele)];
                auto& sc = rev ? rev_strand_counts_[static_cast<size_t>(site_i)]
                               : fwd_strand_counts_[static_cast<size_t>(site_i)];
                ++sc[static_cast<size_t>(allele)];
                read_obs_[rid].push_back(GraphProfileObservation{site_i, allele});
            }
        }
    }

    // The counts onto each candidate.
    void set_site_counts() {
        for (size_t site_i = 0; site_i < out_.chunk.candidates.size(); ++site_i) {
            CandidateVariant& candidate = out_.chunk.candidates[site_i];
            const std::vector<int>& counts = allele_counts_[site_i];
            const std::vector<int>& fwd = fwd_strand_counts_[site_i];
            const std::vector<int>& rev = rev_strand_counts_[site_i];
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
    }

    // 4. Phase 1 of the candidate filtering (three phases mirroring BAM §13
    //    het classification): drop individual alt alleles below min_alt_depth
    //    (noise walks). allele_remap_[site][old allele] = new allele (or -1
    //    if dropped); per-site counts recomputed from surviving alleles only;
    //    allele_orig_idx_ maps each surviving allele back to its allele walk
    //    (callers read GraphSite.alts[orig - 1]).
    void drop_weak_alleles() {
        n_cands_ = out_.chunk.candidates.size();
        allele_remap_.assign(n_cands_, {});
        for (size_t i = 0; i < n_cands_; ++i) {
            CandidateVariant& cand = out_.chunk.candidates[i];
            std::vector<int>& ac = allele_counts_[i];
            std::vector<int>& fc = fwd_strand_counts_[i];
            std::vector<int>& rc = rev_strand_counts_[i];
            allele_remap_[i].assign(ac.size(), -1);
            allele_remap_[i][0] = 0;  // ref walk always kept at index 0
            std::vector<int> new_ac = {ac[0]};
            std::vector<int> new_fc = {fc[0]};
            std::vector<int> new_rc = {rc[0]};
            int next_allele = 1;
            for (size_t a = 1; a < ac.size(); ++a) {
                if (ac[a] >= opts_.min_alt_depth) {
                    allele_remap_[i][a] = next_allele++;
                    new_ac.push_back(ac[a]);
                    new_fc.push_back(a < fc.size() ? fc[a] : 0);
                    new_rc.push_back(a < rc.size() ? rc[a] : 0);
                }
            }
            {
                static const char* dbg = getenv("PGPHASE_DEBUG_SITE");
                static const long dbg_pos = dbg ? atol(dbg) : -1;
                if (dbg_pos >= 0 && i < out_.site_meta.size() &&
                    std::llabs(static_cast<long long>(out_.site_meta[i].pos) - dbg_pos) <= 2) {
                    std::string before, after;
                    for (int x : ac) before += std::to_string(x) + ",";
                    for (int x : new_ac) after += std::to_string(x) + ",";
                    std::fprintf(stderr, "[compact %ld] cate=0x%x min_alt_depth=%d counts before=[%s] after=[%s]\n",
                                 static_cast<long>(out_.site_meta[i].pos), cand.lcd_var_i_to_cate,
                                 opts_.min_alt_depth, before.c_str(), after.c_str());
                }
            }
            allele_counts_[i] = new_ac;
            fwd_strand_counts_[i] = new_fc;
            rev_strand_counts_[i] = new_rc;
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
        allele_orig_idx_.assign(n_cands_, {});
        for (size_t i = 0; i < n_cands_; ++i) {
            allele_orig_idx_[i].push_back(0);  // new index 0 = ref = orig index 0
            for (size_t old_a = 1; old_a < allele_remap_[i].size(); ++old_a) {
                if (allele_remap_[i][old_a] >= 0) {
                    allele_orig_idx_[i].push_back(static_cast<int>(old_a));
                }
            }
        }
    }

    // The candidate tables Phase 2 builds, swapped in when it is done.
    struct NewTables {
        std::vector<CandidateVariant> cands;
        std::vector<std::string> ids;
        std::vector<GraphSiteMeta> meta;
        std::vector<std::vector<int>> allele_counts;
        std::vector<std::vector<int>> fwd_strand;
        std::vector<std::vector<int>> rev_strand;
        std::vector<std::vector<int>> orig_idx;
    };

    // 5. Phase 2: decompose each surviving site into biallelic (ref vs alt_i)
    //    pairs, mirroring the BAM path where every candidate is a single
    //    ref/alt pair. AF and depth filters are applied per pair; filtered
    //    pairs go to out.filtered_sites. A read observing ref contributes
    //    allele 0 to every pair from that site, and a read observing alt_i
    //    allele 1 to its own pair only.
    void decompose_sites() {
        old_site_to_pairs_.assign(n_cands_, {});

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
        canonical_by_fp.reserve(n_cands_ * 2);

        NewTables t;
        // Most sites are biallelic (1 surviving alt); reserve n_cands as a lower bound.
        t.cands.reserve(n_cands_);
        t.ids.reserve(n_cands_);
        t.meta.reserve(n_cands_);
        t.allele_counts.reserve(n_cands_);
        t.fwd_strand.reserve(n_cands_);
        t.rev_strand.reserve(n_cands_);
        t.orig_idx.reserve(n_cands_);

        for (size_t i = 0; i < n_cands_; ++i) {
            const std::vector<int>& ac = allele_counts_[i];
            if (ac.size() < 2) {
                // Chunks overlap, so a site near a boundary is built in two chunks and
                // sees reads in only one of them.  The empty copy is not a filtered
                // site -- it is an artifact of chunking -- and labelling it "ref_only"
                // makes the dump read as though the site had no alt evidence anywhere.
                // Distinguish the two so the dump can be deduplicated honestly.
                const char* reason = (ac[0] == 0) ? "no_reads_in_chunk" : "ref_only";
                out_.filtered_sites.push_back({out_.site_ids[i], out_.site_meta[i].chrom,
                                               out_.site_meta[i].pos, ac[0], 0, ac[0], 0.0, reason});
                continue;
            }
            const std::vector<int>& fc = fwd_strand_counts_[i];
            const std::vector<int>& rc_s = rev_strand_counts_[i];
            const int ref_c = ac[0];
            const std::string& meta_chrom = out_.site_meta[i].chrom;
            const hts_pos_t meta_pos = out_.site_meta[i].pos;
            for (size_t a = 1; a < ac.size(); ++a) {
                const int alt_c   = ac[a];
                const int pair_ref_c = ref_c;
                const int total_c = ref_c + alt_c;
                const double af   = total_c > 0 ? static_cast<double>(alt_c) / total_c : 0.0;
                const int orig_alt = allele_orig_idx_[i][a];

                const std::string pair_id = (ac.size() > 2)
                    ? out_.site_ids[i] + ":" + std::to_string(orig_alt)
                    : out_.site_ids[i];

                std::string drop_reason;
                if (total_c < opts_.min_depth) {
                    drop_reason = "low_depth";
                } else if (af > opts_.max_af) {
                    drop_reason = "high_af";
                } else if (af < opts_.min_af) {
                    drop_reason = "low_af";
                }
                if (!drop_reason.empty()) {
                    if (pgphase_trace::on()) {
                        const GraphSiteMeta& m = out_.site_meta[i];
                        const size_t ai = static_cast<size_t>(orig_alt - 1);
                        pgphase_trace::emit("graph", "drop", m.pos, m.pos + std::max<hts_pos_t>(1, m.ref.size()) - 1,
                                            m.ref + ">" + (ai < m.alts.size() ? m.alts[ai].substr(0, 60) : "?"),
                                            "pair_" + drop_reason,
                                            "pair=" + pair_id + " ref_alt=" + std::to_string(pair_ref_c) + "/" +
                                                std::to_string(alt_c));
                    }
                    out_.filtered_sites.push_back({pair_id, meta_chrom, meta_pos,
                                                   pair_ref_c, alt_c, total_c, af, std::move(drop_reason)});
                    continue;
                }

                const int fwd_ref = fc.empty() ? 0 : fc[0];
                const int rev_ref = rc_s.empty() ? 0 : rc_s[0];
                const int fwd_alt = a < fc.size() ? fc[a] : 0;
                const int rev_alt = a < rc_s.size() ? rc_s[a] : 0;

                // Compute the minimal-VCF identity (pos, ref, alt) of this pair so
                // duplicates from overlapping snarls share a key. Identity views point
                // into meta.ref / meta.alts[alt_idx], which outlive this loop.
                const GraphSiteMeta& meta = out_.site_meta[i];
                const size_t alt_idx = static_cast<size_t>(orig_alt - 1);
                const std::string& id_alt_src = alt_idx < meta.alts.size() ? meta.alts[alt_idx] : meta.ref;
                const MinimalVcfId identity = minimal_vcf_id(meta.pos, meta.ref, id_alt_src);
                const uint64_t fp = minimal_vcf_fingerprint(identity);

                int canon = -1;
                std::vector<CanonEntry>& bucket = canonical_by_fp[fp];
                for (const CanonEntry& e : bucket) {
                    if (e.old_i != i && e.id == identity) { canon = e.new_idx; break; }
                }
                if (canon >= 0) {
                    if (pgphase_trace::on())
                        pgphase_trace::emit("graph", "merge", meta.pos,
                                            meta.pos + std::max<hts_pos_t>(1, meta.ref.size()) - 1,
                                            meta.ref + ">" + id_alt_src.substr(0, 60), "duplicate_snarl_routed",
                                            "pair=" + pair_id);
                    // Duplicate variant from an overlapping snarl: route this
                    // snarl's observations to the canonical anchor so the variant
                    // is not double-counted in k-means, but keep the canonical
                    // candidate's original counts: pooling coverage across
                    // overlapping snarls over-weights these anchors and
                    // net-regresses Hamming. The per-read observation dedup in
                    // Phase 3 still collapses a read seen at both snarls into one
                    // vote.
                    old_site_to_pairs_[i].push_back({canon, static_cast<int>(a)});
                    continue;
                }

                const int new_idx = static_cast<int>(t.cands.size());
                bucket.push_back({identity, new_idx, i});
                old_site_to_pairs_[i].push_back({new_idx, static_cast<int>(a)});

                CandidateVariant pair_cand = out_.chunk.candidates[i];
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
                t.cands.push_back(pair_cand);
                t.ids.push_back(pair_id);
                t.meta.push_back(out_.site_meta[i]);
                t.allele_counts.push_back({ref_c, alt_c});
                t.fwd_strand.push_back({fwd_ref, fwd_alt});
                t.rev_strand.push_back({rev_ref, rev_alt});
                t.orig_idx.push_back({0, orig_alt});
            }
        }

        out_.chunk.candidates = std::move(t.cands);
        out_.site_ids         = std::move(t.ids);
        out_.site_meta        = std::move(t.meta);
        allele_counts_        = std::move(t.allele_counts);
        fwd_strand_counts_    = std::move(t.fwd_strand);
        rev_strand_counts_    = std::move(t.rev_strand);
        allele_orig_idx_      = std::move(t.orig_idx);
    }

    // 6. Phase 3: remap read observations to the new biallelic pair space.
    //    Allele 0 (ref) fans out to all new pairs from its original site
    //    (allele 0 in each). Allele j (alt) maps to allele 1 in the one pair
    //    that holds alt j. Observations at pruned candidates are dropped.
    void remap_observations() {
        const size_t n_final_cands = out_.chunk.candidates.size();
        std::vector<bool> cand_pruned(n_final_cands, false);
        for (size_t ci = 0; ci < n_final_cands; ++ci) {
            const VariantCategory cat = out_.chunk.candidates[ci].counts.category;
            if (cat == VariantCategory::LowCoverage || cat == VariantCategory::StrandBias ||
                cat == VariantCategory::NonVariant) {
                cand_pruned[ci] = true;
            }
        }
        for (auto& obs : read_obs_) {
            std::vector<GraphProfileObservation> new_obs;
            new_obs.reserve(obs.size());
            for (const auto& o : obs) {
                if (o.site_index < 0 || static_cast<size_t>(o.site_index) >= n_cands_) continue;
                const size_t old_si = static_cast<size_t>(o.site_index);
                const auto& pairs = old_site_to_pairs_[old_si];
                if (pairs.empty()) continue;
                // Fast path: biallelic site (one surviving pair) -- no loop needed.
                if (pairs.size() == 1) {
                    if (cand_pruned[static_cast<size_t>(pairs[0].new_idx)]) continue;
                    if (o.allele == 0) {
                        new_obs.push_back({pairs[0].new_idx, 0});
                    } else if (o.allele > 0 && static_cast<size_t>(o.allele) < allele_remap_[old_si].size()) {
                        if (allele_remap_[old_si][static_cast<size_t>(o.allele)] == pairs[0].old_alt_phase1)
                            new_obs.push_back({pairs[0].new_idx, 1});
                    }
                    continue;
                }
                // General path: multiallelic site with multiple surviving pairs.
                if (o.allele == 0) {
                    for (const NewPairEntry& pe : pairs)
                        if (!cand_pruned[static_cast<size_t>(pe.new_idx)]) new_obs.push_back({pe.new_idx, 0});
                } else if (o.allele > 0 && static_cast<size_t>(o.allele) < allele_remap_[old_si].size()) {
                    const int phase1_alt = allele_remap_[old_si][static_cast<size_t>(o.allele)];
                    if (phase1_alt >= 0) {
                        for (const NewPairEntry& pe : pairs) {
                            if (cand_pruned[static_cast<size_t>(pe.new_idx)]) continue;
                            if (pe.old_alt_phase1 == phase1_alt) new_obs.push_back({pe.new_idx, 1});
                        }
                    }
                }
            }
            obs = std::move(new_obs);
        }
    }

    // 7. One read profile per read with observations, in read-name order for
    //    deterministic output; a read seen at one candidate with two alleles
    //    loses that call.
    void add_read_profiles() {
        const uint32_t n_reads = static_cast<uint32_t>(read_id_to_name_.size());
        std::vector<uint32_t> read_order(n_reads);
        std::iota(read_order.begin(), read_order.end(), 0);
        std::sort(read_order.begin(), read_order.end(),
                  [&](uint32_t a, uint32_t b) { return read_id_to_name_[a] < read_id_to_name_[b]; });
        for (uint32_t rid : read_order) {
            auto& observations = read_obs_[rid];
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
            add_read_profile(out_, read_id_to_name_[rid], dedup, read_max_mapq_[rid]);
        }
    }
};

}  // namespace

// Main entry point: convert graph-space allele observations into a PhasingChunk
// (GraphChunkBuilder, one method per step):
//   1. Enumerate eligible catalog sites -> initial multi-allelic candidates
//   2. Single pass over rows: intern read names, collect per-read allele obs
//   3. Observation counting
//   4. Phase 1: drop individual alt alleles below min_alt_depth
//   5. Phase 2: decompose multi-allelic sites into biallelic ref/alt pairs,
//      apply depth + AF filters per pair
//   6. Classify surviving candidates (het/hom/low-cov/strand-bias), then
//      Phase 3: remap read observations to biallelic pair space, drop pruned
//   7. Build ReadVariantProfile + cgranges index for k-means
GraphChunkBuildResult build_graph_chunk(const GraphSiteCatalogView& catalog,
                                        const std::vector<GraphReadAllele>& rows,
                                        const std::string& /*contig*/,
                                        hts_pos_t beg,
                                        hts_pos_t end,
                                        int chunk_id,
                                        const Options& opts) {
    return GraphChunkBuilder(catalog, rows, beg, end, chunk_id, opts).build();
}

} // namespace pgphase_collect
