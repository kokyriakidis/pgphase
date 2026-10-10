// Union gap phasing (--union-gap-phasing): the catalog's sites are phased
// together with the alignment's sample-specific heterozygotes -- clean calls
// and MSA-verified noisy calls the catalog lacks -- by EM over the whole
// read x site matrix, plus local haplotype windows over noisy loci.
// See docs/IMPLEMENTATION.md, "Union gap phasing".

#include "union_phase.hpp"

#include "collect_phase.hpp"
#include "collect_var.hpp"
#include "edlib.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <htslib/sam.h>

namespace pgphase_collect {

namespace {

using CandKey = CandidateIdentityKey;

CandKey cand_key_of(const CandidateVariant& cand) {
    return CandKey{cand.key.sort_pos(), static_cast<int>(cand.key.type), cand.key.ref_len,
                   cand.key.alt};
}

}  // namespace

// A site is admitted only if each of its two alleles is seen by at least this
// many reads mapped at kDefaultMinMapq or better.
constexpr int kUnionMinConfidentAlleleReads = 2;
// A hidden complement site takes calls only from reads the graph left silent.
constexpr bool kComplementSilentReadsOnly = true;

/// VCF-form metadata for an alignment candidate the catalog never held, in the
/// convention the graph writer emits and vcf_to_variant_key reads back. An
/// empty REF makes the writer skip the record while index-parallel arrays stay
/// aligned.
static std::pair<GraphSiteMeta, std::vector<int>> synthesize_alignment_site_meta(
        const CandidateVariant& cand, const PhasingChunk& source, const char* contig_name) {
    const auto ref_at = [&source](hts_pos_t pos, int len) -> std::string {
        if (len <= 0 || source.ref_seq.empty()) return std::string();
        const hts_pos_t off = pos - source.ref_beg;
        if (off < 0 || static_cast<size_t>(off) + static_cast<size_t>(len) > source.ref_seq.size())
            return std::string();
        return source.ref_seq.substr(static_cast<size_t>(off), static_cast<size_t>(len));
    };
    const VariantKey& key = cand.key;
    const auto anchored = [&](const std::string& allele) -> std::string {
        if (key.type == VariantType::Snp) return allele;
        if (key.type == VariantType::Insertion && key.ref_len != 0) return allele;
        return ref_at(key.pos - 1, 1) + allele;
    };
    hts_pos_t vcf_pos = key.pos;
    std::string vcf_ref;
    if (key.type == VariantType::Snp) {
        vcf_ref = ref_at(key.pos, std::max(1, key.ref_len));
    } else if (key.type == VariantType::Insertion && key.ref_len == 0) {
        vcf_pos = key.pos - 1;
        vcf_ref = ref_at(vcf_pos, 1);
    } else if (key.type == VariantType::Insertion) {
        vcf_ref = ref_at(key.pos, key.ref_len);
    } else {
        vcf_pos = key.pos - 1;
        vcf_ref = ref_at(vcf_pos, 1) + ref_at(key.pos, key.ref_len);
    }
    const std::string vcf_alt = anchored(key.alt);
    bool usable = !vcf_ref.empty() && !vcf_alt.empty();
    if (usable) {
        const VariantKey round_trip = vcf_to_variant_key(key.tid, vcf_pos, vcf_ref, vcf_alt);
        usable = round_trip.type == key.type && round_trip.pos == key.pos &&
                 round_trip.ref_len == key.ref_len && round_trip.alt == key.alt;
    }
    std::vector<std::string> alts;
    if (cand.msa_insertion_alts.size() >= 2) {
        for (const std::string& allele : cand.msa_insertion_alts) {
            const std::string a = anchored(allele);
            if (a.empty()) { alts.clear(); break; }
            alts.push_back(a);
        }
    }
    if (alts.empty() && usable) alts.push_back(vcf_alt);
    GraphSiteMeta meta;
    meta.chrom = contig_name != nullptr ? contig_name : std::string();
    if (usable && !alts.empty()) {
        meta.pos = vcf_pos;
        meta.ref = vcf_ref;
        meta.alts = alts;
    }
    std::vector<int> orig;
    for (int i = 0; i <= static_cast<int>(alts.size()); ++i) orig.push_back(i);
    if (orig.size() < 2) orig = std::vector<int>{0, 1};
    return {std::move(meta), std::move(orig)};
}

// Left-aligned identity of a simple insertion or deletion, so placements of
// one allele at different repeat copies compare equal. Other rows keep their key.
static CandKey left_aligned_indel_key(const CandidateVariant& c, const PhasingChunk& source) {
    VariantKey k = c.key;
    const auto base = [&source](hts_pos_t pos) -> char {
        const hts_pos_t off = pos - source.ref_beg;
        if (off < 0 || static_cast<size_t>(off) >= source.ref_seq.size()) return 'N';
        return static_cast<char>(std::toupper(static_cast<unsigned char>(source.ref_seq[static_cast<size_t>(off)])));
    };
    if (!c.msa_insertion_alts.empty()) return cand_key_of(c);
    if (k.type == VariantType::Insertion && k.ref_len == 0 && !k.alt.empty()) {
        // The inserted bases sit after pos - 1; shift while the base before
        // the insertion point equals the last inserted base.
        while (true) {
            const char prev = base(k.pos - 1);
            if (prev == 'N' || std::toupper(static_cast<unsigned char>(k.alt.back())) != prev) break;
            k.alt = std::string(1, prev) + k.alt.substr(0, k.alt.size() - 1);
            --k.pos;
        }
    } else if (k.type == VariantType::Deletion && k.ref_len > 0 && k.alt.empty()) {
        while (true) {
            const char prev = base(k.pos - 1);
            if (prev == 'N' || prev != base(k.pos + k.ref_len - 1)) break;
            --k.pos;
        }
    }
    return CandKey{k.sort_pos(), static_cast<int>(k.type), k.ref_len, k.alt};
}

size_t inject_alignment_sites(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                              const char* contig_name) {
    PhasingChunk& graph = graph_chunk.chunk;

    // Rows that already vote describe this locus. An alignment call overlapping
    // one is another description of it: it joins as a hidden complement site
    // (no VCF record) carrying the alignment's calls, because the graph often
    // leaves reads at such a locus without any allele.
    std::vector<std::tuple<hts_pos_t, hts_pos_t, size_t>> voting_spans;
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
        const CandidateVariant& row = graph.candidates[ci];
        if ((row.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        voting_spans.emplace_back(row.key.pos - 1, row.key.pos + std::max(1, row.key.ref_len), ci);
    }
    std::sort(voting_spans.begin(), voting_spans.end());
    hts_pos_t widest = 0;
    for (const auto& span : voting_spans) widest = std::max(widest, std::get<1>(span) - std::get<0>(span));
    const auto voting_rows_over = [&](hts_pos_t beg, hts_pos_t end) {
        std::vector<size_t> rows;
        auto it = std::lower_bound(voting_spans.begin(), voting_spans.end(),
                                   std::make_tuple(beg - widest, hts_pos_t{0}, size_t{0}));
        for (; it != voting_spans.end() && std::get<0>(*it) < end; ++it)
            if (std::get<1>(*it) > beg) rows.push_back(std::get<2>(*it));
        return rows;
    };
    std::unordered_map<size_t, std::vector<size_t>> complement_rows;  // bam candidate -> graph rows it overlaps

    // Reads admitted at low MAPQ can represent a divergent haplotype but cannot
    // create a site on their own: each of the site's two alleles (a 1/2 site
    // carries no REF reads) needs confidently mapped support.
    std::vector<std::map<int, int>> confident(bam.candidates.size());
    std::vector<std::map<int, int>> observed(bam.candidates.size());
    for (const ReadVariantProfile& prof : bam.read_var_profile) {
        if (prof.read_id < 0 || prof.start_var_idx < 0) continue;
        const ReadRecord& read = bam.reads[static_cast<size_t>(prof.read_id)];
        if (read.is_skipped) continue;
        for (size_t k = 0; k < prof.alleles.size(); ++k) {
            const int a = prof.alleles[k];
            if (a < 0) continue;
            const size_t ci = static_cast<size_t>(prof.start_var_idx) + k;
            ++observed[ci][a];
            if (read.mapq >= kDefaultMinMapq) ++confident[ci][a];
        }
    }
    const auto confidently_supported = [&](size_t bi) {
        std::vector<std::pair<int, int>> top;
        for (const auto& [a, n] : observed[bi]) top.emplace_back(n, a);
        std::sort(top.rbegin(), top.rend());
        if (top.size() < 2) return false;
        for (size_t i = 0; i < 2; ++i) {
            const auto it = confident[bi].find(top[i].second);
            if (it == confident[bi].end() || it->second < kUnionMinConfidentAlleleReads) return false;
        }
        return true;
    };

    // Clean heterozygotes and MSA-verified noisy heterozygotes, keeping the
    // alignment's own category bits.
    std::vector<size_t> selected;
    for (size_t bi = 0; bi < bam.candidates.size(); ++bi) {
        const CandidateVariant& c = bam.candidates[bi];
        const VariantCategory cat = c.counts.category;
        const bool clean = cat == VariantCategory::CleanHetSnp || cat == VariantCategory::CleanHetIndel;
        const bool msa = cat == VariantCategory::NoisyCandHet && c.msa_verified;
        if (!clean && !msa) continue;
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        if (c.key.pos < graph.ref_beg || c.key.pos > graph.ref_end) continue;
        if (!confidently_supported(bi)) continue;
        std::vector<size_t> rows = voting_rows_over(c.key.pos - 1, c.key.pos + std::max(1, c.key.ref_len));
        if (!rows.empty()) complement_rows.emplace(bi, std::move(rows));
        selected.push_back(bi);
    }
    if (selected.empty()) return 0;

    // One site per left-aligned allele: other placements become aliases whose
    // read calls fold into the canonical row.
    std::unordered_map<size_t, size_t> alias_of;
    {
        std::map<CandKey, size_t> canonical;
        std::vector<size_t> kept;
        for (const size_t bi : selected) {
            const auto [it, inserted] = canonical.emplace(left_aligned_indel_key(bam.candidates[bi], bam), bi);
            if (inserted) kept.push_back(bi);
            else alias_of.emplace(bi, it->second);
        }
        selected = std::move(kept);
    }

    // Merge in key order; every index-parallel array is rebuilt with the table.
    struct Slot { CandKey key; long graph_index; long bam_index; };
    std::vector<Slot> slots;
    slots.reserve(graph.candidates.size() + selected.size());
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci)
        slots.push_back(Slot{cand_key_of(graph.candidates[ci]), static_cast<long>(ci), -1});
    for (const size_t bi : selected) {
        CandidateVariant relabeled = bam.candidates[bi];
        relabeled.key.tid = graph.region.tid;
        slots.push_back(Slot{cand_key_of(relabeled), -1, static_cast<long>(bi)});
    }
    std::stable_sort(slots.begin(), slots.end(),
                     [](const Slot& a, const Slot& b) { return a.key < b.key; });
    const bool had_meta = graph_chunk.site_meta.size() == graph.candidates.size();
    CandidateTable merged;
    std::vector<std::string> merged_ids;
    std::vector<GraphSiteMeta> merged_meta;
    std::vector<std::vector<int>> merged_orig;
    std::vector<long> graph_to_new(graph.candidates.size(), -1);
    std::unordered_map<size_t, size_t> bam_to_new;
    for (const Slot& slot : slots) {
        const size_t ni = merged.size();
        if (slot.graph_index >= 0) {
            const size_t ci = static_cast<size_t>(slot.graph_index);
            graph_to_new[ci] = static_cast<long>(ni);
            merged.push_back(graph.candidates[ci]);
            merged_ids.push_back(ci < graph_chunk.site_ids.size() ? graph_chunk.site_ids[ci] : std::string());
            merged_meta.push_back(had_meta ? graph_chunk.site_meta[ci] : GraphSiteMeta{});
            merged_orig.push_back(ci < graph_chunk.site_allele_orig_idx.size()
                                      ? graph_chunk.site_allele_orig_idx[ci] : std::vector<int>{});
            continue;
        }
        const size_t bi = static_cast<size_t>(slot.bam_index);
        bam_to_new.emplace(bi, ni);
        CandidateVariant cand = bam.candidates[bi];
        cand.key.tid = graph.region.tid;
        cand.bam_injected = true;
        cand.alignment_verified = true;
        // Phase comes from the union's own solve, not the alignment's gauge.
        cand.phase_set = kUnsetCandidatePhaseSet;
        cand.hap_to_cons_alle[1] = -1;
        cand.hap_to_cons_alle[2] = -1;
        cand.hap_alt = 0;
        cand.hap_ref = 0;
        if (cand.counts.alle_covs.size() < 2)
            cand.counts.alle_covs = {cand.counts.ref_cov, cand.counts.alt_cov};
        auto [meta, orig] = synthesize_alignment_site_meta(bam.candidates[bi], bam, contig_name);
        if (complement_rows.count(bi) != 0) {
            meta = GraphSiteMeta{};  // empty REF: the writer emits no record
            meta.chrom = contig_name != nullptr ? contig_name : std::string();
        }
        merged.push_back(std::move(cand));
        merged_ids.emplace_back();
        merged_meta.push_back(std::move(meta));
        merged_orig.push_back(std::move(orig));
    }

    // Alignment observations at the injected sites, keyed by read name. Reads
    // the GAF never had are added: they are what carries a catalog gap.
    std::unordered_map<std::string_view, size_t> graph_read_by_qname;
    graph_read_by_qname.reserve(graph.reads.size());
    for (size_t ri = 0; ri < graph.reads.size(); ++ri)
        graph_read_by_qname.try_emplace(graph.reads[ri].qname, ri);
    std::vector<std::vector<std::tuple<size_t, int, int>>> added_obs(graph.reads.size());
    const size_t old_profiles = graph.reads.size();
    std::unordered_map<size_t, const std::vector<size_t>*> complement_of_new;  // merged index -> graph rows
    for (const auto& [bi, rows] : complement_rows) {
        const auto it = bam_to_new.find(bi);
        if (it != bam_to_new.end()) complement_of_new.emplace(it->second, &rows);
    }
    for (const ReadVariantProfile& profile : bam.read_var_profile) {
        if (profile.read_id < 0 || profile.start_var_idx < 0) continue;
        const ReadRecord& source_read = bam.reads[static_cast<size_t>(profile.read_id)];
        if (source_read.is_skipped) continue;
        std::vector<std::tuple<size_t, int, int>> calls;
        for (size_t k = 0; k < profile.alleles.size(); ++k) {
            if (profile.alleles[k] < 0) continue;
            size_t source_i = static_cast<size_t>(profile.start_var_idx) + k;
            const auto alias = alias_of.find(source_i);
            if (alias != alias_of.end()) source_i = alias->second;
            const auto found = bam_to_new.find(source_i);
            if (found == bam_to_new.end()) continue;
            const int qi = k < profile.alt_qi.size() ? profile.alt_qi[k] : 0;
            // An alias and its canonical row describe one allele: ALT wins.
            const auto prior = std::find_if(calls.begin(), calls.end(),
                [&](const auto& call) { return std::get<0>(call) == found->second; });
            if (prior == calls.end())
                calls.emplace_back(found->second, profile.alleles[k], qi);
            else if (profile.alleles[k] > std::get<1>(*prior))
                *prior = std::make_tuple(found->second, profile.alleles[k], qi);
        }
        if (calls.empty()) continue;
        auto read_i = graph_read_by_qname.find(source_read.qname);
        size_t ri;
        if (read_i != graph_read_by_qname.end()) {
            ri = read_i->second;
        } else {
            ReadRecord read;
            read.tid = graph.region.tid;
            read.qname = source_read.qname;
            read.mapq = source_read.mapq;
            read.beg = source_read.beg;
            read.end = source_read.end;
            read.reverse = source_read.reverse;
            ri = graph.reads.size();
            graph.reads.push_back(std::move(read));
            graph.read_var_profile.push_back(ReadVariantProfile{});
            added_obs.emplace_back();
            graph_read_by_qname.try_emplace(graph.reads.back().qname, ri);
        }
        for (auto& call : calls) {
            const auto comp = complement_of_new.find(std::get<0>(call));
            if (kComplementSilentReadsOnly && comp != complement_of_new.end() && ri < old_profiles) {
                // The read speaks through the graph row already.
                const ReadVariantProfile& old = graph.read_var_profile[ri];
                bool graph_call = false;
                for (const size_t row : *comp->second) {
                    if (old.start_var_idx < 0 || static_cast<int>(row) < old.start_var_idx ||
                        static_cast<int>(row) > old.end_var_idx) continue;
                    if (old.alleles[row - static_cast<size_t>(old.start_var_idx)] >= 0) graph_call = true;
                }
                if (graph_call) continue;
            }
            added_obs[ri].push_back(call);
        }
    }

    // Rebuild each profile over the merged index space, every channel included.
    for (size_t ri = 0; ri < graph.reads.size(); ++ri) {
        const ReadVariantProfile& old = graph.read_var_profile[ri];
        std::map<size_t, std::array<int, 5>> by_index;  // allele, qi, graph, bam, bam_qi
        const auto slot = [&by_index](size_t i) -> std::array<int, 5>& {
            return by_index.try_emplace(i, std::array<int, 5>{-1, 0, -1, -1, 0}).first->second;
        };
        std::map<size_t, uint8_t> qualities;
        if (old.start_var_idx >= 0) {
            for (size_t k = 0; k < old.alleles.size(); ++k) {
                const size_t ci = static_cast<size_t>(old.start_var_idx) + k;
                if (ci >= graph_to_new.size() || graph_to_new[ci] < 0) continue;
                const size_t ni = static_cast<size_t>(graph_to_new[ci]);
                const bool any = old.alleles[k] >= 0 ||
                    (k < old.graph_alleles.size() && old.graph_alleles[k] >= 0) ||
                    (k < old.bam_alleles.size() && old.bam_alleles[k] != -1);
                if (!any) continue;
                auto& s = slot(ni);
                s[0] = old.alleles[k];
                s[1] = k < old.alt_qi.size() ? old.alt_qi[k] : 0;
                s[2] = k < old.graph_alleles.size() ? old.graph_alleles[k] : -1;
                s[3] = k < old.bam_alleles.size() ? old.bam_alleles[k] : -1;
                s[4] = k < old.bam_qi.size() ? old.bam_qi[k] : 0;
                if (k < old.bam_base_qualities.size() && old.bam_base_qualities[k] > 0)
                    qualities.emplace(ni, old.bam_base_qualities[k]);
            }
        }
        for (const auto& [ni, allele, qi] : added_obs[ri]) {
            auto& s = slot(ni);
            s[0] = allele;
            s[1] = qi;
            s[3] = allele;
            s[4] = qi;
        }
        ReadVariantProfile prof;
        prof.read_id = static_cast<int>(ri);
        prof.bam_mapq = old.bam_mapq;
        if (!by_index.empty()) {
            prof.start_var_idx = static_cast<int>(by_index.begin()->first);
            prof.end_var_idx = static_cast<int>(by_index.rbegin()->first);
            const size_t span = static_cast<size_t>(prof.end_var_idx - prof.start_var_idx + 1);
            prof.alleles.assign(span, -1);
            prof.alt_qi.assign(span, 0);
            const bool graph_channel = !old.graph_alleles.empty();
            const bool bam_channel = !old.bam_alleles.empty() || !added_obs[ri].empty();
            if (graph_channel) prof.graph_alleles.assign(span, -1);
            if (bam_channel) { prof.bam_alleles.assign(span, -1); prof.bam_qi.assign(span, 0); }
            if (!qualities.empty()) prof.bam_base_qualities.assign(span, 0);
            for (const auto& [ni, s] : by_index) {
                const size_t off = ni - static_cast<size_t>(prof.start_var_idx);
                prof.alleles[off] = s[0];
                prof.alt_qi[off] = s[1];
                if (graph_channel) prof.graph_alleles[off] = s[2];
                if (bam_channel) { prof.bam_alleles[off] = s[3]; prof.bam_qi[off] = s[4]; }
            }
            for (const auto& [ni, q] : qualities)
                prof.bam_base_qualities[ni - static_cast<size_t>(prof.start_var_idx)] = q;
        }
        graph.read_var_profile[ri] = std::move(prof);
    }
    graph.candidates = std::move(merged);
    graph_chunk.site_ids = std::move(merged_ids);
    graph_chunk.site_meta = std::move(merged_meta);
    graph_chunk.site_allele_orig_idx = std::move(merged_orig);
    graph_chunk.candidate_allele_contrasts.assign(graph.candidates.size(), std::nullopt);
    graph.haps.resize(graph.reads.size(), 0);
    graph.phase_sets.resize(graph.reads.size(), kUnphasedReadPhaseSet);
    rebuild_read_var_cr(graph);
    return selected.size();
}
// ── Global EM ───────────────────────────────────────────────────────────────

// Error-rate bounds; starting error for clean and other sites; the learned
// error at which a site is neither reported phased nor used to label reads;
// iteration limits; the read posterior needed for a label; the boundary
// log-likelihood below which blocks are cut; and the MAPQ from which a read
// teaches the model (lower reads are labelled from it but do not shape it).
constexpr double kEmMinError = 0.01;
constexpr double kEmMaxError = 0.45;
constexpr double kEmCleanStartError = 0.05;
constexpr double kEmOtherStartError = 0.2;
constexpr double kEmUnreliableError = 0.3;
constexpr int kEmIterations = 15;
constexpr int kEmSwitchRounds = 30;
constexpr double kEmReadLabelPosterior = 0.8;
constexpr double kEmBlockCutLogLikelihood = 4.0;
// The cut for joins carried only by weak sites (see bridge_weak_sites).
constexpr double kEmWeakBridgeCutLogLikelihood = 10.0;
constexpr int kEmLearnMinMapq = 20;
// A site whose less observed allele carries under this fraction of its reads
// does not take part: a skewed split is the signature of a paralog or error
// call, and in a homozygous stretch such a site labels reads at random.
constexpr double kEmMinMinorAlleleFraction = 0.25;

static double log_sum_exp(double a, double b) {
    const double m = std::max(a, b);
    return m + std::log(std::exp(a - m) + std::exp(b - m));
}

static void set_site_phase(CandidateVariant& c, hts_pos_t ps, int hap1_allele, int hap2_allele) {
    c.phase_set = ps;
    c.hap_to_cons_alle[1] = hap1_allele;
    c.hap_to_cons_alle[2] = hap2_allele;
    const bool h1 = hap1_allele != 0, h2 = hap2_allele != 0;
    c.hap_alt = h1 && h2 ? 3 : h1 ? 1 : h2 ? 2 : 0;
    c.hap_ref = h1 && h2 ? 0 : h1 ? 2 : h2 ? 1 : 0;
}

size_t phase_chunk_by_global_em(PhasingChunk& chunk, const std::vector<LocusWindowSite>* loci,
                                const std::vector<char>* bridge_weak) {
    struct EmSite { long locus; size_t ci; int x; int y; hts_pos_t pos; };  // locus < 0: candidate
    struct EmObs { int s; int b; double e; };
    // Sites in position order with their two-allele contrast; a locus window
    // contrasts its two local haplotype sequences (0 and 1).
    std::vector<std::map<int, int>> allele_counts(chunk.candidates.size());
    for (size_t ri = 0; ri < chunk.read_var_profile.size(); ++ri) {
        const ReadVariantProfile& prof = chunk.read_var_profile[ri];
        if (prof.start_var_idx < 0 || chunk.reads[ri].is_skipped) continue;
        for (size_t k = 0; k < prof.alleles.size(); ++k)
            if (prof.alleles[k] >= 0)
                ++allele_counts[static_cast<size_t>(prof.start_var_idx) + k][prof.alleles[k]];
    }
    std::vector<EmSite> sites;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& c = chunk.candidates[ci];
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        int x = c.hap_to_cons_alle[1], y = c.hap_to_cons_alle[2];
        if (x < 0 || y < 0 || x == y) {
            std::vector<std::pair<int, int>> by_count;
            for (const auto& [a, n] : allele_counts[ci]) by_count.emplace_back(n, a);
            std::sort(by_count.rbegin(), by_count.rend());
            if (by_count.size() < 2) continue;
            x = by_count[0].second;
            y = by_count[1].second;
        }
        const auto count_of = [&](int a) {
            const auto it = allele_counts[ci].find(a);
            return it == allele_counts[ci].end() ? 0 : it->second;
        };
        const int nx = count_of(x), ny = count_of(y);
        if (nx + ny == 0 || std::min(nx, ny) < kEmMinMinorAlleleFraction * (nx + ny)) continue;
        sites.push_back(EmSite{-1, ci, x, y, c.key.sort_pos()});
    }
    const size_t n_locus = loci != nullptr ? loci->size() : 0;
    for (size_t li = 0; li < n_locus; ++li)
        sites.push_back(EmSite{static_cast<long>(li), 0, 0, 1, (*loci)[li].pos});
    std::stable_sort(sites.begin(), sites.end(),
                     [](const EmSite& a, const EmSite& b) { return a.pos < b.pos; });
    const size_t n = sites.size();
    if (n < 2) return 0;
    std::vector<int> site_of(chunk.candidates.size(), -1);
    std::vector<int> locus_site(n_locus, -1);
    for (size_t k = 0; k < n; ++k) {
        if (sites[k].locus >= 0) locus_site[static_cast<size_t>(sites[k].locus)] = static_cast<int>(k);
        else site_of[sites[k].ci] = static_cast<int>(k);
    }
    std::vector<int> phase(n, 0);  // 0: hap1 carries x
    std::vector<double> err(n);
    for (size_t k = 0; k < n; ++k) {
        const bool clean = sites[k].locus < 0 &&
            (chunk.candidates[sites[k].ci].lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0;
        err[k] = clean ? kEmCleanStartError : kEmOtherStartError;
    }
    // Read observations as (site, side, own error) in site order.
    std::vector<std::vector<EmObs>> robs(chunk.reads.size());
    for (size_t ri = 0; ri < chunk.read_var_profile.size(); ++ri) {
        const ReadVariantProfile& prof = chunk.read_var_profile[ri];
        if (prof.start_var_idx < 0 || chunk.reads[ri].is_skipped) continue;
        for (size_t k = 0; k < prof.alleles.size(); ++k) {
            const int s = site_of[static_cast<size_t>(prof.start_var_idx) + k];
            if (s < 0 || prof.alleles[k] < 0) continue;
            const EmSite& site = sites[static_cast<size_t>(s)];
            if (prof.alleles[k] == site.x) robs[ri].push_back(EmObs{s, 0, 0.0});
            else if (prof.alleles[k] == site.y) robs[ri].push_back(EmObs{s, 1, 0.0});
        }
    }
    for (size_t li = 0; li < n_locus; ++li)
        for (const auto& [ri, side, e] : (*loci)[li].observations)
            if (ri < robs.size()) robs[ri].push_back(EmObs{locus_site[li], side, static_cast<double>(e)});
    std::vector<double> post(chunk.reads.size(), 0.5);
    std::vector<char> learns(chunk.reads.size(), 0);
    for (size_t ri = 0; ri < robs.size(); ++ri) {
        std::sort(robs[ri].begin(), robs[ri].end(),
                  [](const EmObs& a, const EmObs& b) { return a.s < b.s; });
        const int hap = ri < chunk.haps.size() ? chunk.haps[ri] : 0;
        post[ri] = hap == 1 ? 0.95 : hap == 2 ? 0.05 : 0.5;
        learns[ri] = chunk.reads[ri].mapq >= kEmLearnMinMapq;
    }
    const auto obs_error = [&](const EmObs& o) {
        return std::min(kEmMaxError, 1.0 - (1.0 - err[static_cast<size_t>(o.s)]) * (1.0 - o.e));
    };
    const auto m_step = [&]() {
        std::vector<std::array<double, 2>> w(n, {0.0, 0.0});  // hap1 weight per side
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            if (!learns[ri]) continue;
            for (const EmObs& o : robs[ri]) {
                w[static_cast<size_t>(o.s)][static_cast<size_t>(o.b)] += post[ri];
                w[static_cast<size_t>(o.s)][static_cast<size_t>(1 - o.b)] += 1.0 - post[ri];
            }
        }
        for (size_t k = 0; k < n; ++k) {
            const double total = w[k][0] + w[k][1];
            if (total < 1e-6) continue;
            phase[k] = w[k][0] >= w[k][1] ? 0 : 1;
            err[k] = std::min(kEmMaxError, std::max(kEmMinError, std::min(w[k][0], w[k][1]) / total));
        }
    };
    const auto read_log_likelihoods = [&](const EmObs& o) -> std::array<double, 2> {
        const double e = obs_error(o);
        const bool match1 = o.b == phase[static_cast<size_t>(o.s)];
        return {std::log(match1 ? 1.0 - e : e), std::log(match1 ? e : 1.0 - e)};
    };
    const auto e_step = [&]() {
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            if (robs[ri].empty()) continue;
            double l1 = 0.0, l2 = 0.0;
            for (const EmObs& o : robs[ri]) {
                const auto ll = read_log_likelihoods(o);
                l1 += ll[0];
                l2 += ll[1];
            }
            const double d = l2 - l1;
            post[ri] = d > 700 ? 0.0 : d < -700 ? 1.0 : 1.0 / (1.0 + std::exp(d));
        }
    };
    // delta[k] = log L(current) - log L(every site right of k flipped), over
    // the reads that teach the model or, for block cuts, over every read: a
    // read too ambiguously mapped to shape the model still shows the molecule
    // continues across the boundary.
    // Weak sites (see bridge_weak_sites) are phased but cannot hold blocks
    // together: the final cuts are judged on the other sites' evidence.
    std::vector<char> site_weak(n, 0);
    if (bridge_weak != nullptr)
        for (size_t k = 0; k < n; ++k)
            if (sites[k].locus < 0 && sites[k].ci < bridge_weak->size()) site_weak[k] = (*bridge_weak)[sites[k].ci];
    const auto boundary_delta = [&](bool every_read, bool strong_only = false) {
        std::vector<double> diff(n + 1, 0.0);
        std::vector<EmObs> strong;
        for (size_t ri = 0; ri < robs.size(); ++ri) {
            const std::vector<EmObs>* items_ptr = &robs[ri];
            if (strong_only) {
                strong.clear();
                for (const EmObs& o : robs[ri])
                    if (!site_weak[static_cast<size_t>(o.s)]) strong.push_back(o);
                items_ptr = &strong;
            }
            const auto& items = *items_ptr;
            if ((!every_read && !learns[ri]) || items.size() < 2 || items.front().s == items.back().s)
                continue;
            std::vector<std::array<double, 2>> ll(items.size());
            double tot1 = 0.0, tot2 = 0.0;
            for (size_t i = 0; i < items.size(); ++i) {
                ll[i] = read_log_likelihoods(items[i]);
                tot1 += ll[i][0];
                tot2 += ll[i][1];
            }
            double c1 = 0.0, c2 = 0.0;
            for (size_t i = 0; i + 1 < items.size(); ++i) {
                c1 += ll[i][0];
                c2 += ll[i][1];
                const int k = items[i].s, k_next = items[i + 1].s;
                if (k_next == k) continue;
                const double gain = log_sum_exp(tot1, tot2) -
                                    log_sum_exp(c1 + (tot2 - c2), c2 + (tot1 - c1));
                diff[static_cast<size_t>(k)] += gain;
                diff[static_cast<size_t>(k_next)] -= gain;
            }
        }
        std::vector<double> delta(n, 0.0);
        double run = 0.0;
        for (size_t k = 0; k < n; ++k) {
            run += diff[k];
            delta[k] = run;
        }
        return delta;
    };
    // EM, then switch moves: flip everything right of the worst boundary while
    // that raises the likelihood.
    std::vector<double> delta;
    for (int round = 0; round < kEmSwitchRounds; ++round) {
        for (int it = 0; it < kEmIterations; ++it) {
            m_step();
            e_step();
        }
        delta = boundary_delta(false);
        size_t worst = n;
        for (size_t k = 0; k + 1 < n; ++k)
            if (delta[k] < -1e-6 && (worst == n || delta[k] < delta[worst])) worst = k;
        if (worst == n) break;
        for (size_t k = worst + 1; k < n; ++k) phase[k] ^= 1;
        for (size_t ri = 0; ri < robs.size(); ++ri)
            if (!robs[ri].empty() && static_cast<size_t>(robs[ri].front().s) > worst)
                post[ri] = 1.0 - post[ri];
    }
    delta = boundary_delta(true);
    // A join that the non-weak sites alone do not carry (its margin comes from
    // weak sites: homopolymer length calls) needs a larger margin. Measured on
    // chr20: wrong joins of this kind sit just above the ordinary cut
    // (median 6), right ones well above it (median 19).
    const std::vector<double> delta_strong =
        bridge_weak != nullptr ? boundary_delta(true, true) : delta;
    // Blocks: cut where flipping the rest would cost less than the threshold.
    std::vector<hts_pos_t> block(n);
    hts_pos_t current = sites[0].pos;
    size_t blocks = 1;
    for (size_t k = 0; k < n; ++k) {
        block[k] = current;
        const bool weak_bridge = delta_strong[k] < kEmBlockCutLogLikelihood;
        if (k + 1 < n && (delta[k] < kEmBlockCutLogLikelihood ||
                          (weak_bridge && delta[k] < kEmWeakBridgeCutLogLikelihood))) {
            current = sites[k + 1].pos;
            ++blocks;
        }
    }
    for (CandidateVariant& c : chunk.candidates)
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0) c.phase_set = kUnsetCandidatePhaseSet;
    for (size_t k = 0; k < n; ++k) {
        if (sites[k].locus >= 0) {
            const LocusWindowSite& w = (*loci)[static_cast<size_t>(sites[k].locus)];
            w.phase = phase[k];
            w.block = block[k];
            w.error = err[k];
            continue;
        }
        CandidateVariant& c = chunk.candidates[sites[k].ci];
        if (err[k] >= kEmUnreliableError) {
            c.hap_to_cons_alle[1] = c.hap_to_cons_alle[2] = -1;
            c.hap_alt = c.hap_ref = 0;
            continue;
        }
        const int hap1 = phase[k] == 0 ? sites[k].x : sites[k].y;
        const int hap2 = phase[k] == 0 ? sites[k].y : sites[k].x;
        set_site_phase(c, block[k], hap1, hap2);
    }
    // Reads: label in the block where their own evidence is strongest.
    for (size_t ri = 0; ri < robs.size() && ri < chunk.haps.size(); ++ri) {
        chunk.haps[ri] = 0;
        chunk.phase_sets[ri] = kUnphasedReadPhaseSet;
        std::map<hts_pos_t, std::array<double, 2>> per_block;
        for (const EmObs& o : robs[ri]) {
            if (err[static_cast<size_t>(o.s)] >= kEmUnreliableError) continue;
            const auto ll = read_log_likelihoods(o);
            auto& l = per_block[block[static_cast<size_t>(o.s)]];
            l[0] += ll[0];
            l[1] += ll[1];
        }
        double best = 0.0;
        for (const auto& [ps, l] : per_block) {
            const double margin = l[0] - l[1];
            const double p1 = margin > 700 ? 1.0 : margin < -700 ? 0.0 : 1.0 / (1.0 + std::exp(-margin));
            if (std::max(p1, 1.0 - p1) < kEmReadLabelPosterior || std::abs(margin) <= best) continue;
            best = std::abs(margin);
            chunk.haps[ri] = p1 >= 0.5 ? 1 : 2;
            chunk.phase_sets[ri] = ps;
        }
    }
    return blocks;
}

// ── Local haplotype windows ─────────────────────────────────────────────────

// Noisy calls closer than this merge into one locus; padding; the longest
// window; the MAPQ a read needs to take part; reads per window; seed sample
// per side; minimum reads per side and minor fraction; the edit-distance margin
// needed to assign a read; repeat scan bounds; the largest per-read error.
constexpr hts_pos_t kLocusMergeGap = 100;
constexpr hts_pos_t kLocusPad = 40;
constexpr hts_pos_t kLocusMaxLength = 2000;
constexpr int kLocusMinMapq = 5;
constexpr size_t kLocusMaxReads = 60;
constexpr size_t kLocusSeedSample = 10;
constexpr int kLocusMinSide = 3;
constexpr double kLocusMinMinorFraction = 0.2;
constexpr int kLocusMinMargin = 2;
constexpr int kLocusSplitRounds = 3;  // refinements of an unseeded two-way split
// Allele windows: flank around the repeat, and reads per window. A read is
// assigned at any edit margin: both sequences are exact alleles, so a one-edit
// difference is the allele itself, not a read's sequencing error.
constexpr hts_pos_t kAlleleWindowFlank = 20;
constexpr size_t kAlleleWindowMaxReads = 100;
// An indel is a candidate for a two-allele pair when at most this fraction of
// the reads at it carry the reference.
constexpr double kPairMaxRefFraction = 0.2;
// Merge two-allele rows into one EM site (measured: more phased variants, more
// flips, fewer correct reads -- off).
constexpr bool kMergeTwoAlleleRows = false;
// Keep the MSA's own per-read calls at injected indel sites.
constexpr bool kRealignKeepsMsaCalls = true;
// A site is a weak bridge when fewer than this fraction of its covering reads call it.
constexpr double kBridgeMinCalledFraction = 0.5;
constexpr bool kBridgeDepthRule = false;
constexpr bool kBridgeHomopolymerRule = true;
// Alignment reads missing from the graph chunk are added from this MAPQ:
// confidently placed, so not a paralog the graph put elsewhere.
constexpr int kAddedReadMinMapq = 20;
// Last-resort labels: an MSA-verified indel orients itself from at least this
// many labelled reads of one phase set, at this concordance.
constexpr int kLastResortMinLabelled = 4;
constexpr double kLastResortMinAgreement = 0.9;
// Haplotype-consensus labels: an event marks a haplotype when each side has at
// least this many reads covering it (with this flank) and it is on at least
// this fraction of one side's reads and at most the complement of the other's.
// A read is labelled when this share of its votes agrees.
constexpr int kConsensusMinReads = 3;
constexpr hts_pos_t kConsensusFlank = 10;
constexpr double kConsensusCarrierFraction = 0.8;
constexpr double kConsensusMinAgreement = 0.9;
constexpr int kEventMinBaseQuality = 10;
constexpr int kLocusMaxRepeatScan = 500;
constexpr int kLocusMaxPeriod = 50;
constexpr double kLocusMaxObsError = 0.45;
// A read's label for seeding: decisive agreement with one haplotype's
// consensus at a phase set's phased sites.
constexpr int kSeedLabelMinMargin = 2;
constexpr double kSeedLabelMinAgreement = 0.9;

static int locus_edit_distance(const std::string& a, const std::string& b) {
    EdlibAlignResult r = edlibAlign(a.data(), static_cast<int>(a.size()), b.data(),
                                    static_cast<int>(b.size()),
                                    edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0));
    const int d = r.editDistance;
    edlibFreeAlignResult(r);
    return d;
}

static const std::string* locus_medoid(const std::vector<const std::string*>& seqs) {
    const std::string* best = nullptr;
    long best_cost = -1;
    for (size_t i = 0; i < seqs.size(); ++i) {
        long cost = 0;
        for (size_t j = 0; j < seqs.size(); ++j)
            if (j != i) cost += locus_edit_distance(*seqs[i], *seqs[j]);
        if (best_cost < 0 || cost < best_cost) { best = seqs[i]; best_cost = cost; }
    }
    return best;
}

// Query index aligned at (or, inside a deletion, after) a 0-based reference position.
static hts_pos_t locus_query_index(const bam1_t* b, hts_pos_t ref_pos) {
    if (ref_pos < b->core.pos || ref_pos >= bam_endpos(b)) return -1;
    hts_pos_t r = b->core.pos, q = 0;
    const uint32_t* cigar = bam_get_cigar(b);
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const hts_pos_t len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (r + len > ref_pos) return q + (ref_pos - r);
            r += len;
            q += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            if (r + len > ref_pos) return q;
            r += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            q += len;
        }
    }
    return -1;
}

// (haplotype, phase set) per read from the chunk's phased sites.
static std::vector<std::pair<int, hts_pos_t>> site_derived_read_labels(const PhasingChunk& chunk) {
    std::vector<std::pair<int, hts_pos_t>> labels(chunk.reads.size(), {0, kUnphasedReadPhaseSet});
    for (size_t ri = 0; ri < chunk.read_var_profile.size() && ri < chunk.reads.size(); ++ri) {
        const ReadVariantProfile& prof = chunk.read_var_profile[ri];
        if (prof.start_var_idx < 0 || chunk.reads[ri].is_skipped) continue;
        std::map<hts_pos_t, std::array<int, 2>> votes;
        for (size_t k = 0; k < prof.alleles.size(); ++k) {
            const CandidateVariant& c = chunk.candidates[static_cast<size_t>(prof.start_var_idx) + k];
            const int a = prof.alleles[k];
            if (a < 0 || c.phase_set <= 0 || c.hap_to_cons_alle[1] < 0 ||
                c.hap_to_cons_alle[1] == c.hap_to_cons_alle[2]) continue;
            if (a == c.hap_to_cons_alle[1]) ++votes[c.phase_set][0];
            else if (a == c.hap_to_cons_alle[2]) ++votes[c.phase_set][1];
        }
        int best_margin = 0;
        for (const auto& [ps, v] : votes) {
            const int margin = std::abs(v[0] - v[1]);
            if (margin < kSeedLabelMinMargin ||
                std::max(v[0], v[1]) < kSeedLabelMinAgreement * (v[0] + v[1])) continue;
            if (margin > best_margin) {
                best_margin = margin;
                labels[ri] = {v[0] > v[1] ? 1 : 2, ps};
            }
        }
    }
    return labels;
}

// The alignment chunk's reference in 0-based coordinates (ref_seq is 1-based
// from ref_beg), with tandem-repeat extents for window bounds.
struct ReferenceView {
    const PhasingChunk& chunk;
    char base(hts_pos_t p0) const {
        const hts_pos_t off = p0 + 1 - chunk.ref_beg;
        if (off < 0 || static_cast<size_t>(off) >= chunk.ref_seq.size()) return 'N';
        return static_cast<char>(std::toupper(static_cast<unsigned char>(chunk.ref_seq[static_cast<size_t>(off)])));
    }
    std::string slice(hts_pos_t p0, hts_pos_t p1) const {
        std::string out;
        for (hts_pos_t p = p0; p < p1; ++p) out += base(p);
        return out;
    }
    // Length of the longest tandem repeat starting at `start` (rightward) or
    // ending just before it (leftward).
    hts_pos_t repeat_right(hts_pos_t start) const {
        hts_pos_t best = 0;
        for (int period = 1; period <= kLocusMaxPeriod; ++period) {
            hts_pos_t k = start + period;
            while (k - start < kLocusMaxRepeatScan && base(k) != 'N' && base(k) == base(k - period)) ++k;
            if (k - start >= 2 * period) best = std::max(best, k - start);
        }
        return best;
    }
    hts_pos_t repeat_left(hts_pos_t end) const {
        hts_pos_t best = 0;
        for (int period = 1; period <= kLocusMaxPeriod; ++period) {
            hts_pos_t k = end - 1 - period;
            while (end - 1 - k < kLocusMaxRepeatScan && base(k) != 'N' && base(k) == base(k + period)) --k;
            if (end - 1 - k >= 2 * period) best = std::max(best, end - 1 - k);
        }
        return best;
    }
};

// Primary alignment-chunk reads the graph chunk also holds, by start.
struct SpanningReadIndex {
    std::vector<const ReadRecord*> reads;
    std::vector<hts_pos_t> starts;
    std::vector<size_t> graph_index;

    SpanningReadIndex(const PhasingChunk& bam, const PhasingChunk& graph) {
        std::unordered_map<std::string_view, size_t> graph_read;
        graph_read.reserve(graph.reads.size());
        for (size_t ri = 0; ri < graph.reads.size(); ++ri) graph_read.try_emplace(graph.reads[ri].qname, ri);
        std::vector<std::pair<const ReadRecord*, size_t>> kept;
        for (const ReadRecord& r : bam.reads) {
            if (!r.alignment) continue;
            if ((r.alignment->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FUNMAP)) != 0)
                continue;
            if (r.alignment->core.qual < kLocusMinMapq) continue;
            const auto g = graph_read.find(r.qname);
            if (g == graph_read.end()) continue;
            kept.emplace_back(&r, g->second);
        }
        std::sort(kept.begin(), kept.end(), [](const auto& a, const auto& b) {
            return a.first->alignment->core.pos < b.first->alignment->core.pos;
        });
        for (const auto& [r, gi] : kept) {
            reads.push_back(r);
            starts.push_back(r->alignment->core.pos);
            graph_index.push_back(gi);
        }
    }

    // Reads covering the whole of [l0, l1), cut to it: (graph read, sequence).
    std::vector<std::pair<size_t, std::string>> segments(hts_pos_t l0, hts_pos_t l1, size_t max_reads) const {
        std::vector<std::pair<size_t, std::string>> segs;
        const size_t last = static_cast<size_t>(std::upper_bound(starts.begin(), starts.end(), l0) - starts.begin());
        for (size_t i = 0; i < last && segs.size() < max_reads; ++i) {
            const bam1_t* aln = reads[i]->alignment.get();
            if (bam_endpos(aln) < l1) continue;
            const hts_pos_t q0 = locus_query_index(aln, l0), q1 = locus_query_index(aln, l1 - 1);
            if (q0 < 0 || q1 < q0) continue;
            std::string seq(static_cast<size_t>(q1 - q0 + 1), 'N');
            const uint8_t* packed = bam_get_seq(aln);
            for (hts_pos_t k = 0; k <= q1 - q0; ++k)
                seq[static_cast<size_t>(k)] = seq_nt16_str[bam_seqi(packed, q0 + k)];
            segs.emplace_back(graph_index[i], std::move(seq));
        }
        return segs;
    }
};

std::vector<LocusWindowSite> build_locus_window_sites(const PhasingChunk& bam,
                                                      const PhasingChunk& graph) {
    std::vector<LocusWindowSite> out;
    if (bam.ref_seq.empty()) return out;
    const std::vector<std::pair<int, hts_pos_t>> seed_labels = site_derived_read_labels(graph);
    const ReferenceView ref{bam};
    // Loci: the alignment's noisy-region calls, merged.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> calls;
    for (const CandidateVariant& c : bam.candidates) {
        if (c.counts.category != VariantCategory::NoisyCandHet &&
            c.counts.category != VariantCategory::NoisyCandHom) continue;
        const hts_pos_t p0 = c.key.pos - 1;
        const hts_pos_t len = std::max<hts_pos_t>(
            1, std::max<hts_pos_t>(c.key.ref_len, static_cast<hts_pos_t>(c.key.alt.size())));
        calls.emplace_back(p0, p0 + len);
    }
    std::sort(calls.begin(), calls.end());
    std::vector<std::pair<hts_pos_t, hts_pos_t>> merged;
    for (const auto& [a, b] : calls) {
        if (!merged.empty() && a - merged.back().second <= kLocusMergeGap)
            merged.back().second = std::max(merged.back().second, b);
        else
            merged.emplace_back(a, b);
    }
    const SpanningReadIndex index(bam, graph);
    for (const auto& [a, b] : merged) {
        const hts_pos_t l0 = std::max<hts_pos_t>(0, a - kLocusPad);
        const hts_pos_t l1 = b + ref.repeat_right(b) + kLocusPad;
        if (l1 - l0 > kLocusMaxLength) continue;
        const std::vector<std::pair<size_t, std::string>> segs = index.segments(l0, l1, kLocusMaxReads);
        if (segs.size() < static_cast<size_t>(2 * kLocusMinSide)) continue;
        // Seeds: reads already labelled in the dominant phase set here.
        std::map<hts_pos_t, int> ps_count;
        for (const auto& [ri, seq] : segs)
            if (seed_labels[ri].first != 0) ++ps_count[seed_labels[ri].second];
        std::vector<const std::string*> side1, side2;
        if (!ps_count.empty()) {
            const hts_pos_t ps = std::max_element(ps_count.begin(), ps_count.end(),
                [](const auto& x, const auto& y) { return x.second < y.second; })->first;
            for (const auto& [ri, seq] : segs) {
                if (seed_labels[ri].second != ps) continue;
                if (seed_labels[ri].first == 1 && side1.size() < kLocusSeedSample) side1.push_back(&seq);
                if (seed_labels[ri].first == 2 && side2.size() < kLocusSeedSample) side2.push_back(&seq);
            }
        }
        const std::string* hap1 = nullptr;
        const std::string* hap2 = nullptr;
        if (side1.size() >= 2 && side2.size() >= 2) {
            hap1 = locus_medoid(side1);
            hap2 = locus_medoid(side2);
        }
        // Seeds can all come from one haplotype (the other's reads unlabelled
        // here); identical seeded medoids then say nothing, and the reads
        // themselves are split instead.
        if (hap1 == nullptr || *hap1 == *hap2) {
            hap1 = hap2 = nullptr;
            // No usable seeds: split the spanning reads themselves. The EM
            // orients the locus through the reads it shares with other sites.
            std::vector<const std::string*> pool;
            for (const auto& [ri, seq] : segs)
                if (pool.size() < 2 * kLocusSeedSample) pool.push_back(&seq);
            hap1 = locus_medoid(pool);
            int far = -1;
            for (const std::string* seq : pool) {
                const int d = locus_edit_distance(*seq, *hap1);
                if (d > far) { far = d; hap2 = seq; }
            }
            for (int round = 0; round < kLocusSplitRounds && hap2 != nullptr; ++round) {
                std::vector<const std::string*> a, b;
                for (const std::string* seq : pool)
                    (locus_edit_distance(*seq, *hap1) <= locus_edit_distance(*seq, *hap2) ? a : b).push_back(seq);
                if (a.size() < 2 || b.size() < 2) break;
                if (a.size() > kLocusSeedSample) a.resize(kLocusSeedSample);
                if (b.size() > kLocusSeedSample) b.resize(kLocusSeedSample);
                hap1 = locus_medoid(a);
                hap2 = locus_medoid(b);
            }
            if (hap2 == nullptr) continue;
        }
        if (*hap1 == *hap2) continue;  // homozygous here
        LocusWindowSite site;
        site.pos = l0 + 1;
        site.end = l1;
        int n0 = 0, n1 = 0;
        for (const auto& [ri, seq] : segs) {
            const int d1 = locus_edit_distance(seq, *hap1), d2 = locus_edit_distance(seq, *hap2);
            const int margin = std::abs(d1 - d2);
            if (margin < kLocusMinMargin) continue;
            const int side = d1 < d2 ? 0 : 1;
            (side == 0 ? n0 : n1)++;
            site.observations.emplace_back(
                ri, side, static_cast<float>(std::min(kLocusMaxObsError, 0.5 * std::exp(-margin))));
        }
        const int called = n0 + n1;
        if (std::min(n0, n1) < kLocusMinSide || std::min(n0, n1) < kLocusMinMinorFraction * called) continue;
        out.push_back(std::move(site));
    }
    return out;
}

std::vector<LocusWindowSite> build_allele_window_sites(const PhasingChunk& bam,
                                                       const GraphChunkBuildResult& graph_chunk,
                                                       const std::vector<LocusWindowSite>& taken) {
    const PhasingChunk& graph = graph_chunk.chunk;
    std::vector<LocusWindowSite> out;
    if (bam.ref_seq.empty()) return out;
    const ReferenceView ref{bam};
    const auto is_acgt = [](const std::string& s) {
        return s.find_first_not_of("ACGTacgt") == std::string::npos;
    };
    // Windows already observed by the EM: haplotype windows, and voting indel
    // rows (their pileup calls are the site's evidence already).
    std::vector<std::pair<hts_pos_t, hts_pos_t>> busy;
    for (const LocusWindowSite& w : taken) busy.emplace_back(w.pos - 1, w.end);
    for (const CandidateVariant& c : graph.candidates) {
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0 || c.key.type == VariantType::Snp) continue;
        busy.emplace_back(c.key.pos - 1, c.key.pos - 1 + std::max(1, c.key.ref_len));
    }
    std::sort(busy.begin(), busy.end());
    hts_pos_t widest = 0;
    for (const auto& b : busy) widest = std::max(widest, b.second - b.first);
    const auto overlaps_busy = [&](hts_pos_t a, hts_pos_t b) {
        auto it = std::lower_bound(busy.begin(), busy.end(), std::make_pair(a - widest, hts_pos_t{0}));
        for (; it != busy.end() && it->first < b; ++it)
            if (it->second > a) return true;
        return false;
    };
    // Heterozygous indels the alignment called (noisy-region MSA or repeat
    // context) that the EM does not otherwise observe: injected ones are now
    // voting rows and fall under `busy`.
    // Each candidate is a reference span [p0, p0 + ref_len) and the two allele
    // sequences that replace it on the two haplotypes.
    struct Allele { hts_pos_t p0; int ref_len; std::string seq0; std::string seq1; };
    std::vector<Allele> alleles;
    for (const CandidateVariant& c : bam.candidates) {
        if (c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        const VariantCategory cat = c.counts.category;
        if (cat != VariantCategory::RepeatHetIndel && cat != VariantCategory::NoisyCandHet) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        alleles.push_back(Allele{c.key.pos - 1, c.key.ref_len, ref.slice(c.key.pos - 1, c.key.pos - 1 + c.key.ref_len),
                                 c.key.alt});
    }
    // Two-allele heterozygotes (1/2): at a repeat, the two haplotypes can carry
    // two different indels and no reference allele. Called against REF one at a
    // time, each looks homozygous; together they are one heterozygous site.
    // Pair the two best-supported REF-less indels within one repeat window and
    // contrast the reference carrying one with the reference carrying the other.
    {
        struct Lone { hts_pos_t p0; int ref_len; std::string alt; int reads; hts_pos_t w0; hts_pos_t w1; };
        std::vector<Lone> lone;
        for (const CandidateVariant& c : bam.candidates) {
            if (c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
            if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
            const int alt_reads = c.counts.alt_cov, ref_reads = c.counts.ref_cov;
            if (alt_reads < kLocusMinSide || ref_reads > kPairMaxRefFraction * (alt_reads + ref_reads)) continue;
            const hts_pos_t p0 = c.key.pos - 1, end0 = p0 + c.key.ref_len;
            lone.push_back(Lone{p0, c.key.ref_len, c.key.alt, alt_reads,
                                p0 - ref.repeat_left(p0), end0 + ref.repeat_right(end0)});
        }
        std::sort(lone.begin(), lone.end(), [](const Lone& a, const Lone& b) { return a.w0 < b.w0; });
        for (size_t i = 0; i < lone.size();) {
            // One group: indels whose repeat spans overlap.
            size_t j = i + 1;
            hts_pos_t group_end = lone[i].w1;
            while (j < lone.size() && lone[j].w0 <= group_end) group_end = std::max(group_end, lone[j++].w1);
            std::vector<const Lone*> group;
            for (size_t k = i; k < j; ++k) group.push_back(&lone[k]);
            i = j;
            if (group.size() < 2) continue;
            std::sort(group.begin(), group.end(), [](const Lone* a, const Lone* b) { return a->reads > b->reads; });
            const Lone& a = *group[0];
            const Lone& b = *group[1];
            if (b.reads < kLocusMinMinorFraction * (a.reads + b.reads)) continue;
            const hts_pos_t lo = std::min(a.p0, b.p0);
            const hts_pos_t hi = std::max(a.p0 + a.ref_len, b.p0 + b.ref_len);
            const auto apply = [&](const Lone& e) {
                return ref.slice(lo, e.p0) + e.alt + ref.slice(e.p0 + e.ref_len, hi);
            };
            std::string seq_a = apply(a), seq_b = apply(b);
            for (char& ch : seq_a) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            for (char& ch : seq_b) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            if (seq_a == seq_b) continue;
            alleles.push_back(Allele{lo, static_cast<int>(hi - lo), std::move(seq_a), std::move(seq_b)});
        }
    }
    // The catalog's non-voting repeat indels, from their VCF-form alleles.
    if (graph_chunk.site_meta.size() == graph.candidates.size() &&
        graph_chunk.site_allele_orig_idx.size() == graph.candidates.size()) {
        for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
            const CandidateVariant& c = graph.candidates[ci];
            if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0) continue;
            if (c.counts.category != VariantCategory::RepeatHetIndel &&
                c.counts.category != VariantCategory::NoisyCandHet) continue;
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
            const std::vector<int>& orig = graph_chunk.site_allele_orig_idx[ci];
            if (meta.pos <= 0 || meta.ref.empty() || meta.non_selected_alt_class || orig.size() < 2) continue;
            const auto allele_seq = [&meta](int k) -> const std::string* {
                if (k == 0) return &meta.ref;
                return k >= 1 && static_cast<size_t>(k) <= meta.alts.size() ? &meta.alts[static_cast<size_t>(k) - 1]
                                                                            : nullptr;
            };
            const std::string* a0 = allele_seq(orig[0]);
            const std::string* a1 = allele_seq(orig[1]);
            if (a0 == nullptr || a1 == nullptr || !is_acgt(*a0) || !is_acgt(*a1) || *a0 == *a1) continue;
            const hts_pos_t p0 = meta.pos - 1;
            const int ref_len = static_cast<int>(meta.ref.size());
            std::string expected = meta.ref;
            for (char& ch : expected) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            if (ref.slice(p0, p0 + ref_len) != expected) continue;
            std::string s0 = *a0, s1 = *a1;
            for (char& ch : s0) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            for (char& ch : s1) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            alleles.push_back(Allele{p0, ref_len, std::move(s0), std::move(s1)});
        }
    }
    std::sort(alleles.begin(), alleles.end(), [](const Allele& a, const Allele& b) { return a.p0 < b.p0; });
    const SpanningReadIndex index(bam, graph);
    hts_pos_t last_end = -1;
    for (const Allele& al : alleles) {
        const hts_pos_t end0 = al.p0 + al.ref_len;
        const hts_pos_t w0 = std::max<hts_pos_t>(0, al.p0 - ref.repeat_left(al.p0) - kAlleleWindowFlank);
        const hts_pos_t w1 = end0 + ref.repeat_right(end0) + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength || w0 < last_end || overlaps_busy(w0, w1)) continue;
        const std::string left = ref.slice(w0, al.p0), right = ref.slice(end0, w1);
        const std::string hap_ref = left + al.seq0 + right;
        const std::string hap_alt = left + al.seq1 + right;
        if (hap_ref == hap_alt || hap_ref.find('N') != std::string::npos) continue;
        LocusWindowSite site;
        site.pos = w0 + 1;
        site.end = w1;
        int n0 = 0, n1 = 0;
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            const int d0 = locus_edit_distance(seq, hap_ref), d1 = locus_edit_distance(seq, hap_alt);
            if (std::abs(d0 - d1) < kLocusMinMargin) continue;
            const int side = d0 < d1 ? 0 : 1;
            (side == 0 ? n0 : n1)++;
            site.observations.emplace_back(
                ri, side, static_cast<float>(std::min(kLocusMaxObsError, 0.5 * std::exp(-std::abs(d0 - d1)))));
        }
        const int called = n0 + n1;
        if (std::min(n0, n1) < kLocusMinSide || std::min(n0, n1) < kLocusMinMinorFraction * called) continue;
        last_end = w1;
        out.push_back(std::move(site));
    }
    return out;
}

// Alleles competing with an indel at its locus: the reference and every other
// called allele whose edit lies inside the same window. At a two-allele
// heterozygote (1/2) the "other" haplotype carries the second allele, not the
// reference; contrasting REF with ALT alone calls its reads ALT.
struct LocusAllele { hts_pos_t p0; hts_pos_t end0; std::string alt; };

static std::vector<std::string> competing_haplotypes(const ReferenceView& ref, hts_pos_t w0, hts_pos_t w1,
                                                     const LocusAllele& self,
                                                     const std::vector<LocusAllele>& sorted_alleles) {
    std::vector<std::string> haps{ref.slice(w0, w1)};
    // Only alleles of the same locus: edits inside this allele's own repeat.
    // A different variant elsewhere in the window would otherwise stand in
    // for the reference, and reads carrying both would be miscalled.
    const hts_pos_t lo = self.p0 - ref.repeat_left(self.p0);
    const hts_pos_t hi = self.end0 + ref.repeat_right(self.end0);
    auto it = std::lower_bound(sorted_alleles.begin(), sorted_alleles.end(), lo,
                               [](const LocusAllele& a, hts_pos_t p) { return a.p0 < p; });
    for (; it != sorted_alleles.end() && it->p0 <= hi; ++it) {
        if (it->end0 > hi || it->p0 < w0 || it->end0 > w1) continue;
        if (it->p0 == self.p0 && it->end0 == self.end0 && it->alt == self.alt) continue;
        std::string h = ref.slice(w0, it->p0) + it->alt + ref.slice(it->end0, w1);
        if (std::find(haps.begin(), haps.end(), h) == haps.end()) haps.push_back(std::move(h));
    }
    return haps;
}

static int closest_allele_call(const std::string& seq, const std::string& hap_alt,
                               const std::vector<std::string>& competitors, int min_margin) {
    const int d1 = locus_edit_distance(seq, hap_alt);
    int d0 = -1;
    for (const std::string& h : competitors) {
        if (h == hap_alt) continue;
        const int d = locus_edit_distance(seq, h);
        if (d0 < 0 || d < d0) d0 = d;
    }
    if (d0 < 0 || std::abs(d0 - d1) < min_margin) return -1;
    return d1 < d0 ? 1 : 0;
}

// Index of the closest of several haplotype sequences, or -1 when the best
// does not beat the runner-up by min_margin edits.
static int closest_haplotype(const std::string& seq, const std::vector<std::string>& haps, int min_margin) {
    int best = -1, best_d = -1, second_d = -1;
    for (size_t i = 0; i < haps.size(); ++i) {
        const int d = locus_edit_distance(seq, haps[i]);
        if (best_d < 0 || d < best_d) {
            second_d = best_d;
            best_d = d;
            best = static_cast<int>(i);
        } else if (second_d < 0 || d < second_d) {
            second_d = d;
        }
    }
    if (second_d >= 0 && second_d - best_d < min_margin) return -1;
    return best;
}

// Haplotype sequences of a multi-allelic MSA insertion row over [w0, w1):
// index 0 the reference, index i the reference with msa_insertion_alts[i-1].
static std::vector<std::string> msa_allele_haplotypes(const ReferenceView& ref, const CandidateVariant& c,
                                                      hts_pos_t w0, hts_pos_t w1) {
    const hts_pos_t p0 = c.key.pos - 1;
    std::vector<std::string> haps{ref.slice(w0, w1)};
    for (const std::string& ins : c.msa_insertion_alts) {
        if (ins.find_first_not_of("ACGTacgt") != std::string::npos) return {};
        haps.push_back(ref.slice(w0, p0) + ins + ref.slice(p0, w1));
    }
    return haps;
}

// Widen a read profile's candidate range to [lo, hi], keeping every channel aligned.
static void widen_profile(ReadVariantProfile& p, int lo, int hi) {
    if (lo == p.start_var_idx && hi == p.end_var_idx) return;
    const size_t span = static_cast<size_t>(hi - lo + 1);
    if (p.start_var_idx < 0) {  // no calls yet: every channel starts empty
        p.alleles.clear();
        p.alt_qi.clear();
        p.graph_alleles.clear();
        p.bam_alleles.clear();
        p.bam_qi.clear();
        p.bam_base_qualities.clear();
    }
    const size_t shift = p.start_var_idx >= 0 ? static_cast<size_t>(p.start_var_idx - lo) : 0;
    const auto widen = [&](auto& v, auto fill, bool required) {
        if (v.empty() && !required) return;
        std::remove_reference_t<decltype(v)> w(span, fill);
        for (size_t k = 0; k < v.size(); ++k) w[shift + k] = v[k];
        v = std::move(w);
    };
    widen(p.alleles, -1, true);
    widen(p.alt_qi, 0, true);
    widen(p.graph_alleles, -1, false);
    widen(p.bam_alleles, -1, false);
    widen(p.bam_qi, 0, false);
    widen(p.bam_base_qualities, static_cast<uint8_t>(0), false);
    p.start_var_idx = lo;
    p.end_var_idx = hi;
}

size_t realign_indel_observations(PhasingChunk& chunk, const PhasingChunk& bam) {
    if (bam.ref_seq.empty()) return 0;
    const ReferenceView ref{bam};
    const SpanningReadIndex index(bam, chunk);
    const auto is_acgt = [](const std::string& s) {
        return s.find_first_not_of("ACGTacgt") == std::string::npos;
    };
    // Read -> its profile, for writing calls at a candidate index.
    std::vector<ReadVariantProfile*> profile_of(chunk.reads.size(), nullptr);
    for (ReadVariantProfile& p : chunk.read_var_profile)
        if (p.read_id >= 0 && static_cast<size_t>(p.read_id) < profile_of.size()) profile_of[p.read_id] = &p;
    std::map<size_t, std::vector<std::pair<size_t, int>>> calls;  // read -> (candidate, allele or -1)
    std::vector<LocusAllele> locus_alleles;
    for (const CandidateVariant& c : chunk.candidates) {
        if (!c.bam_injected || c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        locus_alleles.push_back(LocusAllele{c.key.pos - 1, c.key.pos - 1 + c.key.ref_len, c.key.alt});
    }
    std::sort(locus_alleles.begin(), locus_alleles.end(), [](const LocusAllele& a, const LocusAllele& b) { return a.p0 < b.p0; });
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& c = chunk.candidates[ci];
        if (!c.bam_injected || c.key.type == VariantType::Snp) continue;
        if (!c.msa_insertion_alts.empty()) {
            // One row, several insertion alleles: the closest allele by index.
            if (c.key.type != VariantType::Insertion || (c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
            const hts_pos_t p0 = c.key.pos - 1;
            const hts_pos_t w0 = std::max<hts_pos_t>(0, p0 - ref.repeat_left(p0) - kAlleleWindowFlank);
            const hts_pos_t w1 = p0 + ref.repeat_right(p0) + kAlleleWindowFlank;
            if (w1 - w0 > kLocusMaxLength) continue;
            const std::vector<std::string> haps = msa_allele_haplotypes(ref, c, w0, w1);
            if (haps.size() < 3) continue;
            for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads))
                if (profile_of[ri] != nullptr) calls[ri].emplace_back(ci, closest_haplotype(seq, haps, 1));
            continue;
        }
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        const hts_pos_t p0 = c.key.pos - 1;
        const hts_pos_t end0 = p0 + c.key.ref_len;
        const hts_pos_t w0 = std::max<hts_pos_t>(0, p0 - ref.repeat_left(p0) - kAlleleWindowFlank);
        const hts_pos_t w1 = end0 + ref.repeat_right(end0) + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength) continue;
        const std::string hap_ref = ref.slice(w0, w1);
        const std::string hap_alt = ref.slice(w0, p0) + c.key.alt + ref.slice(end0, w1);
        if (hap_ref == hap_alt || hap_ref.find('N') != std::string::npos) continue;
        // Allele 0 is the closest of REF and every other allele of the locus.
        const std::vector<std::string> competitors =
            competing_haplotypes(ref, w0, w1, LocusAllele{p0, end0, c.key.alt}, locus_alleles);
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            if (profile_of[ri] == nullptr) continue;
            calls[ri].emplace_back(ci, closest_allele_call(seq, hap_alt, competitors, 1));
        }
    }
    // Every read spanning the site's window is re-called there, its range
    // widened if its other calls stop short of the site.
    size_t changed = 0;
    for (const auto& [ri, list] : calls) {
        ReadVariantProfile* prof = profile_of[ri];
        int lo = prof->start_var_idx, hi = prof->end_var_idx;
        for (const auto& [ci, a] : list) {
            if (a < 0 && (lo < 0 || static_cast<int>(ci) < lo || static_cast<int>(ci) > hi)) continue;
            lo = lo < 0 ? static_cast<int>(ci) : std::min(lo, static_cast<int>(ci));
            hi = std::max(hi, static_cast<int>(ci));
        }
        if (lo < 0) continue;
        widen_profile(*prof, lo, hi);
        for (const auto& [ci, a] : list) {
            if (static_cast<int>(ci) < lo || static_cast<int>(ci) > hi) continue;
            int& slot = prof->alleles[ci - static_cast<size_t>(lo)];
            // The MSA already called this read here; realignment only fills reads
            // it did not call.
            if (kRealignKeepsMsaCalls && slot >= 0) continue;
            if (slot != a) { slot = a; ++changed; }
        }
    }
    if (changed > 0) rebuild_read_var_cr(chunk);
    return changed;
}

size_t fill_missing_observations(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                                 const std::vector<char>& only_reads) {
    PhasingChunk& chunk = graph_chunk.chunk;
    if (bam.ref_seq.empty()) return 0;
    const ReferenceView ref{bam};
    const SpanningReadIndex index(bam, chunk);
    const auto upper = [](std::string s) {
        for (char& ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        return s;
    };
    const auto is_acgt = [](const std::string& s) {
        return s.find_first_not_of("ACGTacgt") == std::string::npos;
    };
    const bool have_meta = graph_chunk.site_meta.size() == chunk.candidates.size() &&
                           graph_chunk.site_allele_orig_idx.size() == chunk.candidates.size();
    std::vector<ReadVariantProfile*> profile_of(chunk.reads.size(), nullptr);
    for (ReadVariantProfile& p : chunk.read_var_profile)
        if (p.read_id >= 0 && static_cast<size_t>(p.read_id) < profile_of.size()) profile_of[p.read_id] = &p;
    const auto has_call = [&](size_t ri, size_t ci) {
        const ReadVariantProfile* p = profile_of[ri];
        if (p == nullptr || p->start_var_idx < 0 || static_cast<int>(ci) < p->start_var_idx ||
            static_cast<int>(ci) > p->end_var_idx) return false;
        return p->alleles[ci - static_cast<size_t>(p->start_var_idx)] >= 0;
    };
    std::map<size_t, std::vector<std::pair<size_t, int>>> fills;  // read -> (candidate, allele)
    std::vector<LocusAllele> locus_alleles;
    for (const CandidateVariant& c : chunk.candidates) {
        if (!c.bam_injected || c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        locus_alleles.push_back(LocusAllele{c.key.pos - 1, c.key.pos - 1 + c.key.ref_len, c.key.alt});
    }
    std::sort(locus_alleles.begin(), locus_alleles.end(), [](const LocusAllele& a, const LocusAllele& b) { return a.p0 < b.p0; });
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& c = chunk.candidates[ci];
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        if (!c.msa_insertion_alts.empty()) {
            if (!c.bam_injected || c.key.type != VariantType::Insertion) continue;
            const hts_pos_t p0 = c.key.pos - 1;
            const hts_pos_t w0 = std::max<hts_pos_t>(0, p0 - ref.repeat_left(p0) - kAlleleWindowFlank);
            const hts_pos_t w1 = p0 + ref.repeat_right(p0) + kAlleleWindowFlank;
            if (w1 - w0 > kLocusMaxLength) continue;
            const std::vector<std::string> haps = msa_allele_haplotypes(ref, c, w0, w1);
            if (haps.size() < 3) continue;
            for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
                if (ri >= only_reads.size() || !only_reads[ri] || chunk.reads[ri].is_skipped || has_call(ri, ci)) continue;
                const int a = closest_haplotype(seq, haps, kLocusMinMargin);
                if (a >= 0) fills[ri].emplace_back(ci, a);
            }
            continue;
        }
        // The site's reference span and its two candidate alleles as sequence.
        hts_pos_t p0 = 0;
        int ref_len = 0;
        std::string seq0, seq1;
        std::vector<std::string> decomposed_others;
        if (c.bam_injected) {
            if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
            p0 = c.key.pos - 1;
            ref_len = c.key.ref_len;
            seq0 = ref.slice(p0, p0 + ref_len);
            seq1 = upper(c.key.alt);
        } else {
            if (!have_meta) continue;
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
            const std::vector<int>& orig = graph_chunk.site_allele_orig_idx[ci];
            if (meta.pos <= 0 || meta.ref.empty() || orig.size() < 2) continue;
            const auto allele_seq = [&meta](int k) -> const std::string* {
                if (k == 0) return &meta.ref;
                return k >= 1 && static_cast<size_t>(k) <= meta.alts.size() ? &meta.alts[static_cast<size_t>(k) - 1]
                                                                            : nullptr;
            };
            const std::string* a0 = allele_seq(orig[0]);
            const std::string* a1 = allele_seq(orig[1]);
            if (a0 == nullptr || a1 == nullptr || !is_acgt(*a0) || !is_acgt(*a1)) continue;
            p0 = meta.pos - 1;
            ref_len = static_cast<int>(meta.ref.size());
            if (ref.slice(p0, p0 + ref_len) != upper(meta.ref)) continue;
            seq0 = upper(*a0);
            seq1 = upper(*a1);
            // Allele 0 of a decomposed multi-allelic row is "any other allele of
            // the site": every one of them competes with the selected allele.
            if (meta.non_selected_alt_class) {
                decomposed_others.push_back(upper(meta.ref));
                bool usable = is_acgt(meta.ref);
                for (size_t a = 0; a < meta.alts.size(); ++a) {
                    if (static_cast<int>(a) + 1 == orig[1]) continue;
                    if (!is_acgt(meta.alts[a])) { usable = false; break; }
                    decomposed_others.push_back(upper(meta.alts[a]));
                }
                if (!usable) continue;
            }
        }
        if (seq0 == seq1) continue;
        const hts_pos_t end0 = p0 + ref_len;
        const bool snp = seq0.size() == 1 && seq1.size() == 1 && ref_len == 1;
        const hts_pos_t w0 = std::max<hts_pos_t>(0, p0 - (snp ? 0 : ref.repeat_left(p0)) - kAlleleWindowFlank);
        const hts_pos_t w1 = end0 + (snp ? 0 : ref.repeat_right(end0)) + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength) continue;
        const std::string left = ref.slice(w0, p0), right = ref.slice(end0, w1);
        const std::string hap0 = left + seq0 + right, hap1 = left + seq1 + right;
        if (hap0.find('N') != std::string::npos) continue;
        std::vector<std::string> competitors{hap0};
        if (c.bam_injected && !snp)
            competitors = competing_haplotypes(ref, w0, w1, LocusAllele{p0, end0, seq1}, locus_alleles);
        for (const std::string& o : decomposed_others) competitors.push_back(left + o + right);
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            if (ri >= only_reads.size() || !only_reads[ri] || chunk.reads[ri].is_skipped || has_call(ri, ci)) continue;
            // A one-edit difference decides a SNP; at an indel it is as often a
            // homopolymer length error as the allele.
            const int a = closest_allele_call(seq, hap1, competitors, snp ? 1 : kLocusMinMargin);
            if (a >= 0) fills[ri].emplace_back(ci, a);
        }
    }
    // Write the calls, widening a profile's index range where needed.
    size_t added = 0;
    for (auto& [ri, calls] : fills) {
        ReadVariantProfile* p = profile_of[ri];
        if (p == nullptr) continue;
        int lo = p->start_var_idx, hi = p->end_var_idx;
        for (const auto& [ci, a] : calls) {
            lo = lo < 0 ? static_cast<int>(ci) : std::min(lo, static_cast<int>(ci));
            hi = std::max(hi, static_cast<int>(ci));
        }
        widen_profile(*p, lo, hi);
        for (const auto& [ci, a] : calls) {
            p->alleles[ci - static_cast<size_t>(lo)] = a;
            ++added;
        }
    }
    if (added > 0) rebuild_read_var_cr(chunk);
    return added;
}

std::vector<char> add_alignment_only_reads(PhasingChunk& chunk, const PhasingChunk& bam) {
    std::unordered_set<std::string_view> present;
    present.reserve(chunk.reads.size());
    for (const ReadRecord& r : chunk.reads) present.insert(r.qname);
    size_t added = 0;
    std::vector<char> is_added(chunk.reads.size(), 0);
    for (const ReadRecord& source : bam.reads) {
        if (!source.alignment || source.is_skipped) continue;
        if ((source.alignment->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FUNMAP)) != 0)
            continue;
        if (source.alignment->core.qual < kAddedReadMinMapq || present.count(source.qname) != 0) continue;
        ReadRecord read;
        read.tid = chunk.region.tid;
        read.qname = source.qname;
        read.mapq = source.mapq;
        read.beg = source.beg;
        read.end = source.end;
        read.reverse = source.reverse;
        ReadVariantProfile prof;
        prof.read_id = static_cast<int>(chunk.reads.size());
        chunk.reads.push_back(std::move(read));
        chunk.read_var_profile.push_back(std::move(prof));
        present.insert(chunk.reads.back().qname);
        is_added.push_back(1);
        ++added;
    }
    chunk.haps.resize(chunk.reads.size(), 0);
    chunk.phase_sets.resize(chunk.reads.size(), kUnphasedReadPhaseSet);
    if (added > 0) rebuild_read_var_cr(chunk);
    return is_added;
}

// One alignment event of a read against the reference: a mismatching base, an
// insertion or a deletion, with indels left-normalized so placements of one
// homopolymer or repeat length change compare equal.
struct ReadEvent {
    hts_pos_t pos;    // 0-based: the base (SNP), the first deleted base, or the base an insertion precedes
    std::string what; // "S" + base, "I" + inserted sequence, "D" + deleted length
    hts_pos_t end;    // 0-based exclusive reference end the event touches
    bool operator<(const ReadEvent& o) const { return pos != o.pos ? pos < o.pos : what < o.what; }
    bool operator==(const ReadEvent& o) const { return pos == o.pos && what == o.what; }
};

static std::vector<ReadEvent> read_events(const bam1_t* b, const ReferenceView& ref, hts_pos_t lo, hts_pos_t hi) {
    std::vector<ReadEvent> out;
    const uint32_t* cigar = bam_get_cigar(b);
    const uint8_t* seq = bam_get_seq(b);
    const uint8_t* qual = bam_get_qual(b);
    hts_pos_t r = b->core.pos, q = 0;
    for (uint32_t i = 0; i < b->core.n_cigar && r < hi; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const hts_pos_t len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (hts_pos_t k = 0; k < len; ++k) {
                const hts_pos_t p = r + k;
                if (p < lo || p >= hi) continue;
                const char base = seq_nt16_str[bam_seqi(seq, q + k)];
                if (base == 'N' || qual[q + k] < kEventMinBaseQuality) continue;
                const char rb = ref.base(p);
                if (rb != 'N' && base != rb) out.push_back(ReadEvent{p, std::string("S") + base, p + 1});
            }
            r += len;
            q += len;
        } else if (op == BAM_CINS) {
            if (r >= lo && r < hi) {
                std::string ins(static_cast<size_t>(len), 'N');
                for (hts_pos_t k = 0; k < len; ++k) ins[static_cast<size_t>(k)] = seq_nt16_str[bam_seqi(seq, q + k)];
                hts_pos_t p = r;
                while (p > 0 && ref.base(p - 1) != 'N' && ref.base(p - 1) == ins.back()) {
                    ins = std::string(1, ref.base(p - 1)) + ins.substr(0, ins.size() - 1);
                    --p;
                }
                out.push_back(ReadEvent{p, "I" + ins, p});
            }
            q += len;
        } else if (op == BAM_CDEL) {
            if (r >= lo && r < hi) {
                hts_pos_t p = r;
                while (p > 0 && ref.base(p - 1) != 'N' && ref.base(p - 1) == ref.base(p + len - 1)) --p;
                out.push_back(ReadEvent{p, "D" + std::to_string(len), p + len});
            }
            r += len;
        } else if (op == BAM_CREF_SKIP) {
            r += len;
        } else if (op == BAM_CSOFT_CLIP) {
            q += len;
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

size_t label_reads_from_haplotype_consensus(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam) {
    PhasingChunk& chunk = graph_chunk.chunk;
    if (bam.ref_seq.empty()) return 0;
    const ReferenceView ref{bam};
    const SpanningReadIndex index(bam, chunk);
    struct Span { size_t i; hts_pos_t beg; hts_pos_t end; };  // index entry, aligned span
    std::unordered_set<size_t> decided_before;
    for (const auto& d : graph_chunk.deferred_read_labels) decided_before.insert(std::get<0>(d));
    std::vector<Span> unlabelled, labelled;
    for (size_t i = 0; i < index.reads.size(); ++i) {
        const size_t ri = index.graph_index[i];
        if (ri >= chunk.haps.size() || chunk.reads[ri].is_skipped) continue;
        if (chunk.haps[ri] == 0 && decided_before.count(ri) != 0) continue;
        const bam1_t* aln = index.reads[i]->alignment.get();
        const Span sp{i, aln->core.pos, bam_endpos(aln)};
        if (chunk.haps[ri] == 0) unlabelled.push_back(sp);
        else if (chunk.phase_sets[ri] > 0) labelled.push_back(sp);
    }
    if (unlabelled.empty() || labelled.empty()) return 0;
    // Only the stretches unlabelled reads cover matter.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> regions;
    for (const Span& u : unlabelled) regions.emplace_back(u.beg, u.end);
    std::sort(regions.begin(), regions.end());
    std::vector<std::pair<hts_pos_t, hts_pos_t>> merged;
    for (const auto& r : regions) {
        if (!merged.empty() && r.first <= merged.back().second) merged.back().second = std::max(merged.back().second, r.second);
        else merged.push_back(r);
    }
    // Events of the labelled reads inside those stretches, by (phase set, haplotype).
    std::map<std::pair<hts_pos_t, ReadEvent>, std::array<int, 2>> carriers;  // (ps, event) -> reads per haplotype
    std::map<hts_pos_t, std::vector<std::pair<Span, int>>> by_set;           // ps -> (span, hap index)
    for (const Span& l : labelled) {
        const size_t ri = index.graph_index[l.i];
        const hts_pos_t ps = chunk.phase_sets[ri];
        const int h = chunk.haps[ri] == 1 ? 0 : 1;
        auto it = std::upper_bound(merged.begin(), merged.end(), std::make_pair(l.end, hts_pos_t{0}));
        bool touches = false;
        for (auto m = merged.begin(); m != it; ++m)
            if (m->second > l.beg) { touches = true; break; }
        if (!touches) continue;
        by_set[ps].emplace_back(l, h);
        for (const auto& [a, b] : merged) {
            if (b <= l.beg || a >= l.end) continue;
            for (const ReadEvent& e : read_events(index.reads[l.i]->alignment.get(), ref, std::max(a, l.beg), std::min(b, l.end)))
                ++carriers[{ps, e}][h];
        }
    }
    // Discriminating events: carried by one haplotype's covering reads and
    // absent from the other's.
    struct Marker { hts_pos_t ps; ReadEvent e; int carrier; };
    std::vector<Marker> markers;
    for (const auto& [key, cnt] : carriers) {
        if (std::max(cnt[0], cnt[1]) < kConsensusMinReads) continue;
        const auto& [ps, e] = key;
        std::array<int, 2> cover{0, 0};
        for (const auto& [sp, h] : by_set[ps])
            if (sp.beg <= e.pos - kConsensusFlank && sp.end >= e.end + kConsensusFlank) ++cover[h];
        if (cover[0] < kConsensusMinReads || cover[1] < kConsensusMinReads) continue;
        const double f0 = static_cast<double>(cnt[0]) / cover[0], f1 = static_cast<double>(cnt[1]) / cover[1];
        if (f0 >= kConsensusCarrierFraction && f1 <= 1.0 - kConsensusCarrierFraction) markers.push_back(Marker{ps, e, 0});
        else if (f1 >= kConsensusCarrierFraction && f0 <= 1.0 - kConsensusCarrierFraction) markers.push_back(Marker{ps, e, 1});
    }
    if (markers.empty()) return 0;
    std::sort(markers.begin(), markers.end(), [](const Marker& a, const Marker& b) { return a.e.pos < b.e.pos; });
    std::map<hts_pos_t, size_t> anchor_of;
    for (size_t ri = 0; ri < chunk.haps.size(); ++ri)
        if (chunk.haps[ri] != 0 && chunk.phase_sets[ri] > 0) anchor_of.try_emplace(chunk.phase_sets[ri], ri);
    size_t decided = 0;
    for (const Span& u : unlabelled) {
        const auto first = std::lower_bound(markers.begin(), markers.end(), u.beg + kConsensusFlank,
            [](const Marker& m, hts_pos_t p) { return m.e.pos < p; });
        std::map<hts_pos_t, std::array<int, 2>> votes;  // ps -> votes for haplotype index
        std::vector<ReadEvent> own;
        bool have_own = false;
        for (auto m = first; m != markers.end() && m->e.pos < u.end - kConsensusFlank; ++m) {
            if (m->e.end + kConsensusFlank > u.end) continue;
            if (!have_own) {
                own = read_events(index.reads[u.i]->alignment.get(), ref, u.beg, u.end);
                have_own = true;
            }
            const bool carries = std::binary_search(own.begin(), own.end(), m->e);
            ++votes[m->ps][carries ? m->carrier : 1 - m->carrier];
        }
        if (votes.empty()) continue;
        const auto best = std::max_element(votes.begin(), votes.end(), [](const auto& a, const auto& b) {
            return a.second[0] + a.second[1] < b.second[0] + b.second[1];
        });
        const int n = best->second[0] + best->second[1];
        const int lead = std::max(best->second[0], best->second[1]);
        if (lead < kConsensusMinAgreement * n) continue;
        const int hap = best->second[0] >= best->second[1] ? 1 : 2;
        const auto anchor = anchor_of.find(best->first);
        if (anchor == anchor_of.end()) continue;
        const size_t ri = index.graph_index[u.i];
        graph_chunk.deferred_read_labels.emplace_back(ri, anchor->second, chunk.haps[anchor->second] == hap);
        ++decided;
    }
    return decided;
}

size_t label_reads_from_verified_indels(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam) {
    PhasingChunk& chunk = graph_chunk.chunk;
    if (bam.ref_seq.empty()) return 0;
    const ReferenceView ref{bam};
    const SpanningReadIndex index(bam, chunk);
    const auto is_acgt = [](const std::string& s) {
        return s.find_first_not_of("ACGTacgt") == std::string::npos;
    };
    // (phase set, haplotype) votes per unlabelled read; a read with conflicting
    // votes stays unlabelled.
    std::map<size_t, std::map<std::pair<hts_pos_t, int>, int>> votes;
    std::vector<LocusAllele> locus_alleles;
    for (const CandidateVariant& c : bam.candidates) {
        if (c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        if (c.counts.alt_cov < kLocusMinSide) continue;
        locus_alleles.push_back(LocusAllele{c.key.pos - 1, c.key.pos - 1 + c.key.ref_len, c.key.alt});
    }
    std::sort(locus_alleles.begin(), locus_alleles.end(), [](const LocusAllele& a, const LocusAllele& b) { return a.p0 < b.p0; });
    hts_pos_t last_end = -1;
    for (const CandidateVariant& c : bam.candidates) {
        if (c.counts.category != VariantCategory::NoisyCandHet || !c.msa_verified) continue;
        if (c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        const hts_pos_t p0 = c.key.pos - 1;
        const hts_pos_t end0 = p0 + c.key.ref_len;
        const hts_pos_t w0 = std::max<hts_pos_t>(0, p0 - ref.repeat_left(p0) - kAlleleWindowFlank);
        const hts_pos_t w1 = end0 + ref.repeat_right(end0) + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength || w0 < last_end) continue;  // one call per repeat
        const std::string hap_ref = ref.slice(w0, w1);
        const std::string hap_alt = ref.slice(w0, p0) + c.key.alt + ref.slice(end0, w1);
        if (hap_ref == hap_alt || hap_ref.find('N') != std::string::npos) continue;
        std::vector<std::pair<size_t, int>> calls;  // read, allele
        const std::vector<std::string> competitors =
            competing_haplotypes(ref, w0, w1, LocusAllele{p0, end0, c.key.alt}, locus_alleles);
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            const int a = closest_allele_call(seq, hap_alt, competitors, 1);
            if (a >= 0) calls.emplace_back(ri, a);
        }
        // The site's orientation from the reads already labelled: which allele
        // haplotype 1 carries, in the phase set most of them share.
        std::map<hts_pos_t, std::array<int, 2>> by_set;  // (hap1 carries ALT, hap1 carries REF)
        for (const auto& [ri, a] : calls) {
            if (ri >= chunk.haps.size() || chunk.haps[ri] == 0 || chunk.phase_sets[ri] <= 0) continue;
            ++by_set[chunk.phase_sets[ri]][(chunk.haps[ri] == 1) == (a == 1) ? 0 : 1];
        }
        if (by_set.empty()) continue;
        const auto best = std::max_element(by_set.begin(), by_set.end(), [](const auto& x, const auto& y) {
            return x.second[0] + x.second[1] < y.second[0] + y.second[1];
        });
        const int n = best->second[0] + best->second[1];
        const int agree = std::max(best->second[0], best->second[1]);
        if (n < kLastResortMinLabelled || agree < kLastResortMinAgreement * n) continue;
        last_end = w1;
        const bool hap1_alt = best->second[0] >= best->second[1];
        for (const auto& [ri, a] : calls) {
            if (ri >= chunk.haps.size() || chunk.haps[ri] != 0) continue;
            const int hap = (a == 1) == hap1_alt ? 1 : 2;
            ++votes[ri][{best->first, hap}];
        }
    }
    // Decisions are tied to an anchor read of the same phase set and applied
    // after chunk stitching, so a last-resort label never votes in it and
    // follows any flip the stitch applies to its block.
    std::map<hts_pos_t, size_t> anchor_of;
    for (size_t ri = 0; ri < chunk.haps.size(); ++ri)
        if (chunk.haps[ri] != 0 && chunk.phase_sets[ri] > 0) anchor_of.try_emplace(chunk.phase_sets[ri], ri);
    size_t labelled = 0;
    for (const auto& [ri, v] : votes) {
        if (v.size() != 1) continue;  // the read's verified indels disagree
        const auto& [ps, hap] = v.begin()->first;
        const auto anchor = anchor_of.find(ps);
        if (anchor == anchor_of.end()) continue;
        graph_chunk.deferred_read_labels.emplace_back(ri, anchor->second, chunk.haps[anchor->second] == hap);
        ++labelled;
    }
    return labelled;
}

size_t apply_deferred_read_labels(GraphChunkBuildResult& graph_chunk) {
    PhasingChunk& chunk = graph_chunk.chunk;
    size_t applied = 0;
    for (const auto& [ri, anchor, same] : graph_chunk.deferred_read_labels) {
        if (ri >= chunk.haps.size() || anchor >= chunk.haps.size() || chunk.haps[ri] != 0) continue;
        const int anchor_hap = chunk.haps[anchor];
        if ((anchor_hap != 1 && anchor_hap != 2) || chunk.phase_sets[anchor] <= 0) continue;
        chunk.haps[ri] = same ? anchor_hap : 3 - anchor_hap;
        chunk.phase_sets[ri] = chunk.phase_sets[anchor];
        ++applied;
    }
    graph_chunk.deferred_read_labels.clear();
    return applied;
}

// A two-allele (1/2) heterozygote at a repeat: the alignment solve reports it
// as two REF-versus-ALT rows, each with almost no REF reads. Judged one row at
// a time, a read carrying the other allele looks like ALT to both rows, and its
// calls cancel. The two rows become one EM site instead: every spanning read is
// realigned to the reference carrying allele A and carrying allele B.
struct TwoAlleleLocus { size_t row_a; size_t row_b; size_t window; };

static std::vector<TwoAlleleLocus> two_allele_loci(const PhasingChunk& chunk, const PhasingChunk& bam,
                                                   std::vector<LocusWindowSite>& windows) {
    std::vector<TwoAlleleLocus> out;
    if (bam.ref_seq.empty()) return out;
    const ReferenceView ref{bam};
    const auto is_acgt = [](const std::string& s) {
        return s.find_first_not_of("ACGTacgt") == std::string::npos;
    };
    struct Row { size_t ci; hts_pos_t p0; hts_pos_t end0; std::string alt; int reads; int depth; hts_pos_t lo; hts_pos_t hi; };
    std::vector<Row> rows;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& c = chunk.candidates[ci];
        if (!c.bam_injected || (c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        if (c.key.type == VariantType::Snp || !c.msa_insertion_alts.empty()) continue;
        if (!is_acgt(c.key.alt) || (c.key.ref_len == 0 && c.key.alt.empty())) continue;
        // A row's REF count holds every read without its ALT, the other allele's
        // carriers included, so REF-lessness shows only for the pair together.
        const int alt_reads = c.counts.alt_cov;
        if (alt_reads < kLocusMinSide) continue;
        const hts_pos_t p0 = c.key.pos - 1, end0 = p0 + c.key.ref_len;
        rows.push_back(Row{ci, p0, end0, c.key.alt, alt_reads, alt_reads + c.counts.ref_cov,
                           p0 - ref.repeat_left(p0), end0 + ref.repeat_right(end0)});
    }
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.lo < b.lo; });
    const SpanningReadIndex index(bam, chunk);
    for (size_t i = 0; i < rows.size();) {
        size_t j = i + 1;
        hts_pos_t group_hi = rows[i].hi;
        while (j < rows.size() && rows[j].lo <= group_hi) group_hi = std::max(group_hi, rows[j++].hi);
        std::vector<const Row*> group;
        for (size_t k = i; k < j; ++k) group.push_back(&rows[k]);
        i = j;
        if (group.size() < 2) continue;
        std::sort(group.begin(), group.end(), [](const Row* a, const Row* b) { return a->reads > b->reads; });
        const Row& a = *group[0];
        const Row& b = *group[1];
        // Together the two alleles must hold the locus: little room for REF.
        if (a.reads + b.reads < (1.0 - kPairMaxRefFraction) * std::max(a.depth, b.depth)) continue;
        const hts_pos_t lo = std::min(a.p0, b.p0), hi = std::max(a.end0, b.end0);
        const hts_pos_t w0 = std::max<hts_pos_t>(0, lo - ref.repeat_left(lo) - kAlleleWindowFlank);
        const hts_pos_t w1 = hi + ref.repeat_right(hi) + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength || a.p0 < w0 || b.p0 < w0 || a.end0 > w1 || b.end0 > w1) continue;
        const std::string hap_a = ref.slice(w0, a.p0) + a.alt + ref.slice(a.end0, w1);
        const std::string hap_b = ref.slice(w0, b.p0) + b.alt + ref.slice(b.end0, w1);
        if (hap_a == hap_b || hap_a.find('N') != std::string::npos) continue;
        LocusWindowSite site;
        site.pos = w0 + 1;
        site.end = w1;
        int n0 = 0, n1 = 0;
        for (const auto& [ri, seq] : index.segments(w0, w1, kAlleleWindowMaxReads)) {
            const int da = locus_edit_distance(seq, hap_a), db = locus_edit_distance(seq, hap_b);
            if (da == db) continue;
            const int side = da < db ? 0 : 1;
            (side == 0 ? n0 : n1)++;
            site.observations.emplace_back(
                ri, side, static_cast<float>(std::min(kLocusMaxObsError, 0.5 * std::exp(-std::abs(da - db)))));
        }
        const int called = n0 + n1;
        if (std::min(n0, n1) < kLocusMinSide || std::min(n0, n1) < kLocusMinMinorFraction * called) continue;
        out.push_back(TwoAlleleLocus{a.ci, b.ci, windows.size()});
        windows.push_back(std::move(site));
    }
    return out;
}

std::vector<char> bridge_weak_sites(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam) {
    const PhasingChunk& chunk = graph_chunk.chunk;
    std::vector<char> weak(chunk.candidates.size(), 0);
    // Reads covering each site, and reads calling it.
    std::vector<hts_pos_t> begs, ends;
    for (const ReadRecord& r : chunk.reads) {
        if (r.is_skipped || r.end <= r.beg) continue;
        begs.push_back(r.beg);
        ends.push_back(r.end);
    }
    std::sort(begs.begin(), begs.end());
    std::sort(ends.begin(), ends.end());
    std::vector<int> called(chunk.candidates.size(), 0);
    for (size_t ri = 0; ri < chunk.read_var_profile.size() && ri < chunk.reads.size(); ++ri) {
        const ReadVariantProfile& p = chunk.read_var_profile[ri];
        if (p.start_var_idx < 0 || chunk.reads[ri].is_skipped) continue;
        for (size_t k = 0; k < p.alleles.size(); ++k)
            if (p.alleles[k] >= 0) ++called[static_cast<size_t>(p.start_var_idx) + k];
    }
    const bool have_ref = !bam.ref_seq.empty();
    const ReferenceView ref{bam};
    const bool have_meta = graph_chunk.site_meta.size() == chunk.candidates.size() &&
                           graph_chunk.site_allele_orig_idx.size() == chunk.candidates.size();
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& c = chunk.candidates[ci];
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        const hts_pos_t pos = c.key.pos;
        const long covering = static_cast<long>(std::upper_bound(begs.begin(), begs.end(), pos) - begs.begin()) -
                              static_cast<long>(std::upper_bound(ends.begin(), ends.end(), pos) - ends.begin());
        if (kBridgeDepthRule && covering > 0 && called[ci] < kBridgeMinCalledFraction * covering) {
            weak[ci] = 1;
            continue;
        }
        if (!kBridgeHomopolymerRule || c.key.type == VariantType::Snp) continue;
        // One homopolymer base between the two alleles.
        if (!c.msa_insertion_alts.empty()) {
            for (size_t a = 0; a < c.msa_insertion_alts.size(); ++a)
                for (size_t b = a + 1; b < c.msa_insertion_alts.size(); ++b)
                    if (std::abs(static_cast<long>(c.msa_insertion_alts[a].size()) -
                                 static_cast<long>(c.msa_insertion_alts[b].size())) == 1) weak[ci] = 1;
            continue;
        }
        if (c.bam_injected && have_ref) {
            const hts_pos_t p0 = c.key.pos - 1;
            if (c.key.type == VariantType::Insertion && c.key.ref_len == 0 && c.key.alt.size() == 1) {
                const char b = static_cast<char>(std::toupper(static_cast<unsigned char>(c.key.alt[0])));
                if (ref.base(p0 - 1) == b || ref.base(p0) == b) weak[ci] = 1;
            } else if (c.key.type == VariantType::Deletion && c.key.ref_len == 1) {
                const char b = ref.base(p0);
                if (ref.base(p0 - 1) == b || ref.base(p0 + 1) == b) weak[ci] = 1;
            }
            continue;
        }
        if (have_meta && !c.bam_injected) {
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
            const std::vector<int>& orig = graph_chunk.site_allele_orig_idx[ci];
            if (orig.size() < 2 || orig[1] < 1 || static_cast<size_t>(orig[1]) > meta.alts.size()) continue;
            const long sel = static_cast<long>(meta.alts[static_cast<size_t>(orig[1]) - 1].size());
            if (std::abs(sel - static_cast<long>(meta.ref.size())) == 1 && !meta.non_selected_alt_class) weak[ci] = 1;
            if (meta.non_selected_alt_class)
                for (size_t a = 0; a < meta.alts.size(); ++a)
                    if (static_cast<int>(a) + 1 != orig[1] && std::abs(sel - static_cast<long>(meta.alts[a].size())) == 1)
                        weak[ci] = 1;
        }
    }
    return weak;
}

void phase_chunk_with_alignment_sites(GraphChunkBuildResult& graph_chunk, const PhasingChunk* bam) {
    PhasingChunk& chunk = graph_chunk.chunk;
    // Reads the graph chunk dropped for lack of an informative catalog allele
    // can still carry a variant only the windows see (a two-allele indel).
    if (bam != nullptr) {
        const std::vector<char> added = add_alignment_only_reads(chunk, *bam);
        // They carry no catalog calls; give them calls by realignment.
        fill_missing_observations(graph_chunk, *bam, added);
    }
    // Pileup calls at repeat indels are unreliable; the alignment's indel sites
    // are re-called per read by realignment to their two exact alleles.
    if (bam != nullptr) realign_indel_observations(chunk, *bam);
    // Two-allele loci enter the EM as one site; their two rows sit out the
    // solve and take their phase from it afterwards.
    std::vector<LocusWindowSite> pair_windows;
    std::vector<TwoAlleleLocus> pairs;
    std::vector<std::pair<size_t, uint32_t>> parked;  // row, original category bits
    if (bam != nullptr) {
        if (kMergeTwoAlleleRows) pairs = two_allele_loci(chunk, *bam, pair_windows);
        for (const TwoAlleleLocus& p : pairs)
            for (const size_t row : {p.row_a, p.row_b}) {
                parked.emplace_back(row, chunk.candidates[row].lcd_var_i_to_cate);
                chunk.candidates[row].lcd_var_i_to_cate &= ~kCandGermlineVarCate;
            }
    }
    // The first solve labels reads where the clean stage could not (dense
    // noisy regions); those labels seed the windows the second solve adds.
    const std::vector<char> weak = bam != nullptr ? bridge_weak_sites(graph_chunk, *bam) : std::vector<char>{};
    const std::vector<char>* weak_ptr = bam != nullptr ? &weak : nullptr;
    phase_chunk_by_global_em(chunk, pair_windows.empty() ? nullptr : &pair_windows, weak_ptr);
    std::vector<LocusWindowSite> loci = pair_windows;
    if (bam != nullptr) {
        std::vector<LocusWindowSite> haplotype_windows = build_locus_window_sites(*bam, chunk);
        // A haplotype window over a two-allele locus would count its reads twice.
        for (LocusWindowSite& w : haplotype_windows) {
            bool overlaps = false;
            for (const LocusWindowSite& p : pair_windows)
                if (w.pos <= p.end && p.pos <= w.end) { overlaps = true; break; }
            if (!overlaps) loci.push_back(std::move(w));
        }
        std::vector<LocusWindowSite> alleles = build_allele_window_sites(*bam, graph_chunk, loci);
        loci.insert(loci.end(), std::make_move_iterator(alleles.begin()), std::make_move_iterator(alleles.end()));
    }
    phase_chunk_by_global_em(chunk, &loci, weak_ptr);
    for (const auto& [row, bits] : parked) chunk.candidates[row].lcd_var_i_to_cate = bits;
    for (size_t i = 0; i < pairs.size(); ++i) {
        const LocusWindowSite& w = loci[i];  // pair windows lead the list
        CandidateVariant& a = chunk.candidates[pairs[i].row_a];
        CandidateVariant& b = chunk.candidates[pairs[i].row_b];
        if (w.phase < 0 || w.error >= kEmUnreliableError) {
            for (CandidateVariant* c : {&a, &b}) {
                c->phase_set = kUnsetCandidatePhaseSet;
                c->hap_to_cons_alle[1] = c->hap_to_cons_alle[2] = -1;
                c->hap_alt = c->hap_ref = 0;
            }
            continue;
        }
        // Side 0 is allele A; phase 0 puts side 0 on haplotype 1.
        const bool a_on_hap1 = w.phase == 0;
        set_site_phase(a, w.block, a_on_hap1 ? 1 : 0, a_on_hap1 ? 0 : 1);
        set_site_phase(b, w.block, a_on_hap1 ? 0 : 1, a_on_hap1 ? 1 : 0);
    }
    // Last resort: reads still unlabelled take their haplotype from MSA-verified
    // indels whose phase the labelled reads already establish.
    if (bam != nullptr) {
        // Verified indels first (exact alleles), then the labelled reads' own
        // haplotype markers for reads still undecided.
        label_reads_from_verified_indels(graph_chunk, *bam);
        label_reads_from_haplotype_consensus(graph_chunk, *bam);
    }
}

} // namespace pgphase_collect
