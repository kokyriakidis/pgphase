/**
 * @file graph_post.cpp
 * @brief Post-construction and post-solve graph phasing passes: physical SNP
 *        validation, branch supplementation, joint candidate phasing, chunk
 *        overlap linking, unphased read rescue, and read/BAM/TSV output helpers.
 */

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

} // namespace

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
        const bool low_fraction = low_fraction_sites.count(out.chunk.candidates[i].graph_site_key) != 0 &&
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
        if (alt_i >= site.allele_walks.size() || site.allele_walks.size() <= 2) continue;
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
        graph_chunk.chunk.candidates[candidate_index].counts.n_uniq_alles != 2 ||
        candidate_index >= graph_chunk.site_allele_orig_idx.size() ||
        graph_chunk.site_allele_orig_idx[candidate_index].size() != 2) {
        return nullptr;
    }
    const std::optional<RowAllele> allele = row_allele(graph_chunk, candidate_index);
    if (!allele) return nullptr;
    return allele->alt->empty() || *allele->alt == "*" ? nullptr : allele->alt;
}

ChunkReadLabel read_label(const PhasingChunk& chunk, size_t ri) {
    const int hap = ri < chunk.haps.size() ? chunk.haps[ri] : 0;
    const hts_pos_t ps = ri < chunk.phase_sets.size() ? chunk.phase_sets[ri] : kUnphasedReadPhaseSet;
    if ((hap == 1 || hap == 2) && ps > 0) return ChunkReadLabel{hap, ps, true};
    const int gap = ri < chunk.gap_haps.size() ? chunk.gap_haps[ri] : 0;
    const hts_pos_t gap_ps = ri < chunk.gap_phase_sets.size() ? chunk.gap_phase_sets[ri] : 0;
    if ((gap == 1 || gap == 2) && gap_ps > kGapFillPsOffset)
        return ChunkReadLabel{gap, gap_ps - kGapFillPsOffset, false};
    return ChunkReadLabel{};
}

std::unordered_map<size_t, UnlinkedHetPlacement> place_unlinked_hets(
        const PhasingChunk& chunk, std::unordered_map<size_t, std::string>* trace_votes) {
    std::unordered_map<size_t, UnlinkedHetPlacement> placements;
    const auto phased_het = [](const CandidateVariant& c) {
        const auto& h = c.hap_to_cons_alle;
        return (c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0 && c.phase_set > 0 &&
               h.size() > 2 && h[1] >= 0 && h[2] >= 0 && h[1] != h[2];
    };
    // Hets per phase set by position: a catalog site and the
    // alignment's injected copy of the same SNP are two candidates
    // here and one record after duplicates collapse, so counting
    // candidates skipped every lone het held twice -- chr12:10,999,331
    // -351, four SNPs whose 86 reads split cleanly in one block.
    std::unordered_map<hts_pos_t, int> hets;
    {
        std::unordered_map<hts_pos_t, std::set<hts_pos_t>> at;
        for (const CandidateVariant& c : chunk.candidates)
            if (phased_het(c)) at[c.phase_set].insert(c.key.sort_pos());
        for (const auto& [ps, positions] : at) hets[ps] = static_cast<int>(positions.size());
    }
    std::unordered_map<size_t, std::map<hts_pos_t, std::array<std::map<int, int>, 3>>> votes;
    std::unordered_map<size_t, int> labelled;
    // Two het SNPs at one position with different ALTs, both left
    // unphased: a site with two non-reference alleles (1|2). No read
    // carries REF there, so each record's reads carry its ALT or
    // nothing and the phaser has no two-allele contrast to solve: the
    // catalog site's pairs (REF vs each ALT) look 100% ALT and are
    // dropped (high_af), and the alignment's two records each hear one
    // allele. HLA-B (chr6:31,212,612 G>A,T and 31,212,651) was written
    // unphased where HiPhase phases it. Each is placed below from the
    // final labels of the reads carrying its ALT, so no label, join or
    // block changes. Over chr12/18/20 and the MHC: 8 more truth 1|2
    // sites phased, reads, switches and flips unchanged (the one
    // placement against the truth VCF agrees with its read labels).
    // Indel pairs as well: +2 switches on chr12.
    std::unordered_map<size_t, size_t> partner;
    std::unordered_map<size_t, std::map<hts_pos_t, std::array<int, 3>>> alt_votes;
    {
        const auto unphased_het = [](const CandidateVariant& c) {
            return (c.lcd_var_i_to_cate & (kCandCleanHetSnp | kCandCleanHetIndel)) != 0 && c.phase_set <= 0;
        };
        std::map<hts_pos_t, std::vector<size_t>> at;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const CandidateVariant& c = chunk.candidates[ci];
            // Alignment calls only: a catalog row's key is its site's, not its allele's.
            if (!unphased_het(c) || is_catalog_row(c)) continue;
            if (!(c.key.type == VariantType::Snp && c.key.ref_len == 1 && c.key.alt.size() == 1)) continue;
            at[c.key.pos].push_back(ci);
        }
        for (const auto& [pos, cis] : at)
            if (cis.size() == 2 && chunk.candidates[cis[0]].key.alt != chunk.candidates[cis[1]].key.alt)
                partner[cis[0]] = cis[1], partner[cis[1]] = cis[0];
    }
    for (const ReadVariantProfile& prof : chunk.read_var_profile) {
        if (prof.read_id < 0 || prof.start_var_idx < 0) continue;
        const ChunkReadLabel label = read_label(chunk, static_cast<size_t>(prof.read_id));
        const int hap = label.hap;
        const hts_pos_t ps = label.phase_set;
        if (hap == 0) continue;
        for (size_t k = 0; k < prof.alleles.size(); ++k) {
            const int a = prof.alleles[k];
            if (a < 0) continue;
            const size_t ci = static_cast<size_t>(prof.start_var_idx) + k;
            if (a == 1 && partner.count(ci) != 0) ++alt_votes[ci][ps][static_cast<size_t>(hap)];
            if (ci >= chunk.candidates.size() || !phased_het(chunk.candidates[ci]) ||
                hets[chunk.candidates[ci].phase_set] != 1) continue;
            ++votes[ci][ps][static_cast<size_t>(hap)][a];
            ++labelled[ci];
        }
    }
    for (const auto& [ci, by_ps] : votes) {
        hts_pos_t best = 0;
        int best_n = 0;
        for (const auto& [ps, h] : by_ps) {
            int n = 0;
            for (int hh = 1; hh <= 2; ++hh) for (const auto& [a, c] : h[static_cast<size_t>(hh)]) n += c;
            if (n > best_n) best_n = n, best = ps;
        }
        if (best == chunk.candidates[ci].phase_set || best_n < 3) continue;
        const auto& h = by_ps.at(best);
        int alle[3] = {-1, -1, -1};
        bool ok = true;
        for (int hh = 1; hh <= 2 && ok; ++hh) {
            int tot = 0, top = 0, top_a = -1;
            for (const auto& [a, c] : h[static_cast<size_t>(hh)]) { tot += c; if (c > top) top = c, top_a = a; }
            ok = tot >= 2 && top >= 0.8 * tot;
            alle[hh] = top_a;
        }
        if (!ok || alle[1] == alle[2]) continue;
        placements[ci] = UnlinkedHetPlacement{best, alle[1], alle[2]};
        if (trace_votes != nullptr) {
            std::string v;
            for (int hh = 1; hh <= 2; ++hh) {
                int tot = 0, top = 0;
                for (const auto& [a, c] : h[static_cast<size_t>(hh)]) { tot += c; top = std::max(top, c); }
                v += "votes_h" + std::to_string(hh) + "=" + std::to_string(top) + "/" + std::to_string(tot) + " ";
            }
            v += "ps_n=" + std::to_string(by_ps.size()) + " labelled=" + std::to_string(labelled[ci]);
            (*trace_votes)[ci] = v;
        }
    }
    // Each ALT goes on the haplotype holding 80% of its labelled
    // carriers, in the phase set holding most of them (at least 3);
    // both in one phase set, on opposite haplotypes.
    const auto place = [&](size_t ci, hts_pos_t& ps, int& hap) {
        const auto it = alt_votes.find(ci);
        if (it == alt_votes.end()) return false;
        int best_n = 0;
        for (const auto& [p, h] : it->second)
            if (h[1] + h[2] > best_n) best_n = h[1] + h[2], ps = p;
        if (best_n < 3) return false;
        const auto& h = it->second.at(ps);
        hap = h[1] >= h[2] ? 1 : 2;
        return h[static_cast<size_t>(hap)] >= 0.8 * best_n;
    };
    for (const auto& [a, b] : partner) {
        if (a > b) continue;
        hts_pos_t ps_a = 0, ps_b = 0;
        int hap_a = 0, hap_b = 0;
        if (!place(a, ps_a, hap_a) || !place(b, ps_b, hap_b) || ps_a != ps_b || hap_a == hap_b) continue;
        placements[a] = UnlinkedHetPlacement{ps_a, hap_a == 1 ? 1 : 0, hap_a == 2 ? 1 : 0};
        placements[b] = UnlinkedHetPlacement{ps_b, hap_b == 1 ? 1 : 0, hap_b == 2 ? 1 : 0};
    }
    return placements;
}


const std::string* row_allele_seq(const GraphChunkBuildResult& graph_chunk, size_t ci, size_t k) {
    if (ci >= graph_chunk.site_meta.size() || ci >= graph_chunk.site_allele_orig_idx.size()) return nullptr;
    const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
    const std::vector<int>& orig = graph_chunk.site_allele_orig_idx[ci];
    if (meta.ref.empty() || k >= orig.size()) return nullptr;
    const int walk = orig[k];
    if (walk == 0) return &meta.ref;
    return walk > 0 && static_cast<size_t>(walk) <= meta.alts.size() ? &meta.alts[static_cast<size_t>(walk) - 1]
                                                                     : nullptr;
}

std::optional<RowAllele> row_allele(const GraphChunkBuildResult& graph_chunk, size_t ci) {
    if (ci >= graph_chunk.site_meta.size() || ci >= graph_chunk.site_allele_orig_idx.size()) return std::nullopt;
    const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
    const std::vector<int>& orig = graph_chunk.site_allele_orig_idx[ci];
    if (meta.pos <= 0 || orig.size() < 2 || orig[0] != 0) return std::nullopt;
    const std::string* alt = row_allele_seq(graph_chunk, ci, 1);
    if (alt == nullptr || alt == &meta.ref) return std::nullopt;
    return RowAllele{meta.pos, &meta.ref, alt};
}

std::optional<VcfAllele> row_vcf_allele(const GraphChunkBuildResult& graph_chunk, size_t ci) {
    const std::optional<RowAllele> allele = row_allele(graph_chunk, ci);
    if (!allele) return std::nullopt;
    const auto upper = [](std::string x) {
        for (char& ch : x) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        return x;
    };
    VcfAllele out{allele->pos, upper(*allele->ref), upper(*allele->alt)};
    trim_to_minimal_vcf(out.pos, out.ref, out.alt);
    return out;
}

void rebuild_joint_candidate_loci(GraphChunkBuildResult& graph_chunk,
        const std::function<char(hts_pos_t)>& reference_base) {
    std::map<CandidateIdentityKey, std::vector<size_t>> descriptions;
    std::map<AlleleContrastKey, std::vector<size_t>> contrasts;
    graph_chunk.candidate_allele_contrasts.assign(graph_chunk.chunk.candidates.size(), std::nullopt);
    for (size_t ci = 0; ci < graph_chunk.chunk.candidates.size(); ++ci) {
        if (ci >= graph_chunk.site_meta.size() || ci >= graph_chunk.site_allele_orig_idx.size()) continue;
        const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
        const CandidateVariant& row = graph_chunk.chunk.candidates[ci];
        const auto& original = graph_chunk.site_allele_orig_idx[ci];
        std::array<int, 2> local = {row.hap_to_cons_alle[1], row.hap_to_cons_alle[2]};
        if (local[0] < 0 || local[1] < 0 || local[0] == local[1]) {
            if (row.counts.n_uniq_alles != 2 || original.size() != 2) continue;
            local = {0, 1};
        }
        const auto sequence = [&](int allele) -> const std::string* {
            if (allele < 0 || allele >= row.counts.n_uniq_alles) return nullptr;
            return row_allele_seq(graph_chunk, ci, static_cast<size_t>(allele));
        };
        const std::string* first = sequence(local[0]);
        const std::string* second = sequence(local[1]);
        if (first == nullptr || second == nullptr) continue;
        const auto contrast = normalize_allele_contrast(meta.pos, meta.ref, *first, *second, reference_base);
        if (!contrast) continue;
        if (contrast->reversed) std::swap(local[0], local[1]);
        graph_chunk.candidate_allele_contrasts[ci] = CandidateAlleleContrast{contrast->key, local};
        contrasts[contrast->key].push_back(ci);
        // Keep the established REF/ALT consumers in their original 0/1 gauge.
        if (row.counts.n_uniq_alles == 2 && original.size() == 2 && original[0] == 0 &&
            local == std::array<int, 2>{0, 1} && !contrast->key.alleles[0])
            descriptions[*contrast->key.alleles[1]].push_back(ci);
    }
    graph_chunk.joint_candidate_loci.clear();
    for (auto& [identity, members] : descriptions) {
        if (members.size() > 1)
            graph_chunk.joint_candidate_loci.emplace_back(identity, std::move(members));
    }
    graph_chunk.joint_allele_contrasts.clear();
    for (auto& [contrast, members] : contrasts) {
        if (members.size() > 1)
            graph_chunk.joint_allele_contrasts.emplace_back(contrast, std::move(members));
    }
}

void phase_joint_graph_candidates(GraphChunkBuildResult& graph_chunk, const Options& opts,
        const std::function<char(hts_pos_t)>& reference_base) {
    rebuild_joint_candidate_loci(graph_chunk, reference_base);
    PhasingChunk& chunk = graph_chunk.chunk;
    std::vector<const std::vector<size_t>*> projected;
    for (const auto& [identity, members] : graph_chunk.joint_candidate_loci) {
        (void)identity;
        const CandidateVariant& first = chunk.candidates[members.front()];
        // Identity cannot widen MSA/repeat admission or equate phase-set gauges.
        if ((first.lcd_var_i_to_cate != kCandCleanHetSnp &&
             first.lcd_var_i_to_cate != kCandCleanHetIndel) ||
            !first.msa_insertion_alts.empty() || first.is_homopolymer_indel ||
            first.phase_set > 0) continue;
        const bool compatible = std::all_of(members.begin(), members.end(), [&](size_t ci) {
            const CandidateVariant& row = chunk.candidates[ci];
            return row.lcd_var_i_to_cate == first.lcd_var_i_to_cate &&
                row.counts.category == first.counts.category && row.phase_set <= 0 &&
                row.hap_to_cons_alle[1] < 0 && row.hap_to_cons_alle[2] < 0 &&
                row.msa_insertion_alts.empty() && !row.is_homopolymer_indel &&
                !row.is_homopolymer_indel;
        });
        if (!compatible) continue;
        std::vector<VariantKey> edits;
        for (const size_t ci : members) {
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
            edits.push_back(vcf_to_variant_key(chunk.region.tid, meta.pos, meta.ref,
                *selected_graph_candidate_alt(graph_chunk, ci)));
        }
        const bool covered = std::all_of(chunk.read_var_profile.begin(), chunk.read_var_profile.end(),
            [&](const ReadVariantProfile& profile) {
                for (size_t i = 1; i < members.size(); ++i) {
                    const int offset = static_cast<int>(members[i]) - profile.start_var_idx;
                    if (profile.start_var_idx < 0 || offset < 0 ||
                        static_cast<size_t>(offset) >= profile.alleles.size()) continue;
                    const int allele = profile.alleles[offset];
                    if (allele != 0 && allele != 1) continue;
                    const GraphSiteMeta& source_meta = graph_chunk.site_meta[members[i]];
                    const GraphSiteMeta& parent_meta = graph_chunk.site_meta[members.front()];
                    const bool same_contrast = source_meta.pos == parent_meta.pos &&
                        source_meta.ref == parent_meta.ref &&
                        *selected_graph_candidate_alt(graph_chunk, members[i]) ==
                        *selected_graph_candidate_alt(graph_chunk, members.front());
                    const int parent_offset = static_cast<int>(members.front()) - profile.start_var_idx;
                    const auto& graph_calls = profile.graph_alleles.empty()
                        ? profile.alleles : profile.graph_alleles;
                    const int parent_call = parent_offset >= 0 &&
                        static_cast<size_t>(parent_offset) < graph_calls.size()
                        ? graph_calls[parent_offset] : -1;
                    // A primitive event's REF does not certify the full allele
                    // at a containing snarl, which may have other sample ALTs.
                    if (!same_contrast && parent_call != 0 && parent_call != 1) return false;
                    const ReadRecord& read = chunk.reads[profile.read_id];
                    // A graph walk through a shifted site cannot establish a
                    // complete contrast outside the molecule's physical span.
                    if ((edits[i].pos != edits[0].pos || edits[i].ref_len != edits[0].ref_len) &&
                        !supports_normalized_observation(allele, read.beg, read.end,
                            edits[i].pos, edits[i].ref_len, edits[0].pos, edits[0].ref_len)) return false;
                }
                return true;
            });
        if (covered) projected.push_back(&members);
    }
    if (!opts.phase_matrix_dump_prefix.empty())
        std::fprintf(stderr, "[joint-candidate-loci] duplicates=%zu projected=%zu\n",
                     graph_chunk.joint_candidate_loci.size(), projected.size());
    if (projected.empty()) {
        assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineClean);
        return;
    }
    std::vector<ReadVariantProfile> original_profiles = chunk.read_var_profile;
    std::vector<std::pair<size_t, uint32_t>> original_masks;
    for (const auto* members : projected) {
        for (size_t i = 1; i < members->size(); ++i) {
            CandidateVariant& row = chunk.candidates[(*members)[i]];
            original_masks.emplace_back((*members)[i], row.lcd_var_i_to_cate);
            row.lcd_var_i_to_cate = 0;
        }
        for (ReadVariantProfile& profile : chunk.read_var_profile) {
            MoleculeAlleleEvidence evidence;
            for (const size_t ci : *members) {
                const int offset = static_cast<int>(ci) - profile.start_var_idx;
                if (profile.start_var_idx < 0 || offset < 0 ||
                    static_cast<size_t>(offset) >= profile.alleles.size()) continue;
                evidence.add(ci, profile.alleles[offset],
                    static_cast<size_t>(offset) < profile.alt_qi.size() ? profile.alt_qi[offset] : 0);
            }
            if (evidence.observations().empty()) continue;
            update_read_var_profile_with_allele(static_cast<int>(members->front()),
                evidence.allele() < 0 ? -1 : evidence.allele(), evidence.query_index(), profile);
        }
    }
    const auto restore = [&]() {
        for (const auto& [ci, mask] : original_masks) chunk.candidates[ci].lcd_var_i_to_cate = mask;
        chunk.read_var_profile = std::move(original_profiles);
        rebuild_read_var_cr(chunk);
    };
    try {
        rebuild_read_var_cr(chunk);
        assign_hap_based_on_germline_het_vars_kmeans(chunk, opts, kCandGermlineClean);
        for (const auto* members : projected) {
            const CandidateVariant& phased = chunk.candidates[members->front()];
            for (size_t i = 1; i < members->size(); ++i) {
                CandidateVariant& alias = chunk.candidates[(*members)[i]];
                alias.phase_set = phased.phase_set;
                alias.hap_alt = phased.hap_alt;
                alias.hap_ref = phased.hap_ref;
                alias.hap_to_cons_alle = phased.hap_to_cons_alle;
                alias.hap_to_alle_profile = phased.hap_to_alle_profile;
            }
        }
    } catch (...) {
        restore();
        throw;
    }
    restore();
}

// Overlap reads between adjacent chunks, paired by name. build_graph_chunk
// inserts reads in name order, but the union path appends more after it
// (alignment reads at injected sites, alignment-only reads), and a sorted
// merge missed most pairs involving them: the reads that carry a catalog gap
// across a chunk edge, i.e. the stitch votes that matter there.
static void populate_graph_chunk_pair_overlap_impl(PhasingChunk& pre, PhasingChunk& cur) {
    pre.down_ovlp_read_i.assign(1, {});
    cur.up_ovlp_read_i.assign(1, {});
    pre.n_down_ovlp_skip_reads.assign(1, 0);
    cur.n_up_ovlp_skip_reads.assign(1, 0);

    std::unordered_map<std::string_view, size_t> cur_by_name;
    cur_by_name.reserve(cur.reads.size());
    for (size_t ci = 0; ci < cur.reads.size(); ++ci) cur_by_name.try_emplace(cur.reads[ci].qname, ci);
    for (size_t pi = 0; pi < pre.reads.size(); ++pi) {
        const auto it = cur_by_name.find(pre.reads[pi].qname);
        if (it == cur_by_name.end()) continue;
        pre.down_ovlp_read_i[0].push_back(static_cast<int>(pi));
        cur.up_ovlp_read_i[0].push_back(static_cast<int>(it->second));
    }
}

void populate_graph_chunk_overlaps(std::vector<GraphChunkBuildResult>& graph_chunks) {
    for (size_t i = 1; i < graph_chunks.size(); ++i)
        populate_graph_chunk_pair_overlap_impl(graph_chunks[i - 1].chunk, graph_chunks[i].chunk);
}

namespace {

constexpr double kRescueSiteOrientationPValue = 0.01;
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
        const ChunkReadLabel label = read_label(chunk, read_i);
        const int hap = label.hap;
        const hts_pos_t phase_set = label.phase_set;
        const bool primary_assignment = label.primary;
        if (hap == 0) continue;

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
            if (marker.is_snp) {
                const auto found = site_votes[site_i].find(marker.phase_set);
                marker.singleton_safe = found != site_votes[site_i].end() &&
                    rescue_singleton_has_primary_support(found->second, marker);
            }
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
            supported.singleton_safe =
                rescue_singleton_has_primary_support(entry.second, supported);
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
        if (pgphase_trace::on())
            pgphase_trace::read_trace(chunk, read_i)->fill =
                std::string("rescue obs=") + (use_bam_observations ? "alignment" : "graph") +
                " hap=" + std::to_string(best_hap) + " block=" + std::to_string(best_phase_set) +
                " margin=" + std::to_string(best_margin) + " votes=" + std::to_string(best_total);
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
    const GraphChunkBuildResult& gc) {
    const PhasingChunk& chunk = gc.chunk;
    for (const ReadVariantProfile& profile : chunk.read_var_profile) {
        const size_t read_i = static_cast<size_t>(profile.read_id);
        const ReadRecord& read = chunk.reads[read_i];
        PhaseReadOutputRow& row = rows_by_read[read.qname];
        if (row.read_name.empty()) row.read_name = read.qname;
        if (pgphase_trace::on()) {
            row.beg = row.beg == 0 ? read.beg : std::min(row.beg, read.beg);
            row.end = std::max(row.end, read.end);
        }

        // A fill label is written in its block's gauge (rescue markers, the EM's block, an island's site), under that
        // block's PS: written with kGapFillPsOffset, each block's fill reads formed a phase set of their own, oriented
        // on their own -- over chr12/18/20 it hid 261 arm reads phased against their block.
        const ChunkReadLabel label = read_label(chunk, read_i);
        const int hap = label.hap;
        const hts_pos_t phase_set = label.phase_set;
        const bool is_primary = label.primary;

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
        if (pgphase_trace::on()) {
            row.trace += (row.trace.empty() ? "" : " | ") + std::string("chunk=") +
                         std::to_string(chunk.region.chunk_id) + " hap=" + std::to_string(hap) +
                         " ps=" + std::to_string(phase_set) +
                         (is_primary ? " primary" : hap != 0 ? (fill_only ? " fill_only" : " fill") : " none") +
                         ": " + (read_i < chunk.read_trace.size()
                                     ? pgphase_trace::describe_read(chunk.read_trace[read_i]) : std::string("-"));
        }
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
        if (pgphase_trace::on())
            pgphase_trace::emit_read(kv.second.beg, kv.second.end, kv.first,
                                     kv.second.has_phased_assignment ? kv.second.hap : 0,
                                     kv.second.has_phased_assignment ? kv.second.phase_set : -1, kv.second.copies,
                                     kv.second.trace);
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

} // namespace pgphase_collect
